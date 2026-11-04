import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  FILES, USER_AGENT, PHAROS_QUERY, pharosLevelQuery,
} from '../../src/queries.js';

test('every data file is an https URL on a known host', () => {
  const hosts = new Set([
    'ftp.ncbi.nlm.nih.gov', 'hgdownload.soe.ucsc.edu', 'ftp.ebi.ac.uk', 'storage.googleapis.com', 'search.clinicalgenome.org',
    'ftp.clinicalgenome.org', 'github.com', 'www.proteinatlas.org', 'downloads.wenglab.org', 'rest.uniprot.org', 'api.cpicpgx.org',
    'eutils.ncbi.nlm.nih.gov', 'api.ncbi.nlm.nih.gov', 'rest.ensembl.org', 'www.ebi.ac.uk', 'gtexportal.org', 'pharos-api.ncats.io',
  ]);
  for (const [key, url] of Object.entries(FILES)) {
    const u = new URL(url);
    assert.equal(u.protocol, 'https:', key);
    assert.ok(hosts.has(u.hostname), `${key}: ${u.hostname}`);
  }
  assert.equal(new Set(Object.values(FILES)).size, Object.keys(FILES).length);
});

test('assembly files point at GRCh38.p14 and T2T-CHM13v2.0', () => {
  assert.match(FILES.grch38Report, /GCF_000001405\.40_GRCh38\.p14_assembly_report\.txt$/);
  assert.match(FILES.grch38Stats, /GCF_000001405\.40_GRCh38\.p14_assembly_stats\.txt$/);
  assert.match(FILES.t2tReport, /GCF_009914755\.1_T2T-CHM13v2\.0_assembly_report\.txt$/);
  assert.match(FILES.t2tStats, /GCF_009914755\.1_T2T-CHM13v2\.0_assembly_stats\.txt$/);
  for (const key of ['cytobands', 'gaps', 'centromeres', 'rmsk', 'superDups', 'simpleRepeat']) assert.match(FILES[key], /^https:\/\/hgdownload\.soe\.ucsc\.edu\/goldenPath\/hg38\/database\/\w+\.txt\.gz$/, key);
});

test('sequence files point at the hg38 2bit and the current MANE release folder', () => {
  assert.equal(FILES.hg38TwoBit, 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit');
  assert.equal(FILES.maneDir, 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/');
});

test('gene annotation comes from the latest GENCODE release folder', () => {
  assert.equal(FILES.gencodeDir, 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/');
});

test('gene catalogs point at HGNC and NCBI Gene', () => {
  assert.equal(FILES.hgnc, 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt');
  assert.equal(FILES.geneInfo, 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz');
});

test('protein files ask UniProt for every reviewed human entry with the fields the summary reads', () => {
  const u = new URL(FILES.uniprot);
  assert.equal(`${u.origin}${u.pathname}`, 'https://rest.uniprot.org/uniprotkb/stream');
  assert.equal(u.searchParams.get('query'), 'organism_id:9606 AND reviewed:true');
  assert.equal(u.searchParams.get('format'), 'tsv');
  assert.equal(u.searchParams.get('compressed'), 'true');
  assert.deepEqual(u.searchParams.get('fields').split(','), ['accession', 'gene_primary', 'length', 'protein_existence', 'xref_pdb', 'ft_transmem', 'ft_signal', 'cc_subcellular_location', 'xref_pfam', 'ft_region', 'ft_mod_res', 'ft_carbohyd', 'ft_crosslnk']);
  assert.equal(FILES.interproPfam, 'https://www.ebi.ac.uk/interpro/api/entry/pfam/');
});

test('expression and target files point at HPA, GTEx v10 and Pharos', () => {
  assert.equal(FILES.hpa, 'https://www.proteinatlas.org/download/proteinatlas.tsv.zip');
  assert.match(FILES.gtexMedian, /^https:\/\/storage\.googleapis\.com\/adult-gtex\/bulk-gex\/v10\/rna-seq\/GTEx_Analysis_v10_\S+_gene_median_tpm\.gct\.gz$/);
  assert.equal(FILES.gtexTissues, 'https://gtexportal.org/api/v2/dataset/tissueSiteDetail');
  assert.equal(FILES.pharos, 'https://pharos-api.ncats.io/graphql');
  assert.match(PHAROS_QUERY, /targets\s*\{\s*count\s+facets\s*\{\s*facet\s+values\s*\{\s*name\s+value/);
  const level = pharosLevelQuery('Tclin');
  assert.match(level, /facet: "Target Development Level", values: \["Tclin"\]/);
  assert.match(level, /facets\(include: \["IDG Family"\]\)/);
});

test('requests identify the project', () => {
  assert.match(USER_AGENT, /human-genome/);
  assert.doesNotMatch(USER_AGENT, /@/);
});
