import { test } from 'node:test';
import assert from 'node:assert/strict';
import { FILES, USER_AGENT } from '../../src/queries.js';

test('every data file is an https URL on a known host', () => {
  const hosts = new Set([
    'ftp.ncbi.nlm.nih.gov', 'hgdownload.soe.ucsc.edu', 'ftp.ebi.ac.uk', 'storage.googleapis.com', 'search.clinicalgenome.org',
    'ftp.clinicalgenome.org', 'github.com', 'www.proteinatlas.org', 'downloads.wenglab.org', 'rest.uniprot.org', 'api.cpicpgx.org',
    'eutils.ncbi.nlm.nih.gov', 'api.ncbi.nlm.nih.gov', 'rest.ensembl.org',
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

test('requests identify the project', () => {
  assert.match(USER_AGENT, /human-genome/);
  assert.doesNotMatch(USER_AGENT, /@/);
});
