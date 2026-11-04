import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  FILES, USER_AGENT, PHAROS_QUERY, pharosLevelQuery, stringFiles, LOLLIPOPS,
} from '../../src/queries.js';

test('network, pathway and ontology files point at STRING, Reactome and the GO release', () => {
  assert.equal(FILES.stringVersion, 'https://string-db.org/api/json/version');
  assert.deepEqual(stringFiles('12.5'), {
    detailed: 'https://stringdb-downloads.org/download/protein.links.detailed.v12.5/9606.protein.links.detailed.v12.5.txt.gz',
    physical: 'https://stringdb-downloads.org/download/protein.physical.links.v12.5/9606.protein.physical.links.v12.5.txt.gz',
    info: 'https://stringdb-downloads.org/download/protein.info.v12.5/9606.protein.info.v12.5.txt.gz',
  });
  assert.equal(FILES.reactomeVersion, 'https://reactome.org/ContentService/data/database/version');
  for (const key of ['reactomePathways', 'reactomeRelations', 'reactomeUniprot']) assert.match(FILES[key], /^https:\/\/reactome\.org\/download\/current\/\w+\.txt$/, key);
  assert.match(FILES.reactomeUniprot, /UniProt2Reactome_All_Levels\.txt$/);
  assert.equal(FILES.goaHuman, 'https://current.geneontology.org/annotations/goa_human.gaf.gz');
  assert.equal(FILES.goBasic, 'https://current.geneontology.org/ontology/go-basic.obo');
});

test('every data file is an https URL on a known host', () => {
  const hosts = new Set([
    'ftp.ncbi.nlm.nih.gov', 'hgdownload.soe.ucsc.edu', 'ftp.ebi.ac.uk', 'storage.googleapis.com', 'search.clinicalgenome.org',
    'ftp.clinicalgenome.org', 'github.com', 'www.proteinatlas.org', 'downloads.wenglab.org', 'rest.uniprot.org', 'api.cpicpgx.org',
    'eutils.ncbi.nlm.nih.gov', 'api.ncbi.nlm.nih.gov', 'rest.ensembl.org', 'www.ebi.ac.uk', 'gtexportal.org', 'pharos-api.ncats.io',
    'string-db.org', 'reactome.org', 'current.geneontology.org', 'gtrnadb.ucsc.edu', 'www.encodeproject.org', 's3-us-west-2.amazonaws.com',
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

test('regulation files point at the ENCODE cCRE registry, UCSC CpG islands, GtRNAdb and the ENCODE portal', () => {
  assert.equal(FILES.ccre, 'https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed');
  assert.equal(FILES.cpgIslands, 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz');
  assert.equal(FILES.trnas, 'https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa');
  const e = new URL(FILES.encode);
  assert.equal(`${e.origin}${e.pathname}`, 'https://www.encodeproject.org/search/');
  assert.equal(e.searchParams.get('type'), 'Experiment');
  assert.equal(e.searchParams.get('status'), 'released');
  assert.equal(e.searchParams.get('replicates.library.biosample.donor.organism.scientific_name'), 'Homo sapiens');
  assert.equal(e.searchParams.get('limit'), '0');
  assert.equal(e.searchParams.get('format'), 'json');
});

test('structure files point at the rCRS feature table, gnomAD mtDNA and the T2T annotations', () => {
  const ft = new URL(FILES.mitoFeatures);
  assert.equal(`${ft.origin}${ft.pathname}`, 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi');
  assert.equal(ft.searchParams.get('id'), 'NC_012920.1');
  assert.equal(ft.searchParams.get('rettype'), 'ft');
  assert.equal(FILES.gnomadMito, 'https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv');
  assert.match(FILES.censat, /^https:\/\/s3-us-west-2\.amazonaws\.com\/human-pangenomics\/T2T\/CHM13\/assemblies\/annotation\/chm13v2\.0_censat_v[\d.]+\.bed$/);
  assert.equal(FILES.telomeres, 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_telomere.bed');
});

test('variation and constraint files point at dbSNP, dbVar, gnomAD v4.1, ClinGen and IMPC', () => {
  assert.equal(FILES.dbsnpNotes, 'https://ftp.ncbi.nlm.nih.gov/snp/latest_release/release_notes.txt');
  for (const t of ['deletions', 'duplications', 'insertions']) assert.equal(FILES[`nr_${t}`], `https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/${t}/GRCh38.nr_${t}.tsv.gz`, t);
  assert.equal(FILES.gnomadConstraint, 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv');
  assert.equal(FILES.clingenGenes, 'https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv');
  assert.equal(FILES.impcViability, 'https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/viability.csv.gz');
});

test('clinical files point at ClinVar, its archive, ClinGen validity and regions and InterPro', () => {
  const clinvar = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/';
  assert.equal(FILES.clinvarVcf, `${clinvar}vcf_GRCh38/clinvar.vcf.gz`);
  assert.equal(FILES.clinvarArchive, `${clinvar}vcf_GRCh38/archive_2.0/`);
  assert.equal(FILES.clinvarTable, `${clinvar}tab_delimited/variant_summary.txt.gz`);
  assert.equal(FILES.clinvarGenes, `${clinvar}tab_delimited/gene_specific_summary.txt`);
  assert.equal(FILES.clinvarOrgs, `${clinvar}tab_delimited/organization_summary.txt`);
  assert.equal(FILES.clingenValidity, 'https://search.clinicalgenome.org/kb/gene-validity/download');
  assert.equal(FILES.clingenRegions, 'https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv');
  assert.equal(FILES.interproProtein, 'https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/');
});

test('lollipop proteins name their UniProt accessions', () => {
  assert.deepEqual(LOLLIPOPS, [
    { gene: 'TP53', accession: 'P04637' },
    { gene: 'BRCA1', accession: 'P38398' },
    { gene: 'CFTR', accession: 'P13569' },
    { gene: 'MYH7', accession: 'P12883' },
  ]);
});

test('requests identify the project', () => {
  assert.match(USER_AGENT, /human-genome/);
  assert.doesNotMatch(USER_AGENT, /@/);
});
