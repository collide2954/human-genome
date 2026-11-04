export const USER_AGENT = 'human-genome-snapshot (https://github.com/collide2954/human-genome)';

const NCBI_GENOMES = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF';
const GRCH38 = `${NCBI_GENOMES}/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14`;
const T2T = `${NCBI_GENOMES}/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0`;
const UCSC = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38';
const UNIPROT_FIELDS = ['accession', 'gene_primary', 'length', 'protein_existence', 'xref_pdb', 'ft_transmem', 'ft_signal', 'cc_subcellular_location', 'xref_pfam', 'ft_region', 'ft_mod_res', 'ft_carbohyd', 'ft_crosslnk'];
const UNIPROT_QUERY = new URLSearchParams({
  query: 'organism_id:9606 AND reviewed:true', fields: UNIPROT_FIELDS.join(','), format: 'tsv', compressed: 'true',
});

export const PHAROS_QUERY = '{ targets { count facets { facet values { name value } } } }';

const STRING_DIR = 'https://stringdb-downloads.org/download/';

export const stringFiles = (v) => ({
  detailed: `${STRING_DIR}protein.links.detailed.v${v}/9606.protein.links.detailed.v${v}.txt.gz`,
  physical: `${STRING_DIR}protein.physical.links.v${v}/9606.protein.physical.links.v${v}.txt.gz`,
  info: `${STRING_DIR}protein.info.v${v}/9606.protein.info.v${v}.txt.gz`,
});

const REACTOME = 'https://reactome.org/download/current/';
const T2T_ANNOTATION = 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/';
const DBVAR_NR = 'https://ftp.ncbi.nlm.nih.gov/pub/dbVar/sandbox/sv_datasets/nonredundant/';

const CLINVAR = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/';

export const LOLLIPOPS = [
  { gene: 'TP53', accession: 'P04637' },
  { gene: 'BRCA1', accession: 'P38398' },
  { gene: 'CFTR', accession: 'P13569' },
  { gene: 'MYH7', accession: 'P12883' },
];

const HPO = 'https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/';
const ORPHADATA = 'https://www.orphadata.com/data/xml/';
const TRIALS = 'https://clinicaltrials.gov/api/v2/studies';

export const TRIAL_MODALITIES = [
  { key: 'aav', label: 'AAV vectors', query: 'adeno-associated OR AAV' },
  { key: 'lentiviral', label: 'Lentiviral vectors', query: 'lentiviral OR lentivirus' },
  { key: 'retroviral', label: 'Gammaretroviral vectors', query: 'retroviral OR gammaretroviral' },
  { key: 'cart', label: 'CAR T cells', query: '"chimeric antigen receptor" OR CAR-T OR "CAR T"' },
  { key: 'editing', label: 'Gene editing', query: 'CRISPR OR Cas9 OR "base editing" OR "prime editing" OR "zinc finger nuclease" OR TALEN' },
  { key: 'antisense', label: 'Antisense oligonucleotides', query: 'antisense OR oligonucleotide' },
  { key: 'sirna', label: 'siRNA', query: 'siRNA OR "small interfering RNA"' },
];

const ASSEMBLIES = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/9606/dataset_report';

export const assemblyPage = (token = null) => `${ASSEMBLIES}?${new URLSearchParams({ page_size: '1000', ...(token ? { page_token: token } : {}) })}`;

export const trialPage = (query, token = null) => `${TRIALS}?${new URLSearchParams({
  'query.intr': query, 'filter.advanced': 'AREA[StudyType]INTERVENTIONAL', pageSize: '1000', fields: 'NCTId,StartDate,Phase,OverallStatus', ...(token ? { pageToken: token } : {}),
})}`;

export const pharosLevelQuery = (level) => `{ targets(filter: { facets: [{ facet: "Target Development Level", values: ["${level}"] }] }) { count facets(include: ["IDG Family"]) { facet values { name value } } } }`;

export const FILES = {
  grch38Report: `${GRCH38}_assembly_report.txt`,
  grch38Stats: `${GRCH38}_assembly_stats.txt`,
  t2tReport: `${T2T}_assembly_report.txt`,
  t2tStats: `${T2T}_assembly_stats.txt`,
  cytobands: `${UCSC}/database/cytoBandIdeo.txt.gz`,
  gaps: `${UCSC}/database/gap.txt.gz`,
  centromeres: `${UCSC}/database/centromeres.txt.gz`,
  rmsk: `${UCSC}/database/rmsk.txt.gz`,
  superDups: `${UCSC}/database/genomicSuperDups.txt.gz`,
  simpleRepeat: `${UCSC}/database/simpleRepeat.txt.gz`,
  hg38TwoBit: `${UCSC}/bigZips/hg38.2bit`,
  maneDir: 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/',
  gencodeDir: 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/',
  hgnc: 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt',
  geneInfo: 'https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz',
  uniprot: `https://rest.uniprot.org/uniprotkb/stream?${UNIPROT_QUERY}`,
  interproPfam: 'https://www.ebi.ac.uk/interpro/api/entry/pfam/',
  hpa: 'https://www.proteinatlas.org/download/proteinatlas.tsv.zip',
  gtexMedian: 'https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz',
  gtexTissues: 'https://gtexportal.org/api/v2/dataset/tissueSiteDetail',
  pharos: 'https://pharos-api.ncats.io/graphql',
  stringVersion: 'https://string-db.org/api/json/version',
  reactomeVersion: 'https://reactome.org/ContentService/data/database/version',
  reactomePathways: `${REACTOME}ReactomePathways.txt`,
  reactomeRelations: `${REACTOME}ReactomePathwaysRelation.txt`,
  reactomeUniprot: `${REACTOME}UniProt2Reactome_All_Levels.txt`,
  goaHuman: 'https://current.geneontology.org/annotations/goa_human.gaf.gz',
  goBasic: 'https://current.geneontology.org/ontology/go-basic.obo',
  mitoFeatures: 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=ft&retmode=text',
  gnomadMito: 'https://storage.googleapis.com/gcp-public-data--gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv',
  censat: `${T2T_ANNOTATION}chm13v2.0_censat_v2.1.bed`,
  telomeres: `${T2T_ANNOTATION}chm13v2.0_telomere.bed`,
  dbsnpNotes: 'https://ftp.ncbi.nlm.nih.gov/snp/latest_release/release_notes.txt',
  nr_deletions: `${DBVAR_NR}deletions/GRCh38.nr_deletions.tsv.gz`,
  nr_duplications: `${DBVAR_NR}duplications/GRCh38.nr_duplications.tsv.gz`,
  nr_insertions: `${DBVAR_NR}insertions/GRCh38.nr_insertions.tsv.gz`,
  gnomadConstraint: 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv',
  clingenGenes: 'https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv',
  impcViability: 'https://ftp.ebi.ac.uk/pub/databases/impc/all-data-releases/latest/results/viability.csv.gz',
  clinvarVcf: `${CLINVAR}vcf_GRCh38/clinvar.vcf.gz`,
  clinvarArchive: `${CLINVAR}vcf_GRCh38/archive_2.0/`,
  clinvarTable: `${CLINVAR}tab_delimited/variant_summary.txt.gz`,
  clinvarGenes: `${CLINVAR}tab_delimited/gene_specific_summary.txt`,
  clinvarOrgs: `${CLINVAR}tab_delimited/organization_summary.txt`,
  clingenValidity: 'https://search.clinicalgenome.org/kb/gene-validity/download',
  clingenRegions: 'https://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv',
  interproProtein: 'https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/',
  strchive: 'https://raw.githubusercontent.com/dashnowlab/STRchive/main/data/STRchive-loci.json',
  hpoAnnotations: `${HPO}phenotype.hpoa`,
  hpoOntology: `${HPO}hp.obo`,
  hpoGenes: `${HPO}genes_to_disease.txt`,
  orphaGenes: `${ORPHADATA}en_product6.xml`,
  orphaPrevalence: `${ORPHADATA}en_product9_prev.xml`,
  orphaAges: `${ORPHADATA}en_product9_ages.xml`,
  gtrTests: 'https://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_version.gz',
  fdaCgt: 'https://www.fda.gov/vaccines-blood-biologics/cellular-gene-therapy-products/approved-cellular-and-gene-therapy-products',
  trials: TRIALS,
  gwasAncestry: 'https://www.ebi.ac.uk/gwas/api/search/downloads/ancestry',
  aadrDataset: `https://dataverse.harvard.edu/api/datasets/:persistentId/?${new URLSearchParams({ persistentId: 'doi:10.7910/DVN/FFIDCW' })}`,
  dataverseFile: 'https://dataverse.harvard.edu/api/access/datafile/',
  igsrSamples: 'https://www.internationalgenome.org/api/beta/sample/_search',
  strchiveCommits: `https://api.github.com/repos/dashnowlab/STRchive/commits?${new URLSearchParams({ path: 'data/STRchive-loci.json', per_page: '1' })}`,
  ccre: 'https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed',
  cpgIslands: `${UCSC}/database/cpgIslandExt.txt.gz`,
  trnas: 'https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa',
  encode: `https://www.encodeproject.org/search/?${new URLSearchParams({
    type: 'Experiment', status: 'released', 'replicates.library.biosample.donor.organism.scientific_name': 'Homo sapiens', limit: '0', format: 'json',
  })}`,
};
