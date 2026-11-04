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
  ccre: 'https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed',
  cpgIslands: `${UCSC}/database/cpgIslandExt.txt.gz`,
  trnas: 'https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-mature-tRNAs.fa',
  encode: `https://www.encodeproject.org/search/?${new URLSearchParams({
    type: 'Experiment', status: 'released', 'replicates.library.biosample.donor.organism.scientific_name': 'Homo sapiens', limit: '0', format: 'json',
  })}`,
};
