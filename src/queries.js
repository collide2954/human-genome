export const USER_AGENT = 'human-genome-snapshot (https://github.com/collide2954/human-genome)';

const NCBI_GENOMES = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF';
const GRCH38 = `${NCBI_GENOMES}/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14`;
const T2T = `${NCBI_GENOMES}/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0`;
const UCSC = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38';

export const FILES = {
  grch38Report: `${GRCH38}_assembly_report.txt`,
  grch38Stats: `${GRCH38}_assembly_stats.txt`,
  t2tReport: `${T2T}_assembly_report.txt`,
  t2tStats: `${T2T}_assembly_stats.txt`,
  cytobands: `${UCSC}/database/cytoBandIdeo.txt.gz`,
  gaps: `${UCSC}/database/gap.txt.gz`,
  centromeres: `${UCSC}/database/centromeres.txt.gz`,
};
