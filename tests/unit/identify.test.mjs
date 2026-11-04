import { test } from 'node:test';
import assert from 'node:assert/strict';
import { identify, bandInterval, hgvsChange } from '../../src/lib.js';

const part = (hit, key) => hit?.parts.find(([k]) => k === key)?.[1];

test('Ensembl stable IDs decode by feature type and keep their version', () => {
  const gene = identify('ENSG00000141510');
  assert.equal(gene.type, 'Ensembl gene');
  assert.equal(gene.url, 'https://www.ensembl.org/id/ENSG00000141510');
  const tx = identify('enst00000269305.9');
  assert.equal(tx.type, 'Ensembl transcript');
  assert.equal(tx.id, 'ENST00000269305.9');
  assert.equal(part(tx, 'Version'), '9');
  assert.equal(tx.url, 'https://www.ensembl.org/id/ENST00000269305');
  assert.equal(identify('ENSP00000269305').type, 'Ensembl protein');
  assert.equal(identify('ENSE00001146308').type, 'Ensembl exon');
  assert.equal(identify('ENSR00000082383').type, 'Ensembl regulatory feature');
  assert.equal(identify('ENSG0000014151'), null);
});

test('RefSeq accessions decode by prefix', () => {
  const nm = identify('NM_000546.6');
  assert.equal(nm.type, 'RefSeq mRNA');
  assert.equal(part(nm, 'Prefix'), 'NM_: curated mRNA');
  assert.equal(part(nm, 'Version'), '6');
  assert.equal(nm.url, 'https://www.ncbi.nlm.nih.gov/nuccore/NM_000546.6');
  const np = identify('NP_000537.3');
  assert.equal(np.type, 'RefSeq protein');
  assert.equal(np.url, 'https://www.ncbi.nlm.nih.gov/protein/NP_000537.3');
  assert.equal(identify('XM_011523961.3').type, 'RefSeq model mRNA');
  assert.equal(identify('NR_003051.4').type, 'RefSeq noncoding RNA');
  assert.equal(identify('NG_017013.2').type, 'RefSeqGene region');
});

test('RefSeq chromosome accessions name the chromosome and assembly', () => {
  const chr17 = identify('NC_000017.11');
  assert.equal(chr17.type, 'RefSeq chromosome');
  assert.equal(part(chr17, 'Molecule'), 'Chromosome 17');
  assert.equal(part(chr17, 'Assembly'), 'GRCh38');
  assert.equal(part(identify('NC_000017.10'), 'Assembly'), 'GRCh37');
  assert.equal(part(identify('NC_000023.11'), 'Molecule'), 'Chromosome X');
  assert.equal(part(identify('NC_000024.10'), 'Molecule'), 'Chromosome Y');
  assert.equal(part(identify('NC_012920.1'), 'Molecule'), 'Mitochondrion (rCRS)');
  assert.equal(part(identify('NC_060941.1'), 'Assembly'), 'T2T-CHM13v2.0');
  assert.equal(part(identify('NC_060941.1'), 'Molecule'), 'Chromosome 17');
  assert.equal(part(identify('NC_000017.9'), 'Assembly'), undefined);
});

test('variant IDs from dbSNP, ClinVar, dbVar, ClinGen and COSMIC resolve', () => {
  assert.equal(identify('RS334').id, 'rs334');
  assert.equal(identify('rs334').url, 'https://www.ncbi.nlm.nih.gov/snp/rs334');
  const vcv = identify('VCV000017661.80');
  assert.equal(vcv.type, 'ClinVar variation');
  assert.equal(vcv.url, 'https://www.ncbi.nlm.nih.gov/clinvar/variation/17661/');
  assert.equal(identify('RCV000013601.32').url, 'https://www.ncbi.nlm.nih.gov/clinvar/RCV000013601/');
  assert.equal(identify('SCV000255967').type, 'ClinVar submission');
  assert.equal(identify('nsv1067916').type, 'dbVar variant region');
  assert.equal(identify('nssv13638').type, 'dbVar variant call');
  assert.equal(identify('nstd102').url, 'https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd102/');
  assert.equal(identify('CA000330').type, 'ClinGen canonical allele');
  assert.equal(identify('COSV52661282').type, 'COSMIC genomic variant');
});

test('gene, disease and phenotype IDs resolve with links', () => {
  assert.equal(identify('HGNC:11998').url, 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:11998');
  assert.equal(identify('GeneID:7157').url, 'https://www.ncbi.nlm.nih.gov/gene/7157');
  assert.equal(identify('NCBIGene:7157').type, 'NCBI Gene');
  const mim = identify('MIM:191170');
  assert.equal(mim.type, 'OMIM entry');
  assert.equal(mim.url, 'https://www.omim.org/entry/191170');
  assert.equal(part(mim, 'Number'), 'Autosomal, entered before May 1994');
  assert.equal(part(identify('MIM:602421'), 'Number'), 'Autosomal, entered since May 1994');
  const cf = identify('#219700');
  assert.equal(part(cf, 'Prefix'), '#: phenotype, molecular basis known');
  assert.equal(part(cf, 'Number'), 'Autosomal, entered before May 1994');
  assert.equal(part(identify('OMIM:300377'), 'Number'), 'X-linked');
  assert.equal(identify('ORPHA:558').url, 'https://www.orpha.net/en/disease/detail/558');
  assert.equal(identify('MONDO:0007254').type, 'MONDO disease');
  assert.equal(identify('HP:0001250').url, 'https://hpo.jax.org/browse/term/HP:0001250');
  assert.equal(identify('GO:0006915').type, 'Gene Ontology term');
  assert.equal(identify('C0010674').type, 'MedGen concept');
  assert.equal(identify('NBK1247').url, 'https://www.ncbi.nlm.nih.gov/books/NBK1247/');
});

test('protein IDs from UniProt and PDB resolve, years do not', () => {
  const p53 = identify('P04637');
  assert.equal(p53.type, 'UniProtKB accession');
  assert.equal(p53.url, 'https://www.uniprot.org/uniprotkb/P04637/entry');
  assert.equal(part(identify('P04637-2'), 'Isoform'), '2');
  assert.equal(identify('A0A024R1R8').type, 'UniProtKB accession');
  assert.equal(identify('P53_HUMAN').type, 'UniProtKB entry name');
  assert.equal(identify('1tup').id, '1TUP');
  assert.equal(identify('1TUP').url, 'https://www.rcsb.org/structure/1TUP');
  assert.equal(identify('pdb_00001tup').type, 'PDB entry');
  assert.equal(identify('1990'), null);
});

test('study, literature, pathway and assembly IDs resolve', () => {
  assert.equal(identify('GCST000001').url, 'https://www.ebi.ac.uk/gwas/studies/GCST000001');
  assert.equal(identify('EFO_0000270').id, 'EFO:0000270');
  assert.equal(identify('PMID:11237011').url, 'https://pubmed.ncbi.nlm.nih.gov/11237011/');
  assert.equal(identify('ENCSR000AKB').type, 'ENCODE experiment');
  assert.equal(identify('EH38E1516972').type, 'ENCODE candidate cis-regulatory element');
  assert.equal(identify('R-HSA-109582').url, 'https://reactome.org/content/detail/R-HSA-109582');
  assert.equal(identify('CCDS11118.1').type, 'Consensus CDS');
  assert.equal(identify('LRG_292').type, 'Locus Reference Genomic');
  assert.equal(identify('DB00316').type, 'DrugBank drug');
  const asm = identify('GCF_000001405.40');
  assert.equal(asm.type, 'NCBI assembly');
  assert.equal(part(asm, 'Name'), 'GRCh38.p14');
  assert.equal(part(identify('GCA_009914755.4'), 'Name'), 'T2T-CHM13v2.0');
});

test('star alleles and HLA alleles split into their fields', () => {
  const star = identify('cyp2d6*4');
  assert.equal(star.type, 'PharmVar star allele');
  assert.equal(star.id, 'CYP2D6*4');
  assert.equal(part(star, 'Gene'), 'CYP2D6');
  assert.equal(star.url, 'https://www.pharmvar.org/gene/CYP2D6');
  assert.equal(part(identify('CYP2D6*4.001'), 'Suballele'), '001');
  assert.equal(part(identify('DPYD*2A'), 'Allele'), '*2A');
  const hla = identify('HLA-B*57:01');
  assert.equal(hla.type, 'HLA allele');
  assert.equal(part(hla, 'Gene'), 'HLA-B');
  assert.equal(part(hla, 'Allele group'), '57');
  assert.equal(part(hla, 'Protein'), '01');
  assert.equal(part(identify('HLA-A*24:09N'), 'Expression'), 'N: null, not expressed');
  assert.equal(part(identify('HLA-DRB1*15:01:01:02'), 'Noncoding'), '02');
});

test('HGVS names decode reference, coordinate system and change', () => {
  const sub = identify('NM_000546.6:c.215C>G');
  assert.equal(sub.type, 'HGVS variant');
  assert.equal(part(sub, 'Reference'), 'NM_000546.6, RefSeq mRNA');
  assert.equal(part(sub, 'Coordinates'), 'c. coding DNA, numbered from the A of the ATG start codon');
  assert.equal(part(sub, 'Change'), 'substitution');
  assert.equal(sub.url, 'https://mutalyzer.nl/normalizer/NM_000546.6%3Ac.215C%3EG');
  assert.equal(part(identify('NC_012920.1:m.3243A>G'), 'Coordinates'), 'm. mitochondrial DNA');
  assert.equal(part(identify('NM_000492.4:c.1521_1523del'), 'Change'), 'deletion');
  assert.equal(part(identify('NP_000537.3:p.(Pro72Arg)'), 'Coordinates'), 'p. protein, predicted');
  assert.equal(part(identify('NC_000017.11:g.7676154G>C'), 'Reference'), 'NC_000017.11, chromosome 17, GRCh38');
});

test('hgvsChange names the kind of change and where it falls', () => {
  assert.deepEqual(hgvsChange('c', '215C>G'), { kind: 'substitution', region: 'coding sequence' });
  assert.deepEqual(hgvsChange('c', '-14C>T'), { kind: 'substitution', region: '5′ UTR, before the start codon' });
  assert.deepEqual(hgvsChange('c', '*32G>A'), { kind: 'substitution', region: '3′ UTR, after the stop codon' });
  assert.deepEqual(hgvsChange('c', '123+1G>A'), { kind: 'substitution', region: 'intron, next to an exon' });
  assert.deepEqual(hgvsChange('c', '68_69del'), { kind: 'deletion', region: 'coding sequence' });
  assert.deepEqual(hgvsChange('c', '5266dup'), { kind: 'duplication', region: 'coding sequence' });
  assert.deepEqual(hgvsChange('c', '1A>G'), { kind: 'substitution', region: 'start codon' });
  assert.equal(hgvsChange('c', '76_78delinsTT').kind, 'deletion-insertion');
  assert.equal(hgvsChange('p', 'Glu23ValfsTer17').kind, 'frameshift');
  assert.equal(hgvsChange('p', 'Trp53Ter').kind, 'nonsense');
  assert.equal(hgvsChange('p', 'W53*').kind, 'nonsense');
  assert.equal(hgvsChange('p', 'Leu54=').kind, 'synonymous');
  assert.equal(hgvsChange('p', 'Arg72Pro').kind, 'missense');
});

test('coordinates, bands and gnomAD and SPDI variant IDs decode', () => {
  const span = identify('chr17:7,661,779-7,687,538');
  assert.equal(span.type, 'GRCh38 interval');
  assert.equal(part(span, 'Length'), '25.8 kb');
  assert.equal(span.url, 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A7661779-7687538');
  assert.equal(part(identify('7:117559590'), 'Length'), '1 bp');
  assert.equal(identify('chr17:200-100'), null);
  const band = identify('17p13.1');
  assert.equal(band.type, 'Cytogenetic band');
  assert.equal(part(band, 'Arm'), 'p, short arm');
  assert.equal(part(band, 'Region'), '1');
  assert.equal(part(band, 'Band'), '3');
  assert.equal(part(band, 'Sub-band'), '1');
  const gnomad = identify('17-7676154-G-C');
  assert.equal(gnomad.type, 'gnomAD variant');
  assert.equal(gnomad.url, 'https://gnomad.broadinstitute.org/variant/17-7676154-G-C?dataset=gnomad_r4');
  const spdi = identify('NC_000017.11:7676153:G:C');
  assert.equal(spdi.type, 'SPDI variant');
  assert.equal(part(spdi, 'Position'), '7,676,153 counted from 0, so 7,676,154 counted from 1');
});

test('bandInterval spans every sub-band under a band name', () => {
  const bands = [
    { chrom: '17', start: 0, end: 3400000, name: 'p13.3', stain: 'gneg' },
    { chrom: '17', start: 3400000, end: 6500000, name: 'p13.2', stain: 'gpos50' },
    { chrom: '17', start: 6500000, end: 10800000, name: 'p13.1', stain: 'gneg' },
    { chrom: '17', start: 10800000, end: 16100000, name: 'p12', stain: 'gpos75' },
  ];
  assert.deepEqual(bandInterval(bands, '17', 'p13'), { chrom: '17', start: 0, end: 10800000 });
  assert.deepEqual(bandInterval(bands, '17', 'p13.1'), { chrom: '17', start: 6500000, end: 10800000 });
  assert.equal(bandInterval(bands, '17', 'q21'), null);
  const hit = identify('17p13.1', { bands });
  assert.equal(part(hit, 'GRCh38'), 'chr17:6,500,001-10,800,000');
  assert.equal(hit.url, 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=chr17%3A6500001-10800000');
});

test('empty and unknown input gives no match', () => {
  assert.equal(identify(''), null);
  assert.equal(identify('   '), null);
  assert.equal(identify('not an id'), null);
});
