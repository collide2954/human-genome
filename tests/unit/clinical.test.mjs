import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  germlineClass, consequenceOf, starsOf, CHANNELS, channelOf, clinvarVcf, vcfClassCounts, earliestFile,
  geneSummary, organizations, proteinChange, variantTable, pfamDomains, validitySummary, supportedGenes, regionSummary, repeatLoci,
} from '../../scripts/clinical.mjs';

const locus = (fields) => ({
  id: 'X', disease_id: 'X', gene: 'X', disease: 'x', evidence: ['Definitive'], inheritance: ['AD'], location_in_gene: null, mechanism: null, year: null, benign_min: null, benign_max: null, intermediate_min: null, intermediate_max: null, pathogenic_min: null, pathogenic_max: null, reference_motif_reference_orientation: [], pathogenic_motif_reference_orientation: [], motif_len: 3, typ_age_onset_min: null, typ_age_onset_max: null, genereviews: [], ...fields,
});

test('repeatLoci reads STRchive ranges, motifs, gene locations, discovery years and onset ages', () => {
  const r = repeatLoci([
    locus({
      id: 'HD_HTT', disease_id: 'HD', gene: 'HTT', disease: 'Huntington disease', location_in_gene: 'Coding Exon 1', mechanism: 'GoF', year: '1993 [@pmid:8458085]', benign_min: 6, benign_max: 26, intermediate_min: 27, intermediate_max: 35, pathogenic_min: 36, pathogenic_max: 250, reference_motif_reference_orientation: ['CAG'], pathogenic_motif_reference_orientation: ['CAG'], typ_age_onset_min: 35, typ_age_onset_max: 44, genereviews: ['NBK1305'],
    }),
    locus({
      id: 'CANVAS_RFC1', disease_id: 'CANVAS', gene: 'RFC1', disease: 'CANVAS', evidence: ['Definitive'], inheritance: ['AR'], location_in_gene: 'Intron 2', year: 'Causation identified in 2019 [@pmid:31230722]; clinical association found in 2004', benign_min: 0, benign_max: 11, pathogenic_min: 400, pathogenic_max: 2750, reference_motif_reference_orientation: ['AAAAG'], pathogenic_motif_reference_orientation: ['AAGGG', 'ACAGG'], motif_len: 5,
    }),
    locus({
      id: 'MODY8_CEL', disease_id: 'MODY8', gene: 'CEL', disease: 'Maturity-onset diabetes', evidence: ['Provisional'], location_in_gene: 'Exon 11', year: '2005', motif_len: 33,
    }),
    locus({
      id: 'OPDM2_GIPC1', gene: 'GIPC1', location_in_gene: "5' UTR", year: null,
    }),
  ]);
  assert.deepEqual(r.evidence, { Definitive: 3, Provisional: 1 });
  assert.deepEqual(r.loci[0], {
    id: 'HD_HTT', gene: 'HTT', disease: 'HD', name: 'Huntington disease', evidence: 'Definitive', inheritance: 'AD', motif: 'CAG', motifLength: 3, location: 'coding exon', mechanism: 'GoF', year: 1993, benign: [6, 26], intermediate: [27, 35], pathogenic: [36, 250], onset: [35, 44], genereviews: 'NBK1305',
  });
  assert.equal(r.loci[1].motif, 'AAGGG');
  assert.equal(r.loci[1].year, 2019);
  assert.equal(r.loci[1].location, 'intron');
  assert.equal(r.loci[1].intermediate, null);
  assert.equal(r.loci[1].onset, null);
  assert.equal(r.loci[1].genereviews, null);
  assert.equal(r.loci[2].location, 'exon');
  assert.equal(r.loci[2].pathogenic, null);
  assert.equal(r.loci[2].year, 2005);
  assert.equal(r.loci[3].location, "5' UTR");
  assert.equal(r.loci[3].year, null);
});

test('germlineClass folds VCF and table spellings into nine classes', () => {
  const cases = {
    Pathogenic: 'Pathogenic',
    'Likely_pathogenic,_low_penetrance': 'Likely pathogenic',
    'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance': 'Pathogenic/Likely pathogenic',
    'Pathogenic; risk factor': 'Pathogenic',
    'Benign/Likely_benign|other': 'Benign/Likely benign',
    'Likely benign': 'Likely benign',
    'Uncertain_significance/VUS-high': 'Uncertain significance',
    'VUS-low': 'Uncertain significance',
    'Conflicting classifications of pathogenicity': 'Conflicting',
    'Conflicting_interpretations_of_pathogenicity|risk_factor': 'Conflicting',
    'conflicting data from submitters': 'Conflicting',
    'Uncertain_significance/Uncertain_risk_allele': 'Other',
    drug_response: 'Other',
    'not provided': 'Other',
  };
  for (const [input, expected] of Object.entries(cases)) assert.equal(germlineClass(input), expected, input);
  for (const none of ['no_classification_for_the_single_variant', 'no classifications from unflagged records', '-', '', undefined]) assert.equal(germlineClass(none), null, String(none));
});

test('consequenceOf picks the most severe molecular consequence in VEP order', () => {
  assert.equal(consequenceOf('SO:0001583|missense_variant,SO:0001627|intron_variant'), 'missense');
  assert.equal(consequenceOf('SO:0001587|nonsense,SO:0001575|splice_donor_variant'), 'splice site');
  assert.equal(consequenceOf('SO:0001589|frameshift_variant,SO:0001587|nonsense'), 'nonsense');
  assert.equal(consequenceOf('SO:0001822|inframe_deletion'), 'in-frame indel');
  assert.equal(consequenceOf('SO:0001582|initiator_codon_variant'), 'start lost');
  assert.equal(consequenceOf('SO:0001623|5_prime_UTR_variant,SO:0001627|intron_variant'), 'UTR');
  assert.equal(consequenceOf('SO:0001619|non-coding_transcript_variant'), 'non-coding RNA');
  assert.equal(consequenceOf('SO:0001986|genic_upstream_transcript_variant'), 'near gene');
  assert.equal(consequenceOf('SO:0002073|no_sequence_alteration'), 'none');
  assert.equal(consequenceOf(undefined), 'none');
});

test('starsOf turns review status into ClinVar stars', () => {
  assert.equal(starsOf('practice_guideline'), 4);
  assert.equal(starsOf('reviewed by expert panel'), 3);
  assert.equal(starsOf('criteria_provided,_multiple_submitters,_no_conflicts'), 2);
  assert.equal(starsOf('criteria provided, single submitter'), 1);
  assert.equal(starsOf('criteria_provided,_conflicting_classifications'), 1);
  assert.equal(starsOf('no_assertion_criteria_provided'), 0);
  assert.equal(starsOf('no classification provided'), 0);
  assert.equal(starsOf('-'), null);
  assert.equal(starsOf(undefined), null);
});

test('CHANNELS lists the 96 pyrimidine-centred trinucleotide substitutions and channelOf folds purines', () => {
  assert.equal(CHANNELS.length, 96);
  assert.equal(CHANNELS[0], 'A[C>A]A');
  assert.equal(CHANNELS[34], 'A[C>T]G');
  assert.equal(CHANNELS[95], 'T[T>G]T');
  assert.equal(channelOf('A', 'C', 'G', 'T'), 34);
  assert.equal(channelOf('C', 'G', 'T', 'A'), 34);
  assert.equal(channelOf('A', 'T', 'A', 'C'), 64 + 0 + 0);
  assert.equal(channelOf('N', 'C', 'G', 'T'), -1);
  assert.equal(channelOf('A', 'C', 'G', 'C'), -1);
});

const vcf = (chrom, pos, ref, alt, info) => [chrom, pos, '1', ref, alt, '.', '.', info].join('\t');
const GENOME = { 1: 'AACGTTGCA' };
const base = (chrom, pos) => GENOME[chrom]?.[pos - 1] ?? null;

test('clinvarVcf counts classes, consequences, stars, types, origins, somatic calls and the SNV spectrum', () => {
  const acc = clinvarVcf({ base });
  [
    '##fileformat=VCFv4.1',
    '##fileDate=2026-11-01',
    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    vcf('1', 3, 'C', 'T', 'ALLELEID=1;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;MC=SO:0001587|nonsense;ORIGIN=1'),
    vcf('1', 4, 'G', 'A', 'ALLELEID=2;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;MC=SO:0001819|synonymous_variant;ORIGIN=1'),
    vcf('1', 5, 'T', 'G', 'ALLELEID=3;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;MC=SO:0001583|missense_variant,SO:0001627|intron_variant;ORIGIN=0'),
    vcf('1', 6, 'TG', 'T', 'ALLELEID=4;CLNREVSTAT=reviewed_by_expert_panel;CLNSIG=Pathogenic/Likely_pathogenic;CLNVC=Deletion;MC=SO:0001589|frameshift_variant;ORIGIN=33'),
    vcf('1', 7, 'A', 'T', 'ALLELEID=5;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Pathogenic;CLNVC=single_nucleotide_variant;MC=SO:0001583|missense_variant;ORIGIN=3'),
    vcf('1', 8, 'C', 'T', 'ALLELEID=6;CLNVC=single_nucleotide_variant;ONC=Oncogenic;ONCREVSTAT=criteria_provided,_single_submitter;SCI=Tier_I_-_Strong;ORIGIN=2'),
    vcf('1', 2, 'A', 'AT', 'ALLELEID=7;CLNREVSTAT=criteria_provided,_conflicting_classifications;CLNSIG=Conflicting_classifications_of_pathogenicity|risk_factor;CLNVC=Duplication;ORIGIN=1'),
    vcf('1', 1, 'A', 'G', 'ALLELEID=8;CLNREVSTAT=no_assertion_criteria_provided;CLNSIG=Benign;CLNVC=single_nucleotide_variant;MC=SO:0001623|5_prime_UTR_variant;ORIGIN=1'),
  ].forEach((l) => acc.add(l));
  const r = acc.result();
  assert.equal(r.release, '2026-11-01');
  assert.equal(r.records, 8);
  assert.deepEqual(r.classes, {
    Pathogenic: 2,
    'Pathogenic/Likely pathogenic': 1,
    'Likely pathogenic': 0,
    'Uncertain significance': 1,
    'Likely benign': 1,
    'Benign/Likely benign': 0,
    Benign: 1,
    Conflicting: 1,
    Other: 0,
    none: 1,
  });
  assert.equal(r.consequences.nonsense.Pathogenic, 1);
  assert.equal(r.consequences.missense.Pathogenic, 1);
  assert.equal(r.consequences.missense['Uncertain significance'], 1);
  assert.equal(r.consequences.frameshift['Pathogenic/Likely pathogenic'], 1);
  assert.equal(r.consequences.synonymous['Likely benign'], 1);
  assert.equal(r.consequences.UTR.Benign, 1);
  assert.equal(r.consequences.none.Conflicting, 1);
  assert.equal(Object.keys(r.consequences).length, 13);
  assert.equal(Object.values(r.consequences).reduce((n, row) => n + Object.values(row).reduce((a, b) => a + b, 0), 0), 7);
  assert.deepEqual(r.stars.Pathogenic, [0, 1, 1, 0, 0]);
  assert.deepEqual(r.stars['Pathogenic/Likely pathogenic'], [0, 0, 0, 1, 0]);
  assert.deepEqual(r.stars.Benign, [1, 0, 0, 0, 0]);
  assert.deepEqual(r.stars.Conflicting, [0, 1, 0, 0, 0]);
  assert.deepEqual(r.types, { 'single nucleotide variant': 6, Deletion: 1, Duplication: 1 });
  assert.deepEqual(r.origin, {
    germline: 5, somatic: 1, both: 1, unknown: 1,
  });
  assert.deepEqual(r.oncogenicity, { Oncogenic: 1 });
  assert.deepEqual(r.clinicalImpact, { 'Tier I - Strong': 1 });
  assert.equal(r.spectrum.plp.length, 96);
  assert.equal(r.spectrum.plp[34], 1);
  assert.equal(r.spectrum.blb[34], 1);
  assert.equal(r.spectrum.plp.reduce((a, b) => a + b, 0), 1);
  assert.equal(r.spectrum.blb.reduce((a, b) => a + b, 0), 1);
  assert.equal(r.spectrum.mismatches, 1);
  assert.equal(r.spectrum.skipped, 1);
});

test('vcfClassCounts reads an archived VCF date and its class counts', () => {
  const acc = vcfClassCounts();
  [
    '##fileDate=2018-01-28',
    vcf('1', 1, 'A', 'G', 'ALLELEID=1;CLNSIG=Conflicting_interpretations_of_pathogenicity;CLNVC=single_nucleotide_variant'),
    vcf('1', 2, 'A', 'G', 'CLNSIG=Pathogenic;ALLELEID=2'),
    vcf('1', 3, 'A', 'G', 'ALLELEID=3;CLNSIG=Likely_benign'),
    vcf('1', 4, 'A', 'G', 'ALLELEID=4;CLNSIG=Uncertain_significance'),
    vcf('1', 5, 'A', 'G', 'ALLELEID=5;ONC=Oncogenic'),
  ].forEach((l) => acc.add(l));
  const r = acc.result();
  assert.equal(r.release, '2018-01-28');
  assert.equal(r.records, 5);
  assert.equal(r.classes.Conflicting, 1);
  assert.equal(r.classes.Pathogenic, 1);
  assert.equal(r.classes['Likely benign'], 1);
  assert.equal(r.classes['Uncertain significance'], 1);
  assert.equal(r.classes.none, 1);
});

test('earliestFile picks the first dated release in a folder listing', () => {
  const listing = ['clinvar_20180225.vcf.gz', 'clinvar_20180128.vcf.gz.md5', 'clinvar_20180128_papu.vcf.gz', 'clinvar_20180128.vcf.gz', 'clinvar_20180107_papu.vcf.gz'].map((f) => `<a href="${f}">${f}</a>`).join('\n');
  assert.deepEqual(earliestFile(listing, /clinvar_(\d{8})\.vcf\.gz/), { file: 'clinvar_20180128.vcf.gz', version: '20180128' });
  assert.equal(earliestFile('<a href="x">x</a>', /clinvar_(\d{8})\.vcf\.gz/), null);
});

const GSS_HEAD = '#Symbol\tGeneID\tTotal_submissions\tTotal_alleles\tSubmissions_reporting_this_gene\tAlleles_reported_Pathogenic_Likely_pathogenic\tGene_MIM_number\tNumber_uncertain\tNumber_with_conflicts';

test('geneSummary ranks protein-coding genes by P/LP alleles and traces how concentrated they are', () => {
  const text = [
    '#Overview of data in ClinVar by gene, dated October 31, 2026', GSS_HEAD,
    'TTN\t7273\t97618\t40949\t96593\t6183\t188840\t12270\t4459',
    'TTN-AS1\t100506866\t56352\t23410\t130\t4506\t-\t7289\t6',
    'BRCA2\t675\t77762\t22068\t77593\t6017\t600185\t4356\t5744',
    'A1BG\t1\t107\t106\t87\t-\t138670\t85\t0',
  ].join('\n');
  const g = geneSummary(text, new Set(['TTN', 'BRCA2', 'A1BG']));
  assert.equal(g.date, 'October 31, 2026');
  assert.equal(g.genes, 3);
  assert.equal(g.withPlp, 2);
  assert.equal(g.plp, 12200);
  assert.equal(g.vus, 16711);
  assert.deepEqual(g.top[0], {
    gene: 'TTN', plp: 6183, vus: 12270, conflicts: 4459, alleles: 40949,
  });
  assert.deepEqual(g.top.map((t) => t.gene), ['TTN', 'BRCA2']);
  assert.equal(g.lorenz.length, 101);
  assert.equal(g.lorenz[0], 0);
  assert.equal(g.lorenz[100], 1);
  assert.equal(g.lorenz[50], 0.5068);
});

const ORG_HEAD = ['#organization', 'organization ID', 'institution type', 'street address', 'city', 'country', 'number of ClinVar submissions', 'date last submitted', 'maximum review status', 'collection methods', 'novel and updates', 'clinical significance categories submitted', 'number of submissions from clinical testing', 'number of submissions from research', 'number of submissions from literature only', 'number of submissions from curation', 'number of submissions from phenotyping', 'somatic clinical impact values submitted', 'somatic oncogenicity values submitted'].join('\t');
const org = (name, type, country, n, methods) => [name, '1', type, 'street', 'city', country, n, 'Oct 1, 2026', 'criteria provided, single submitter', 'clinical testing', 'novel and updates', 'Pathogenic', ...methods, '', ''].join('\t');

test('organizations totals ClinVar submitters by type, country and collection method', () => {
  const text = [
    ORG_HEAD,
    org('Labcorp Genetics (formerly Invitae); LabCorp', 'lab', 'United States', 2034787, [2034000, 787, 0, 0, 0]),
    org('OMIM; Johns Hopkins University', 'resource', 'United States', 37918, [0, 0, 37915, 0, 0]),
    org('Some Clinic', 'clinic', 'China', 12, [12, 0, 0, 0, 0]),
    org('Blank Country Lab', 'lab', '', 5, [5, 0, 0, 0, 0]),
  ].join('\n');
  const o = organizations(text);
  assert.equal(o.organizations, 4);
  assert.equal(o.submissions, 2072722);
  assert.deepEqual(o.types, [{ type: 'lab', orgs: 2, submissions: 2034792 }, { type: 'resource', orgs: 1, submissions: 37918 }, { type: 'clinic', orgs: 1, submissions: 12 }]);
  assert.deepEqual(o.countries, [{ country: 'United States', orgs: 2, submissions: 2072705 }, { country: 'China', orgs: 1, submissions: 12 }, { country: 'not stated', orgs: 1, submissions: 5 }]);
  assert.deepEqual(o.top[0], {
    name: 'Labcorp Genetics (formerly Invitae); LabCorp', type: 'lab', country: 'United States', submissions: 2034787,
  });
  assert.deepEqual(o.methods, {
    'clinical testing': 2034017, research: 787, 'literature only': 37915, curation: 0, phenotyping: 0,
  });
});

test('proteinChange reads the residue and kind from a ClinVar preferred name', () => {
  const tp53 = (c) => `NM_000546.6(TP53):${c}`;
  assert.deepEqual(proteinChange(tp53('c.818G>A (p.Arg273His)')), {
    transcript: 'NM_000546.6', position: 273, kind: 'missense', aa: 'Arg',
  });
  assert.equal(proteinChange(tp53('c.637C>T (p.Arg213Ter)')).kind, 'truncating');
  assert.equal(proteinChange(tp53('c.637C>T (p.Arg213*)')).kind, 'truncating');
  assert.deepEqual(proteinChange(tp53('c.215del (p.Pro72fs)')), {
    transcript: 'NM_000546.6', position: 72, kind: 'truncating', aa: 'Pro',
  });
  assert.equal(proteinChange(tp53('c.215del (p.Pro72LeufsTer3)')).kind, 'truncating');
  assert.equal(proteinChange(tp53('c.639A>G (p.Arg213=)')).kind, 'synonymous');
  assert.deepEqual(proteinChange('NM_000492.4(CFTR):c.1521_1523del (p.Phe508del)'), {
    transcript: 'NM_000492.4', position: 508, kind: 'inframe', aa: 'Phe',
  });
  assert.equal(proteinChange('NM_000492.4(CFTR):c.1519_1521del (p.Ile507_Phe508del)').position, 507);
  assert.equal(proteinChange('NM_014855.3(AP5Z1):c.80_83delinsTGCTGTAAACTGTAACTGTAAA (p.Arg27_Ile28delinsLeuLeuTer)').kind, 'truncating');
  assert.equal(proteinChange(tp53('c.1A>G (p.Met1?)')).kind, 'other');
  assert.equal(proteinChange(tp53('c.1181A>C (p.Ter394Ser)')).kind, 'other');
  assert.deepEqual(proteinChange(tp53('c.672+1G>A')), {
    transcript: 'NM_000546.6', position: 224, kind: 'splice', aa: '',
  });
  assert.equal(proteinChange(tp53('c.673-2A>G')).position, 225);
  for (const none of [tp53('c.672+31A>G'), tp53('c.-29G>A'), tp53('c.818G>A (p.?)'), 'GRCh38/hg38 17p13.1(chr17:7668402-7687550)x1']) assert.equal(proteinChange(none), null, none);
});

const TABLE_COLS = ['AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinicalSignificance', 'ClinSigSimple', 'LastEvaluated', 'RS# (dbSNP)', 'nsv/esv (dbVar)', 'RCVaccession', 'PhenotypeIDS', 'PhenotypeList', 'Origin', 'OriginSimple', 'Assembly', 'ChromosomeAccession', 'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele', 'Cytogenetic', 'ReviewStatus', 'NumberSubmitters'];
const tableRow = (fields) => TABLE_COLS.map((c) => fields[c] ?? '-').join('\t');
const tp53Row = (name, sig, extra = {}) => tableRow({
  Type: 'single nucleotide variant', Name: `NM_000546.6(TP53):${name}`, GeneSymbol: 'TP53', ClinicalSignificance: sig, Assembly: 'GRCh38', NumberSubmitters: '1', ...extra,
});

test('variantTable counts GRCh38 rows by type, submitters and evaluation year and builds mirrored lollipops', () => {
  const acc = variantTable({ genes: ['TP53', 'BRCA1', 'CFTR', 'MYH7'] });
  [
    `#${TABLE_COLS.join('\t')}`,
    tp53Row('c.818G>A (p.Arg273His)', 'Pathogenic', { LastEvaluated: 'Dec 17, 2024', NumberSubmitters: '12' }),
    tp53Row('c.818G>A (p.Arg273His)', 'Pathogenic', { Assembly: 'GRCh37' }),
    tp53Row('c.817C>T (p.Arg273Cys)', 'Likely pathogenic', { LastEvaluated: 'Jan 3, 2025', NumberSubmitters: '3' }),
    tp53Row('c.637C>T (p.Arg213Ter)', 'Pathogenic/Likely pathogenic'),
    tp53Row('c.215del (p.Pro72fs)', 'Pathogenic', { Type: 'Deletion', LastEvaluated: 'Feb 1, 2023' }),
    tp53Row('c.672+1G>A', 'Pathogenic'),
    tp53Row('c.639A>G (p.Arg213=)', 'Benign', { NumberSubmitters: '25' }),
    tp53Row('c.672+31A>G', 'Likely benign'),
    tp53Row('c.100A>G (p.Thr34Ala)', 'Uncertain significance'),
    tableRow({
      Type: 'single nucleotide variant', Name: 'NM_001126112.3(TP53):c.818G>A (p.Arg273His)', GeneSymbol: 'TP53', ClinicalSignificance: 'Pathogenic', Assembly: 'GRCh38', NumberSubmitters: '1',
    }),
    tableRow({
      Type: 'Deletion', Name: 'NM_000492.4(CFTR):c.1521_1523del (p.Phe508del)', GeneSymbol: 'CFTR', ClinicalSignificance: 'Pathogenic', Assembly: 'GRCh38', NumberSubmitters: '40',
    }),
    tableRow({
      Type: 'copy number loss', Name: 'GRCh38/hg38 7q31.2(chr7:117480025-117668665)x1', GeneSymbol: 'CFTR', ClinicalSignificance: 'Pathogenic', Assembly: 'GRCh38', NumberSubmitters: '1',
    }),
    tableRow({
      Type: 'single nucleotide variant', Name: 'NM_000059.4(BRCA2):c.100G>A (p.Ala34Thr)', GeneSymbol: 'BRCA2', ClinicalSignificance: 'Conflicting classifications of pathogenicity', Assembly: 'GRCh38', NumberSubmitters: '4',
    }),
  ].forEach((l) => acc.add(l));
  const r = acc.result();
  assert.equal(r.records, 12);
  assert.deepEqual(r.types, { 'single nucleotide variant': 9, Deletion: 2, 'copy number loss': 1 });
  assert.equal(r.submitters.bins.length, 11);
  assert.deepEqual(r.submitters['P/LP'], [5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2]);
  assert.deepEqual(r.submitters['B/LB'], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);
  assert.deepEqual(r.submitters.VUS, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
  assert.deepEqual(r.submitters.Conflicting, [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]);
  assert.deepEqual(r.evaluated['P/LP'], { 2023: 1, 2024: 1, 2025: 1, none: 5 });
  assert.deepEqual(r.lollipops.TP53, {
    transcript: 'NM_000546.6',
    up: [[72, 'truncating', 1, 'Pro'], [213, 'truncating', 1, 'Arg'], [224, 'splice', 1, ''], [273, 'missense', 2, 'Arg']],
    down: [[213, 'synonymous', 1, 'Arg']],
    unplaced: { up: 1, down: 1 },
  });
  assert.deepEqual(r.lollipops.CFTR.up, [[508, 'inframe', 1, 'Phe']]);
  assert.deepEqual(r.lollipops.CFTR.unplaced, { up: 1, down: 0 });
  assert.equal(r.lollipops.BRCA1, undefined);
});

test('pfamDomains lists each Pfam match on a protein in order', () => {
  const json = {
    results: [
      { metadata: { accession: 'PF00870', name: 'P53 DNA-binding domain' }, proteins: [{ entry_protein_locations: [{ fragments: [{ start: 100, end: 288 }] }] }] },
      { metadata: { accession: 'PF07710', name: 'P53 tetramerisation motif' }, proteins: [{ entry_protein_locations: [{ fragments: [{ start: 319, end: 357 }] }] }] },
      { metadata: { accession: 'PF08563', name: 'P53 transactivation motif' }, proteins: [{ entry_protein_locations: [{ fragments: [{ start: 6, end: 30 }] }, { fragments: [{ start: 40, end: 50 }] }] }] },
    ],
  };
  assert.deepEqual(pfamDomains(json), [
    {
      id: 'PF08563', name: 'P53 transactivation motif', start: 6, end: 30,
    },
    {
      id: 'PF08563', name: 'P53 transactivation motif', start: 40, end: 50,
    },
    {
      id: 'PF00870', name: 'P53 DNA-binding domain', start: 100, end: 288,
    },
    {
      id: 'PF07710', name: 'P53 tetramerisation motif', start: 319, end: 357,
    },
  ]);
});

const quoted = (cells) => cells.map((c) => `"${c}"`).join(',');
const RULE = quoted(Array(10).fill('+++'));
const curation = (gene, moi, cls, date, gcep, mondo) => quoted([gene, 'HGNC:1', 'a disease', mondo, moi, 'SOP10', cls, 'https://search.clinicalgenome.org/kb/gene-validity/x', `${date}T16:00:00.000Z`, gcep]);

test('validitySummary counts ClinGen curations by class, inheritance, year and panel', () => {
  const text = [
    quoted(['CLINGEN GENE DISEASE VALIDITY CURATIONS', '', '', '', '', '', '', '', '', '']),
    quoted(['FILE CREATED: 2026-11-04', '', '', '', '', '', '', '', '', '']),
    quoted(['WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity', '', '', '', '', '', '', '', '', '']),
    RULE,
    quoted(['GENE SYMBOL', 'GENE ID (HGNC)', 'DISEASE LABEL', 'DISEASE ID (MONDO)', 'MOI', 'SOP', 'CLASSIFICATION', 'ONLINE REPORT', 'CLASSIFICATION DATE', 'GCEP']),
    RULE,
    curation('AARS1', 'AD', 'Definitive', '2024-03-14', 'Charcot-Marie-Tooth Disease Gene Curation Expert Panel', 'MONDO:1'),
    curation('AARS2', 'AR', 'Definitive', '2022-04-18', 'Mitochondrial Diseases Gene Curation Expert Panel', 'MONDO:2'),
    curation('AARS2', 'AD', 'Disputed', '2023-01-02', 'Mitochondrial Diseases Gene Curation Expert Panel', 'MONDO:3'),
    curation('SCN2A', 'AD', 'Definitive', '2019-05-01', 'Epilepsy Gene Curation Expert Panel', 'MONDO:4'),
    curation('SCN2A', 'AR', 'Moderate', '2023-05-01', 'Epilepsy Gene Curation Expert Panel', 'MONDO:4'),
    curation('MT-X', 'MT', 'Strong', '2023-06-01', 'Mitochondrial Diseases Gene Curation Expert Panel', 'MONDO:5'),
    curation('GENE5', 'UD', 'Limited', '2018-01-01', 'Hearing Loss Gene Curation Expert Panel', 'MONDO:6'),
  ].join('\n');
  const v = validitySummary(text);
  assert.equal(v.date, '2026-11-04');
  assert.equal(v.curations, 7);
  assert.equal(v.genes, 5);
  assert.equal(v.diseases, 6);
  assert.deepEqual(v.classes, {
    Definitive: 3, Strong: 1, Moderate: 1, Limited: 1, Disputed: 1, Refuted: 0, 'No Known Disease Relationship': 0,
  });
  assert.deepEqual(v.byMoi, {
    AD: { Definitive: 2, Disputed: 1 }, AR: { Definitive: 1, Moderate: 1 }, MT: { Strong: 1 }, UD: { Limited: 1 },
  });
  assert.deepEqual(v.byYear, {
    2018: { Limited: 1 }, 2019: { Definitive: 1 }, 2022: { Definitive: 1 }, 2023: { Disputed: 1, Moderate: 1, Strong: 1 }, 2024: { Definitive: 1 },
  });
  assert.deepEqual(v.gceps[0], { name: 'Mitochondrial Diseases Gene Curation Expert Panel', n: 3 });
  assert.equal(v.gceps.length, 4);
  assert.deepEqual([...v.moiGenes], [['AARS1', 'AD'], ['AARS2', 'AR'], ['SCN2A', 'both'], ['MT-X', 'MT']]);
});

test('supportedGenes lists the HGNC IDs of genes with a definitive, strong or moderate ClinGen curation', () => {
  const row = (gene, hgnc, cls) => quoted([gene, hgnc, 'a disease', 'MONDO:1', 'AD', 'SOP10', cls, 'https://search.clinicalgenome.org/kb/gene-validity/x', '2024-01-01T16:00:00.000Z', 'Panel']);
  const text = [
    quoted(['FILE CREATED: 2026-11-04', '', '', '', '', '', '', '', '', '']),
    RULE,
    quoted(['GENE SYMBOL', 'GENE ID (HGNC)', 'DISEASE LABEL', 'DISEASE ID (MONDO)', 'MOI', 'SOP', 'CLASSIFICATION', 'ONLINE REPORT', 'CLASSIFICATION DATE', 'GCEP']),
    RULE,
    row('AARS1', 'HGNC:20', 'Definitive'),
    row('AARS2', 'HGNC:21022', 'Limited'),
    row('AARS2', 'HGNC:21022', 'Moderate'),
    row('GENE5', 'HGNC:5', 'Disputed'),
    row('MT-X', 'HGNC:7', 'Strong'),
  ].join('\n');
  assert.deepEqual([...supportedGenes(text)], ['HGNC:20', 'HGNC:21022', 'HGNC:7']);
});

const REGION_HEAD = ['#ISCA ID', 'ISCA Region Name', 'cytoBand', 'Genomic Location', 'Haploinsufficiency Score', 'Haploinsufficiency Description', 'Haploinsufficiency PMID1', 'Triplosensitivity Score', 'Triplosensitivity Description', 'Date Last Evaluated'].join('\t');
const region = (id, name, band, loc, hi, ts) => [id, name, band, loc, hi, 'x', '1', ts, 'y', '2025-10-27'].join('\t');

test('regionSummary scores ClinGen dosage regions and lists the sensitive ones with their spans', () => {
  const text = [
    '#ClinGen Region Curation Results', '#28 Apr,2026', '#Genomic Locations are reported on GRCh38 (hg38): GCF_000001405.36', REGION_HEAD,
    region('ISCA-37446', '22q11.2 recurrent (DGS/VCFS) region (proximal, A-D) (includes TBX1)', '22q11.21', 'chr22:18924718-21111383', '3', '3'),
    region('ISCA-1', 'Xq25 region (includes STAG2)', 'Xq25', 'chrX:123900469-124102669', '3', '1'),
    region('ISCA-2', '16p12.3 recurrent region (BP3-BP4)', '16p12.3', 'chr16:16762939-18073544', '0', '0'),
    region('ISCA-3', '21q21.3 region (includes APP)', '21q21.3', 'chr21:25734570-26573286', '0', '3'),
    region('ISCA-4', 'A region unlikely to matter', '1p36.33', 'chr1:100-200', '40', '40'),
    region('ISCA-5', 'Not yet placed', '1p36', 'tbd', '', ''),
  ].join('\n');
  const r = regionSummary(text);
  assert.equal(r.date, '28 Apr,2026');
  assert.equal(r.regions, 6);
  assert.deepEqual(r.hi, { 3: 2, 0: 2, 40: 1 });
  assert.deepEqual(r.ts, { 3: 2, 1: 1, 0: 1, 40: 1 });
  assert.deepEqual(r.sensitive.map((s) => s.id), ['ISCA-3', 'ISCA-37446', 'ISCA-1']);
  assert.deepEqual(r.sensitive[1], {
    id: 'ISCA-37446', name: '22q11.2 recurrent (DGS/VCFS) region (proximal, A-D) (includes TBX1)', band: '22q11.21', chrom: '22', start: 18924718, end: 21111383, hi: 3, ts: 3, recurrent: true,
  });
  assert.equal(r.sensitive[2].recurrent, false);
});
