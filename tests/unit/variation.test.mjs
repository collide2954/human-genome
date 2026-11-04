import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  dbsnpNotes, svSummary, constraintSummary, clingenDosage, impcViability, hgncMouse, groupedValues,
} from '../../scripts/variation.mjs';

const NOTES = [
  'dbSNP build 157 release notes', '', 'Organism name: Homo sapiens', '', '1. RefSnp (RS)', '', 'Total RS: 1206053617', '',
  '1.1 RS counts by location:', 'chr1\t93459666', 'chrX\t45803018', 'chrY\t3048013', 'chrM\t9229', 'Alt Only\t35320', 'PAR\t1427972', '',
  '1.2 RS counts by type:', 'Live RS\t1172689405', 'Unsupported RS\t4922330', 'Withdrawn RS\t6854924', 'Locationless RS\t1286', 'Merged RS\t21585672', '',
  '2. SubSnp (SS)', '', 'Total SS: 4849775973', 'Unmapped SS: 153387',
].join('\n');

test('dbsnpNotes reads the build, totals and counts by location and status', () => {
  assert.deepEqual(dbsnpNotes(NOTES), {
    build: 157,
    total: 1206053617,
    byChrom: {
      1: 93459666, X: 45803018, Y: 3048013, MT: 9229,
    },
    other: { 'Alt Only': 35320, PAR: 1427972 },
    status: {
      'Live RS': 1172689405, 'Unsupported RS': 4922330, 'Withdrawn RS': 6854924, 'Locationless RS': 1286, 'Merged RS': 21585672,
    },
    ss: { total: 4849775973, unmapped: 153387 },
  });
});

const sv = (chr, start, stop, type, study, clin = '') => [chr, start, stop, 1, type, 'Merging', 'Merging', 'x', study, 'nssv1', clin, '', 'medium'].join('\t');

test('svSummary counts non-redundant SVs by type, size, chromosome, study and clinical assertion', () => {
  const acc = svSummary();
  [
    '#NR_SVs GRCh38',
    '#chr\touterm_start\touterm_stop\tvariant_count\tvariant_type\tmethod\tanalysis\tplatform\tstudy\tvariant\tclinical_assertion\tclinvar_accession\tbin_size',
    sv('1', 1001, 1100, 'deletion', 'Ebert2021'),
    sv('1', 5001, 15000, 'deletion', 'gnomAD', 'Pathogenic'),
    sv('X', 1, 1000000, 'duplication', 'Ebert2021', 'Benign;Likely benign'),
    sv('2', 500, 501, 'insertion', 'Ebert2021', 'Pathogenic;Benign'),
    sv('Un_KI270742v1', 1, 100, 'deletion', 'x'),
  ].forEach((l) => acc.add(l));
  const s = acc.result();
  assert.deepEqual(s.classes, { deletion: 2, duplication: 1, insertion: 1 });
  assert.deepEqual(s.byChrom, { 1: { deletion: 2 }, X: { duplication: 1 }, 2: { insertion: 1 } });
  assert.equal(s.sizes.deletion.counts[4], 1);
  assert.equal(s.sizes.deletion.counts[8], 1);
  assert.equal(s.sizes.duplication.counts[12], 1);
  assert.equal(s.sizes.insertion, undefined);
  assert.deepEqual(s.studies[0], { study: 'Ebert2021', n: 3 });
  assert.deepEqual(s.clinical, {
    pathogenic: 1, benign: 1, conflicting: 1, uncertain: 0, other: 0, none: 1,
  });
});

const CONSTRAINT_HEAD = ['gene', 'gene_id', 'transcript', 'canonical', 'mane_select', 'lof.exp', 'lof.pLI', 'lof.oe_ci.upper', 'mis.z_score', 'constraint_flags'].join('\t');
const row = (gene, mane, exp, pli, loeuf, misZ, flags = '[]') => [gene, 'ENSG1', 'ENST1', 'true', mane, exp, pli, loeuf, misZ, flags].join('\t');

test('svSummary folds dbVar compound types into deletion, duplication, CNV and insertion classes and counts mobile elements', () => {
  const acc = svSummary();
  [
    sv('1', 1, 300, 'copy_number_loss;deletion', 's'),
    sv('1', 1, 300, 'alu_deletion', 's'),
    sv('1', 1, 3000, 'copy_number_gain;copy_number_variation', 's'),
    sv('1', 1, 3000, 'copy_number_variation', 's'),
    sv('1', 1, 2, 'alu_insertion;insertion', 's'),
    sv('1', 1, 2, 'line1_insertion', 's'),
    sv('1', 1, 2, 'sva_insertion;mobile_element_insertion', 's'),
  ].forEach((l) => acc.add(l));
  const s = acc.result();
  assert.deepEqual(s.classes, {
    deletion: 2, duplication: 1, cnv: 1, insertion: 3,
  });
  assert.deepEqual(s.mobile, {
    Alu: { insertion: 1, deletion: 1 }, L1: { insertion: 1, deletion: 0 }, SVA: { insertion: 1, deletion: 0 }, HERV: { insertion: 0, deletion: 0 },
  });
  assert.deepEqual(Object.keys(s.sizes).sort(), ['cnv', 'deletion', 'duplication']);
  assert.equal(s.sizes.deletion.counts.reduce((a, b) => a + b, 0), 2);
  assert.deepEqual(s.byChrom, {
    1: {
      deletion: 2, duplication: 1, cnv: 1, insertion: 3,
    },
  });
});

test('constraintSummary counts a MANE transcript listed under Ensembl and RefSeq IDs once', () => {
  const acc = constraintSummary();
  [
    CONSTRAINT_HEAD,
    row('XPO6', 'true', '119', '1.00', '0.04', '5'),
    row('XPO6', 'true', '119', '1.00', '0.04', '5'),
    row('TTN', 'true', '400', '0.00', '0.85', '5.1'),
  ].forEach((l) => acc.add(l));
  const c = acc.result();
  assert.equal(c.genes, 2);
  assert.deepEqual(c.top.map((t) => t.gene), ['XPO6', 'TTN']);
});

test('constraintSummary reads one MANE Select transcript per gene', () => {
  const acc = constraintSummary();
  [
    CONSTRAINT_HEAD,
    row('TTN', 'true', '400', '0.00', '0.85', '5.1'),
    row('SCN2A', 'true', '80', '1.00', '0.12', '6.2'),
    row('SCN2A', 'false', '70', '1.00', '0.20', '6.0'),
    row('OR4F5', 'true', '2.1', '0.30', '1.95', '-1.2'),
    row('NEW1', 'true', 'NA', 'NA', 'NA', 'NA'),
  ].forEach((l) => acc.add(l));
  const c = acc.result();
  assert.equal(c.genes, 3);
  assert.equal(c.loeuf.summary.median, 0.85);
  assert.equal(c.below06, 1);
  assert.equal(c.pli09, 1);
  assert.deepEqual(c.top[0], {
    gene: 'SCN2A', loeuf: 0.12, pli: 1, expected: 80,
  });
  assert.equal(c.misZ.counts.reduce((a, b) => a + b, 0), 3);
  assert.equal(acc.values().get('SCN2A'), 0.12);
});

test('clingenDosage counts haploinsufficiency and triplosensitivity scores', () => {
  const text = [
    '#ClinGen Gene Curation Results', '#04 Nov,2026',
    '#Gene Symbol\tGene ID\tcytoBand\tGenomic Location\tHaploinsufficiency Score\tHaploinsufficiency Description\tTriplosensitivity Score',
    'AAGAB\t79719\t15q23\tchr15:1-2\t3\tSufficient evidence\t0',
    'A4GALT\t53947\t22q13.2\tchr22:1-2\t30\tAutosomal recessive\t0',
    'MECP2\t4204\tXq28\tchrX:1-2\t2\tSome evidence\t3',
    'X1\t1\t1p\tchr1:1-2\tNot yet evaluated\t\tNot yet evaluated',
  ].join('\n');
  const d = clingenDosage(text);
  assert.equal(d.date, '04 Nov,2026');
  assert.deepEqual(d.hi, { 3: 1, 30: 1, 2: 1 });
  assert.deepEqual(d.ts, { 0: 2, 3: 1 });
  assert.deepEqual([...d.hiGenes], ['AAGAB']);
});

test('impcViability gives each mouse gene its most severe call', () => {
  const head = '"Gene Symbol","Gene Accession Id","Zygosity","Viability Phenotype HOMs/HEMIs","Comment"';
  const v = impcViability([head, '"Abc","MGI:1","homozygote","viable",""', '"Abc","MGI:1","homozygote","lethal","Conflicting data"', '"Def","MGI:2","homozygote","subviable",""', '"Ghi","MGI:3","hemizygote","viable",""'].join('\n'));
  assert.deepEqual(v.calls, { lethal: 1, subviable: 1, viable: 1 });
  assert.equal(v.genes, 3);
  assert.equal(v.byMgi.get('MGI:1'), 'lethal');
});

test('hgncMouse maps MGI IDs to human symbols and groupedValues splits a measure by group', () => {
  const text = ['hgnc_id\tsymbol\tmgd_id', 'HGNC:1\tSCN2A\tMGI:1', 'HGNC:2\tTTN\tMGI:2|MGI:9', 'HGNC:3\tOR4F5\t'].join('\n');
  const map = hgncMouse(text);
  assert.equal(map.get('MGI:1'), 'SCN2A');
  assert.equal(map.get('MGI:9'), 'TTN');
  const values = new Map([['SCN2A', 0.12], ['TTN', 0.85], ['OR4F5', 1.95]]);
  const groups = new Map([['SCN2A', 'lethal'], ['TTN', 'viable']]);
  const g = groupedValues(values, groups, [0, 0.5, 1, 1.5, 2]);
  assert.deepEqual(g.lethal, { n: 1, median: 0.12, counts: [1, 0, 0, 0] });
  assert.deepEqual(g.viable, { n: 1, median: 0.85, counts: [0, 1, 0, 0] });
  assert.deepEqual(g.other, { n: 1, median: 1.95, counts: [0, 0, 0, 1] });
});
