import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  stringNames, stringSummary, physicalSummary, reactomeSummary, oboTerms, goSummary,
} from '../../scripts/networks.mjs';

const INFO = [
  '#string_protein_id\tpreferred_name\tprotein_size\tannotation',
  '9606.A\tALPHA\t100\tFirst protein',
  '9606.B\tBETA\t200\tSecond protein',
].join('\n');

const DETAILED = [
  'protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score',
  '9606.A 9606.B 0 0 0 0 800 900 500 990',
  '9606.B 9606.A 0 0 0 0 800 900 500 990',
  '9606.A 9606.C 0 0 0 0 0 0 750 750',
  '9606.C 9606.A 0 0 0 0 0 0 750 750',
  '9606.A 9606.D 0 0 0 300 0 0 300 450',
  '9606.D 9606.A 0 0 0 300 0 0 300 450',
  '9606.B 9606.C 0 0 0 0 0 0 200 200',
  '9606.C 9606.B 0 0 0 0 0 0 200 200',
];

test('stringNames maps STRING protein IDs to preferred names', () => {
  assert.deepEqual(stringNames(INFO), { '9606.A': 'ALPHA', '9606.B': 'BETA' });
});

test('stringSummary counts each pair once at four confidence cut-offs', () => {
  const acc = stringSummary({ names: stringNames(INFO) });
  DETAILED.forEach((l) => acc.add(l));
  const s = acc.result();
  assert.deepEqual(s.edges, {
    all: 4, medium: 3, high: 2, highest: 1,
  });
  assert.deepEqual(s.proteins, { all: 4, high: 3 });
  assert.deepEqual(s.hubs[0], { id: '9606.A', name: 'ALPHA', degree: 2 });
  assert.equal(s.hubs[1].degree, 1);
  assert.equal(s.degree.counts[0], 2);
  assert.equal(s.degree.counts[1], 1);
  assert.equal(s.scores.counts[16], 1);
  assert.equal(s.scores.counts[12], 1);
  assert.equal(s.scores.counts[6], 1);
  assert.equal(s.scores.counts[1], 1);
});

test('stringSummary counts the evidence channels behind high-confidence links', () => {
  const acc = stringSummary({ names: {} });
  DETAILED.forEach((l) => acc.add(l));
  const { channels, only } = acc.result();
  assert.deepEqual(channels, {
    neighborhood: 0, fusion: 0, cooccurence: 0, coexpression: 0, experimental: 1, database: 1, textmining: 2,
  });
  assert.equal(only.textmining, 1);
  assert.equal(only.experimental, 0);
});

test('physicalSummary counts physical links once per pair', () => {
  const acc = physicalSummary();
  ['protein1 protein2 combined_score', '9606.A 9606.B 950', '9606.B 9606.A 950', '9606.A 9606.C 500', '9606.C 9606.A 500', '9606.B 9606.C 300', '9606.C 9606.B 300']
    .forEach((l) => acc.add(l));
  assert.deepEqual(acc.result(), {
    edges: {
      all: 3, medium: 2, high: 1, highest: 1,
    },
    proteins: { all: 3, high: 2 },
  });
});

const PATHWAYS = [
  'R-HSA-1\tSignal Transduction\tHomo sapiens',
  'R-HSA-2\tSignaling by EGFR\tHomo sapiens',
  'R-HSA-3\tMetabolism\tHomo sapiens',
  'R-HSA-4\tGlycolysis\tHomo sapiens',
  'R-HSA-5\tGlucose metabolism\tHomo sapiens',
  'R-MMU-1\tSignal Transduction\tMus musculus',
].join('\n');

const RELATIONS = ['R-HSA-1\tR-HSA-2', 'R-HSA-3\tR-HSA-5', 'R-HSA-5\tR-HSA-4', 'R-MMU-1\tR-MMU-2'].join('\n');

const map = (acc, id, name, evidence = 'TAS', species = 'Homo sapiens') => [acc, id, `https://reactome.org/PathwayBrowser/#/${id}`, name, evidence, species].join('\t');

const MAPPINGS = [
  map('P1', 'R-HSA-1', 'Signal Transduction'),
  map('P1', 'R-HSA-2', 'Signaling by EGFR'),
  map('P2', 'R-HSA-1', 'Signal Transduction'),
  map('P2', 'R-HSA-3', 'Metabolism', 'IEA'),
  map('P2-2', 'R-HSA-3', 'Metabolism'),
  map('P3', 'R-HSA-3', 'Metabolism'),
  map('P3', 'R-HSA-5', 'Glucose metabolism'),
  map('P3', 'R-HSA-4', 'Glycolysis'),
  map('P4', 'R-HSA-1', 'Signal Transduction'),
  map('Q9', 'R-MMU-1', 'Signal Transduction', 'IEA', 'Mus musculus'),
];

const reactome = () => {
  const acc = reactomeSummary({ pathways: PATHWAYS, relations: RELATIONS });
  MAPPINGS.forEach((l) => acc.add(l));
  return acc.result();
};

test('reactomeSummary sizes the human top-level pathways and their hierarchy', () => {
  const r = reactome();
  assert.equal(r.pathways, 5);
  assert.equal(r.proteins, 4);
  assert.deepEqual(r.topLevel, [
    {
      id: 'R-HSA-1', name: 'Signal Transduction', proteins: 3, pathways: 2,
    },
    {
      id: 'R-HSA-3', name: 'Metabolism', proteins: 2, pathways: 3,
    },
  ]);
  assert.deepEqual(r.depth, [2, 2, 1]);
  assert.deepEqual(r.breadth.slice(0, 3), [0, 3, 1]);
  assert.equal(r.largest[0].name, 'Glucose metabolism');
  assert.deepEqual(r.evidence, { TAS: 8, IEA: 1 });
});

test('reactomeSummary measures how much the top-level pathways share', () => {
  const r = reactome();
  assert.equal(r.overlap.matrix[0][1], 0.25);
  assert.equal(r.overlap.matrix[1][0], 0.25);
  assert.equal(r.overlap.matrix[0][0], 1);
  assert.deepEqual([...r.overlap.order].sort(), [0, 1]);
  assert.equal(r.sizes.counts.reduce((a, b) => a + b, 0), 5);
});

const OBO = [
  'format-version: 1.2',
  'data-version: releases/2026-08-08',
  '',
  '[Term]',
  'id: GO:0005515',
  'name: protein binding',
  'namespace: molecular_function',
  '',
  '[Term]',
  'id: GO:0006915',
  'name: apoptotic process',
  'namespace: biological_process',
  '',
  '[Term]',
  'id: GO:0000001',
  'name: mitochondrion inheritance',
  'namespace: biological_process',
  'is_obsolete: true',
  '',
  '[Typedef]',
  'id: part_of',
  'name: part of',
].join('\n');

test('oboTerms counts live terms by namespace and keeps their names', () => {
  const o = oboTerms(OBO);
  assert.equal(o.version, 'releases/2026-08-08');
  assert.deepEqual(o.counts, { biological_process: 1, molecular_function: 1, cellular_component: 0 });
  assert.equal(o.obsolete, 1);
  assert.equal(o.terms.get('GO:0005515').name, 'protein binding');
  assert.equal(o.terms.has('part_of'), false);
});

const gaf = (acc, symbol, qualifier, term, evidence, aspect, by, db = 'UniProtKB', type = 'protein') => [db, acc, symbol, qualifier, term, 'PMID:1', evidence, '', aspect, 'name', '', type, 'taxon:9606', '20240101', by, '', ''].join('\t');

test('goSummary keeps RNAs and complexes out of the protein counts but in the annotation totals', () => {
  const acc = goSummary(oboTerms(OBO).terms);
  [
    gaf('P1', 'TP53', 'involved_in', 'GO:0006915', 'IDA', 'P', 'UniProt'),
    gaf('URS0001', 'MIR21', 'involved_in', 'GO:0035195', 'IDA', 'P', 'RNAcentral', 'RNAcentral', 'miRNA'),
    gaf('URS0001', 'MIR21', 'involved_in', 'GO:0035195', 'IMP', 'P', 'RNAcentral', 'RNAcentral', 'miRNA'),
    gaf('CPX-1', 'Complex', 'part_of', 'GO:0005634', 'IDA', 'C', 'ComplexPortal', 'ComplexPortal', 'protein_complex'),
  ].forEach((l) => acc.add(l));
  const g = acc.result();
  assert.equal(g.annotations, 4);
  assert.equal(g.proteins, 1);
  assert.deepEqual(g.types, { protein: 1, miRNA: 1, protein_complex: 1 });
  assert.equal(g.byAspect.P.experimental, 3);
  assert.deepEqual(g.topTerms.P.map((t) => t.id), ['GO:0006915']);
  assert.deepEqual(g.topTerms.C, []);
  assert.deepEqual(g.topProteins.map((p) => p.acc), ['P1']);
});

test('goSummary counts annotations by aspect and evidence and ranks proteins and terms', () => {
  const acc = goSummary(oboTerms(OBO).terms);
  [
    '!gaf-version: 2.2',
    gaf('P1', 'TP53', 'enables', 'GO:0005515', 'IPI', 'F', 'IntAct'),
    gaf('P1', 'TP53', 'involved_in', 'GO:0006915', 'IDA', 'P', 'UniProt'),
    gaf('P1', 'TP53', 'involved_in', 'GO:0006915', 'IEA', 'P', 'InterPro'),
    gaf('P1', 'TP53', 'located_in', 'GO:0005634', 'IBA', 'C', 'GO_Central'),
    gaf('P2', 'MDM2', 'enables', 'GO:0005515', 'IPI', 'F', 'IntAct'),
    gaf('P2', 'MDM2', 'NOT|enables', 'GO:0003677', 'IDA', 'F', 'UniProt'),
    gaf('P3', 'BRCA1', 'involved_in', 'GO:0006281', 'TAS', 'P', 'UniProt'),
    gaf('P3', 'BRCA1', 'involved_in', 'GO:0006281', 'ISS', 'P', 'UniProt'),
  ].forEach((l) => acc.add(l));
  const g = acc.result();
  assert.equal(g.annotations, 8);
  assert.equal(g.not, 1);
  assert.equal(g.proteins, 3);
  assert.deepEqual(g.byAspect.F, {
    experimental: 2, phylogenetic: 0, computational: 0, author: 0, electronic: 0,
  });
  assert.deepEqual(g.byAspect.P, {
    experimental: 1, phylogenetic: 0, computational: 1, author: 1, electronic: 1,
  });
  assert.deepEqual(g.byAspect.C, {
    experimental: 0, phylogenetic: 1, computational: 0, author: 0, electronic: 0,
  });
  assert.deepEqual(g.topProteins[0], { acc: 'P1', symbol: 'TP53', terms: 3 });
  assert.deepEqual(g.topTerms.F[0], { id: 'GO:0005515', name: 'protein binding', n: 2 });
  assert.deepEqual(g.topTerms.P.map((t) => t.id), ['GO:0006281', 'GO:0006915']);
  assert.deepEqual(g.assignedBy.slice(0, 2), [{ group: 'UniProt', n: 4 }, { group: 'IntAct', n: 2 }]);
  assert.equal(g.perProtein.counts[0], 2);
  assert.equal(g.perProtein.counts[2], 1);
});
