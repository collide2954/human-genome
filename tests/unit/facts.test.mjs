import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  SOURCES, SOURCE_TITLES, SOURCE_GROUPS, EXAMPLES, PARS, GENE_THERAPIES, GENE_THERAPIES_CHECKED, DIAGNOSTIC_YIELD, CFDNA_SCREENING, CARRIER_TIER3,
} from '../../src/facts.js';

test('chromosome structure sources are cited', () => {
  for (const key of ['anderson1981', 'andrews1999', 'laricchia2022', 'gnomadMito', 'altemose2022', 'rhie2023', 't2tAnnotation']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Assemblies and sequence', key);
  }
});

test('variation and constraint sources are cited', () => {
  for (const key of ['dbsnp', 'sherry2001', 'dbvar', 'gnomad', 'karczewski2020']) assert.equal(SOURCE_TITLES[key]?.[0], 'Variation and populations', key);
  for (const key of ['clingen', 'rehm2015', 'riggs2020']) assert.equal(SOURCE_TITLES[key]?.[0], 'Clinical genetics', key);
  for (const key of ['impc', 'dickinson2016']) assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
});

test('the GRCh38 pseudoautosomal regions match on X and Y', () => {
  assert.deepEqual(PARS.map((p) => p.name), ['PAR1', 'PAR2']);
  for (const p of PARS) assert.equal(p.x[1] - p.x[0], p.y[1] - p.y[0], p.name);
  assert.deepEqual(PARS[0].x, [10001, 2781479]);
  assert.deepEqual(PARS[1].y, [56887903, 57217415]);
});
import { identify } from '../../src/lib.js';

test('clinical variant sources are cited', () => {
  for (const key of ['clinvar', 'landrum2018', 'landrum2025', 'richards2015', 'strande2017', 'clingenValidity', 'dosageMap', 'riggs2012', 'strchive', 'hiatt2025']) assert.equal(SOURCE_TITLES[key]?.[0], 'Clinical genetics', key);
  for (const key of ['mclaren2016', 'alexandrov2020']) assert.equal(SOURCE_TITLES[key]?.[0], 'Variation and populations', key);
});

test('rare disease and testing sources are cited', () => {
  for (const key of ['hpo', 'gargano2024', 'orphanet', 'nguengang2020', 'gtr', 'rubinstein2013', 'clinicaltrials', 'fdaCgt', 'miller2010', 'clark2018', 'srivastava2019', 'gil2017', 'gregg2021']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Clinical genetics', key);
  }
});

test('gene therapies list each FDA product with its first approval year, vector and group', () => {
  assert.equal(GENE_THERAPIES.length, 34);
  assert.equal(new Set(GENE_THERAPIES.map((g) => g.trade)).size, 34);
  const groups = new Set(['AAV', 'Other in vivo vectors', 'CAR-T and TCR-T cells', 'Blood stem cells', 'Other gene-modified cells']);
  for (const g of GENE_THERAPIES) {
    assert.ok(groups.has(g.group), g.trade);
    assert.ok(g.year >= 2015 && g.year <= 2026, g.trade);
    assert.ok(g.vector && g.condition && g.proper, g.trade);
  }
  assert.deepEqual(GENE_THERAPIES.find((g) => g.trade === 'LUXTURNA'), {
    trade: 'LUXTURNA', proper: 'voretigene neparvovec-rzyl', year: 2017, vector: 'AAV2', group: 'AAV', condition: 'biallelic RPE65 retinal dystrophy',
  });
  assert.equal(GENE_THERAPIES.find((g) => g.trade === 'IMLYGIC').year, 2015);
  assert.equal(GENE_THERAPIES.filter((g) => g.group === 'AAV').length, 11);
  assert.equal(GENE_THERAPIES_CHECKED, '2026-11-04');
});

test('diagnostic yields, cfDNA screening and carrier screening come from cited meta-analyses and guidance', () => {
  for (const y of DIAGNOSTIC_YIELD) {
    assert.ok(['miller2010', 'clark2018', 'srivastava2019'].includes(y.source), y.key);
    assert.ok(y.value !== null || (y.lo !== null && y.hi !== null), y.key);
    if (y.lo !== null) assert.ok(y.lo <= (y.value ?? y.lo) && (y.value ?? y.hi) <= y.hi, y.key);
  }
  assert.deepEqual(DIAGNOSTIC_YIELD.find((y) => y.key === 'wgs'), {
    key: 'wgs', test: 'Genome sequencing', measure: 'diagnostic', value: 41, lo: 34, hi: 48, source: 'clark2018',
  });
  assert.equal(CFDNA_SCREENING.length, 5);
  assert.deepEqual(CFDNA_SCREENING[0], {
    condition: 'Trisomy 21', detection: 99.7, lo: 99.1, hi: 99.9, fpr: 0.04, fprLo: 0.02, fprHi: 0.07, cases: 1963, unaffected: 223932,
  });
  assert.deepEqual(CARRIER_TIER3, { genes: 97, definitive: 39, frequency: '1 in 200' });
});

test('every source is an https URL with a title in a known group', () => {
  assert.ok(Object.keys(SOURCES).length >= 6);
  for (const [key, url] of Object.entries(SOURCES)) {
    assert.match(url, /^https:\/\//, key);
    const [group, title] = SOURCE_TITLES[key] ?? [];
    assert.ok(SOURCE_GROUPS.includes(group), `${key}: ${group}`);
    assert.ok(title?.length > 5 && !/^https?:|^www\./.test(title), key);
  }
  for (const key of Object.keys(SOURCE_TITLES)) assert.ok(SOURCES[key], `title without source: ${key}`);
});

test('protein, target and expression sources are cited in their groups', () => {
  const proteins = ['uniprot', 'uniprot2025', 'interpro', 'pfam2021', 'interpro2025', 'pdb', 'burley2019', 'pharos', 'pharos2023', 'oprea2018'];
  const expression = ['hpa', 'uhlen2015', 'karlsson2021', 'uhlen2019', 'gtex', 'gtex2020'];
  for (const key of proteins) assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
  for (const key of expression) assert.equal(SOURCE_TITLES[key]?.[0], 'Expression and regulation', key);
  for (const key of [...proteins, ...expression]) assert.ok(SOURCES[key], key);
});

test('network, pathway and ontology sources are cited', () => {
  for (const key of ['string', 'szklarczyk2023', 'reactome', 'milacic2024', 'geneOntology', 'go2023', 'ashburner2000']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
    assert.ok(SOURCES[key], key);
  }
});

test('regulation and RNA gene sources are cited', () => {
  for (const key of ['encode', 'encode2020', 'encode2012', 'gardinerGarden1987']) assert.equal(SOURCE_TITLES[key]?.[0], 'Expression and regulation', key);
  for (const key of ['gtrnadb', 'chan2016']) assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
});

test('every decoder example is recognised, unique and described', () => {
  assert.ok(EXAMPLES.length >= 18);
  assert.equal(new Set(EXAMPLES.map((e) => e.id)).size, EXAMPLES.length);
  for (const e of EXAMPLES) {
    assert.ok(identify(e.id), e.id);
    assert.ok(e.is.length > 3, e.id);
  }
  const types = new Set(EXAMPLES.map((e) => identify(e.id).type));
  assert.ok(types.size >= 15, [...types].join(', '));
});
