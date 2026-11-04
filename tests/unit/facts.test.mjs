import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  SOURCES, SOURCE_TITLES, SOURCE_GROUPS, EXAMPLES, PARS,
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
