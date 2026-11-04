import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  hpoTerms, hpoAnnotations, geneDiseases, orphaGenes, orphaPrevalence, orphaAges, gtrTests, fdaProducts, trialYears,
} from '../../scripts/rare.mjs';

const OBO = [
  'format-version: 1.2', 'data-version: hp/releases/2026-09-01', '',
  '[Term]', 'id: HP:0000001', 'name: All', '',
  '[Term]', 'id: HP:0000118', 'name: Phenotypic abnormality', 'is_a: HP:0000001 ! All', '',
  '[Term]', 'id: HP:0000707', 'name: Abnormality of the nervous system', 'is_a: HP:0000118 ! Phenotypic abnormality', '',
  '[Term]', 'id: HP:0012638', 'name: Abnormal nervous system physiology', 'is_a: HP:0000707 ! Abnormality of the nervous system', '',
  '[Term]', 'id: HP:0001250', 'name: Seizure', 'is_a: HP:0012638 ! Abnormal nervous system physiology', '',
  '[Term]', 'id: HP:0000478', 'name: Abnormality of the eye', 'is_a: HP:0000118 ! Phenotypic abnormality', '',
  '[Term]', 'id: HP:0000505', 'name: Visual impairment', 'is_a: HP:0000478 ! Abnormality of the eye', 'is_a: HP:0012638 ! Abnormal nervous system physiology', '',
  '[Term]', 'id: HP:0000006', 'name: Autosomal dominant inheritance', 'is_a: HP:0000005 ! Mode of inheritance', '',
  '[Term]', 'id: HP:9999999', 'name: obsolete thing', 'is_obsolete: true', '',
  '[Typedef]', 'id: part_of', 'name: part of',
].join('\n');

test('hpoTerms reads the HPO release, live terms with their parents and the organ systems', () => {
  const o = hpoTerms(OBO);
  assert.equal(o.version, '2026-09-01');
  assert.equal(o.terms.size, 8);
  assert.deepEqual(o.terms.get('HP:0000505'), { name: 'Visual impairment', parents: ['HP:0000478', 'HP:0012638'] });
  assert.deepEqual(o.systems, ['HP:0000707', 'HP:0000478']);
});

const HPOA = [
  '#description: "HPO annotations for rare diseases [2: OMIM; 0: DECIPHER; 1 ORPHANET]"', '#version: 2026-09-02',
  'database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration',
  'OMIM:1\tA\t\tHP:0001250\tPMID:1\tPCS\t\t1/2\t\t\tP\tx',
  'OMIM:1\tA\t\tHP:0000505\tPMID:1\tTAS\t\t\t\t\tP\tx',
  'OMIM:1\tA\t\tHP:0000006\tPMID:1\tTAS\t\t\t\t\tI\tx',
  'OMIM:2\tB\t\tHP:0001250\tPMID:2\tIEA\t\t\t\t\tP\tx',
  'OMIM:2\tB\tNOT\tHP:0000505\tPMID:2\tIEA\t\t\t\t\tP\tx',
  'ORPHA:3\tC\t\tHP:0000505\tORPHA:3\tTAS\t\tHP:0040283\t\t\tP\tx',
].join('\n');

test('hpoAnnotations counts diseases, annotations, organ systems, frequent terms and inheritance', () => {
  const a = hpoAnnotations(HPOA, hpoTerms(OBO));
  assert.equal(a.version, '2026-09-02');
  assert.equal(a.diseases, 3);
  assert.deepEqual(a.sources, { OMIM: 2, ORPHA: 1 });
  assert.equal(a.annotations, 4);
  assert.deepEqual(a.evidence, { PCS: 1, TAS: 2, IEA: 1 });
  assert.deepEqual(a.perDisease.counts.slice(0, 3), [2, 1, 0]);
  assert.equal(a.perDisease.edges[0], 1);
  assert.deepEqual(a.systems, [
    { id: 'HP:0000707', name: 'Abnormality of the nervous system', diseases: 3 },
    { id: 'HP:0000478', name: 'Abnormality of the eye', diseases: 2 },
  ]);
  assert.deepEqual(a.topTerms, [
    { id: 'HP:0001250', name: 'Seizure', diseases: 2 },
    { id: 'HP:0000505', name: 'Visual impairment', diseases: 2 },
  ]);
  assert.deepEqual(a.inheritance, [{ id: 'HP:0000006', name: 'Autosomal dominant inheritance', diseases: 1 }]);
});

test('geneDiseases counts Mendelian genes, diseases and diseases per gene', () => {
  const g = geneDiseases([
    'ncbi_gene_id\tgene_symbol\tassociation_type\tdisease_id\tsource',
    'NCBIGene:1\tA\tMENDELIAN\tOMIM:1\tx', 'NCBIGene:1\tA\tMENDELIAN\tOMIM:2\tx', 'NCBIGene:2\tB\tMENDELIAN\tOMIM:1\tx',
    'NCBIGene:3\tC\tPOLYGENIC\tOMIM:3\tx', 'NCBIGene:4\tD\tUNKNOWN\tORPHA:4\tx',
  ].join('\n'));
  assert.deepEqual(g.types, { MENDELIAN: 3, POLYGENIC: 1, UNKNOWN: 1 });
  assert.equal(g.genes, 2);
  assert.equal(g.diseases, 2);
  assert.deepEqual(g.perGene, { bins: ['1', '2', '3-5', '6-10', '>10'], counts: [1, 1, 0, 0, 0] });
  assert.deepEqual(g.top, [{ gene: 'A', diseases: 2 }, { gene: 'B', diseases: 1 }]);
});

const named = (tag, name) => `<${tag} id="1">\n  <Name lang="en">${name}</Name>\n</${tag}>`;
const disorder = (code, group, inner) => `<Disorder id="${code}"><OrphaCode>${code}</OrphaCode><Name lang="en">D${code}</Name>${named('DisorderType', 'Disease')}${named('DisorderGroup', group)}${inner}</Disorder>`;
const jdbor = (...disorders) => `<?xml version="1.0"?>\n<JDBOR date="2026-06-23 07:57:31" version="1.3.42"><DisorderList count="${disorders.length}">${disorders.join('\n')}</DisorderList></JDBOR>`;
const assoc = (symbol, type) => `<DisorderGeneAssociation><Gene id="1"><Name lang="en">a gene</Name><Symbol>${symbol}</Symbol></Gene>${named('DisorderGeneAssociationType', type)}${named('DisorderGeneAssociationStatus', 'Assessed')}</DisorderGeneAssociation>`;

test('orphaGenes counts Orphanet disorders by gene association type and genes per disorder', () => {
  const o = orphaGenes(jdbor(
    disorder(1, 'Disorder', assoc('KIF7', 'Disease-causing germline mutation(s) in') + assoc('ABC', 'Candidate gene tested in')),
    disorder(2, 'Disorder', assoc('KIF7', 'Disease-causing germline mutation(s) (loss of function) in')),
    disorder(3, 'Subtype of disorder', assoc('XYZ', 'Role in the phenotype of')),
  ));
  assert.equal(o.date, '2026-06-23');
  assert.equal(o.disorders, 2);
  assert.equal(o.genes, 3);
  assert.deepEqual(o.associations, {
    'Disease-causing germline mutation(s) in': 1, 'Candidate gene tested in': 1, 'Disease-causing germline mutation(s) (loss of function) in': 1, 'Role in the phenotype of': 1,
  });
  assert.deepEqual(o.perDisorder, { bins: ['1', '2', '3-5', '6-10', '>10'], counts: [1, 1, 0, 0, 0] });
});

const prevalence = (type, cls, geo, status) => `<Prevalence id="1">${named('PrevalenceType', type)}${cls ? named('PrevalenceClass', cls) : '<PrevalenceClass/>'}${named('PrevalenceGeographic', geo)}${named('PrevalenceValidationStatus', status)}</Prevalence>`;

test('orphaPrevalence gives each disorder one point prevalence class, worldwide and validated first', () => {
  const p = orphaPrevalence(jdbor(
    disorder(1, 'Disorder', prevalence('Cases/families', null, 'Worldwide', 'Validated') + prevalence('Point prevalence', '&lt;1 / 1 000 000', 'Worldwide', 'Validated')),
    disorder(2, 'Disorder', prevalence('Point prevalence', '1-9 / 100 000', 'Europe', 'Validated') + prevalence('Point prevalence', '1-5 / 10 000', 'Worldwide', 'Not yet validated')),
    disorder(3, 'Disorder', prevalence('Annual incidence', '1-9 / 1 000 000', 'Worldwide', 'Validated')),
    disorder(4, 'Disorder', prevalence('Cases/families', null, 'Worldwide', 'Validated')),
    disorder(5, 'Subtype of disorder', prevalence('Point prevalence', '&gt;1 / 1000', 'Worldwide', 'Validated')),
  ));
  assert.equal(p.disorders, 4);
  assert.deepEqual(p.classes, {
    '<1 / 1 000 000': 1, '1-5 / 10 000': 1, 'No point prevalence': 1, 'Case reports only': 1,
  });
});

const listOf = (tag, item, names) => `<${tag}List count="${names.length}">${names.map((n) => named(item, n)).join('')}</${tag}List>`;

test('orphaAges counts disorders by age of onset and type of inheritance', () => {
  const a = orphaAges(jdbor(
    disorder(1, 'Disorder', listOf('AverageAgeOfOnset', 'AverageAgeOfOnset', ['Infancy', 'Neonatal']) + listOf('TypeOfInheritance', 'TypeOfInheritance', ['Autosomal recessive'])),
    disorder(2, 'Disorder', listOf('AverageAgeOfOnset', 'AverageAgeOfOnset', ['Childhood']) + listOf('TypeOfInheritance', 'TypeOfInheritance', [])),
    disorder(3, 'Disorder', listOf('AverageAgeOfOnset', 'AverageAgeOfOnset', []) + listOf('TypeOfInheritance', 'TypeOfInheritance', ['Autosomal dominant', 'Not applicable'])),
    disorder(4, 'Group of disorders', listOf('AverageAgeOfOnset', 'AverageAgeOfOnset', ['Adult']) + listOf('TypeOfInheritance', 'TypeOfInheritance', ['Unknown'])),
  ));
  assert.equal(a.disorders, 3);
  assert.deepEqual(a.onset, { Infancy: 1, Neonatal: 1, Childhood: 1 });
  assert.deepEqual(a.inheritance, { 'Autosomal recessive': 1, 'Autosomal dominant': 1, 'Not applicable': 1 });
  assert.equal(a.noOnset, 1);
  assert.equal(a.noInheritance, 1);
});

const GTR_COLS = ['test_accession_ver', 'name_of_laboratory', 'facility_country', 'lab_id', 'test_development', 'method_categories', 'genes', 'now_current', 'test_create_date'];
const gtr = (fields) => GTR_COLS.map((c) => fields[c] ?? '').join('\t');

test('gtrTests counts current GTR tests by lab, country, method, development and genes', () => {
  const g = gtrTests([
    GTR_COLS.join('\t'),
    gtr({
      test_accession_ver: 'GTR000000001.2', facility_country: '"United States"', lab_id: '1', test_development: '"Test developed by laboratory (no manufacturer test name)"', method_categories: '"Sequence analysis of the entire coding region|Deletion/duplication analysis"', genes: '"BRCA1|BRCA2"', now_current: '1', test_create_date: '2016-03-01',
    }),
    gtr({
      test_accession_ver: 'GTR000000002.1', facility_country: 'Turkey', lab_id: '2', method_categories: '"Targeted variant analysis"', genes: 'FMR1', now_current: '1', test_create_date: '2020-05-05',
    }),
    gtr({
      test_accession_ver: 'GTR000000003.1', facility_country: '"United States"', lab_id: '1', test_development: '"FDA-reviewed (has FDA test name)"', method_categories: 'Karyotyping', now_current: '1', test_create_date: '2020-07-07',
    }),
    gtr({
      test_accession_ver: 'GTR000000001.1', facility_country: 'Spain', lab_id: '3', method_categories: 'Karyotyping', genes: 'TP53', now_current: '0', test_create_date: '2014-01-01',
    }),
  ].join('\n'));
  assert.equal(g.tests, 3);
  assert.equal(g.labs, 2);
  assert.deepEqual(g.countries, [{ country: 'United States', tests: 2, labs: 1 }, { country: 'Turkey', tests: 1, labs: 1 }]);
  assert.deepEqual(g.methods.map((m) => m.method), ['Deletion/duplication analysis', 'Karyotyping', 'Sequence analysis of the entire coding region', 'Targeted variant analysis']);
  assert.deepEqual(g.genesPerTest, { bins: ['0', '1', '2-10', '11-100', '101-1000', '>1000'], counts: [1, 1, 1, 0, 0, 0] });
  assert.deepEqual(g.topGenes, [{ gene: 'BRCA1', tests: 1 }, { gene: 'BRCA2', tests: 1 }, { gene: 'FMR1', tests: 1 }]);
  assert.deepEqual(g.development, { 'Test developed by laboratory (no manufacturer test name)': 1, 'FDA-reviewed (has FDA test name)': 1, 'not stated': 1 });
  assert.deepEqual(g.created, { 2016: 1, 2020: 2 });
});

test('fdaProducts lists the products on FDA’s approved cellular and gene therapy page', () => {
  const html = [
    '<header><a href="/home">Home</a></header><main>',
    '<a href="/vaccines-blood-biologics/cellular-gene-therapy-products">Cellular &amp; Gene Therapy Products</a>',
    '<a href="/vaccines-blood-biologics/cellular-gene-therapy-products/approved-cellular-and-gene-therapy-products">Approved Cellular and Gene Therapy Products</a>',
    '<ul><li><a href="/vaccines-blood-biologics/abecma">ABECMA (idecabtagene vicleucel)</a></li>',
    '<li><a href="/vaccines-blood-biologics/cellular-gene-therapy-products/clevecord-hpc-cord-blood">CLEVECORD&nbsp;(HPC Cord Blood)</a></li>',
    '<li><a href="/vaccines-blood-biologics/cellular-gene-therapy-products/hpc-cord-blood-lifesouth">HPC, Cord Blood - LifeSouth</a></li>',
    '<li><a href="/vaccines-blood-biologics/otarmeni">OTARMENI™ (lunsotogene parvec cwha)</a></li>',
    '<li><a href="https://www.fda.gov/vaccines-blood-biologics/tudriqev">TUDRIQEV (vusolimogene oderparepvec-wtpg</a></li></ul>',
    '<a href="/about-fda/contact">Contact FDA</a></main>',
  ].join('\n');
  const p = fdaProducts(html);
  assert.deepEqual(p.map((x) => x.trade), ['ABECMA', 'CLEVECORD', 'HPC, Cord Blood - LifeSouth', 'OTARMENI', 'TUDRIQEV']);
  assert.deepEqual(p[0], { trade: 'ABECMA', name: 'ABECMA (idecabtagene vicleucel)', url: 'https://www.fda.gov/vaccines-blood-biologics/abecma' });
  assert.equal(p[1].name, 'CLEVECORD (HPC Cord Blood)');
  assert.equal(p[4].url, 'https://www.fda.gov/vaccines-blood-biologics/tudriqev');
});

const study = (date, phases, status) => ({
  protocolSection: {
    identificationModule: { nctId: 'NCT1' },
    statusModule: { overallStatus: status, ...(date ? { startDateStruct: { date } } : {}) },
    designModule: { phases },
  },
});

test('trialYears counts studies by start year, highest phase and recruitment', () => {
  const t = trialYears([
    study('2025-08', ['PHASE1', 'PHASE2'], 'RECRUITING'),
    study('2020-09-08', ['PHASE3'], 'COMPLETED'),
    study(null, ['NA'], 'WITHDRAWN'),
    study('2020-01', ['EARLY_PHASE1'], 'RECRUITING'),
  ]);
  assert.equal(t.total, 4);
  assert.deepEqual(t.byYear, { 2020: 2, 2025: 1 });
  assert.equal(t.undated, 1);
  assert.deepEqual(t.phases, {
    'Phase 2': 1, 'Phase 3': 1, 'Not applicable': 1, 'Early phase 1': 1,
  });
  assert.equal(t.recruiting, 2);
});
