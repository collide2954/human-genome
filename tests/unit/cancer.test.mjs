import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  latestIntogen, intogenDrivers, intogenCohorts, gdcRelease, gdcPrograms, facetCounts, casesByProject,
} from '../../scripts/cancer.mjs';

test('latestIntogen picks the newest drivers and cohorts files from the download page', () => {
  const html = [
    '<a href="?file=IntOGen-Drivers-20230531.zip">a</a>',
    '<a href="?file=IntOGen-Cohorts-20240920.zip">b</a>',
    '<a href="?file=IntOGen-Drivers-20240920.zip">c</a>',
    '<a href="?file=IntOGen-Cohorts-20230531.zip">d</a>',
    '<a href="?file=intogen_driver_mutations_catalog-2016.5.zip">e</a>',
  ].join('');
  assert.deepEqual(latestIntogen(html), { release: '2024-09-20', drivers: 'IntOGen-Drivers-20240920.zip', cohorts: 'IntOGen-Cohorts-20240920.zip' });
  assert.equal(latestIntogen('<a href="?file=other.zip">x</a>'), null);
});

const DRIVER_HEAD = ['SYMBOL', 'TRANSCRIPT', 'COHORT', 'CANCER_TYPE', 'METHODS', 'MUTATIONS', 'SAMPLES', '%_SAMPLES_COHORT', 'QVALUE_COMBINATION', 'ROLE', 'DOMAINS', '2D_CLUSTERS', '3D_CLUSTERS', 'EXCESS_MIS', 'EXCESS_NON', 'EXCESS_SPL', 'IS_DRIVER', 'TOTAL_SAMPLES'];
const driver = (gene, cohort, type, role, samples, total) => [gene, 'ENST1', cohort, type, 'dndscv', samples, samples, samples / total, '0.001', role, '', '', '', '0', '0', '0', 'True', total].join('\t');

test('intogenDrivers counts driver genes per cancer type and ranks genes by the cancer types they drive, by role', () => {
  const d = intogenDrivers([
    DRIVER_HEAD.join('\t'),
    driver('TP53', 'C1', 'BRCA', 'LoF', 290, 1000),
    driver('TP53', 'C2', 'BRCA', 'LoF', 10, 40),
    driver('TP53', 'C3', 'LUAD', 'LoF', 180, 400),
    driver('KRAS', 'C3', 'LUAD', 'Act', 118, 400),
    driver('KRAS', 'C4', 'PAAD', 'Act', 150, 167),
    driver('NOTCH1', 'C5', 'CLL', 'Act', 20, 200),
    driver('NOTCH1', 'C6', 'HNSC', 'LoF', 30, 200),
    driver('NOTCH1', 'C7', 'HNSC', 'LoF', 5, 50),
  ].join('\n'));
  assert.equal(d.rows, 8);
  assert.equal(d.genes, 3);
  assert.equal(d.cancerTypes, 5);
  assert.equal(d.cohorts, 7);
  assert.deepEqual(d.roles, { LoF: 2, Act: 1 });
  assert.deepEqual(d.byType, [
    { type: 'LUAD', genes: 2 }, { type: 'BRCA', genes: 1 }, { type: 'CLL', genes: 1 }, { type: 'HNSC', genes: 1 }, { type: 'PAAD', genes: 1 },
  ]);
  assert.deepEqual(d.top, [
    {
      gene: 'KRAS', types: 2, Act: 2, LoF: 0, ambiguous: 0,
    },
    {
      gene: 'NOTCH1', types: 2, Act: 1, LoF: 1, ambiguous: 0,
    },
    {
      gene: 'TP53', types: 2, Act: 0, LoF: 2, ambiguous: 0,
    },
  ]);
});

const COHORT_HEAD = ['COHORT', 'CANCER', 'CANCER_NAME', 'SOURCE', 'PLATFORM', 'REFERENCE', 'TYPE', 'TREATED', 'AGE', 'SAMPLES', 'MUTATIONS', 'COHORT_NICK', 'COHORT_NAME', 'DRIVERS'];
const cohort = (id, cancer, name, platform, type, age, samples, mutations) => [id, cancer, name, 'TCGA/PanCatAtlas, phs000178', platform, 'PMID:1', type, 'Untreated', age, samples, mutations, id, `${name} cohort`, '3.0'].join('\t');

test('intogenCohorts sums samples and mutations and ranks cancer types by median mutations per sample', () => {
  const c = intogenCohorts([
    COHORT_HEAD.join('\t'),
    cohort('TCGA_WXS_ACC', 'ACC', 'Adrenocortical Carcinoma', 'WXS', 'Primary', 'Adult', 91, 11660),
    cohort('STJUDE_WGS_D_ACC_2018', 'ACC', 'Adrenocortical Carcinoma', 'WGS', 'Primary', 'Adult', 20, 17560),
    cohort('PEDCBIOP_WXS_ACC_PRY', 'ACC', 'Adrenocortical Carcinoma', 'WXS', 'Primary', 'Pediatric', 8, 239),
    cohort('TCGA_WXS_SKCM', 'MEL', 'Cutaneous Melanoma', 'WXS', 'Metastasis', 'Adult', 100, 120000),
  ].join('\n'));
  assert.equal(c.cohorts, 4);
  assert.equal(c.samples, 219);
  assert.equal(c.mutations, 149459);
  assert.deepEqual(c.platforms, { WXS: 3, WGS: 1 });
  assert.deepEqual(c.kinds, { Primary: 3, Metastasis: 1 });
  assert.deepEqual(c.ages, { Adult: 3, Pediatric: 1 });
  assert.deepEqual(c.names, { ACC: 'Adrenocortical Carcinoma', MEL: 'Cutaneous Melanoma' });
  assert.deepEqual(c.burden, [
    {
      cancer: 'ACC',
      name: 'Adrenocortical Carcinoma',
      median: 128.1,
      cohorts: [
        {
          cohort: 'PEDCBIOP_WXS_ACC_PRY', platform: 'WXS', age: 'Pediatric', samples: 8, perSample: 29.9,
        },
        {
          cohort: 'TCGA_WXS_ACC', platform: 'WXS', age: 'Adult', samples: 91, perSample: 128.1,
        },
        {
          cohort: 'STJUDE_WGS_D_ACC_2018', platform: 'WGS', age: 'Adult', samples: 20, perSample: 878,
        },
      ],
    },
    {
      cancer: 'MEL',
      name: 'Cutaneous Melanoma',
      median: 1200,
      cohorts: [{
        cohort: 'TCGA_WXS_SKCM', platform: 'WXS', age: 'Adult', samples: 100, perSample: 1200,
      }],
    },
  ]);
});

test('gdcRelease reads the data release from the status endpoint', () => {
  assert.deepEqual(gdcRelease({ data_release: 'Data Release 46.0 - August 10, 2026', data_release_version: { major: 46, minor: 0, release_date: '2026-08-10' } }), { release: '46.0', date: '2026-08-10' });
  assert.deepEqual(gdcRelease({}), { release: null, date: null });
});

test('gdcPrograms sums cases and projects by program and lists the TCGA projects', () => {
  const project = (id, name, program, cases) => ({
    project_id: id, name, program: { name: program }, summary: { case_count: cases },
  });
  assert.deepEqual(gdcPrograms({
    data: {
      hits: [
        project('TCGA-LUAD', 'Lung Adenocarcinoma', 'TCGA', 585),
        project('TCGA-BRCA', 'Breast Invasive Carcinoma', 'TCGA', 1098),
        project('FM-AD', 'Foundation Medicine Adult Cancer Clinical Dataset (FM-AD)', 'FM', 18004),
      ],
    },
  }), {
    projects: 3,
    cases: 19687,
    programs: [{ program: 'FM', cases: 18004, projects: 1 }, { program: 'TCGA', cases: 1683, projects: 2 }],
    tcga: [{ project: 'TCGA-BRCA', name: 'Breast Invasive Carcinoma', cases: 1098 }, { project: 'TCGA-LUAD', name: 'Lung Adenocarcinoma', cases: 585 }],
  });
});

test('facetCounts lists the buckets of one GDC facet, largest first', () => {
  const json = { data: { aggregations: { primary_site: { buckets: [{ doc_count: 5, key: 'breast' }, { doc_count: 9, key: 'bronchus and lung' }, { doc_count: 5, key: 'brain' }] } } } };
  assert.deepEqual(facetCounts(json, 'primary_site'), [{ key: 'bronchus and lung', n: 9 }, { key: 'brain', n: 5 }, { key: 'breast', n: 5 }]);
  assert.deepEqual(facetCounts({ data: {} }, 'primary_site'), []);
});

test('casesByProject counts each case once per project', () => {
  const hit = (id, project) => ({ case: { case_id: id, project: { project_id: project } } });
  assert.deepEqual(casesByProject([hit('a', 'TCGA-PAAD'), hit('a', 'TCGA-PAAD'), hit('b', 'TCGA-PAAD'), hit('c', 'TCGA-LUAD')]), { 'TCGA-PAAD': 2, 'TCGA-LUAD': 1 });
});
