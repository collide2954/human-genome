import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  CPIC_GROUPS, SCORE_BINS, cpicPairs, cpicAlleles, cpicPhenotypes, cyp2d6, fdaBiomarkers, clinpgxAnnotations, clinpgxLabels, createdDate,
} from '../../scripts/pgx.mjs';

const pair = (gene, drug, level, clinpgx, testing, removed = false) => ({
  genesymbol: gene, drugid: drug, cpiclevel: level, clinpgxlevel: clinpgx, pgxtesting: testing, removed,
});

test('cpicPairs counts current gene-drug pairs by CPIC level, gene, ClinPGx level and label testing level', () => {
  const p = cpicPairs([
    pair('CYP2D6', 'RxNorm:1', 'A', '1A', 'Actionable PGx'),
    pair('CYP2D6', 'RxNorm:2', 'A', '1A', null),
    pair('CYP2D6', 'RxNorm:3', 'B', '2A', null),
    pair('CYP2C19', 'RxNorm:1', 'A', '1A', 'Testing Recommended'),
    pair('ABCB1', 'RxNorm:4', 'D', '3', null),
    pair('ABCB1', 'RxNorm:5', 'C', null, null, true),
  ]);
  assert.equal(p.pairs, 5);
  assert.equal(p.removed, 1);
  assert.equal(p.genes, 3);
  assert.equal(p.drugs, 4);
  assert.deepEqual(p.levels, { A: 3, B: 1, D: 1 });
  assert.deepEqual(p.byGene, [
    { gene: 'CYP2D6', total: 3, levels: { A: 2, B: 1 } },
    { gene: 'CYP2C19', total: 1, levels: { A: 1 } },
    { gene: 'ABCB1', total: 1, levels: { D: 1 } },
  ]);
  assert.deepEqual(p.clinpgx, { A: { '1A': 3 }, B: { '2A': 1 }, D: { 3: 1 } });
  assert.deepEqual(p.testing, { 'Actionable PGx': 1, 'Testing Recommended': 1, none: 3 });
});

test('cpicAlleles counts curated alleles per gene by clinical function and keeps gene-specific classes', () => {
  const a = cpicAlleles([
    { genesymbol: 'CYP2D6', clinicalfunctionalstatus: 'Normal function' },
    { genesymbol: 'CYP2D6', clinicalfunctionalstatus: 'No function' },
    { genesymbol: 'CYP2D6', clinicalfunctionalstatus: 'Uncertain function' },
    { genesymbol: 'CYP2D6', clinicalfunctionalstatus: null },
    { genesymbol: 'G6PD', clinicalfunctionalstatus: 'II/Deficient' },
    { genesymbol: 'G6PD', clinicalfunctionalstatus: 'IV/Normal' },
    { genesymbol: 'MT-RNR1', clinicalfunctionalstatus: 'normal risk of aminoglycoside-induced hearing loss' },
  ]);
  assert.equal(a.alleles, 7);
  assert.deepEqual(a.genes, [
    {
      gene: 'CYP2D6', total: 4, functions: { 'Normal function': 1, 'No function': 1, 'Uncertain function': 1, 'Not assigned': 1 }, classes: [],
    },
    {
      gene: 'G6PD', total: 2, functions: { 'Gene-specific class': 2 }, classes: ['II/Deficient', 'IV/Normal'],
    },
    {
      gene: 'MT-RNR1', total: 1, functions: { 'Gene-specific class': 1 }, classes: ['normal risk of aminoglycoside-induced hearing loss'],
    },
  ]);
});

const result = (gene, name, priority, frequency, score = null) => ({
  genesymbol: gene, result: name, ehrpriority: priority, activityscore: score, frequency,
});

test('cpicPhenotypes sums the frequency of high-risk phenotypes per gene and biogeographic group', () => {
  const p = cpicPhenotypes([
    result('CYP2C19', 'Poor Metabolizer', 'Abnormal/Priority/High Risk', { European: 0.02, 'East Asian': 0.13 }),
    result('CYP2C19', 'Intermediate Metabolizer', 'Abnormal/Priority/High Risk', { European: 0.26, 'East Asian': 0.46 }),
    result('CYP2C19', 'Normal Metabolizer', 'Normal/Routine/Low Risk', { European: 0.4, 'East Asian': 0.38 }),
    result('CYP2C19', 'Indeterminate', 'none', { European: 0, 'East Asian': null }),
    result('NAT2', 'Poor Acetylator', 'Priority/High Risk', { European: 0.5, Latino: null }),
    result('DPYD', 'Poor Metabolizer', 'Abnormal/Priority/High Risk', null),
  ]);
  assert.deepEqual(p.groups, CPIC_GROUPS);
  assert.deepEqual(p.genes, [
    {
      gene: 'CYP2C19', results: 4, highRisk: ['Poor Metabolizer', 'Intermediate Metabolizer'], shares: { European: 0.28, 'East Asian': 0.59 },
    },
    {
      gene: 'NAT2', results: 1, highRisk: ['Poor Acetylator'], shares: { European: 0.5 },
    },
  ]);
});

const allele = (name, fn, activity, frequency) => ({
  name, clinicalfunctionalstatus: fn, activityvalue: activity, frequency,
});

test('cyp2d6 ranks alleles by frequency and sums phenotypes, activity scores and allele function by group', () => {
  const c = cyp2d6([
    allele('*1', 'Normal function', '1.0', { European: 0.3, 'East Asian': 0.25 }),
    allele('*4', 'No function', '0.0', { European: 0.19, 'East Asian': 0.005 }),
    allele('*10', 'Decreased function', '0.25', { European: 0.02, 'East Asian': 0.43 }),
    allele('*1x2', 'Increased function', '2.0', { European: 0.01, 'East Asian': null }),
    allele('*99', 'Unknown function', 'n/a', { European: 0, 'East Asian': 0 }),
    allele('*1.002', null, null, null),
  ], [
    result('CYP2D6', 'Poor Metabolizer', 'Abnormal/Priority/High Risk', { European: 0.06, 'East Asian': 0.01 }, '0.0'),
    result('CYP2D6', 'Intermediate Metabolizer', 'Abnormal/Priority/High Risk', { European: 0.07, 'East Asian': 0.08 }, '0.25'),
    result('CYP2D6', 'Intermediate Metabolizer', 'Abnormal/Priority/High Risk', { European: 0.28, 'East Asian': 0.07 }, '1.0'),
    result('CYP2D6', 'Normal Metabolizer', 'Normal/Routine/ Low Risk', { European: 0.32, 'East Asian': 0.15 }, '2.0'),
    result('CYP2D6', 'Ultrarapid Metabolizer', 'Abnormal/Priority/High Risk', { European: 0.02, 'East Asian': 0.006 }, '3.0'),
    result('CYP2D6', 'Ultrarapid Metabolizer', 'Abnormal/Priority/High Risk', { European: 0.0005, 'East Asian': 0.0009 }, '≥3.25'),
    result('CYP2D6', 'Indeterminate', 'none', { European: 0.04, 'East Asian': 0.07 }, 'n/a'),
  ]);
  assert.equal(c.alleles, 6);
  assert.deepEqual(c.functions, {
    'Normal function': 1, 'No function': 1, 'Decreased function': 1, 'Increased function': 1, 'Unknown function': 1, 'Not assigned': 1,
  });
  assert.deepEqual(c.top.map((a) => a.name), ['*10', '*1', '*4', '*1x2']);
  assert.deepEqual(c.top[3], {
    name: '*1x2', function: 'Increased function', activity: '2.0', freqs: { European: 0.01, 'East Asian': null },
  });
  assert.deepEqual(c.phenotypes.European, {
    'Poor Metabolizer': 0.06, 'Intermediate Metabolizer': 0.35, 'Normal Metabolizer': 0.32, 'Ultrarapid Metabolizer': 0.0205, Indeterminate: 0.04,
  });
  assert.equal(c.scores.bins, SCORE_BINS);
  assert.deepEqual(SCORE_BINS, ['0', '0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', '2', '2.25', '2.5', '2.75', '≥3', 'n/a']);
  assert.deepEqual(c.scores.groups['East Asian'], [0.01, 0.08, 0, 0, 0.07, 0, 0, 0, 0.15, 0, 0, 0, 0.0069, 0.07]);
  assert.deepEqual(c.mix.European, {
    'Normal function': 0.3, 'No function': 0.19, 'Decreased function': 0.02, 'Increased function': 0.01, 'Unknown function': 0,
  });
  assert.deepEqual(c.mix['East Asian'], {
    'Normal function': 0.25, 'No function': 0.005, 'Decreased function': 0.43, 'Unknown function': 0,
  });
});

const FDA_PAGE = [
  '<li><h2>Content current as of:</h2><p><time datetime="2026-08-12T00:00:00Z">08/12/2026</time></p></li>',
  '<table class="table"><thead><tr><th scope="col">Drug</th><th>Therapeutic Area*</th><th>Biomarker<sup>†</sup></th><th>Labeling Sections</th></tr></thead><tbody>',
  '<tr><td><a href="x"><strong>Abacavir</strong></a></td><td>Infectious Diseases</td><td>HLA-B</td><td>Boxed Warning, Dosage and Administration, Contraindications, Warnings and Precautions</td></tr>',
  '<tr><td><a href="x"><strong>Abemaciclib (1)</strong></a></td><td>Oncology</td><td>ESR<br>(Hormone Receptor)</td><td>Indications and Usage, Adverse Reactions, Clinical Studies</td></tr>',
  '<tr><td><a href="x"><strong>Abemaciclib (2)</strong></a></td><td>Onoclogy</td><td>ERBB2<br>(HER2)</td><td>Indications and Usage, Clinical Studies</td></tr>',
  '<tr><td><a href="x"><strong>Codeine</strong></a></td><td>Anesthesiology</td><td>CYP2D6</td><td>Boxed Warning, Warnings and Precautions</td></tr>',
  '</tbody></table>',
].join('\n');

test('fdaBiomarkers counts drug-biomarker rows by therapeutic area, biomarker and label section', () => {
  const f = fdaBiomarkers(FDA_PAGE);
  assert.equal(f.rows, 4);
  assert.equal(f.drugs, 3);
  assert.equal(f.biomarkers, 4);
  assert.equal(f.updated, '2026-08-12');
  assert.deepEqual(f.areas, [{ area: 'Oncology', rows: 2 }, { area: 'Anesthesiology', rows: 1 }, { area: 'Infectious Diseases', rows: 1 }]);
  assert.equal(f.folded, 1);
  assert.deepEqual(f.top, [
    { biomarker: 'CYP2D6', rows: 1 }, { biomarker: 'ERBB2 (HER2)', rows: 1 }, { biomarker: 'ESR (Hormone Receptor)', rows: 1 }, { biomarker: 'HLA-B', rows: 1 },
  ]);
  assert.deepEqual(f.sections, [
    { section: 'Boxed Warning', rows: 2 }, { section: 'Clinical Studies', rows: 2 }, { section: 'Indications and Usage', rows: 2 }, { section: 'Warnings and Precautions', rows: 2 },
    { section: 'Adverse Reactions', rows: 1 }, { section: 'Contraindications', rows: 1 }, { section: 'Dosage and Administration', rows: 1 },
  ]);
  assert.equal(fdaBiomarkers('<p>no table</p>'), null);
});

const ANNOTATION_HEAD = ['Clinical Annotation ID', 'Variant/Haplotypes', 'Gene', 'Level of Evidence', 'Level Override', 'Level Modifiers', 'Score', 'Phenotype Category', 'PMID Count', 'Evidence Count', 'Drug(s)', 'Phenotype(s)', 'Latest History Date (YYYY-MM-DD)', 'URL', 'Specialty Population'];
const annotation = (id, gene, level, category, population = '') => [id, 'rs1', gene, level, '', '', '1', category, '1', '1', 'drug', '', '2024-01-01', 'u', population].join('\t');

test('clinpgxAnnotations counts clinical annotations by level of evidence, phenotype category and gene', () => {
  const a = clinpgxAnnotations([
    ANNOTATION_HEAD.join('\t'),
    annotation('1', 'CYP2D6', '1A', 'Metabolism/PK'),
    annotation('2', 'CYP2D6', '3', 'Efficacy;Toxicity', 'Pediatric'),
    annotation('3', 'ABCB1', '3', 'Efficacy'),
    annotation('4', 'CYP2C9;VKORC1', '1A', 'Dosage'),
  ].join('\n'));
  assert.equal(a.annotations, 4);
  assert.deepEqual(a.levels, { '1A': 2, 3: 2 });
  assert.deepEqual(a.categories, {
    'Metabolism/PK': { '1A': 1 }, Efficacy: { 3: 2 }, Toxicity: { 3: 1 }, Dosage: { '1A': 1 },
  });
  assert.deepEqual(a.genes, [
    { gene: 'CYP2D6', total: 2, levels: { '1A': 1, 3: 1 } },
    { gene: 'ABCB1', total: 1, levels: { 3: 1 } },
    { gene: 'CYP2C9', total: 1, levels: { '1A': 1 } },
    { gene: 'VKORC1', total: 1, levels: { '1A': 1 } },
  ]);
  assert.equal(a.pediatric, 1);
});

const LABEL_HEAD = ['PharmGKB ID', 'Name', 'Source', 'Biomarker Flag', 'Testing Level', 'Has Prescribing Info', 'Has Dosing Info', 'Has Alternate Drug', 'Has Other Prescribing Guidance', 'Cancer Genome', 'Prescribing', 'Chemicals', 'Genes', 'Variants/Haplotypes', 'Latest History Date (YYYY-MM-DD)'];
const label = (id, source, flag, level, cancer, genes) => [id, `Annotation of ${source} Label`, source, flag, level, '', '', '', '', cancer, '', 'drug', genes, '', '2025-01-01'].join('\t');

test('clinpgxLabels counts drug-label annotations by agency and testing level', () => {
  const l = clinpgxLabels([
    LABEL_HEAD.join('\t'),
    label('PA1', 'FDA', 'On FDA Biomarker List', 'Testing Required', 'Cancer Genome', 'ERBB2'),
    label('PA2', 'EMA', '', 'Actionable PGx', '', 'CYP2D6; CYP2C19'),
    label('PA3', 'FDA', '', '', '', 'CYP2D6'),
    label('PA4', 'PMDA', '', 'Informative PGx', '', 'G6PD'),
  ].join('\n'));
  assert.equal(l.labels, 4);
  assert.deepEqual(l.sources, {
    FDA: { 'Testing Required': 1, 'No level given': 1 }, EMA: { 'Actionable PGx': 1 }, PMDA: { 'Informative PGx': 1 },
  });
  assert.equal(l.genes, 4);
  assert.deepEqual(l.top, [{ gene: 'CYP2D6', labels: 2 }, { gene: 'CYP2C19', labels: 1 }, { gene: 'ERBB2', labels: 1 }, { gene: 'G6PD', labels: 1 }]);
  assert.equal(l.cancer, 1);
  assert.equal(l.onFdaList, 1);
});

test('createdDate reads the date from the CREATED file inside a ClinPGx download', () => {
  assert.equal(createdDate(['LICENSE.txt', 'CREATED_2025-07-05.txt', 'clinical_annotations.tsv']), '2025-07-05');
  assert.equal(createdDate(['LICENSE.txt']), null);
});
