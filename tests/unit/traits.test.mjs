import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  traitParents, gwasAssociations, studySizes, methodFamily, pgsScores, pgsAncestry, pgsCategories, pgsPerformance, heritability,
} from '../../scripts/traits.mjs';

const MAPPINGS = [
  'Disease trait\tEFO term\tEFO URI\tParent term\tParent URI',
  'Height\tbody height\thttp://www.ebi.ac.uk/efo/EFO_0004339\tBody measurement\thttp://www.ebi.ac.uk/efo/EFO_0004324',
  'Type 1 diabetes\ttype 1 diabetes mellitus\thttp://www.ebi.ac.uk/efo/EFO_0001359\tMetabolic disorder\thttp://www.ebi.ac.uk/efo/EFO_0000589',
].join('\n');

test('traitParents maps EFO URIs and term labels to GWAS Catalog parent categories', () => {
  const m = traitParents(MAPPINGS);
  assert.equal(m.uris.get('http://www.ebi.ac.uk/efo/EFO_0004339'), 'Body measurement');
  assert.equal(m.labels.get('type 1 diabetes mellitus'), 'Metabolic disorder');
});

const ASSOC_HEAD = ['DATE ADDED TO CATALOG', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'JOURNAL', 'LINK', 'STUDY', 'DISEASE/TRAIT', 'INITIAL SAMPLE SIZE', 'REPLICATION SAMPLE SIZE', 'REGION', 'CHR_ID', 'CHR_POS', 'REPORTED GENE(S)', 'MAPPED_GENE', 'UPSTREAM_GENE_ID', 'DOWNSTREAM_GENE_ID', 'SNP_GENE_IDS', 'UPSTREAM_GENE_DISTANCE', 'DOWNSTREAM_GENE_DISTANCE', 'STRONGEST SNP-RISK ALLELE', 'SNPS', 'MERGED', 'SNP_ID_CURRENT', 'CONTEXT', 'INTERGENIC', 'RISK ALLELE FREQUENCY', 'P-VALUE', 'PVALUE_MLOG', 'P-VALUE (TEXT)', 'OR or BETA', '95% CI (TEXT)', 'PLATFORM [SNPS PASSING QC]', 'CNV', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'STUDY ACCESSION', 'GENOTYPING TECHNOLOGY'];
const assoc = (a) => {
  const cells = {
    'DATE ADDED TO CATALOG': a.added ?? '2010-06-01', DATE: a.date, 'CHR_ID': a.chrom, 'CHR_POS': a.pos, MAPPED_GENE: a.gene, SNPS: a.snp, CONTEXT: a.context ?? '', 'RISK ALLELE FREQUENCY': a.raf ?? 'NR', PVALUE_MLOG: a.mlog, 'OR or BETA': a.effect ?? '', '95% CI (TEXT)': a.ci ?? '', MAPPED_TRAIT: a.trait, MAPPED_TRAIT_URI: a.uri, 'STUDY ACCESSION': a.study,
  };
  return ASSOC_HEAD.map((h) => cells[h] ?? 'x').join('\t');
};
const EFO = 'http://www.ebi.ac.uk/efo/';

test('gwasAssociations bins the strongest hit per megabase and counts years, categories, traits, genes, context and disease odds ratios', () => {
  const acc = gwasAssociations({ parents: traitParents(MAPPINGS) });
  for (const line of [
    ASSOC_HEAD.join('\t'),
    assoc({
      study: 'GCST1', date: '2010-05-01', chrom: '1', pos: '1500000', mlog: '9.5', trait: 'body height', uri: `${EFO}EFO_0004339`, gene: 'GENEA', context: 'intron_variant', raf: '0.3', effect: '0.05', ci: '[0.03-0.07] unit increase', snp: 'rs1',
    }),
    assoc({
      study: 'GCST1', date: '2010-05-01', chrom: '1', pos: '1900000', mlog: '20', trait: 'body height', uri: `${EFO}EFO_0004339`, gene: 'GENEB - GENEC', context: 'intergenic_variant', snp: 'rs2',
    }),
    assoc({
      study: 'GCST2', date: '2015-01-02', chrom: '6', pos: '32000000', mlog: '350', trait: 'type 1 diabetes mellitus', uri: 'http://purl.obolibrary.org/obo/MONDO_0005147', gene: 'HLA-DQA1', context: 'missense_variant', raf: '0.15', effect: '3.5', ci: '[3.1-3.9]', snp: 'rs3',
    }),
    assoc({
      study: 'GCST3', date: '2015-03-03', chrom: '6 x 6', pos: '1 x 2', mlog: '6', trait: 'asthma', uri: `${EFO}EFO_0000270`, gene: 'GENED', context: 'intron_variant x intron_variant', snp: 'rs4',
    }),
    assoc({
      study: 'GCST3', date: '2015-03-03', chrom: 'X', pos: '100', mlog: '8', trait: 'asthma', uri: `${EFO}EFO_0000270`, gene: 'GENEE', effect: '1.2', snp: 'rs5',
    }),
    assoc({
      study: 'GCST2', date: '2015-01-02', added: '2016-01-01', chrom: '6', pos: '32500000', mlog: '50', trait: 'type 1 diabetes mellitus', uri: 'http://purl.obolibrary.org/obo/MONDO_0005147', gene: 'HLA-DQA1', context: 'missense_variant; intron_variant', raf: '0.04', effect: '2.5', ci: '[2.0-3.0]', snp: 'rs6',
    }),
  ]) acc.add(line);
  const g = acc.result();
  assert.equal(g.associations, 6);
  assert.equal(g.studies, 3);
  assert.equal(g.significant, 5);
  assert.equal(g.placed, 5);
  assert.equal(g.latestAdded, '2016-01-01');
  assert.deepEqual(g.byYear, { 2010: { associations: 2, studies: 1 }, 2015: { associations: 4, studies: 2 } });
  assert.deepEqual(g.manhattan, {
    1: [[1, 20, 2, 'body height', 'GENEB - GENEC', 'rs2']],
    6: [[32, 350, 2, 'type 1 diabetes mellitus', 'HLA-DQA1', 'rs3']],
    X: [[0, 8, 1, 'asthma', 'GENEE', 'rs5']],
  });
  assert.deepEqual(g.categories.totals, { 'Body measurement': 2, 'Metabolic disorder': 2, 'Not mapped': 2 });
  assert.deepEqual(g.categories.byYear, { 2010: { 'Body measurement': 2 }, 2015: { 'Metabolic disorder': 2, 'Not mapped': 2 } });
  assert.deepEqual(g.topTraits, [
    { trait: 'asthma', n: 2, parent: 'Not mapped' }, { trait: 'body height', n: 2, parent: 'Body measurement' }, { trait: 'type 1 diabetes mellitus', n: 2, parent: 'Metabolic disorder' },
  ]);
  assert.deepEqual(g.topGenes, [
    { gene: 'HLA-DQA1', n: 2 }, { gene: 'GENEA', n: 1 }, { gene: 'GENEB', n: 1 }, { gene: 'GENEC', n: 1 }, { gene: 'GENED', n: 1 }, { gene: 'GENEE', n: 1 },
  ]);
  assert.deepEqual(g.context, {
    intron_variant: 1, intergenic_variant: 1, missense_variant: 1, 'Several variants': 1, Interaction: 1, 'Not given': 1,
  });
  assert.deepEqual(g.effects.orEdges, [1, 1.1, 1.25, 1.5, 2, 3, 5, 10]);
  assert.equal(g.effects.rafEdges.length, 11);
  assert.equal(g.effects.n, 2);
  assert.equal(g.effects.median, 3);
  assert.equal(g.effects.counts[5][1], 1);
  assert.equal(g.effects.counts[4][0], 1);
  assert.equal(g.effects.counts.flat().reduce((a, b) => a + b, 0), 2);
  assert.deepEqual([...acc.perStudy], [['GCST1', 2], ['GCST2', 2], ['GCST3', 2]]);
});

const ANCESTRY_HEAD = ['STUDY ACCESSION', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'INITIAL SAMPLE DESCRIPTION', 'REPLICATION SAMPLE DESCRIPTION', 'STAGE', 'NUMBER OF INDIVDUALS', 'BROAD ANCESTRAL CATEGORY', 'COUNTRY OF ORIGIN', 'COUNTRY OF RECRUITMENT', 'ADDITONAL ANCESTRY DESCRIPTION'];
const ancestryRow = (study, stage, n) => [study, '1', 'A', '2015-01-01', 'x', 'x', stage, n, 'European', 'NR', 'NR', ''].join('\t');

test('studySizes crosses each study’s discovery sample with its number of associations', () => {
  const s = studySizes(new Map([['GCST1', 2], ['GCST2', 2], ['GCST3', 2]]), [
    ANCESTRY_HEAD.join('\t'),
    ancestryRow('GCST1', 'initial', '5000'),
    ancestryRow('GCST1', 'replication', '2000'),
    ancestryRow('GCST2', 'initial', '120000'),
    ancestryRow('GCST2', 'initial', '30000'),
    ancestryRow('GCST3', 'initial', '800'),
    ancestryRow('GCST4', 'initial', '2000000'),
  ].join('\n'));
  assert.deepEqual(s.sizeLabels, ['<1K', '1K–3K', '3K–10K', '10K–30K', '30K–100K', '100K–300K', '300K–1M', '≥1M']);
  assert.deepEqual(s.assocLabels, ['0', '1', '2–3', '4–10', '11–30', '31–100', '101–300', '301–1,000', '>1,000']);
  assert.equal(s.studies, 4);
  assert.equal(s.counts[2][2], 1);
  assert.equal(s.counts[2][5], 1);
  assert.equal(s.counts[2][0], 1);
  assert.equal(s.counts[0][7], 1);
  assert.equal(s.counts.flat().reduce((a, b) => a + b, 0), 4);
});

test('methodFamily groups PGS development methods into families', () => {
  assert.equal(methodFamily('snpnet'), 'Penalized regression');
  assert.equal(methodFamily('Lassosum2'), 'Penalized regression');
  assert.equal(methodFamily('Penalized regression (bigstatsr)'), 'Penalized regression');
  assert.equal(methodFamily('LDpred2-auto'), 'LDpred');
  assert.equal(methodFamily('PRS-CSx'), 'PRS-CS');
  assert.equal(methodFamily('SBayesR'), 'SBayes');
  assert.equal(methodFamily('MegaPRS'), 'MegaPRS and LDAK');
  assert.equal(methodFamily('Pruning and Thresholding (P+T)'), 'Pruning and thresholding');
  assert.equal(methodFamily('PRSice-2'), 'Pruning and thresholding');
  assert.equal(methodFamily('Genome-wide significant variants'), 'Genome-wide significant variants');
  assert.equal(methodFamily('GWAS Hits'), 'Genome-wide significant variants');
  assert.equal(methodFamily('GenoBoost'), 'Other');
});

const LICENSE = '"PGS obtained from the Catalog should be cited appropriately, and used in accordance with any licensing restrictions"';
const SCORE_HEAD = 'Polygenic Score (PGS) ID,PGS Name,Reported Trait,Mapped Trait(s) (EFO label),Mapped Trait(s) (EFO ID),PGS Development Method,PGS Development Details/Relevant Parameters,Original Genome Build,Number of Variants,Number of Interaction Terms,Type of Variant Weight,PGS Publication (PGP) ID,Publication (PMID),Publication (doi),Score and results match the original publication,Ancestry Distribution (%) - Source of Variant Associations (GWAS),Ancestry Distribution (%) - Score Development/Training,Ancestry Distribution (%) - PGS Evaluation,FTP link,Release Date,License/Terms of Use';
const score = (id, label, efo, method, variants, release) => [id, 'name', 'trait', label, efo, method, 'P<5x10-8', 'NR', variants, '0', 'NR', 'PGP1', '1', 'doi', 'True', 'European:100', '', 'European:72.7|East Asian:27.3', 'ftp', release, LICENSE].join(',');
const SCORES = [
  SCORE_HEAD,
  score('PGS1', 'breast carcinoma', 'MONDO_0004989', 'PRS-CS', '1000000', '2021-03-01'),
  score('PGS2', 'breast carcinoma', 'MONDO_0004989', 'LDpred2', '900', '2021-05-01'),
  score('PGS3', 'type 2 diabetes mellitus|body mass index', 'MONDO_0005148|EFO_0004340', 'snpnet', '12000', '2026-02-01'),
  score('PGS4', 'type 2 diabetes mellitus', 'MONDO_0005148', 'Genome-wide significant SNPs', '50', '2026-03-01'),
  score('PGS5', 'body height', 'EFO_0004339', 'GenoBoost', '7', '2019-10-14'),
].join('\n');

test('pgsScores counts scores by release year and method family, variants per score and traits', () => {
  const p = pgsScores(SCORES);
  assert.equal(p.scores, 5);
  assert.equal(p.methods, 5);
  assert.deepEqual(p.families, {
    'PRS-CS': 1, LDpred: 1, 'Penalized regression': 1, 'Genome-wide significant variants': 1, Other: 1,
  });
  assert.deepEqual(p.byYear, { 2019: { Other: 1 }, 2021: { 'PRS-CS': 1, LDpred: 1 }, 2026: { 'Penalized regression': 1, 'Genome-wide significant variants': 1 } });
  assert.equal(p.variants.edges[0], 1);
  assert.equal(p.variants.counts.reduce((a, b) => a + b, 0), 5);
  assert.equal(p.variants.median, 900);
  assert.deepEqual(p.traits, [
    { trait: 'breast carcinoma', n: 2 }, { trait: 'type 2 diabetes mellitus', n: 2 }, { trait: 'body height', n: 1 }, { trait: 'body mass index', n: 1 },
  ]);
});

const DEV_HEAD = 'Polygenic Score (PGS) ID,Stage of PGS Development,Number of Individuals,Number of Cases,Number of Controls,Percent of Participants Who are Male,Sample Age,Broad Ancestry Category,"Ancestry (e.g. French, Chinese)",Country of Recruitment,Additional Ancestry Description,Phenotype Definitions and Methods,Followup Time,GWAS Catalog Study ID (GCST...),Source PubMed ID (PMID),Source DOI,Cohort(s),Additional Sample/Cohort Information';
const EVAL_HEAD = 'PGS Sample Set (PSS),Polygenic Score (PGS) ID,Number of Individuals,Number of Cases,Number of Controls,Percent of Participants Who are Male,Sample Age,Broad Ancestry Category,"Ancestry (e.g. French, Chinese)",Country of Recruitment,Additional Ancestry Description,Phenotype Definitions and Methods,Followup Time,GWAS Catalog Study ID (GCST...),Source PubMed ID (PMID),Source DOI,Cohort(s),Additional Sample/Cohort Information';

test('pgsAncestry sums individuals by broad ancestry for GWAS, development and evaluation samples', () => {
  const a = pgsAncestry([
    DEV_HEAD,
    'PGS1,Source of Variant Associations (GWAS),22627.0,,,,,European,,"Finland, Sweden",,,,GCST1,1,,,',
    'PGS1,Score Development/Training,1000,,,,,East Asian,,,,,,,,,,',
    'PGS2,Source of Variant Associations (GWAS),5000,,,,,"European, East Asian",,,,,,,,,,',
    'PGS3,Source of Variant Associations (GWAS),,,,,,African unspecified,,,,,,,,,,',
  ].join('\n'), [
    EVAL_HEAD,
    'PSS1,PGS1,67054,33673.0,33381.0,,,European,,"Australia, Belarus",,All breast cancer,,,,,ABCFS|MCCS,iCOGS',
    'PSS2,PGS2,2000,,,,,Not reported,,,,,,,,,,',
  ].join('\n'));
  assert.deepEqual(a, {
    gwas: { European: 22627, 'Multiple categories': 5000 },
    development: { 'East Asian': 1000 },
    evaluation: { European: 67054, 'Not reported': 2000 },
  });
});

test('pgsCategories maps each EFO trait to its PGS Catalog trait category', () => {
  const c = pgsCategories({ results: [{ label: 'Cancer', efotraits: [{ id: 'MONDO_0004989' }] }, { label: 'Metabolic disorder', efotraits: [{ id: 'MONDO_0005148' }] }] });
  assert.equal(c.get('MONDO_0004989'), 'Cancer');
  assert.equal(c.get('MONDO_0005148'), 'Metabolic disorder');
});

const METRIC_HEAD = 'PGS Performance Metric (PPM) ID,Evaluated Score,PGS Sample Set (PSS),PGS Publication (PGP) ID,Reported Trait,Covariates Included in the Model,PGS Performance: Other Relevant Information,Publication (PMID),Publication (doi),Hazard Ratio (HR),Odds Ratio (OR),Beta,Area Under the Receiver-Operating Characteristic Curve (AUROC),Concordance Statistic (C-index),Other Metric(s)';
const metric = (id, pgs, auroc, cindex, other, beta = '') => [id, pgs, 'PSS1', 'PGP1', 'trait', '"age, sex"', '', '1', 'doi', '', '', beta, auroc, cindex, other].join(',');

test('pgsPerformance summarizes reported AUROC, C-index and R² overall and by trait category', () => {
  const p = pgsPerformance([
    METRIC_HEAD,
    metric('PPM1', 'PGS1', '"0.622 [0.619,0.627]"', '', ''),
    metric('PPM2', 'PGS1', '0.702', '', ''),
    metric('PPM3', 'PGS3', '0.58', '0.61', ''),
    metric('PPM4', 'PGS4', '', '', 'R²: 0.12'),
    metric('PPM5', 'PGS5', '', '', '', '0.3'),
  ].join('\n'), SCORES, pgsCategories({ results: [{ label: 'Cancer', efotraits: [{ id: 'MONDO_0004989' }] }, { label: 'Metabolic disorder', efotraits: [{ id: 'MONDO_0005148' }] }] }));
  assert.equal(p.metrics, 5);
  assert.equal(p.auroc.n, 3);
  assert.equal(p.auroc.median, 0.622);
  assert.equal(p.auroc.edges[0], 0.4);
  assert.equal(p.auroc.counts.reduce((a, b) => a + b, 0), 3);
  assert.deepEqual(p.auroc.byCategory, [
    {
      category: 'Cancer', n: 2, q1: 0.642, median: 0.662, q3: 0.682,
    },
    {
      category: 'Metabolic disorder', n: 1, q1: 0.58, median: 0.58, q3: 0.58,
    },
  ]);
  assert.deepEqual(p.cindex, { n: 1, median: 0.61 });
  assert.deepEqual(p.r2, { n: 1, median: 0.12 });
});

const H2_HEAD = ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier', 'pop', 'estimates.final.h2_liability', 'estimates.final.h2_liability_se', 'qcflags.pass_all'].join('\t');
const PHENO_HEAD = ['trait_type', 'phenocode', 'pheno_sex', 'coding', 'modifier', 'description', 'category'].join('\t');

test('heritability joins Pan-UKB estimates to phenotype names and summarizes the European-ancestry traits that pass QC', () => {
  const h = heritability([
    H2_HEAD,
    'continuous\t50\tboth_sexes\t\tirnt\tEUR\t0.5876\t0.02\ttrue',
    'continuous\t50\tboth_sexes\t\tirnt\tAFR\t0.3\t0.1\tfalse',
    'biomarkers\t30600\tboth_sexes\t\tirnt\tEUR\t0.14\t0.01\ttrue',
    'icd10\tE11\tboth_sexes\t\t\tEUR\t0.05\t0.01\ttrue',
    'phecode\t250.2\tboth_sexes\t\t\tCSA\t0.2\t0.05\ttrue',
    'categorical\t1234\tboth_sexes\t1\t\tEUR\tNA\tNA\tfalse',
  ].join('\n'), [
    PHENO_HEAD,
    'continuous\t50\tboth_sexes\t\tirnt\tStanding height\tUK Biobank Assessment Centre > Physical measures > Anthropometry',
    'biomarkers\t30600\tboth_sexes\t\tirnt\tAlbumin\tBiological samples > Assay results',
    'icd10\tE11\tboth_sexes\t\t\tE11 Non-insulin-dependent diabetes mellitus\tHealth-related outcomes',
    'phecode\t250.2\tboth_sexes\t\t\tType 2 diabetes\tendocrine/metabolic',
    'categorical\t1234\tboth_sexes\t1\t\tSomething\tOnline follow-up',
  ].join('\n'));
  assert.equal(h.phenotypes, 5);
  assert.deepEqual(h.byPop, { EUR: { estimates: 3, passing: 3 }, AFR: { estimates: 1, passing: 0 }, CSA: { estimates: 1, passing: 1 } });
  assert.equal(h.eur.passing, 3);
  assert.equal(h.eur.median, 0.14);
  assert.equal(h.eur.counts.reduce((a, b) => a + b, 0), 3);
  assert.deepEqual(h.eur.byType.map((t) => [t.type, t.n, t.median]), [['continuous', 1, 0.5876], ['biomarkers', 1, 0.14], ['icd10', 1, 0.05]]);
  assert.deepEqual(h.eur.top[0], {
    description: 'Standing height', type: 'continuous', h2: 0.5876, se: 0.02,
  });
  assert.equal(h.eur.top.length, 3);
});
