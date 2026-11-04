import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  ccreSummary, cpgIslands, tssInIslands, specificityByGene, trnaSummary, encodeSummary,
} from '../../scripts/regulation.mjs';

test('ccreSummary counts cCREs and bases by class and chromosome', () => {
  const acc = ccreSummary();
  [
    'chr1\t100\t300\tEH38D1\tEH38E1\tPLS',
    'chr1\t500\t850\tEH38D2\tEH38E2\tdELS',
    'chr1\t900\t1050\tEH38D3\tEH38E3\tdELS',
    'chrX\t100\t250\tEH38D4\tEH38E4\tCA-CTCF',
    'chrUn_KI270742v1\t1\t200\tEH38D5\tEH38E5\tCA',
  ].forEach((l) => acc.add(l));
  const c = acc.result();
  assert.equal(c.count, 4);
  assert.deepEqual(c.classes, {
    PLS: { count: 1, bp: 200 }, dELS: { count: 2, bp: 500 }, 'CA-CTCF': { count: 1, bp: 150 },
  });
  assert.deepEqual(c.byChrom, { 1: { PLS: 1, dELS: 2 }, X: { 'CA-CTCF': 1 } });
  assert.equal(c.widths.summary.median, 175);
  assert.equal(c.widths.counts[6], 2);
  assert.equal(c.widths.counts[8], 1);
  assert.equal(c.widths.counts[14], 1);
});

const CPG = [
  '585\tchr1\t28735\t29737\tCpG: 111\t1002\t111\t731\t22.2\t73\t0.85',
  '586\tchr1\t135124\t135563\tCpG: 30\t439\t30\t295\t13.7\t67.2\t0.64',
  '1\tchr2\t1000\t3000\tCpG: 200\t2000\t200\t1400\t20\t70\t0.9',
  '1\tchr1_KI270706v1_random\t1\t300\tCpG: 30\t300\t30\t200\t20\t66\t0.7',
].join('\n');

test('cpgIslands summarizes island lengths, CpG ratios and GC content', () => {
  const c = cpgIslands(CPG);
  assert.equal(c.count, 3);
  assert.equal(c.bp, 3441);
  assert.deepEqual(c.byChrom, { 1: 2, 2: 1 });
  assert.equal(c.lengths.summary.median, 1002);
  assert.equal(c.obsExp.counts[5], 1);
  assert.equal(c.obsExp.counts[0], 1);
  assert.equal(c.obsExp.counts[6], 1);
  assert.equal(c.gc.counts[11], 1);
  assert.equal(c.gc.counts[8], 1);
  assert.equal(c.gc.counts[10], 1);
  assert.deepEqual(c.intervals.get('1'), [[28735, 29737], [135124, 135563]]);
});

test('tssInIslands checks each MANE Select start against the islands, by tissue specificity', () => {
  const { intervals } = cpgIslands(CPG);
  const transcripts = [
    {
      gene: 'A', chrom: '1', strand: '+', tag: 'MANE_Select', exons: [[28800, 29000], [30000, 30100]],
    },
    {
      gene: 'B', chrom: '1', strand: '-', tag: 'MANE_Select', exons: [[135500, 135700], [135000, 135100]],
    },
    {
      gene: 'C', chrom: '2', strand: '+', tag: 'MANE_Select', exons: [[2500, 2600]],
    },
    {
      gene: 'D', chrom: '1', strand: '+', tag: 'MANE_Plus_Clinical', exons: [[28800, 29000]],
    },
  ];
  const spec = new Map([['A', 'Low tissue specificity'], ['B', 'Tissue enriched'], ['C', 'Low tissue specificity']]);
  assert.deepEqual(tssInIslands(transcripts, intervals, spec), {
    tss: 3,
    inIsland: 2,
    bySpecificity: {
      'Low tissue specificity': { genes: 2, inIsland: 2 },
      'Tissue enriched': { genes: 1, inIsland: 0 },
    },
  });
});

test('specificityByGene reads the HPA tissue specificity call for each gene symbol', () => {
  const text = ['Gene\tEnsembl\t"RNA tissue specificity"', 'A\tENSG1\t"Low tissue specificity"', 'B\tENSG2\t"Tissue enriched"'].join('\n');
  assert.deepEqual([...specificityByGene(text)], [['A', 'Low tissue specificity'], ['B', 'Tissue enriched']]);
});

test('trnaSummary counts high-confidence tRNA genes by isotype, anticodon and chromosome', () => {
  const fasta = [
    '>Homo_sapiens_tRNA-Ala-AGC-1-1 (tRNAscan-SE ID: chr6.trna116) Ala (AGC) 72 bp mature sequence Sc: 84.9 chr6:28795964-28796035 (-)',
    'GGGGGTGTAGCTCAGTGG',
    '>Homo_sapiens_tRNA-Ala-AGC-10-1 (tRNAscan-SE ID: chr6.trna20) Ala (AGC) 73 bp mature sequence Sc: 60.1 chr6:26687257-26687329 (+)',
    'GGGGAATTAGCTCAAGCGG',
    '>Homo_sapiens_tRNA-iMet-CAT-1-1 (tRNAscan-SE ID: chr1.trna5) iMet (CAT) 72 bp mature sequence Sc: 70.1 chr1:1000-1072 (+)',
    'AGCAGAGTGG',
    '>Homo_sapiens_tRNA-SeC-TCA-1-1 (tRNAscan-SE ID: chr19.trna1) SeC (TCA) 87 bp mature sequence Sc: 120.2 chr19:45478-45565 (+)',
    'GCCCGGATGA',
    '>Homo_sapiens_tRNA-Gly-GCC-1-1 (tRNAscan-SE ID: chr1.trna9) Gly (GCC) 71 bp mature sequence Sc: 70.3 chr1:16861-16932 (-)',
    'GCATTGGTGG',
  ].join('\n');
  assert.deepEqual(trnaSummary(fasta), {
    genes: 5,
    isotypes: {
      Ala: 2, iMet: 1, SeC: 1, Gly: 1,
    },
    anticodons: {
      'Ala-AGC': 2, 'iMet-CAT': 1, 'SeC-TCA': 1, 'Gly-GCC': 1,
    },
    byChrom: { 6: 2, 1: 2, 19: 1 },
  });
});

test('encodeSummary reads ENCODE experiment counts by assay, sample type, target class and phase', () => {
  const e = encodeSummary({
    total: 20450,
    facets: [
      { field: 'assay_title', terms: [{ key: 'TF ChIP-seq', doc_count: 3502 }, { key: 'DNase-seq', doc_count: 3144 }] },
      { field: 'biosample_ontology.classification', terms: [{ key: 'cell line', doc_count: 9781 }] },
      { field: 'award.rfa', terms: [{ key: 'ENCODE4', doc_count: 10367 }] },
      { field: 'target.investigated_as', terms: [{ key: 'histone', doc_count: 4092 }] },
    ],
  });
  assert.deepEqual(e, {
    total: 20450,
    assays: [{ name: 'TF ChIP-seq', n: 3502 }, { name: 'DNase-seq', n: 3144 }],
    biosamples: [{ name: 'cell line', n: 9781 }],
    phases: [{ name: 'ENCODE4', n: 10367 }],
    targets: [{ name: 'histone', n: 4092 }],
  });
});
