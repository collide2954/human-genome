import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  mitoFeatures, mitoVariants, censatSummary, telomeres,
} from '../../scripts/structure.mjs';

const FT = [
  '>Feature ref|NC_012920.1|',
  '576\t1\tD-loop',
  '16569\t16024',
  '577\t647\tgene',
  '\t\t\tgene\tTRNF',
  '577\t647\ttRNA',
  '\t\t\tproduct\ttRNA-Phe',
  '648\t1601\tgene',
  '\t\t\tgene\tRNR1',
  '648\t1601\trRNA',
  '\t\t\tproduct\ts-rRNA',
  '14673\t14149\tgene',
  '\t\t\tgene\tND6',
  '14673\t14149\tCDS',
  '\t\t\tproduct\tNADH dehydrogenase subunit 6',
  '14747\t15887\tgene',
  '\t\t\tgene\tCYTB',
  '14747\t15887\tCDS',
  '\t\t\tproduct\tcytochrome b',
].join('\n');

test('mitoFeatures reads genes, their kind and strand, and the D-loop across the origin', () => {
  const m = mitoFeatures(FT);
  assert.deepEqual(m.genes, [
    {
      name: 'TRNF', start: 577, end: 647, strand: '+', kind: 'tRNA', product: 'tRNA-Phe',
    },
    {
      name: 'RNR1', start: 648, end: 1601, strand: '+', kind: 'rRNA', product: 's-rRNA',
    },
    {
      name: 'ND6', start: 14149, end: 14673, strand: '-', kind: 'protein', product: 'NADH dehydrogenase subunit 6',
    },
    {
      name: 'CYTB', start: 14747, end: 15887, strand: '+', kind: 'protein', product: 'cytochrome b',
    },
  ]);
  assert.deepEqual(m.dloop, [[16024, 16569], [1, 576]]);
});

const GNOMAD = [
  'chromosome\tposition\tref\talt\tfilters\tAC_hom\tAC_het\tAF_hom\tAF_het\tAN\tmax_observed_heteroplasmy',
  'chrM\t3\tT\tC\tPASS\t19\t1\t3.3667646E-4\t1.7719814E-5\t56434\t9.9700e-01',
  'chrM\t6\tC\tCCTCAA\tnpg\t0\t0\t0.0\t0.0\t56433\t0.0000e+00',
  'chrM\t73\tA\tG\tPASS\t43000\t50\t7.6E-1\t8.8E-4\t56434\t1.0000e+00',
  'chrM\t150\tC\tT\tPASS\t5000\t20\t8.9E-2\t3.5E-4\t56434\t1.0000e+00',
  'chrM\t150\tC\tA\tPASS\t0\t3\t0.0\t5.3E-5\t56434\t2.1000e-01',
  'chrM\t16519\tT\tC\tPASS\t36000\t100\t6.4E-1\t1.8E-3\t56434\t1.0000e+00',
];

test('mitoVariants counts PASS variants, homoplasmy and heteroplasmy along the genome', () => {
  const acc = mitoVariants({ length: 16569 });
  GNOMAD.forEach((l) => acc.add(l));
  const v = acc.result();
  assert.equal(v.variants, 5);
  assert.equal(v.positions, 4);
  assert.equal(v.homoplasmic, 4);
  assert.equal(v.heteroplasmicOnly, 1);
  assert.equal(v.samples, 56434);
  assert.equal(v.bins.length, 166);
  assert.deepEqual(v.bins[0], { hom: 2, het: 0 });
  assert.deepEqual(v.bins[1], { hom: 1, het: 1 });
  assert.deepEqual(v.bins[165], { hom: 1, het: 0 });
  assert.deepEqual(v.top[0], {
    pos: 73, ref: 'A', alt: 'G', afHom: 0.76, afHet: 0.00088,
  });
  assert.equal(v.top.length, 4);
  assert.equal(v.heteroplasmy.counts.reduce((a, b) => a + b, 0), 5);
  assert.equal(v.heteroplasmy.counts.at(-1), 4);
});

const CENSAT = [
  'chr1\t116796047\t121405145\tct_1_1(p_arm)\t100\t.\t116796047\t121405145\t224,224,224',
  'chr1\t121405145\t121406286\tcensat_1_1(rnd-6_family-4384)\t100\t.\t121405145\t121406286\t0,204,204',
  'chr1\t121619169\t121625213\thor_1_1(S3C1H2-A,B,C)\t100\t.\t121619169\t121625213\t255,102,0',
  'chr1\t121625213\t121667941\thor_1_2(S3C1H2-A,B)\t100\t.\t121625213\t121667941\t255,102,0',
  'chrY\t10000000\t10100000\thsat3_Y_1\t100\t.\t10000000\t10100000\t0,0,0',
  'chrY\t20000000\t20030000\thsat1B_Y_2\t100\t.\t20000000\t20030000\t0,0,0',
  'chr13\t1000\t46000\trDNA_13_1\t100\t.\t1000\t46000\t0,0,0',
].join('\n');

test('censatSummary totals satellite bases by class and chromosome and keeps transitions apart', () => {
  const c = censatSummary(CENSAT);
  assert.deepEqual(c.classes, {
    censat: 1141, hor: 48772, hsat3: 100000, hsat1B: 30000, rDNA: 45000,
  });
  assert.equal(c.transition, 4609098);
  assert.equal(c.total, 224913);
  assert.deepEqual(c.byChrom, { 1: { censat: 1141, hor: 48772 }, Y: { hsat3: 100000, hsat1B: 30000 }, 13: { rDNA: 45000 } });
  assert.deepEqual(c.arrays.hor, [{ chrom: '1', bp: 48772 }]);
});

test('telomeres names each telomere by chromosome arm', () => {
  const t = telomeres(['chr1\t0\t3000', 'chr1\t248384000\t248387328', 'chrX\t0\t5000', 'chrX\t154257000\t154259566'].join('\n'));
  assert.deepEqual(t, [
    { chrom: '1', p: 3000, q: 3328 },
    { chrom: 'X', p: 5000, q: 2566 },
  ]);
});
