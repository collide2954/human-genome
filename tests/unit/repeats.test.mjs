import { test } from 'node:test';
import assert from 'node:assert/strict';
import { rmskSummary, superDupsSummary, tandemSummary } from '../../scripts/parsers.mjs';

const rm = (chrom, start, end, div, name, cls, family) => ['585', '1000', String(div), '0', '0', chrom, String(start), String(end), '-0', '+', name, cls, family, '1', '1', '0', '1'].join('\t');

const RMSK = [
  rm('chr1', 10000, 10468, 13, '(TAACCC)n', 'Simple_repeat', 'Simple_repeat'),
  rm('chr1', 20000, 20300, 30, 'AluY', 'SINE', 'Alu'),
  rm('chr2', 1000, 1310, 120, 'AluSx', 'SINE', 'Alu'),
  rm('chr3', 5000, 11100, 5, 'L1HS', 'LINE', 'L1'),
  rm('chr3', 20000, 20500, 250, 'L1M5', 'LINE', 'L1'),
  rm('chrX', 100, 300, 200, 'MER5A', 'DNA?', 'hAT-Charlie?'),
  rm('chr1_KI270706v1_random', 1, 50, 100, 'AluY', 'SINE', 'Alu'),
];

test('rmskSummary totals classes and families on the 24 chromosomes', () => {
  const acc = rmskSummary();
  RMSK.forEach((l) => acc.add(l));
  const r = acc.result();
  assert.deepEqual(r.total, { count: 6, bp: 7878 });
  assert.deepEqual(r.classes, [
    { cls: 'LINE', count: 2, bp: 6600 },
    { cls: 'SINE', count: 2, bp: 610 },
    { cls: 'Simple_repeat', count: 1, bp: 468 },
    { cls: 'DNA', count: 1, bp: 200 },
  ]);
  assert.deepEqual(r.families[0], { cls: 'LINE', family: 'L1', count: 2, bp: 6600 });
  assert.deepEqual(r.families.at(-1), { cls: 'DNA', family: 'hAT-Charlie', count: 1, bp: 200 });
  assert.deepEqual(r.perChrom, { 1: 768, 2: 310, 3: 6600, X: 200 });
});

test('rmskSummary bins transposon bases by divergence from their consensus', () => {
  const acc = rmskSummary();
  RMSK.forEach((l) => acc.add(l));
  const { landscape } = acc.result();
  assert.deepEqual(landscape.classes, ['LINE', 'SINE', 'LTR', 'DNA']);
  assert.equal(landscape.edges.length, 51);
  assert.equal(landscape.bp.LINE[0], 6100);
  assert.equal(landscape.bp.LINE[25], 500);
  assert.equal(landscape.bp.SINE[3], 300);
  assert.equal(landscape.bp.SINE[12], 310);
  assert.equal(landscape.bp.DNA[20], 200);
  assert.equal(landscape.bp.LTR.reduce((a, b) => a + b, 0), 0);
});

test('rmskSummary counts Alu subfamilies and full-length L1HS', () => {
  const acc = rmskSummary();
  RMSK.forEach((l) => acc.add(l));
  const r = acc.result();
  assert.deepEqual(r.alu, { AluY: { count: 1, bp: 300 }, AluS: { count: 1, bp: 310 }, AluJ: { count: 0, bp: 0 } });
  assert.deepEqual(r.l1hs, { count: 1, fullLength: 1 });
  assert.deepEqual(r.longest, { name: 'L1HS', cls: 'LINE', chrom: '3', start: 5001, length: 6100 });
});

const dup = (chrom, s, e, other, os, oe, uid, frac) => {
  const c = new Array(30).fill('0');
  Object.assign(c, { 1: chrom, 2: String(s), 3: String(e), 7: other, 8: String(os), 9: String(oe), 11: String(uid), 26: String(frac) });
  return c.join('\t');
};

test('superDupsSummary merges covered bases and bins identity once per pair', () => {
  const acc = superDupsSummary();
  [
    dup('chr1', 100, 1100, 'chr1', 5000, 6000, 1, 0.995),
    dup('chr1', 5000, 6000, 'chr1', 100, 1100, 1, 0.995),
    dup('chr1', 900, 2000, 'chr2', 100, 1200, 2, 0.93),
    dup('chr2', 100, 1200, 'chr1', 900, 2000, 2, 0.93),
  ].forEach((l) => acc.add(l));
  const r = acc.result();
  assert.equal(r.pairs, 2);
  assert.deepEqual(r.kinds, { intra: 1, inter: 1 });
  assert.deepEqual(r.covered, { 1: 2900, 2: 1100 });
  assert.equal(r.identity.edges.length, 21);
  assert.equal(r.identity.counts[19], 1);
  assert.equal(r.identity.counts[6], 1);
  assert.deepEqual(r.longest, { chrom: '2', start: 101, length: 1100, other: 'chr1:901-2,000', identity: 0.93 });
});

const trf = (chrom, s, e, size, seq) => {
  const c = new Array(17).fill('0');
  Object.assign(c, { 1: chrom, 2: String(s), 3: String(e), 7: String(size), 16: seq });
  return c.join('\t');
};

test('tandemSummary groups repeats by unit length and canonical motif', () => {
  const acc = tandemSummary();
  [
    trf('chr1', 10000, 10468, 6, 'TAACCC'),
    trf('chr2', 0, 50, 2, 'AC'),
    trf('chr2', 100, 140, 2, 'TG'),
    trf('chr3', 0, 30, 1, 'A'),
    trf('chr3', 40, 100, 25, 'ACGTACGTACGTACGTACGTACGTA'),
  ].forEach((l) => acc.add(l));
  const r = acc.result();
  assert.deepEqual(r.units.map((u) => u.unit), ['1', '2', '3', '4', '5', '6', '7-12', '13-50', '51+']);
  assert.deepEqual(r.units[1], { unit: '2', count: 2, bp: 90 });
  assert.deepEqual(r.units[7], { unit: '13-50', count: 1, bp: 60 });
  assert.deepEqual(r.motifs[0], { motif: 'AC', count: 2, bp: 90 });
  assert.ok(r.motifs.some((m) => m.motif === 'AACCCT'));
  assert.deepEqual(r.longest, { chrom: '1', start: 10001, length: 468, motif: 'TAACCC' });
});
