import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  exact, cell, compact, pct, usd, decimal, wordy, ordinal, bases, sci, mlog, coords, density, signed, toCsv, perLitre,
} from '../../src/lib.js';

test('perLitre writes a mass concentration in the unit that keeps it between 1 and 1,000', () => {
  assert.equal(perLitre(4e13), '40 g/L');
  assert.equal(perLitre(2.5e10), '25 mg/L');
  assert.equal(perLitre(1.5e6), '1.5 µg/L');
  assert.equal(perLitre(123456), '123 ng/L');
  assert.equal(perLitre(5), '5 pg/L');
  assert.equal(perLitre(1), '1 pg/L');
  assert.equal(perLitre(1e12), '1 g/L');
  assert.equal(perLitre(0.4), '0.4 pg/L');
});

test('toCsv quotes fields that hold commas, quotes or line breaks', () => {
  assert.equal(toCsv(['chrom', 'bases'], [['1', 248956422], ['X', null]]), 'chrom,bases\r\n1,248956422\r\nX,\r\n');
  assert.equal(toCsv(['note'], [['a, b'], ['say "hi"'], ['two\nlines']]), 'note\r\n"a, b"\r\n"say ""hi"""\r\n"two\nlines"\r\n');
});

test('exact groups thousands with commas', () => {
  assert.equal(exact(3117275501), '3,117,275,501');
  assert.equal(exact(0), '0');
  assert.equal(exact(999), '999');
  assert.equal(exact('16,569'), '16,569');
});

test('exact returns a dash for missing values', () => {
  assert.equal(exact(undefined), '—');
  assert.equal(exact(null), '—');
  assert.equal(exact('n/a'), '—');
});

test('compact abbreviates with one decimal', () => {
  assert.equal(compact(41379710), '41.4M');
  assert.equal(compact(1500), '1.5K');
  assert.equal(compact(2_300_000_000), '2.3B');
  assert.equal(compact(4_100_000_000_000), '4.1T');
  assert.equal(compact(950), '950');
  assert.equal(compact(1_000_000), '1M');
});

test('compact rolls up to the next unit instead of printing 1000K', () => {
  assert.equal(compact(999_960), '1M');
  assert.equal(compact(999.6), '1K');
  assert.equal(compact(999_950_000), '1B');
  assert.equal(compact(999_400), '999.4K');
});

test('pct formats a ratio as a percentage', () => {
  assert.equal(pct(0.1234), '12.3%');
  assert.equal(pct(1), '100%');
  assert.equal(pct(0), '0%');
  assert.equal(pct(0.0004), '0.04%');
});

test('usd formats dollars compactly', () => {
  assert.equal(usd(95_263_072), '$95.3M');
  assert.equal(usd(525), '$525');
});

test('cell writes table numbers with thousands grouped and their decimals kept', () => {
  assert.equal(cell(1234567), '1,234,567');
  assert.equal(cell(0.279), '0.279');
  assert.equal(cell(0.15000000000000002), '0.15');
  assert.equal(cell(1234.5678), '1,234.5678');
  assert.equal(cell(-0.5), '-0.5');
  assert.equal(cell(99.7), '99.7');
  assert.equal(cell(0.0000012), '0.0000012');
  assert.equal(cell(null), '—');
});

test('decimal keeps a fixed number of places', () => {
  assert.equal(decimal(0.5106, 2), '0.51');
  assert.equal(decimal(2.866), '2.9');
  assert.equal(decimal(undefined), '—');
});

test('wordy spells out large numbers', () => {
  assert.equal(wordy(1_240_000_000), '1.24 billion');
  assert.equal(wordy(3_117_275_501), '3.12 billion');
  assert.equal(wordy(41_106_423), '41.1 million');
  assert.equal(wordy(9_999), '9,999');
});

test('ordinal adds the English suffix', () => {
  assert.deepEqual([1, 2, 3, 4, 11, 12, 13, 21, 22, 101, 111].map(ordinal), ['1st', '2nd', '3rd', '4th', '11th', '12th', '13th', '21st', '22nd', '101st', '111th']);
});

test('bases scales sequence lengths to bp, kb, Mb and Gb', () => {
  assert.equal(bases(57), '57 bp');
  assert.equal(bases(9_999), '9,999 bp');
  assert.equal(bases(16_569), '16.6 kb');
  assert.equal(bases(248_956_422), '249 Mb');
  assert.equal(bases(46_709_983), '46.7 Mb');
  assert.equal(bases(3_117_275_501), '3.12 Gb');
  assert.equal(bases(999_960), '1 Mb');
  assert.equal(bases(null), '—');
});

test('sci writes small p-values with a superscript power of ten', () => {
  assert.equal(sci(5e-8), '5×10⁻⁸');
  assert.equal(sci(3.24e-12), '3.2×10⁻¹²');
  assert.equal(sci(1e-300), '1×10⁻³⁰⁰');
  assert.equal(sci(0.0123), '0.012');
  assert.equal(sci(0.5), '0.5');
  assert.equal(sci(0), '0');
});

test('mlog gives minus log10 of a p-value', () => {
  assert.equal(mlog(1e-8), 8);
  assert.ok(Math.abs(mlog(5e-8) - 7.30103) < 1e-5);
  assert.equal(mlog(0), Infinity);
});

test('coords writes a 1-based genomic interval', () => {
  assert.equal(coords('17', 7661779, 7687538), 'chr17:7,661,779-7,687,538');
  assert.equal(coords('chrX', 31097677, 33339441), 'chrX:31,097,677-33,339,441');
  assert.equal(coords('MT', 1, 16569), 'chrM:1-16,569');
});

test('density gives a count per megabase with one decimal', () => {
  assert.equal(density(2058, 248_956_422), '8.3 per Mb');
  assert.equal(density(0, 1_000_000), '0 per Mb');
});

test('signed puts a sign in front of changes', () => {
  assert.equal(signed(1200), '+1,200');
  assert.equal(signed(-35), '−35');
  assert.equal(signed(0), '0');
  assert.equal(signed(2.5, (v) => v.toFixed(1)), '+2.5');
});
