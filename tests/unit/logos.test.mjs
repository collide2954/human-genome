import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  information, GENETIC_CODE, codonOf, heat, AMINO_ACIDS,
} from '../../src/lib.js';

const close = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps;

test('information gives 2 bits for a fixed base and 0 for an even mix', () => {
  const [fixed, even] = information([[0, 0, 10, 0], [5, 5, 5, 5]]);
  assert.equal(fixed.bits, 2);
  assert.deepEqual(fixed.letters.map((l) => l.base), ['A', 'C', 'T', 'G']);
  assert.equal(fixed.letters.at(-1).height, 2);
  assert.equal(even.bits, 0);
  assert.ok(even.letters.every((l) => l.height === 0));
});

test('information stacks letters smallest first with heights summing to the bits', () => {
  const [p] = information([[6, 2, 2, 0]]);
  const h = -(0.6 * Math.log2(0.6) + 2 * 0.2 * Math.log2(0.2));
  assert.ok(close(p.bits, 2 - h));
  assert.deepEqual(p.letters.map((l) => l.base), ['T', 'C', 'G', 'A']);
  assert.ok(close(p.letters.reduce((n, l) => n + l.height, 0), p.bits));
  assert.equal(p.n, 10);
});

test('the genetic code translates all 64 codons', () => {
  assert.equal(Object.keys(GENETIC_CODE).length, 64);
  assert.equal(GENETIC_CODE.ATG, 'M');
  assert.equal(GENETIC_CODE.TGG, 'W');
  assert.equal(GENETIC_CODE.TGA, '*');
  assert.equal(GENETIC_CODE.TAA, '*');
  assert.equal(GENETIC_CODE.TAG, '*');
  assert.equal(GENETIC_CODE.AGA, 'R');
  assert.equal(GENETIC_CODE.CTG, 'L');
  assert.equal(Object.values(GENETIC_CODE).filter((a) => a === 'L').length, 6);
  assert.equal(Object.values(GENETIC_CODE).filter((a) => a === 'S').length, 6);
  assert.equal(new Set(Object.values(GENETIC_CODE)).size, 21);
  assert.equal(AMINO_ACIDS.M, 'Methionine');
  assert.equal(Object.keys(AMINO_ACIDS).length, 20);
});

test('codonOf reads the index used by the snapshot, A C G T order', () => {
  assert.equal(codonOf(0), 'AAA');
  assert.equal(codonOf(14), 'ATG');
  assert.equal(codonOf(56), 'TGA');
  assert.equal(codonOf(63), 'TTT');
});

test('heat maps a value to an opacity between a floor and 1', () => {
  assert.equal(heat(0, 0, 10), 0.06);
  assert.equal(heat(10, 0, 10), 1);
  assert.ok(close(heat(5, 0, 10), 0.53));
  assert.equal(heat(20, 0, 10), 1);
  assert.equal(heat(3, 3, 3), 1);
});
