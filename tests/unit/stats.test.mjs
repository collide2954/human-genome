import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  quantile, summary, logEdges, linearEdges, histogram, lorenz, gini, mergeIntervals, coverage, canonicalMotif, topBy, binIndex, pearson, linearFit, clusterOrder,
} from '../../src/stats.js';

test('clusterOrder puts the most similar items next to each other', () => {
  const m = [
    [1, 0.1, 0.9, 0.2],
    [0.1, 1, 0.2, 0.8],
    [0.9, 0.2, 1, 0.1],
    [0.2, 0.8, 0.1, 1],
  ];
  const order = clusterOrder(m);
  assert.deepEqual([...order].sort(), [0, 1, 2, 3]);
  const at = (i) => order.indexOf(i);
  assert.equal(Math.abs(at(0) - at(2)), 1);
  assert.equal(Math.abs(at(1) - at(3)), 1);
  assert.deepEqual(clusterOrder([[1]]), [0]);
  assert.deepEqual(clusterOrder(m), order);
});

test('pearson and linearFit describe a straight-line relation', () => {
  assert.equal(pearson([1, 2, 3, 4], [2, 4, 6, 8]), 1);
  assert.equal(pearson([1, 2, 3, 4], [8, 6, 4, 2]), -1);
  assert.ok(Math.abs(pearson([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]) - 0.8) < 1e-12);
  assert.deepEqual(linearFit([0, 1, 2], [1, 3, 5]), { slope: 2, intercept: 1 });
  assert.ok(Number.isNaN(pearson([1, 1, 1], [1, 2, 3])));
});

const close = (a, b, eps = 1e-9) => Math.abs(a - b) <= eps;

test('quantile interpolates between sorted values', () => {
  const s = [1, 2, 3, 4];
  assert.equal(quantile(s, 0), 1);
  assert.equal(quantile(s, 1), 4);
  assert.equal(quantile(s, 0.5), 2.5);
  assert.equal(quantile([7], 0.9), 7);
});

test('summary gives the five numbers, mean and count without changing its input', () => {
  const values = [9, 1, 5, 3, 7];
  assert.deepEqual(summary(values), { n: 5, min: 1, q1: 3, median: 5, q3: 7, max: 9, mean: 5 });
  assert.deepEqual(values, [9, 1, 5, 3, 7]);
  assert.equal(summary([]), null);
});

test('logEdges spaces bin edges evenly per decade', () => {
  const e = logEdges(0, 2, 2);
  assert.equal(e.length, 5);
  assert.equal(e[0], 1);
  assert.ok(close(e[1], 3.16228, 1e-5));
  assert.equal(e[2], 10);
  assert.equal(e[4], 100);
  assert.deepEqual(linearEdges(0, 1, 0.25), [0, 0.25, 0.5, 0.75, 1]);
});

test('histogram counts into half-open bins and clamps the tails', () => {
  const edges = [0, 10, 20, 30];
  assert.deepEqual(histogram([0, 5, 10, 19.9, 20, 30, 45, -3], edges), [3, 2, 3]);
  assert.equal(binIndex(edges, 25), 2);
  assert.equal(binIndex(edges, 30), 2);
  assert.equal(binIndex(edges, -1), 0);
});

test('lorenz runs from 0,0 to 1,1 along the sorted cumulative share', () => {
  const curve = lorenz([1, 1, 1, 7], 4);
  assert.deepEqual(curve[0], [0, 0]);
  assert.deepEqual(curve.at(-1), [1, 1]);
  assert.deepEqual(curve[2], [0.5, 0.2]);
  assert.equal(curve.length, 5);
});

test('gini is 0 for equal shares and near 1 when one holds everything', () => {
  assert.equal(gini([5, 5, 5, 5]), 0);
  assert.ok(close(gini([0, 0, 0, 10]), 0.75));
  assert.ok(close(gini([1, 2, 3, 4]), 0.25));
});

test('mergeIntervals joins overlapping and touching half-open intervals', () => {
  assert.deepEqual(mergeIntervals([[10, 20], [0, 5], [15, 30], [30, 35], [40, 41]]), [[0, 5], [10, 35], [40, 41]]);
  assert.equal(coverage([[10, 20], [0, 5], [15, 30], [30, 35], [40, 41]]), 31);
  assert.deepEqual(mergeIntervals([]), []);
});

test('canonicalMotif picks the smallest rotation over both strands', () => {
  assert.equal(canonicalMotif('CA'), 'AC');
  assert.equal(canonicalMotif('TG'), 'AC');
  assert.equal(canonicalMotif('GT'), 'AC');
  assert.equal(canonicalMotif('CAG'), 'AGC');
  assert.equal(canonicalMotif('CTG'), 'AGC');
  assert.equal(canonicalMotif('T'), 'A');
  assert.equal(canonicalMotif('GGGGCC'), 'CCCCGG');
});

test('topBy keeps the n largest rows by a key, ties in first-seen order', () => {
  const rows = [{ k: 'a', v: 2 }, { k: 'b', v: 9 }, { k: 'c', v: 9 }, { k: 'd', v: 1 }];
  assert.deepEqual(topBy(rows, 'v', 2).map((r) => r.k), ['b', 'c']);
  assert.deepEqual(topBy(rows, 'v', 10).map((r) => r.k), ['b', 'c', 'a', 'd']);
});
