import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  waffle, scaleLinear, scaleLog, logTicks, stack, zoomDomain, panDomain, easeOut, spread, squarify,
} from '../../src/lib.js';

const close = (a, b, eps = 1e-6) => Math.abs(a - b) <= eps;

test('waffle splits 100 cells by largest remainder', () => {
  const cells = waffle([{ key: 'yes', value: 24.1 }, { key: 'no', value: 75.9 }]);
  assert.equal(cells.length, 100);
  assert.equal(cells.filter((c) => c === 'yes').length, 24);
  assert.equal(cells[0], 'yes');
  const thirds = waffle([{ key: 'a', value: 33.4 }, { key: 'b', value: 33.3 }, { key: 'c', value: 33.3 }]);
  assert.deepEqual(['a', 'b', 'c'].map((k) => thirds.filter((c) => c === k).length), [34, 33, 33]);
});

test('waffle always fills every cell when several parts round down', () => {
  const cells = waffle([1, 1, 1, 1, 1, 1].map((value, key) => ({ key, value })));
  assert.equal(cells.length, 100);
  assert.deepEqual([0, 1, 2, 3, 4, 5].map((k) => cells.filter((c) => c === k).length), [17, 17, 17, 17, 16, 16]);
});

test('scales map linearly and logarithmically', () => {
  assert.equal(scaleLinear([0, 10], [0, 100])(5), 50);
  assert.equal(scaleLinear([0, 10], [100, 0])(10), 0);
  assert.ok(close(scaleLog([1, 1000], [0, 300])(10), 100));
});

test('logTicks lists powers of ten inside the range', () => {
  assert.deepEqual(logTicks(56, 41381383), [100, 1000, 10000, 100000, 1000000, 10000000]);
});

test('stack accumulates series per row', () => {
  const rows = [{ year: 1, a: 1, b: 2 }, { year: 2, a: 3, b: 1 }];
  assert.deepEqual(stack(rows, ['a', 'b']), [
    { key: 'a', points: [[0, 1], [0, 3]] },
    { key: 'b', points: [[1, 3], [3, 4]] },
  ]);
});

test('zoomDomain zooms around a point and stays inside limits', () => {
  assert.deepEqual(zoomDomain([1800, 2000], 1900, 0.5, [1800, 2026]), [1850, 1950]);
  assert.deepEqual(zoomDomain([1950, 2026], 2026, 4, [1800, 2026]), [1800, 2026]);
  assert.deepEqual(zoomDomain([1900, 2000], 2000, 2, [1800, 2026]), [1800, 2000]);
  assert.deepEqual(zoomDomain([1900, 1902], 1901, 0.1, [1800, 2026], 4), [1899, 1903]);
});

test('panDomain shifts and clamps', () => {
  assert.deepEqual(panDomain([1900, 1950], 100, [1800, 2026]), [1976, 2026]);
  assert.deepEqual(panDomain([1900, 1950], -200, [1800, 2026]), [1800, 1850]);
});

test('easeOut eases toward one', () => {
  assert.equal(easeOut(0), 0);
  assert.equal(easeOut(1), 1);
  assert.equal(easeOut(0.5), 0.875);
});

test('spread keeps labels a minimum gap apart and inside the range', () => {
  const out = spread([100, 102, 104, 200], 14, 0, 300);
  for (let i = 1; i < out.length; i++) assert.ok(out[i] - out[i - 1] >= 14 - 1e-9, JSON.stringify(out));
  assert.ok(Math.abs(out.slice(0, 3).reduce((a, b) => a + b) / 3 - 102) < 1e-9);
  assert.equal(out[3], 200);
});

test('spread preserves input order and leaves separated labels alone', () => {
  assert.deepEqual(spread([50, 10, 90], 12, 0, 100), [50, 10, 90]);
  const out = spread([295, 298, 299], 12, 0, 300);
  assert.ok(Math.max(...out) <= 300 && Math.min(...out) >= 0);
  assert.ok(out[0] < out[1] && out[1] < out[2]);
});

test('squarify tiles the whole box with areas in proportion to values', () => {
  const tiles = squarify([{ key: 'a', value: 6 }, { key: 'b', value: 3 }, { key: 'c', value: 1 }], 100, 50);
  assert.equal(tiles.length, 3);
  const area = (t) => t.w * t.h;
  assert.ok(close(tiles.reduce((n, t) => n + area(t), 0), 5000));
  assert.ok(close(area(tiles.find((t) => t.key === 'a')), 3000));
  assert.ok(close(area(tiles.find((t) => t.key === 'c')), 500));
  for (const t of tiles) assert.ok(t.x >= -1e-9 && t.y >= -1e-9 && t.x + t.w <= 100 + 1e-9 && t.y + t.h <= 50 + 1e-9);
});
