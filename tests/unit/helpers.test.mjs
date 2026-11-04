import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  matches, filterTimeline, sortRows, niceTicks, yearOf, longDate, labelStep, LENSES, readLens, inLens, lensSearch,
} from '../../src/lib.js';

test('matches needs every term, ignoring case and accents', () => {
  assert.equal(matches('Francis Collins, NHGRI director', 'collins nhgri'), true);
  assert.equal(matches('Ångström resolution', 'angstrom'), true);
  assert.equal(matches('BRCA1 and BRCA2', 'brca1 tp53'), false);
  assert.equal(matches('anything', ''), true);
  assert.equal(matches('anything', '   '), true);
});

const ITEMS = [
  { date: '1977-12', org: 'Methods', tag: 'sequencing', text: 'Sanger publishes chain-termination sequencing' },
  { date: '1990-10-01', org: 'Projects', tag: 'projects', text: 'Human Genome Project begins' },
  { date: '2022-03-31', org: 'Assemblies', tag: 'assemblies', text: 'T2T-CHM13 published' },
];

test('filterTimeline filters by tag, lane and query', () => {
  assert.equal(filterTimeline(ITEMS, {}).length, 3);
  assert.deepEqual(filterTimeline(ITEMS, { tags: ['projects'] }).map((i) => i.date), ['1990-10-01']);
  assert.deepEqual(filterTimeline(ITEMS, { orgs: ['Methods', 'Assemblies'] }).map((i) => i.date), ['1977-12', '2022-03-31']);
  assert.deepEqual(filterTimeline(ITEMS, { query: 'sanger' }).map((i) => i.date), ['1977-12']);
  assert.deepEqual(filterTimeline(ITEMS, { query: '2022' }).map((i) => i.date), ['2022-03-31']);
});

test('sortRows sorts numbers and strings with nulls last', () => {
  const rows = [{ n: 2, s: 'b' }, { n: null, s: 'c' }, { n: 10, s: 'a' }];
  assert.deepEqual(sortRows(rows, 'n', 'desc').map((r) => r.n), [10, 2, null]);
  assert.deepEqual(sortRows(rows, 'n', 'asc').map((r) => r.n), [2, 10, null]);
  assert.deepEqual(sortRows(rows, 's', 'asc').map((r) => r.s), ['a', 'b', 'c']);
  assert.deepEqual(rows.map((r) => r.n), [2, null, 10]);
});

test('niceTicks returns round steps up to max', () => {
  assert.deepEqual(niceTicks(1_882_376, 4), [0, 500000, 1000000, 1500000, 2000000]);
  assert.deepEqual(niceTicks(9, 3), [0, 5, 10]);
  assert.deepEqual(niceTicks(0.62, 3), [0, 0.5, 1]);
  assert.deepEqual(niceTicks(0.29, 3), [0, 0.1, 0.2, 0.3]);
  assert.deepEqual(niceTicks(0, 3), [0]);
  assert.deepEqual(niceTicks(100.00000000000003, 4), [0, 50, 100]);
});

test('yearOf and longDate read partial dates', () => {
  assert.equal(yearOf('2001-02-15'), 2001);
  assert.equal(yearOf('1953'), 1953);
  assert.equal(longDate('2001-02-15'), '15 February 2001');
  assert.equal(longDate('1977-12'), 'December 1977');
  assert.equal(longDate('1953'), '1953');
});

test('labelStep spaces axis labels to fit the width', () => {
  assert.equal(labelStep(82, 770), 5);
  assert.equal(labelStep(82, 290), 10);
  assert.equal(labelStep(42, 290), 5);
  assert.equal(labelStep(5, 600), 1);
  assert.equal(labelStep(300, 200), 50);
});

test('lenses read from the URL and filter sections', () => {
  assert.deepEqual(LENSES.map((l) => l.id), ['all', 'sequence', 'genes', 'variation', 'clinic']);
  assert.equal(readLens('?lens=clinic'), 'clinic');
  assert.equal(readLens('?lens=bogus'), 'all');
  assert.equal(readLens(''), 'all');
  assert.equal(inLens('sequence genes', 'genes'), true);
  assert.equal(inLens('sequence', 'clinic'), false);
  assert.equal(inLens('', 'clinic'), true);
  assert.equal(inLens('clinic', 'all'), true);
  assert.equal(lensSearch('?lens=genes&x=1', 'all'), '?x=1');
  assert.equal(lensSearch('', 'variation'), '?lens=variation');
});
