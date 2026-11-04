import { test } from 'node:test';
import assert from 'node:assert/strict';
import { runSections, diffCounts } from '../../scripts/sections.mjs';

const NOW = new Date('2026-11-08T06:17:00Z');
const quiet = () => {};

test('every section builds and is stamped with the run time', async () => {
  const out = await runSections({ a: async () => ({ n: 1 }), b: async () => [2] }, { now: NOW, log: quiet });
  assert.equal(out.generated, '2026-11-08T06:17:00.000Z');
  assert.deepEqual(out.a, { n: 1 });
  assert.deepEqual(out.b, [2]);
  assert.deepEqual(out.asOf, { a: '2026-11-08T06:17:00.000Z', b: '2026-11-08T06:17:00.000Z' });
  assert.deepEqual(out.failed, []);
});

test('a failing section keeps last week\'s values and date and is listed', async () => {
  const previous = { generated: '2026-11-01T06:17:00.000Z', asOf: { a: '2026-11-01T06:17:00.000Z' }, a: { n: 0 } };
  const out = await runSections({
    a: async () => { throw new Error('HTTP 503'); },
    b: async () => 'ok',
  }, { now: NOW, previous, log: quiet });
  assert.deepEqual(out.a, { n: 0 });
  assert.equal(out.asOf.a, '2026-11-01T06:17:00.000Z');
  assert.equal(out.b, 'ok');
  assert.deepEqual(out.failed, [{ section: 'a', error: 'HTTP 503' }]);
});

test('a failing section with nothing to fall back on throws', async () => {
  await assert.rejects(runSections({ a: async () => { throw new Error('gone'); } }, { now: NOW, log: quiet }), /a: gone/);
});

test('named sections rebuild alone and the rest carry over', async () => {
  const previous = { generated: '2026-11-01T06:17:00.000Z', asOf: { a: '2026-11-01T06:17:00.000Z', b: '2026-10-25T06:17:00.000Z' }, a: 1, b: 2 };
  let builtA = false;
  const out = await runSections({ a: async () => { builtA = true; return 10; }, b: async () => 20 }, { now: NOW, previous, only: ['b'], log: quiet });
  assert.equal(builtA, false);
  assert.equal(out.a, 1);
  assert.equal(out.b, 20);
  assert.deepEqual(out.asOf, { a: '2026-11-01T06:17:00.000Z', b: '2026-11-08T06:17:00.000Z' });
});

test('diffCounts lists the counts that moved, largest relative change first', () => {
  const before = { clinvar: 1000, gwas: 500, dbsnp: 1197210835, papers: 10 };
  const after = { clinvar: 1100, gwas: 500, dbsnp: 1197210836, papers: 12, pgs: 5 };
  assert.deepEqual(diffCounts(before, after), [
    { key: 'papers', before: 10, after: 12, delta: 2 },
    { key: 'clinvar', before: 1000, after: 1100, delta: 100 },
    { key: 'dbsnp', before: 1197210835, after: 1197210836, delta: 1 },
  ]);
  assert.deepEqual(diffCounts(undefined, after), []);
});
