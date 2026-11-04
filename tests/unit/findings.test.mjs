import { test } from 'node:test';
import assert from 'node:assert/strict';
import { findings, assemblyNumbers } from '../../src/findings.js';
import snapshot from '../../src/snapshot.js';

const chrom = (c, length, ungapped) => ({ chrom: c, length, ungapped, gc: 41, spanned: 1, unspanned: 1, refseq: null, genbank: null });

const fixture = {
  generated: '2026-11-08T06:17:00.000Z',
  assembly: {
    grch38: {
      name: 'GRCh38.p14',
      altLoci: 35,
      roles: [
        { role: 'assembled-molecule', count: 4, length: 400016569 },
        { role: 'unlocalized-scaffold', count: 42, length: 6978808 },
        { role: 'unplaced-scaffold', count: 127, length: 4485509 },
        { role: 'alt-scaffold', count: 261, length: 109535387 },
        { role: 'fix-patch', count: 164, length: 65946498 },
        { role: 'novel-patch', count: 90, length: 23679459 },
      ],
      chromosomes: [chrom('1', 248956422, 230481121), chrom('21', 46709983, 40088623), chrom('Y', 57227415, 26415048), chrom('MT', 16569, 16569)],
    },
    t2t: {
      name: 'T2T-CHM13v2.0',
      chromosomes: [chrom('1', 248387328, 248387328), chrom('21', 45090682, 45090682), chrom('Y', 62460029, 62460029)],
    },
    centromeres: { 1: 2905865, 21: 2050048, Y: 227095 },
    gaps: { 1: { contig: { count: 2, bp: 100 } } },
  },
};

test('assemblyNumbers totals the nuclear chromosomes of both assemblies', () => {
  const n = assemblyNumbers(fixture);
  assert.equal(n.grch38Length, 352893820);
  assert.equal(n.grch38Ungapped, 296984792);
  assert.equal(n.grch38N, 55909028);
  assert.equal(n.t2tLength, 355938039);
  assert.equal(n.added, 58953247);
  assert.equal(n.patches, 254);
  assert.deepEqual(n.mostAdded, { chrom: 'Y', bases: 36044981 });
});

test('findings state each assembly section as a number', () => {
  const f = findings(fixture);
  assert.equal(f.chromosomes, 'Chromosome 1 spans 248,956,422 bp in GRCh38; chromosome 21, the shortest, 46,709,983');
  assert.equal(f.gaps, 'GRCh38 leaves 55,909,028 bases of its chromosomes as N, 15.8% of their length');
  assert.equal(f.scaffolds, 'Besides its chromosomes, GRCh38.p14 holds 261 alternate loci in 35 regions and 254 patches');
  assert.equal(f.t2t, 'T2T-CHM13 sequences 58,953,247 more bases than GRCh38 does; 36,044,981 of them are on chromosome Y');
});

test('the real snapshot gives every finding as plain text without filler words', () => {
  const f = findings(snapshot);
  assert.ok(Object.keys(f).length >= 4);
  for (const [key, text] of Object.entries(f)) {
    assert.equal(typeof text, 'string', key);
    assert.ok(text.length > 10 && text.length < 140, `${key}: ${text}`);
    assert.doesNotMatch(text, /NaN|undefined|null|Infinity/, `${key}: ${text}`);
    assert.doesNotMatch(text, /\b(only|just|mere|merely|soar\w*|plung\w*|stagger\w*|shock\w*|remarkabl\w*|astonish\w*)\b/i, `${key}: ${text}`);
  }
});
