import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readTwoBit } from '../../scripts/twobit.mjs';
import {
  encode, decodeText, composition, homopolymers, g4Motifs,
} from '../../scripts/sequence.mjs';

const PACK = { T: 0, C: 1, A: 2, G: 3, N: 0 };

const twoBitOf = (seqs) => {
  const names = Object.keys(seqs);
  const records = names.map((name) => {
    const dna = seqs[name];
    const upper = dna.toUpperCase();
    const runs = (test) => {
      const out = [];
      for (let i = 0; i < dna.length; i += 1) {
        if (test(i) && (i === 0 || !test(i - 1))) out.push([i, 0]);
        if (test(i)) out.at(-1)[1] += 1;
      }
      return out;
    };
    const nBlocks = runs((i) => upper[i] === 'N');
    const maskBlocks = runs((i) => dna[i] !== upper[i] && upper[i] !== 'N');
    const packed = Buffer.alloc(Math.ceil(dna.length / 4));
    for (let i = 0; i < dna.length; i += 1) packed[i >> 2] |= PACK[upper[i]] << (6 - 2 * (i & 3));
    const head = Buffer.alloc(4 * (2 + nBlocks.length * 2 + 1 + maskBlocks.length * 2 + 1));
    let o = 0;
    const put = (v) => { head.writeUInt32LE(v, o); o += 4; };
    put(dna.length);
    put(nBlocks.length);
    nBlocks.forEach(([s]) => put(s));
    nBlocks.forEach(([, n]) => put(n));
    put(maskBlocks.length);
    maskBlocks.forEach(([s]) => put(s));
    maskBlocks.forEach(([, n]) => put(n));
    put(0);
    return Buffer.concat([head, packed]);
  });
  const indexSize = names.reduce((n, name) => n + 1 + name.length + 4, 0);
  const header = Buffer.alloc(16);
  header.writeUInt32LE(0x1a412743, 0);
  header.writeUInt32LE(names.length, 8);
  const index = Buffer.alloc(indexSize);
  let at = 16 + indexSize;
  let o = 0;
  names.forEach((name, k) => {
    index.writeUInt8(name.length, o);
    index.write(name, o + 1, 'latin1');
    index.writeUInt32LE(at, o + 1 + name.length);
    o += 1 + name.length + 4;
    at += records[k].length;
  });
  return Buffer.concat([header, index, ...records]);
};

test('readTwoBit lists sequences and decodes bases, Ns and soft masking', () => {
  const file = readTwoBit(twoBitOf({ chrA: 'ACGTNNNNacgtGG', chrB: 'TTTTTTTTT' }));
  assert.deepEqual(file.names, ['chrA', 'chrB']);
  const a = file.sequence('chrA');
  assert.equal(a.size, 14);
  assert.equal(a.masked, 4);
  assert.deepEqual(a.nBlocks, [[4, 4]]);
  assert.equal(decodeText(a.codes()), 'ACGTNNNNACGTGG');
  assert.equal(decodeText(file.sequence('chrB').codes()), 'TTTTTTTTT');
});

test('encode and decodeText round-trip', () => {
  assert.equal(decodeText(encode('ACGTN')), 'ACGTN');
});

test('composition counts bases, CpGs and dinucleotides, skipping pairs that touch N', () => {
  const c = composition(encode('ACGCGTNNCGA'), 4);
  assert.deepEqual([c.A, c.C, c.G, c.T, c.N], [2, 3, 3, 1, 2]);
  assert.equal(c.cpg, 3);
  assert.equal(c.dinucs.reduce((a, b) => a + b, 0), 7);
  assert.equal(c.dinucs[0 * 4 + 1], 1);
  assert.deepEqual(c.bins, [0.75, 0.5, 2 / 3]);
});

test('composition leaves a bin empty when most of it is N', () => {
  assert.deepEqual(composition(encode('GGNNNNNA'), 4).bins, [1, null]);
});

test('homopolymers finds the longest run of each base outside N', () => {
  const h = homopolymers(encode('AAAACGGGTTTTTNNNNNNAAAAAAC'));
  assert.deepEqual(h.A, { length: 6, start: 19 });
  assert.deepEqual(h.T, { length: 5, start: 8 });
  assert.deepEqual(h.G, { length: 3, start: 5 });
  assert.deepEqual(h.C, { length: 1, start: 4 });
});

test('g4Motifs counts G3+N1-7 quadruplex motifs on each strand without overlap', () => {
  const motif = 'GGGAGGGTTGGGAGGG';
  assert.deepEqual(g4Motifs(encode(`AA${motif}AA`)), { plus: 1, minus: 0 });
  assert.deepEqual(g4Motifs(encode('CCCTCCCACCCTCCC')), { plus: 0, minus: 1 });
  assert.deepEqual(g4Motifs(encode('GGGAAAAAAAAGGGAGGGAGGG')), { plus: 0, minus: 0 });
  assert.deepEqual(g4Motifs(encode(`${motif}A${motif}`)), { plus: 2, minus: 0 });
  assert.deepEqual(g4Motifs(encode('GGGNGGGAGGGAGGG')), { plus: 0, minus: 0 });
});
