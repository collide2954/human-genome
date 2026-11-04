import { test } from 'node:test';
import assert from 'node:assert/strict';
import { readTwoBit } from '../../scripts/twobit.mjs';
import {
  encode, decodeText, composition, homopolymers, g4Motifs, sequenceSection,
} from '../../scripts/sequence.mjs';
import { latestFile } from '../../scripts/parsers.mjs';

const PACK = { T: 0, C: 1, A: 2, G: 3, N: 0 };

const twoBitOf = (seqs) => {
  const names = Object.keys(seqs);
  const records = names.map((name) => {
    const dna = seqs[name];
    const upper = dna.toUpperCase();
    const runs = (hit) => {
      const out = [];
      for (let i = 0; i < dna.length; i += 1) {
        if (hit(i) && (i === 0 || !hit(i - 1))) out.push([i, 0]);
        if (hit(i)) out.at(-1)[1] += 1;
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

test('sequenceSection summarizes each chromosome, the genome and the MANE transcripts', () => {
  const chr1 = `GCCACCATGGCCGTAAGT${'T'.repeat(22)}AGGATTGAACGTACGTA`;
  const file = readTwoBit(twoBitOf({ chr1, chrM: 'GATCACAGGT', chr1_alt: 'AAAA' }));
  const tx = {
    id: 'T1', gene: 'G1', chrom: '1', strand: '+', tag: 'MANE_Select', exons: [[1, 12], [43, 57]], cds: [[7, 12], [43, 45]], start: [[7, 9]], stop: [[46, 48]],
  };
  const out = sequenceSection(file, [tx], { binSize: 20, chroms: ['1', 'MT'] });
  assert.deepEqual(out.chromosomes.map((c) => c.chrom), ['1', 'MT']);
  const c1 = out.chromosomes[0];
  assert.equal(c1.length, 57);
  assert.equal(c1.A + c1.C + c1.G + c1.T + c1.N, 57);
  assert.equal(out.bins['1'].length, 3);
  assert.equal(out.dinucs.reduce((a, b) => a + b, 0), 56);
  assert.equal(out.homopolymers.T.chrom, '1');
  assert.equal(out.homopolymers.T.length, 23);
  assert.equal(out.homopolymers.T.start, 18);
  assert.equal(out.mane.transcripts, 1);
  assert.equal(out.mane.introns, 1);
  assert.deepEqual(out.mane.classes, { 'GT-AG': 1, 'GC-AG': 0, 'AT-AC': 0, other: 0 });
  assert.deepEqual(out.mane.donor[3], [0, 0, 1, 0]);
  assert.equal(out.mane.kozak.length, 10);
  assert.deepEqual(out.mane.stops, { TAA: 0, TAG: 0, TGA: 1 });
  assert.equal(out.mane.coding, 1);
});

test('latestFile picks the newest versioned file name from a directory listing', () => {
  const listing = '<a href="gencode.v49.annotation.gtf.gz">x</a> <a href="gencode.v50.annotation.gtf.gz">y</a> <a href="gencode.v50.basic.annotation.gtf.gz">z</a>';
  assert.deepEqual(latestFile(listing, /gencode\.v(\d+)\.annotation\.gtf\.gz/), { file: 'gencode.v50.annotation.gtf.gz', version: '50' });
  assert.deepEqual(latestFile('<a href="MANE.GRCh38.v1.5.summary.txt.gz">', /MANE\.GRCh38\.v([\d.]+)\.summary\.txt\.gz/), { file: 'MANE.GRCh38.v1.5.summary.txt.gz', version: '1.5' });
  assert.equal(latestFile('nothing here', /x(\d+)/), null);
});

test('g4Motifs counts G3+N1-7 quadruplex motifs on each strand without overlap', () => {
  const motif = 'GGGAGGGTTGGGAGGG';
  assert.deepEqual(g4Motifs(encode(`AA${motif}AA`)), { plus: 1, minus: 0 });
  assert.deepEqual(g4Motifs(encode('CCCTCCCACCCTCCC')), { plus: 0, minus: 1 });
  assert.deepEqual(g4Motifs(encode('GGGAAAAAAAAGGGAGGGAGGG')), { plus: 0, minus: 0 });
  assert.deepEqual(g4Motifs(encode(`${motif}A${motif}`)), { plus: 2, minus: 0 });
  assert.deepEqual(g4Motifs(encode('GGGNGGGAGGGAGGG')), { plus: 0, minus: 0 });
});
