import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  encode, decodeText, gtfRecord, parseManeGtf, introns, spliceSites, intronClass, kozakWindow, cdsSequence, pwm, codonCounts,
} from '../../scripts/sequence.mjs';
import { reverseComplement } from '../../src/stats.js';

const EXON1 = 'GCCACCATGGCC';
const INTRON = `GTAAGT${'T'.repeat(22)}AG`;
const EXON2 = 'GATTGAACGTACGTA';
const PLUS = EXON1 + INTRON + EXON2;
const MINUS = reverseComplement(PLUS);
const flip = ([s, e]) => [PLUS.length + 1 - e, PLUS.length + 1 - s];

const TX_PLUS = {
  id: 'T1', gene: 'G1', chrom: '1', strand: '+', exons: [[1, 12], [43, 57]], cds: [[7, 12], [43, 45]], start: [[7, 9]], stop: [[46, 48]],
};
const TX_MINUS = {
  id: 'T2', gene: 'G2', chrom: '2', strand: '-', exons: [[46, 57], [1, 15]].map((x) => x), cds: [[7, 12], [43, 45]].map(flip), start: [[7, 9]].map(flip), stop: [[46, 48]].map(flip),
};

test('gtfRecord reads coordinates, strand and attributes with repeated tags', () => {
  const r = gtfRecord('chr1\thavana\tCDS\t14598998\t14599223\t.\t+\t0\tgene_id "ENSG00000189337.17"; transcript_id "ENST00000376030.7"; gene_name "KAZN"; tag "MANE_Select"; tag "basic"; exon_number 1;');
  assert.equal(r.chrom, '1');
  assert.equal(r.feature, 'CDS');
  assert.equal(r.start, 14598998);
  assert.equal(r.end, 14599223);
  assert.equal(r.strand, '+');
  assert.equal(r.attrs.transcript_id, 'ENST00000376030.7');
  assert.deepEqual(r.attrs.tag, ['MANE_Select', 'basic']);
  assert.equal(r.attrs.exon_number, '1');
});

test('parseManeGtf groups exons, CDS, start and stop codons by transcript in transcript order', () => {
  const attrs = (id, gene, tag) => `gene_id "X"; transcript_id "${id}"; gene_name "${gene}"; tag "${tag}";`;
  const lines = [
    '##comment',
    `chr2\tensembl\ttranscript\t1\t57\t.\t-\t.\t${attrs('T2', 'G2', 'MANE_Select')}`,
    `chr2\tensembl\texon\t1\t15\t.\t-\t.\t${attrs('T2', 'G2', 'MANE_Select')}`,
    `chr2\tensembl\texon\t46\t57\t.\t-\t.\t${attrs('T2', 'G2', 'MANE_Select')}`,
    `chr2\tensembl\tCDS\t46\t51\t.\t-\t0\t${attrs('T2', 'G2', 'MANE_Select')}`,
    `chr2\tensembl\tCDS\t13\t15\t.\t-\t0\t${attrs('T2', 'G2', 'MANE_Select')}`,
    `chr2\tensembl\tstart_codon\t49\t51\t.\t-\t0\t${attrs('T2', 'G2', 'MANE_Select')}`,
    `chr2\tensembl\tstop_codon\t10\t12\t.\t-\t0\t${attrs('T2', 'G2', 'MANE_Select')}`,
    `chr1\tensembl\texon\t1\t12\t.\t+\t.\t${attrs('T1', 'G1', 'MANE_Plus_Clinical')}`,
  ];
  const txs = parseManeGtf(lines);
  assert.equal(txs.length, 2);
  const t2 = txs.find((t) => t.id === 'T2');
  assert.deepEqual(t2.exons, [[46, 57], [1, 15]]);
  assert.deepEqual(t2.cds, [[46, 51], [13, 15]]);
  assert.deepEqual(t2.start, [[49, 51]]);
  assert.deepEqual(t2.stop, [[10, 12]]);
  assert.equal(t2.strand, '-');
  assert.equal(t2.chrom, '2');
  assert.equal(t2.gene, 'G2');
  assert.equal(t2.tag, 'MANE_Select');
  assert.equal(txs.find((t) => t.id === 'T1').tag, 'MANE_Plus_Clinical');
});

test('introns lie between exons in transcript order on either strand', () => {
  assert.deepEqual(introns(TX_PLUS), [[13, 42]]);
  assert.deepEqual(introns(TX_MINUS), [[16, 45]]);
});

for (const [name, seq, tx] of [['plus', PLUS, TX_PLUS], ['minus', MINUS, TX_MINUS]]) {
  test(`splice sites, Kozak context and CDS read the same on the ${name} strand`, () => {
    const codes = encode(seq);
    const { donors, acceptors } = spliceSites(codes, tx);
    assert.deepEqual(donors.map(decodeText), ['GCCGTAAGT']);
    assert.deepEqual(acceptors.map(decodeText), [`${'T'.repeat(18)}AGGAT`]);
    assert.equal(intronClass(donors[0], acceptors[0]), 'GT-AG');
    assert.equal(decodeText(kozakWindow(codes, tx)), 'GCCACCATGG');
    assert.equal(decodeText(cdsSequence(codes, tx)), 'ATGGCCGATTGA');
  });
}

test('intronClass names the terminal dinucleotides', () => {
  const donor = (d) => encode(`CAG${d}AGT`);
  const acceptor = (a) => encode(`${'C'.repeat(18)}${a}GTT`);
  assert.equal(intronClass(donor('GC'), acceptor('AG')), 'GC-AG');
  assert.equal(intronClass(donor('AT'), acceptor('AC')), 'AT-AC');
  assert.equal(intronClass(donor('CT'), acceptor('AC')), 'other');
});

test('pwm counts each base at each position and skips N', () => {
  assert.deepEqual(pwm([encode('AC'), encode('AG'), encode('NT')], 2), [[2, 0, 0, 0], [0, 1, 1, 1]]);
});

test('codonCounts tallies sense codons, stop codons and non-ATG starts', () => {
  const out = codonCounts([encode('ATGGCCGATTGA'), encode('ATGGCCTAA'), encode('CTGAAATAG'), encode('ATGAA')]);
  const at = (codon) => out.codons[[...codon].reduce((n, ch) => n * 4 + 'ACGT'.indexOf(ch), 0)];
  assert.equal(at('ATG'), 2);
  assert.equal(at('GCC'), 2);
  assert.equal(at('GAT'), 1);
  assert.equal(at('CTG'), 1);
  assert.equal(at('AAA'), 1);
  assert.deepEqual(out.stops, { TAA: 1, TAG: 1, TGA: 1 });
  assert.equal(out.nonAtg, 1);
  assert.equal(out.broken, 1);
  assert.equal(out.coding, 3);
});
