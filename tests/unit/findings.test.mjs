import { test } from 'node:test';
import assert from 'node:assert/strict';
import { findings, assemblyNumbers, sequenceNumbers } from '../../src/findings.js';
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
  sequence: {
    chromosomes: [
      {
        chrom: '1', length: 1000, A: 300, C: 200, G: 200, T: 300, N: 0, cpg: 10, masked: 500, g4: { plus: 3, minus: 4 },
      },
      {
        chrom: '19', length: 1000, A: 250, C: 250, G: 250, T: 250, N: 0, cpg: 20, masked: 600, g4: { plus: 5, minus: 5 },
      },
      {
        chrom: 'MT', length: 100, A: 25, C: 25, G: 25, T: 25, N: 0, cpg: 4, masked: 0, g4: { plus: 0, minus: 0 },
      },
    ],
    homopolymers: {
      A: { length: 72, chrom: '4', start: 1 }, C: { length: 29, chrom: '7', start: 1 }, G: { length: 81, chrom: '2', start: 1 }, T: { length: 90, chrom: '12', start: 66057593 },
    },
  },
};

fixture.repeats = {
  rmsk: {
    total: { count: 50, bp: 1000 },
    families: [{ cls: 'LINE', family: 'L1', count: 10, bp: 400 }, { cls: 'SINE', family: 'Alu', count: 30, bp: 300 }],
    alu: { AluY: { count: 5, bp: 50 }, AluS: { count: 20, bp: 200 }, AluJ: { count: 5, bp: 50 } },
    l1hs: { count: 8, fullLength: 3 },
  },
  superDups: { pairs: 4, kinds: { intra: 1, inter: 3 }, covered: { 1: 60, 19: 40 } },
  tandem: {
    units: [{ unit: '1', count: 7, bp: 70 }, { unit: '2', count: 5, bp: 50 }],
    motifs: [{ motif: 'A', count: 7, bp: 70 }, { motif: 'AC', count: 4, bp: 40 }],
  },
};

fixture.genes = {
  release: { gencode: '50', ensembl: '116', date: '2026-04-08' },
  genes: {
    total: 100, byType: { protein_coding: 40, lncRNA: 30, processed_pseudogene: 20, miRNA: 10 }, byChrom: { 1: { protein_coding: 30 }, 19: { protein_coding: 10 } },
  },
  transcripts: { total: 500, byType: { protein_coding: 200, nonsense_mediated_decay: 50 } },
  spans: { pc: { summary: { median: 29403 } }, longest: [{ name: 'RBFOX1', span: 2473539 }] },
  isoforms: { summary: { median: 15 }, top: [{ name: 'ZBTB20', n: 361 }] },
  mane: { transcripts: 38, singleExon: 3, exons: { summary: { median: 8 }, top: [{ name: 'TTN', exons: 363 }] } },
  territory: {
    cds: 20, utr: 40, ncExon: 60, intron: 1000,
  },
  immune: {
    IGH: {
      V: 61, D: 37, J: 6, C: 9, pseudo: 95,
    },
  },
};
fixture.catalogs = {
  hgnc: {
    approved: 45289, locusGroups: { 'protein-coding gene': 19299, pseudogene: 14671 }, spreadsheet: { count: 26 }, pseudoParents: [{ parent: 'RN7SK', n: 293 }, { parent: 'RPL21', n: 137 }],
  },
  geneInfo: { types: { 'protein-coding': 20600 } },
  mane: { select: 19363 },
};

const codonFixture = Array.from({ length: 64 }, (_, i) => ([48, 50, 56].includes(i) ? 0 : 10));
codonFixture[34] = 400;
codonFixture[25] = 2;
fixture.sequence.mane = {
  classes: {
    'GT-AG': 90, 'GC-AG': 8, 'AT-AC': 1, other: 1,
  },
  kozak: Array.from({ length: 10 }, (_, i) => (i === 3 ? [6, 1, 2, 1] : [1, 1, 1, 1])),
  stops: { TAA: 3, TAG: 2, TGA: 5 },
  coding: 10,
  codons: codonFixture,
};

test('findings state the splicing, translation and codon sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.splicing, '90% of MANE Select introns start GT and end AG; 8 start GC and 1 runs from AT to AC');
  assert.equal(f.translation, 'A purine sits three bases before the start codon in 80% of MANE coding sequences; 50% end in TGA');
  assert.equal(f.codons, 'GAG is the most used codon, 403.2 per thousand; CGC the least used sense codon, 2.0');
});

test('findings state the gene sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.genecount, 'Protein-coding gene counts run from 40 (GENCODE) to 20,600 (NCBI Gene)');
  assert.equal(f.density, 'Chromosome 1 carries 30,000 protein-coding genes per sequenced megabase; chromosome 19, 10,000');
  assert.equal(f.genesize, 'The median protein-coding gene spans 29,403 bp; RBFOX1, the longest, 2,473,539');
  assert.equal(f.exons, 'A MANE Select transcript has a median of 8 exons; TTN has 363; 3 have one');
  assert.equal(f.territory, 'Coding sequence covers 1% of GRCh38, exons of all genes 6% and introns 50%');
  assert.equal(f.isoforms, 'Protein-coding genes have a median of 15 GENCODE transcripts; ZBTB20 has 361');
  assert.equal(f.names, 'HGNC has approved 45,289 gene symbols; 26 were changed so spreadsheets would stop reading them as dates');
  assert.equal(f.pseudogenes, '14,671 named pseudogenes; RN7SK has 293 copies, RPL21 137');
  assert.equal(f.immune, '61 IGHV, 37 IGHD and 6 IGHJ segments make 13,542 heavy-chain V-D-J combinations before junctional diversity');
});

test('findings state the repeat sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.transposons, 'RepeatMasker marks 50% of GRCh38’s sequenced bases as repeats; L1 elements alone cover 20%');
  assert.equal(f.alu, '30 Alu copies, 5 in the youngest subfamily, AluY; 8 L1HS copies, 3 of them at least 6 kb long');
  assert.equal(f.segdups, 'Segmental duplications cover 100 bases, 5% of GRCh38; 3 of 4 pairs join different chromosomes');
  assert.equal(f.tandem, 'Tandem Repeats Finder calls 12 tandem repeats; the most frequent units are A, then AC');
});

test('sequenceNumbers totals the nuclear chromosomes and the CpG ratio', () => {
  const n = sequenceNumbers(fixture);
  assert.equal(n.acgt, 2000);
  assert.equal(n.gc, 0.45);
  assert.equal(n.cpgOe, 30 / ((450 * 450) / 2000));
  assert.equal(n.g4, 17);
  assert.deepEqual(n.richest, { chrom: '19', gc: 0.5 });
  assert.deepEqual(n.poorest, { chrom: '1', gc: 0.4 });
});

test('findings state the composition and motif sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.composition, 'GRCh38 is 45% G and C; chromosome 19 is the richest at 50%, chromosome 1 the poorest at 40%');
  assert.equal(f.cpg, 'CpG turns up at 29.6% of the rate the C and G content predicts');
  assert.equal(f.motifs, '17 G-quadruplex motifs; the longest run of one letter is 90 Ts on chromosome 12');
});

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
