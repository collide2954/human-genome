import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gencodeSummary } from '../../scripts/genes.mjs';

const attrs = (o) => Object.entries(o).flatMap(([k, v]) => (Array.isArray(v) ? v.map((x) => `${k} "${x}";`) : [`${k} "${v}";`])).join(' ');
const row = (chrom, feature, start, end, strand, a) => [`chr${chrom}`, 'HAVANA', feature, start, end, '.', strand, '.', attrs(a)].join('\t');

const gene = (chrom, start, end, strand, id, name, type) => row(chrom, 'gene', start, end, strand, { gene_id: id, gene_type: type, gene_name: name });
const tx = (chrom, feature, start, end, strand, g, t, type, ttype, tags = []) => row(chrom, feature, start, end, strand, {
  gene_id: g, transcript_id: t, gene_type: type, gene_name: g, transcript_type: ttype, tag: tags,
});

const LINES = [
  '##description: evidence-based annotation of the human genome (GRCh38), version 50 (Ensembl 116)',
  '##date: 2026-04-08',
  gene('1', 1001, 5000, '+', 'A', 'GENEA', 'protein_coding'),
  tx('1', 'transcript', 1001, 5000, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['basic', 'MANE_Select']),
  tx('1', 'exon', 1001, 1100, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['basic', 'MANE_Select']),
  tx('1', 'exon', 2001, 2200, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['basic', 'MANE_Select']),
  tx('1', 'exon', 4001, 5000, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['basic', 'MANE_Select']),
  tx('1', 'CDS', 1051, 1100, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  tx('1', 'CDS', 2001, 2200, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  tx('1', 'CDS', 4001, 4100, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  tx('1', 'UTR', 1001, 1050, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  tx('1', 'UTR', 4101, 5000, '+', 'A', 'A1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  tx('1', 'transcript', 1001, 4500, '+', 'A', 'A2', 'protein_coding', 'nonsense_mediated_decay'),
  tx('1', 'exon', 1001, 1100, '+', 'A', 'A2', 'protein_coding', 'nonsense_mediated_decay'),
  tx('1', 'exon', 4001, 4500, '+', 'A', 'A2', 'protein_coding', 'nonsense_mediated_decay'),
  gene('1', 10001, 10500, '-', 'B', 'LNCB', 'lncRNA'),
  tx('1', 'transcript', 10001, 10500, '-', 'B', 'B1', 'lncRNA', 'lncRNA'),
  tx('1', 'exon', 10001, 10500, '-', 'B', 'B1', 'lncRNA', 'lncRNA'),
  gene('2', 1, 300, '+', 'C', 'GENEC', 'protein_coding'),
  tx('2', 'transcript', 1, 300, '+', 'C', 'C1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  tx('2', 'exon', 1, 300, '+', 'C', 'C1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  tx('2', 'CDS', 1, 300, '+', 'C', 'C1', 'protein_coding', 'protein_coding', ['MANE_Select']),
  ...[['IGHV1-2', 'IG_V_gene', 100, 400], ['IGHD1-1', 'IG_D_gene', 500, 520], ['IGHJ1', 'IG_J_gene', 600, 650], ['IGHM', 'IG_C_gene', 700, 900], ['IGHV1-12', 'IG_V_pseudogene', 1000, 1290]]
    .flatMap(([name, type, s, e]) => [
      row('14', 'gene', s, e, '-', { gene_id: name, gene_type: type, gene_name: name }),
      row('14', 'exon', s, e, '-', { gene_id: name, transcript_id: `${name}-T`, gene_type: type, gene_name: name, transcript_type: type }),
    ]),
  row('7', 'gene', 100, 300, '+', { gene_id: 'TRBV1', gene_type: 'TR_V_gene', gene_name: 'TRBV1' }),
  row('7', 'exon', 100, 300, '+', { gene_id: 'TRBV1', transcript_id: 'TRBV1-T', gene_type: 'TR_V_gene', gene_name: 'TRBV1', transcript_type: 'TR_V_gene' }),
];

const run = () => {
  const acc = gencodeSummary({ binSize: 10000, lengths: { 1: 20000, 2: 1000, 7: 1000, 14: 2000 } });
  LINES.forEach((l) => acc.add(l));
  return acc.result();
};

test('gencodeSummary reads the release and counts genes and transcripts by type', () => {
  const r = run();
  assert.deepEqual(r.release, { gencode: '50', ensembl: '116', date: '2026-04-08' });
  assert.equal(r.genes.total, 9);
  assert.equal(r.genes.byType.protein_coding, 2);
  assert.equal(r.genes.byType.IG_V_pseudogene, 1);
  assert.deepEqual(r.genes.byChrom['1'], { protein_coding: 1, lncRNA: 1 });
  assert.equal(r.transcripts.byType.nonsense_mediated_decay, 1);
  assert.equal(r.transcripts.total, 10);
  assert.deepEqual(r.density.bins['1'], [1, 0]);
  assert.deepEqual(r.density.bins['2'], [1]);
});

test('gencodeSummary measures gene spans and isoforms', () => {
  const r = run();
  assert.equal(r.spans.pc.summary.median, 2150);
  assert.deepEqual(r.spans.longest.map((g) => [g.name, g.span]), [['GENEA', 4000], ['GENEC', 300]]);
  assert.equal(r.spans.lnc.summary.max, 500);
  assert.equal(r.spans.pc.counts.reduce((a, b) => a + b, 0), 2);
  assert.deepEqual(r.isoforms.top.map((g) => [g.name, g.n]), [['GENEA', 2], ['GENEC', 1]]);
});

test('gencodeSummary measures MANE Select structure and UTRs on the right side of the CDS', () => {
  const { mane } = run();
  assert.equal(mane.transcripts, 2);
  assert.equal(mane.singleExon, 1);
  assert.deepEqual(mane.exons.top.map((g) => [g.name, g.exons]), [['GENEA', 3], ['GENEC', 1]]);
  assert.equal(mane.introns.summary.max, 1800);
  assert.deepEqual(mane.introns.longest, { name: 'GENEA', chrom: '1', start: 2201, length: 1800 });
  assert.equal(mane.exonLengths.summary.n, 4);
  assert.equal(mane.cds.summary.max, 350);
  assert.equal(mane.mrna.summary.max, 1300);
  assert.equal(mane.utr5.summary.max, 50);
  assert.equal(mane.utr3.summary.max, 900);
});

test('gencodeSummary splits annotated bases into CDS, UTR, other exons and introns once each', () => {
  assert.deepEqual(run().territory, {
    cds: 650, utr: 950, ncExon: 1566, intron: 2700,
  });
});

test('gencodeSummary counts immune receptor gene segments by locus', () => {
  const { immune } = run();
  assert.deepEqual(immune.IGH, {
    V: 1, D: 1, J: 1, C: 1, pseudo: 1,
  });
  assert.deepEqual(immune.TRB, {
    V: 1, D: 0, J: 0, C: 0, pseudo: 0,
  });
});
