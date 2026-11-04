import { CHROMS } from './parsers.mjs';
import { gtfRecord } from './sequence.mjs';
import {
  logEdges, histogram, summary, coverage,
} from '../src/stats.js';

const STANDARD = new Set([...CHROMS, 'MT']);
const LOCI = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD'];
const SEGMENTS = {
  IG_V_gene: 'V', IG_D_gene: 'D', IG_J_gene: 'J', IG_C_gene: 'C', TR_V_gene: 'V', TR_D_gene: 'D', TR_J_gene: 'J', TR_C_gene: 'C',
};

const own = (s) => (typeof s === 'string' ? Buffer.from(s, 'utf8').toString('utf8') : s);
const count = (map, key) => map.set(key, (map.get(key) ?? 0) + 1);
const push = (store, chrom, iv) => (store[chrom] ??= []).push(iv);
const lengthOf = ([s, e]) => e - s + 1;
const sumOf = (list) => list.reduce((n, iv) => n + lengthOf(iv), 0);
const halfOpen = (list) => list.map(([s, e]) => [s - 1, e]);
const tailHistogram = (values, cap) => {
  const out = new Array(cap).fill(0);
  for (const v of values) out[Math.min(v, cap) - 1] += 1;
  return out;
};
const dist = (values, edges) => ({ edges, counts: histogram(values, edges), summary: summary(values) });

export const gencodeSummary = ({ binSize = 1e6, lengths = {} } = {}) => {
  const release = {};
  const genes = new Map();
  const transcriptTypes = new Map();
  const perGene = new Map();
  let transcriptCount = 0;
  const mane = new Map();
  const cdsIv = {};
  const pcExonIv = {};
  const otherExonIv = {};
  const geneIv = {};
  const territory = {
    cds: 0, utr: 0, ncExon: 0, intron: 0,
  };
  const done = new Set();
  let current = null;
  const flush = (chrom) => {
    const cds = halfOpen(cdsIv[chrom] ?? []);
    const pcExons = halfOpen(pcExonIv[chrom] ?? []);
    const other = halfOpen(otherExonIv[chrom] ?? []);
    const spans = halfOpen(geneIv[chrom] ?? []);
    const c1 = coverage(cds);
    const c2 = coverage([...cds, ...pcExons]);
    const c3 = coverage([...cds, ...pcExons, ...other]);
    const c4 = coverage([...cds, ...pcExons, ...other, ...spans]);
    territory.cds += c1;
    territory.utr += c2 - c1;
    territory.ncExon += c3 - c2;
    territory.intron += c4 - c3;
    for (const store of [cdsIv, pcExonIv, otherExonIv, geneIv]) delete store[chrom];
    done.add(chrom);
  };
  return {
    add(line) {
      if (line.startsWith('##')) {
        const v = line.match(/version (\d+) \(Ensembl (\d+)\)/);
        if (v) [, release.gencode, release.ensembl] = v;
        const d = line.match(/^##date: (\S+)/);
        if (d) [, release.date] = d;
        return;
      }
      if (!line || line.startsWith('#')) return;
      const r = gtfRecord(line);
      if (!STANDARD.has(r.chrom)) return;
      if (r.chrom !== current) {
        if (done.has(r.chrom)) throw new Error(`GTF lines for chr${r.chrom} are not grouped`);
        if (current) flush(current);
        current = own(r.chrom);
      }
      const a = r.attrs;
      if (r.feature === 'gene') {
        genes.set(own(a.gene_id), {
          name: own(a.gene_name), type: own(a.gene_type), chrom: current, start: r.start, end: r.end, strand: own(r.strand),
        });
        push(geneIv, current, [r.start, r.end]);
        return;
      }
      if (r.feature === 'transcript') {
        transcriptCount += 1;
        count(transcriptTypes, own(a.transcript_type));
        count(perGene, own(a.gene_id));
      }
      const t = a.transcript_id;
      if (!t) return;
      const isMane = (a.tag ?? []).includes('MANE_Select');
      if (isMane && !mane.has(t)) {
        mane.set(own(t), {
          gene: own(a.gene_id), name: own(a.gene_name), chrom: current, strand: own(r.strand), exons: [], cds: [], utr: [],
        });
      }
      const iv = [r.start, r.end];
      if (r.feature === 'exon') {
        push(a.gene_type === 'protein_coding' ? pcExonIv : otherExonIv, r.chrom, iv);
        if (isMane) mane.get(t).exons.push(iv);
      } else if (r.feature === 'CDS') {
        push(cdsIv, r.chrom, iv);
        if (isMane) mane.get(t).cds.push(iv);
      } else if (r.feature === 'UTR' && isMane) mane.get(t).utr.push(iv);
    },
    result() {
      if (current) {
        flush(current);
        current = null;
      }
      const byType = new Map();
      const byChrom = {};
      const density = {};
      for (const g of genes.values()) {
        count(byType, g.type);
        byChrom[g.chrom] ??= {};
        byChrom[g.chrom][g.type] = (byChrom[g.chrom][g.type] ?? 0) + 1;
        if (g.type === 'protein_coding') {
          density[g.chrom] ??= new Array(Math.ceil((lengths[g.chrom] ?? g.end) / binSize)).fill(0);
          const bin = Math.floor((g.start - 1) / binSize);
          if (bin < density[g.chrom].length) density[g.chrom][bin] += 1;
        }
      }
      const pc = [...genes.entries()].filter(([, g]) => g.type === 'protein_coding');
      const spanOf = (g) => g.end - g.start + 1;
      const ranked = pc.map(([, g]) => ({
        name: g.name, chrom: g.chrom, start: g.start, end: g.end, span: spanOf(g),
      })).sort((x, y) => y.span - x.span);
      const lnc = [...genes.values()].filter((g) => g.type === 'lncRNA').map(spanOf);
      const isoforms = pc.map(([id, g]) => ({ name: g.name, n: perGene.get(id) ?? 0 })).filter((g) => g.n > 0);
      const structure = { exons: [], introns: [], exonLengths: [], cds: [], mrna: [], utr5: [], utr3: [] };
      const exonTop = [];
      let longestIntron = null;
      let singleExon = 0;
      for (const entry of mane.values()) {
        const m = { ...entry, name: genes.get(entry.gene)?.name ?? entry.name };
        const exons = [...m.exons].sort((x, y) => x[0] - y[0]);
        structure.exons.push(exons.length);
        exonTop.push({ name: m.name, exons: exons.length });
        if (exons.length === 1) singleExon += 1;
        exons.forEach((iv) => structure.exonLengths.push(lengthOf(iv)));
        for (let k = 1; k < exons.length; k += 1) {
          const len = exons[k][0] - exons[k - 1][1] - 1;
          structure.introns.push(len);
          if (!longestIntron || len > longestIntron.length) longestIntron = { name: m.name, chrom: m.chrom, start: exons[k - 1][1] + 1, length: len };
        }
        structure.mrna.push(sumOf(exons));
        if (!m.cds.length) continue;
        structure.cds.push(sumOf(m.cds));
        const lo = Math.min(...m.cds.map((iv) => iv[0]));
        const hi = Math.max(...m.cds.map((iv) => iv[1]));
        const before = sumOf(m.utr.filter((iv) => iv[1] < lo));
        const after = sumOf(m.utr.filter((iv) => iv[0] > hi));
        structure.utr5.push(m.strand === '-' ? after : before);
        structure.utr3.push(m.strand === '-' ? before : after);
      }
      const immune = Object.fromEntries(LOCI.map((l) => [l, {
        V: 0, D: 0, J: 0, C: 0, pseudo: 0,
      }]));
      for (const g of genes.values()) {
        if (!/^(IG|TR)_/.test(g.type)) continue;
        const locus = LOCI.find((l) => g.name.startsWith(l));
        if (!locus) continue;
        const seg = /pseudogene$/.test(g.type) ? 'pseudo' : SEGMENTS[g.type];
        if (seg) immune[locus][seg] += 1;
      }
      return {
        release,
        genes: { total: genes.size, byType: Object.fromEntries([...byType].sort((x, y) => y[1] - x[1])), byChrom },
        transcripts: { total: transcriptCount, byType: Object.fromEntries([...transcriptTypes].sort((x, y) => y[1] - x[1])) },
        density: { binSize, bins: density },
        spans: {
          pc: dist(ranked.map((g) => g.span), logEdges(1, 7, 10)),
          lnc: dist(lnc, logEdges(1, 7, 10)),
          longest: ranked.slice(0, 15),
          shortest: ranked.slice(-10).reverse(),
        },
        isoforms: {
          counts: tailHistogram(isoforms.map((g) => g.n), 50),
          top: [...isoforms].sort((x, y) => y.n - x.n).slice(0, 15),
          summary: summary(isoforms.map((g) => g.n)),
        },
        mane: {
          transcripts: mane.size,
          singleExon,
          exons: { counts: tailHistogram(structure.exons, 100), top: exonTop.sort((x, y) => y.exons - x.exons).slice(0, 15), summary: summary(structure.exons) },
          introns: { ...dist(structure.introns, logEdges(1, 7, 10)), longest: longestIntron },
          exonLengths: dist(structure.exonLengths, logEdges(0, 5, 10)),
          cds: dist(structure.cds, logEdges(1, 6, 10)),
          mrna: dist(structure.mrna, logEdges(2, 6, 10)),
          utr5: dist(structure.utr5, logEdges(0, 5, 10)),
          utr3: dist(structure.utr3, logEdges(0, 5, 10)),
        },
        territory,
        immune,
      };
    },
  };
};
