import { isStandard, bare } from './parsers.mjs';
import {
  linearEdges, logEdges, histogram, summary,
} from '../src/stats.js';

const own = (s) => Buffer.from(s, 'utf8').toString('utf8');
const unquote = (s) => s.replace(/^"|"$/g, '');
const binned = (values, edges) => ({ edges, counts: histogram(values, edges), summary: summary(values) });

const WIDTH_EDGES = linearEdges(0, 500, 25);
const LENGTH_EDGES = logEdges(2, 5, 10);
const OE_EDGES = linearEdges(0.6, 1.4, 0.05);
const GC_EDGES = linearEdges(50, 80, 2);

export const ccreSummary = () => {
  const classes = {};
  const byChrom = {};
  const widths = [];
  return {
    add(line) {
      if (!line) return;
      const c = line.split('\t');
      if (!isStandard(c[0])) return;
      const chrom = bare(c[0]);
      const width = Number(c[2]) - Number(c[1]);
      const cls = c[5];
      if (!classes[cls]) classes[own(cls)] = { count: 0, bp: 0 };
      classes[cls].count += 1;
      classes[cls].bp += width;
      if (!byChrom[chrom]) byChrom[own(chrom)] = {};
      byChrom[chrom][cls] = (byChrom[chrom][cls] ?? 0) + 1;
      widths.push(width);
    },
    result: () => ({
      count: widths.length, classes, byChrom, widths: binned(widths, WIDTH_EDGES),
    }),
  };
};

export const cpgIslands = (text) => {
  const intervals = new Map();
  const byChrom = {};
  const lengths = [];
  const oe = [];
  const gc = [];
  let bp = 0;
  for (const line of text.split('\n')) {
    if (!line) continue;
    const c = line.split('\t');
    if (!isStandard(c[1])) continue;
    const chrom = bare(c[1]);
    const start = Number(c[2]);
    const end = Number(c[3]);
    if (!intervals.has(chrom)) intervals.set(chrom, []);
    intervals.get(chrom).push([start, end]);
    byChrom[chrom] = (byChrom[chrom] ?? 0) + 1;
    lengths.push(end - start);
    bp += end - start;
    oe.push(Number(c[10]));
    gc.push(Number(c[9]));
  }
  for (const list of intervals.values()) list.sort((a, b) => a[0] - b[0]);
  return {
    count: lengths.length,
    bp,
    byChrom,
    lengths: binned(lengths, LENGTH_EDGES),
    obsExp: binned(oe, OE_EDGES),
    gc: binned(gc, GC_EDGES),
    intervals,
  };
};

const inside = (list, pos) => {
  let lo = 0;
  let hi = list.length - 1;
  let found = -1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    if (list[mid][0] <= pos) {
      found = mid;
      lo = mid + 1;
    } else hi = mid - 1;
  }
  return found >= 0 && pos < list[found][1];
};

export const tssInIslands = (transcripts, intervals, specificity = new Map()) => {
  let tss = 0;
  let inIsland = 0;
  const bySpecificity = {};
  for (const t of transcripts) {
    if (t.tag !== 'MANE_Select' || !t.exons.length) continue;
    const pos = (t.strand === '-' ? t.exons[0][1] : t.exons[0][0]) - 1;
    const hit = inside(intervals.get(t.chrom) ?? [], pos);
    tss += 1;
    if (hit) inIsland += 1;
    const cat = specificity.get(t.gene);
    if (!cat) continue;
    bySpecificity[cat] ??= { genes: 0, inIsland: 0 };
    bySpecificity[cat].genes += 1;
    if (hit) bySpecificity[cat].inIsland += 1;
  }
  return { tss, inIsland, bySpecificity };
};

export const specificityByGene = (text) => {
  const lines = text.split('\n').filter(Boolean);
  const head = lines[0].split('\t').map(unquote);
  const g = head.indexOf('Gene');
  const s = head.indexOf('RNA tissue specificity');
  const out = new Map();
  for (const line of lines.slice(1)) {
    const c = line.replace(/\r$/, '').split('\t').map(unquote);
    if (c[g] && c[s]) out.set(own(c[g]), own(c[s]));
  }
  return out;
};

const bumpKey = (obj, key) => {
  obj[key] = (obj[key] ?? 0) + 1;
};

export const trnaSummary = (fasta) => {
  const isotypes = {};
  const anticodons = {};
  const byChrom = {};
  let genes = 0;
  for (const line of fasta.split('\n')) {
    if (!line.startsWith('>')) continue;
    const m = /tRNA-(\w+)-([ACGTN]{3})-\d+-\d+/.exec(line);
    if (!m) continue;
    genes += 1;
    bumpKey(isotypes, m[1]);
    bumpKey(anticodons, `${m[1]}-${m[2]}`);
    const chrom = /\schr([0-9XYM]+):\d+-\d+/.exec(line)?.[1];
    if (chrom) bumpKey(byChrom, chrom);
  }
  return {
    genes, isotypes, anticodons, byChrom,
  };
};

const FACETS = {
  assays: 'assay_title', biosamples: 'biosample_ontology.classification', phases: 'award.rfa', targets: 'target.investigated_as',
};

export const encodeSummary = (json) => ({
  total: json.total,
  ...Object.fromEntries(Object.entries(FACETS).map(([key, field]) => [key, (json.facets.find((f) => f.field === field)?.terms ?? [])
    .filter((t) => t.doc_count > 0)
    .map((t) => ({ name: t.key, n: t.doc_count }))])),
});
