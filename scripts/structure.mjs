import { isStandard, bare } from './parsers.mjs';
import { linearEdges, histogram } from '../src/stats.js';

const KINDS = { CDS: 'protein', tRNA: 'tRNA', rRNA: 'rRNA' };
const position = (s) => Number(String(s).replace(/\D/g, ''));
const span = ([a, b]) => [Math.min(a, b), Math.max(a, b)];

export const mitoFeatures = (text) => {
  const features = [];
  let current = null;
  for (const line of text.split('\n')) {
    if (!line || line.startsWith('>')) continue;
    if (line.startsWith('\t\t\t')) {
      const [, , , key, value] = line.split('\t');
      if (current && !(key in current.quals)) current.quals[key] = value;
      continue;
    }
    const [a, b, key] = line.split('\t');
    if (key) {
      current = { key, ranges: [[position(a), position(b)]], quals: {} };
      features.push(current);
    } else if (current) current.ranges.push([position(a), position(b)]);
  }
  const genes = [];
  features.forEach((f, i) => {
    if (f.key !== 'gene') return;
    const [a, b] = f.ranges[0];
    const product = features.slice(i + 1).find((g) => KINDS[g.key] && g.ranges[0][0] === a && g.ranges[0][1] === b);
    const [start, end] = span([a, b]);
    genes.push({
      name: f.quals.gene, start, end, strand: a > b ? '-' : '+', kind: product ? KINDS[product.key] : null, product: product?.quals.product ?? null,
    });
  });
  const loop = features.find((f) => f.key === 'D-loop');
  const ranges = loop ? (loop.ranges[0][0] > loop.ranges[0][1] ? [...loop.ranges].reverse() : loop.ranges) : [];
  return { genes, dloop: ranges.map(span) };
};

const HETEROPLASMY_EDGES = linearEdges(0, 1, 0.05);

export const mitoVariants = ({ length = 16569 } = {}) => {
  let col = null;
  let variants = 0;
  let homoplasmic = 0;
  let heteroplasmicOnly = 0;
  let samples = 0;
  const positions = new Set();
  const levels = [];
  const list = [];
  const bins = Array.from({ length: Math.ceil(length / 100) }, () => ({ hom: 0, het: 0 }));
  return {
    add(line) {
      if (!line) return;
      const c = line.split('\t');
      if (!col) {
        col = Object.fromEntries(c.map((name, i) => [name, i]));
        return;
      }
      if (c[col.filters] !== 'PASS') return;
      const pos = Number(c[col.position]);
      const hom = Number(c[col.AC_hom]);
      const het = Number(c[col.AC_het]);
      variants += 1;
      positions.add(pos);
      if (hom > 0) homoplasmic += 1;
      else if (het > 0) heteroplasmicOnly += 1;
      bins[Math.min(bins.length - 1, Math.floor((pos - 1) / 100))][hom > 0 ? 'hom' : 'het'] += 1;
      samples = Math.max(samples, Number(c[col.AN]));
      levels.push(Number(c[col.max_observed_heteroplasmy]));
      list.push({
        pos, ref: c[col.ref], alt: c[col.alt], afHom: Number(c[col.AF_hom]), afHet: Number(c[col.AF_het]),
      });
    },
    result: () => ({
      variants,
      positions: positions.size,
      homoplasmic,
      heteroplasmicOnly,
      samples,
      bins,
      top: list.filter((v) => v.afHom > 0).sort((a, b) => b.afHom - a.afHom).slice(0, 15),
      heteroplasmy: { edges: HETEROPLASMY_EDGES, counts: histogram(levels, HETEROPLASMY_EDGES) },
    }),
  };
};

export const censatSummary = (text) => {
  const classes = {};
  const byChrom = {};
  let transition = 0;
  for (const line of text.split('\n')) {
    if (!line) continue;
    const c = line.split('\t');
    if (!isStandard(c[0])) continue;
    const chrom = bare(c[0]);
    const bp = Number(c[2]) - Number(c[1]);
    const cls = c[3].split('(')[0].split('_')[0];
    if (cls === 'ct') {
      transition += bp;
      continue;
    }
    classes[cls] = (classes[cls] ?? 0) + bp;
    byChrom[chrom] ??= {};
    byChrom[chrom][cls] = (byChrom[chrom][cls] ?? 0) + bp;
  }
  return {
    classes,
    transition,
    total: Object.values(classes).reduce((a, b) => a + b, 0),
    byChrom,
    arrays: { hor: Object.entries(byChrom).filter(([, v]) => v.hor).map(([chrom, v]) => ({ chrom, bp: v.hor })) },
  };
};

export const telomeres = (text) => {
  const out = new Map();
  for (const line of text.split('\n')) {
    if (!line) continue;
    const [chr, s, e] = line.split('\t');
    const chrom = bare(chr);
    if (!out.has(chrom)) out.set(chrom, { chrom, p: null, q: null });
    out.get(chrom)[Number(s) === 0 ? 'p' : 'q'] = Number(e) - Number(s);
  }
  return [...out.values()];
};
