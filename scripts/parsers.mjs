import { exact } from '../src/lib.js';
import {
  linearEdges, binIndex, coverage, canonicalMotif,
} from '../src/stats.js';

export const CHROMS = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];
const STANDARD = new Set(CHROMS);

export const bare = (name) => String(name).replace(/^chr/, '');
export const isStandard = (name) => STANDARD.has(bare(name));

const rows = (text) => text.split('\n').filter((l) => l && !l.startsWith('#')).map((l) => l.replace(/\r$/, '').split('\t'));

const header = (text, key) => text.match(new RegExp(`^# ${key}:\\s*(.+?)\\s*$`, 'm'))?.[1] ?? null;

const isPlain = (v) => v === null || typeof v !== 'object';
const scalar = (v) => JSON.stringify(v) ?? 'null';

const layout = (value, indent) => {
  if (value && typeof value.toJSON === 'function') return layout(value.toJSON(), indent);
  if (isPlain(value)) return scalar(value);
  const pad = `${indent}  `;
  if (Array.isArray(value)) {
    if (value.every(isPlain)) return `[${value.map(scalar).join(', ')}]`;
    return `[\n${value.map((v) => `${pad}${layout(v, pad)}`).join(',\n')}\n${indent}]`;
  }
  const entries = Object.entries(value).filter(([, v]) => v !== undefined && typeof v !== 'function');
  if (!entries.length) return '{}';
  return `{\n${entries.map(([k, v]) => `${pad}${JSON.stringify(k)}: ${layout(v, pad)}`).join(',\n')}\n${indent}}`;
};

export const toModule = (value) => `export default ${layout(value, '')};\n`;

const newer = (a, b) => {
  const x = a.split('.').map(Number);
  const y = b.split('.').map(Number);
  for (let i = 0; i < Math.max(x.length, y.length); i += 1) {
    if ((x[i] ?? 0) !== (y[i] ?? 0)) return (x[i] ?? 0) > (y[i] ?? 0);
  }
  return false;
};

export const latestFile = (listing, re) => {
  const hits = [...listing.matchAll(new RegExp(re.source, 'g'))].map((m) => ({ file: m[0], version: m[1] }));
  return hits.length ? hits.reduce((best, h) => (newer(h.version, best.version) ? h : best)) : null;
};

export const parseAssemblyReport = (text) => ({
  name: header(text, 'Assembly name'),
  date: header(text, 'Date'),
  genbank: header(text, 'GenBank assembly accession'),
  refseq: header(text, 'RefSeq assembly accession'),
  altLoci: text.split('\n').filter((l) => /^## ALT_REF_LOCI_\d+/.test(l)).length,
  sequences: rows(text).filter((c) => c.length >= 10).map((c) => ({
    name: c[0], role: c[1], molecule: c[2], genbank: c[4], refseq: c[6], unit: c[7], length: Number(c[8]), ucsc: c[9],
  })),
});

const ROLES = ['assembled-molecule', 'unlocalized-scaffold', 'unplaced-scaffold', 'alt-scaffold', 'fix-patch', 'novel-patch'];

export const summarizeRoles = (sequences) => ROLES.map((role) => {
  const hits = sequences.filter((s) => s.role === role);
  return { role, count: hits.length, length: hits.reduce((n, s) => n + s.length, 0) };
});

const camel = (s) => s.replace(/-([a-z0-9])/gi, (_, c) => c.toUpperCase());

export const parseAssemblyStats = (text) => {
  const total = {};
  const chromosomes = new Map();
  for (const [unit, molecule, , kind, stat, value] of rows(text)) {
    if (unit === 'all' && molecule === 'all' && kind === 'all') total[camel(stat)] = Number(value);
    else if ((unit === 'Primary Assembly' || unit === 'non-nuclear') && kind === 'assembled-molecule') {
      if (!chromosomes.has(molecule)) chromosomes.set(molecule, { chrom: molecule });
      chromosomes.get(molecule)[camel(stat)] = Number(value);
    }
  }
  return { total, chromosomes: [...chromosomes.values()] };
};

export const parseCytobands = (text) => rows(text)
  .filter((c) => c.length >= 5 && isStandard(c[0]) && c[3])
  .map((c) => ({ chrom: bare(c[0]), start: Number(c[1]), end: Number(c[2]), name: c[3], stain: c[4] }));

export const parseGaps = (text) => {
  const out = {};
  for (const c of rows(text)) {
    const chrom = bare(c[1]);
    if (!STANDARD.has(chrom)) continue;
    out[chrom] ??= {};
    out[chrom][c[7]] ??= { count: 0, bp: 0 };
    out[chrom][c[7]].count += 1;
    out[chrom][c[7]].bp += Number(c[3]) - Number(c[2]);
  }
  return out;
};

const joinChromosomes = (report, stats, nuclearOnly) => {
  const byChrom = new Map(stats.chromosomes.map((c) => [c.chrom, c]));
  return report.sequences
    .filter((s) => s.role === 'assembled-molecule' && !(nuclearOnly && s.unit === 'non-nuclear'))
    .map((s) => {
      const st = byChrom.get(s.molecule) ?? {};
      return {
        chrom: s.molecule,
        length: s.length,
        ungapped: st.ungappedLength ?? null,
        gc: st.gcPerc ?? null,
        spanned: st.spannedGaps ?? null,
        unspanned: st.unspannedGaps ?? null,
        refseq: s.refseq === 'na' ? null : s.refseq,
        genbank: s.genbank,
      };
    });
};

export const buildAssembly = ({ report38, stats38, reportT2t, statsT2t, cyto, gaps, cen }) => {
  const r38 = parseAssemblyReport(report38);
  const s38 = parseAssemblyStats(stats38);
  const rT2t = parseAssemblyReport(reportT2t);
  const sT2t = parseAssemblyStats(statsT2t);
  return {
    grch38: {
      name: r38.name, date: r38.date, genbank: r38.genbank, refseq: r38.refseq, altLoci: r38.altLoci,
      roles: summarizeRoles(r38.sequences), total: s38.total, chromosomes: joinChromosomes(r38, s38, false),
    },
    t2t: {
      name: rT2t.name, date: rT2t.date, genbank: rT2t.genbank, refseq: rT2t.refseq, total: sT2t.total, chromosomes: joinChromosomes(rT2t, sT2t, true),
    },
    bands: parseCytobands(cyto),
    gaps: parseGaps(gaps),
    centromeres: parseCentromeres(cen),
  };
};

const LANDSCAPE = ['LINE', 'SINE', 'LTR', 'DNA'];
const ALU = ['AluY', 'AluS', 'AluJ'];

const tally = (map, key, make, len) => {
  if (!map.has(key)) map.set(key, { ...make, count: 0, bp: 0 });
  const t = map.get(key);
  t.count += 1;
  t.bp += len;
};

export const rmskSummary = () => {
  const classes = new Map();
  const families = new Map();
  const perChrom = {};
  const edges = linearEdges(0, 50, 1);
  const bp = Object.fromEntries(LANDSCAPE.map((c) => [c, new Array(50).fill(0)]));
  const alu = Object.fromEntries(ALU.map((a) => [a, { count: 0, bp: 0 }]));
  const l1hs = { count: 0, fullLength: 0 };
  const total = { count: 0, bp: 0 };
  let longest = null;
  return {
    add(line) {
      const c = line.split('\t');
      const chrom = bare(c[5]);
      if (!STANDARD.has(chrom)) return;
      const start = Number(c[6]);
      const len = Number(c[7]) - start;
      const cls = c[11].replace(/\?$/, '');
      const family = c[12].replace(/\?$/, '');
      const name = c[10];
      total.count += 1;
      total.bp += len;
      tally(classes, cls, { cls }, len);
      tally(families, `${cls}\t${family}`, { cls, family }, len);
      perChrom[chrom] = (perChrom[chrom] ?? 0) + len;
      if (bp[cls]) bp[cls][binIndex(edges, Number(c[2]) / 10)] += len;
      if (family === 'Alu' && alu[name.slice(0, 4)]) {
        alu[name.slice(0, 4)].count += 1;
        alu[name.slice(0, 4)].bp += len;
      }
      if (name === 'L1HS') {
        l1hs.count += 1;
        if (len >= 6000) l1hs.fullLength += 1;
      }
      if (!longest || len > longest.length) longest = { name, cls, chrom, start: start + 1, length: len };
    },
    result() {
      const byBp = (a, b) => b.bp - a.bp;
      return {
        total,
        classes: [...classes.values()].sort(byBp),
        families: [...families.values()].sort(byBp),
        perChrom,
        landscape: { classes: LANDSCAPE, edges, bp },
        alu,
        l1hs,
        longest,
      };
    },
  };
};

export const superDupsSummary = () => {
  const intervals = {};
  const pairs = new Map();
  const edges = linearEdges(0.9, 1, 0.005);
  let longest = null;
  return {
    add(line) {
      const c = line.split('\t');
      const chrom = bare(c[1]);
      if (!STANDARD.has(chrom)) return;
      const s = Number(c[2]);
      const e = Number(c[3]);
      (intervals[chrom] ??= []).push([s, e]);
      if (!pairs.has(c[11])) pairs.set(c[11], { identity: Number(c[26]), intra: c[1] === c[7] });
      if (!longest || e - s > longest.length) {
        longest = { chrom, start: s + 1, length: e - s, other: `${c[7]}:${exact(Number(c[8]) + 1)}-${exact(Number(c[9]))}`, identity: Number(c[26]) };
      }
    },
    result() {
      const counts = new Array(edges.length - 1).fill(0);
      const kinds = { intra: 0, inter: 0 };
      for (const p of pairs.values()) {
        counts[binIndex(edges, p.identity)] += 1;
        kinds[p.intra ? 'intra' : 'inter'] += 1;
      }
      return {
        pairs: pairs.size,
        kinds,
        covered: Object.fromEntries(Object.entries(intervals).map(([chrom, list]) => [chrom, coverage(list)])),
        identity: { edges, counts },
        longest,
      };
    },
  };
};

const UNITS = ['1', '2', '3', '4', '5', '6', '7-12', '13-50', '51+'];
const unitOf = (n) => (n <= 6 ? String(n) : n <= 12 ? '7-12' : n <= 50 ? '13-50' : '51+');

export const tandemSummary = () => {
  const units = new Map(UNITS.map((u) => [u, { unit: u, count: 0, bp: 0 }]));
  const motifs = new Map();
  let longest = null;
  return {
    add(line) {
      const c = line.split('\t');
      const chrom = bare(c[1]);
      if (!STANDARD.has(chrom)) return;
      const start = Number(c[2]);
      const len = Number(c[3]) - start;
      const size = Number(c[7]);
      const u = units.get(unitOf(size));
      u.count += 1;
      u.bp += len;
      if (size <= 6) tally(motifs, canonicalMotif(c[16]), { motif: canonicalMotif(c[16]) }, len);
      if (!longest || len > longest.length) longest = { chrom, start: start + 1, length: len, motif: c[16] };
    },
    result() {
      return {
        units: [...units.values()],
        motifs: [...motifs.values()].sort((a, b) => b.count - a.count || b.bp - a.bp).slice(0, 30),
        longest,
      };
    },
  };
};

export const parseCentromeres = (text) => {
  const out = {};
  for (const c of rows(text)) {
    const chrom = bare(c[1]);
    if (!STANDARD.has(chrom)) continue;
    out[chrom] = (out[chrom] ?? 0) + Number(c[3]) - Number(c[2]);
  }
  return out;
};
