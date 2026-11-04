import { isStandard, bare } from './parsers.mjs';
import {
  logEdges, linearEdges, histogram, summary,
} from '../src/stats.js';

const own = (s) => Buffer.from(s, 'utf8').toString('utf8');
const num = (s) => Number(String(s ?? '').replace(/[^\d]/g, ''));
const binned = (values, edges) => ({ edges, counts: histogram(values, edges), summary: summary(values) });

export const dbsnpNotes = (text) => {
  const byChrom = {};
  const other = {};
  const status = {};
  let section = null;
  for (const line of text.split('\n')) {
    if (/^1\.1 /.test(line)) section = 'location';
    else if (/^1\.2 /.test(line)) section = 'type';
    else if (/^\d\. /.test(line)) section = null;
    const m = /^(.+?)\t([\d,]+)\s*$/.exec(line);
    if (!m || !section) continue;
    if (section === 'type') status[m[1]] = num(m[2]);
    else {
      const chr = /^chr(\w+)$/.exec(m[1]);
      if (chr) byChrom[chr[1] === 'M' ? 'MT' : chr[1]] = num(m[2]);
      else other[m[1]] = num(m[2]);
    }
  }
  return {
    build: num(/dbSNP build (\d+)/.exec(text)?.[1]),
    total: num(/Total RS:\s*([\d,]+)/.exec(text)?.[1]),
    byChrom,
    other,
    status,
    ss: { total: num(/Total SS:\s*([\d,]+)/.exec(text)?.[1]), unmapped: num(/Unmapped SS:\s*([\d,]+)/.exec(text)?.[1]) },
  };
};

const SV_EDGES = logEdges(0, 8, 2);

const clinicalOf = (s) => {
  if (!s) return 'none';
  const p = /pathogenic/i.test(s.replace(/conflicting interpretations of pathogenicity/gi, ''));
  const b = /benign/i.test(s);
  if (/conflicting/i.test(s) || (p && b)) return 'conflicting';
  if (p) return 'pathogenic';
  if (b) return 'benign';
  if (/uncertain/i.test(s)) return 'uncertain';
  return 'other';
};

export const svSummary = () => {
  const types = {};
  const byChrom = {};
  const sizes = {};
  const studies = new Map();
  const clinical = {
    pathogenic: 0, benign: 0, conflicting: 0, uncertain: 0, other: 0, none: 0,
  };
  return {
    add(line) {
      if (!line || line.startsWith('#')) return;
      const c = line.split('\t');
      if (!isStandard(c[0])) return;
      const chrom = bare(c[0]);
      const type = c[4];
      if (!(type in types)) types[own(type)] = 0;
      types[type] += 1;
      byChrom[chrom] ??= {};
      byChrom[chrom][type] = (byChrom[chrom][type] ?? 0) + 1;
      if (!/insertion/i.test(type)) (sizes[type] ??= []).push(Number(c[2]) - Number(c[1]) + 1);
      const study = studies.get(c[8]);
      studies.set(study === undefined ? own(c[8]) : c[8], (study ?? 0) + 1);
      clinical[clinicalOf(c[10])] += 1;
    },
    result: () => ({
      types,
      byChrom,
      sizes: Object.fromEntries(Object.entries(sizes).map(([t, list]) => [t, binned(list, SV_EDGES)])),
      studies: [...studies].map(([study, n]) => ({ study, n })).sort((a, b) => b.n - a.n || a.study.localeCompare(b.study)).slice(0, 15),
      clinical,
    }),
  };
};

const LOEUF_EDGES = linearEdges(0, 2, 0.05);
const MISZ_EDGES = linearEdges(-5, 10, 0.5);

export const constraintSummary = () => {
  let col = null;
  const values = new Map();
  const rows = [];
  return {
    add(line) {
      if (!line) return;
      const c = line.split('\t');
      if (!col) {
        col = Object.fromEntries(c.map((name, i) => [name, i]));
        return;
      }
      if (c[col.mane_select] !== 'true') return;
      const loeuf = Number(c[col['lof.oe_ci.upper']]);
      if (!Number.isFinite(loeuf)) return;
      const gene = own(c[col.gene]);
      rows.push({
        gene, loeuf, pli: Number(c[col['lof.pLI']]), expected: Number(c[col['lof.exp']]), misZ: Number(c[col['mis.z_score']]),
      });
      values.set(gene, loeuf);
    },
    values: () => values,
    result() {
      const loeufs = rows.map((r) => r.loeuf);
      return {
        genes: rows.length,
        loeuf: binned(loeufs, LOEUF_EDGES),
        below06: rows.filter((r) => r.loeuf < 0.6).length,
        pli09: rows.filter((r) => r.pli >= 0.9).length,
        misZ: binned(rows.map((r) => r.misZ).filter(Number.isFinite), MISZ_EDGES),
        top: rows.filter((r) => r.expected >= 10).sort((a, b) => a.loeuf - b.loeuf || a.gene.localeCompare(b.gene)).slice(0, 20)
          .map(({
            gene, loeuf, pli, expected,
          }) => ({
            gene, loeuf, pli, expected,
          })),
      };
    },
  };
};

export const clingenDosage = (text) => {
  const lines = text.split('\n');
  const headLine = lines.find((l) => l.startsWith('#Gene Symbol'));
  const col = Object.fromEntries(headLine.slice(1).split('\t').map((h, i) => [h, i]));
  const hi = {};
  const ts = {};
  const hiGenes = new Set();
  for (const line of lines) {
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    const h = c[col['Haploinsufficiency Score']];
    const t = c[col['Triplosensitivity Score']];
    if (/^\d+$/.test(h)) {
      hi[h] = (hi[h] ?? 0) + 1;
      if (h === '3') hiGenes.add(own(c[col['Gene Symbol']]));
    }
    if (/^\d+$/.test(t)) ts[t] = (ts[t] ?? 0) + 1;
  }
  return {
    date: lines[1]?.startsWith('#') ? lines[1].slice(1).trim() : null, hi, ts, hiGenes,
  };
};

const csvLine = (line) => {
  const out = [];
  let cur = '';
  let quoted = false;
  for (let i = 0; i < line.length; i += 1) {
    const ch = line[i];
    if (quoted) {
      if (ch === '"' && line[i + 1] === '"') {
        cur += '"';
        i += 1;
      } else if (ch === '"') quoted = false;
      else cur += ch;
    } else if (ch === '"') quoted = true;
    else if (ch === ',') {
      out.push(cur);
      cur = '';
    } else cur += ch;
  }
  out.push(cur);
  return out;
};

const SEVERITY = { lethal: 3, subviable: 2, viable: 1 };

export const impcViability = (text) => {
  const lines = text.split('\n').filter(Boolean);
  const head = csvLine(lines[0]);
  const g = head.indexOf('Gene Accession Id');
  const v = head.indexOf('Viability Phenotype HOMs/HEMIs');
  const byMgi = new Map();
  for (const line of lines.slice(1)) {
    const c = csvLine(line);
    const call = c[v];
    if (!SEVERITY[call]) continue;
    const prev = byMgi.get(c[g]);
    if (!prev || SEVERITY[call] > SEVERITY[prev]) byMgi.set(own(c[g]), call);
  }
  const calls = { lethal: 0, subviable: 0, viable: 0 };
  for (const call of byMgi.values()) calls[call] += 1;
  return { genes: byMgi.size, calls, byMgi };
};

export const hgncMouse = (text) => {
  const lines = text.split('\n');
  const head = lines[0].split('\t');
  const s = head.indexOf('symbol');
  const m = head.indexOf('mgd_id');
  const out = new Map();
  for (const line of lines.slice(1)) {
    if (!line) continue;
    const c = line.split('\t');
    for (const id of (c[m] ?? '').replace(/"/g, '').split('|').map((x) => x.trim()).filter(Boolean)) out.set(own(id), own(c[s]));
  }
  return out;
};

export const groupedValues = (values, groups, edges) => {
  const lists = {};
  for (const [gene, v] of values) (lists[groups.get(gene) ?? 'other'] ??= []).push(v);
  return Object.fromEntries(Object.entries(lists).map(([g, list]) => [g, { n: list.length, median: summary(list).median, counts: histogram(list, edges) }]));
};
