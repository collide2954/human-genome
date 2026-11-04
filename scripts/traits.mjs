import { csvLine } from './variation.mjs';
import { tsvRecords, bump, round } from './tables.mjs';
import {
  summary, quantile, logEdges, linearEdges, histogram,
} from '../src/stats.js';

const own = (s) => Buffer.from(String(s), 'utf8').toString('utf8');
const MLOG_5E8 = -Math.log10(5e-8);
const CHROMS = new Set([...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y', 'MT']);
const DISEASES = new Set(['Cancer', 'Cardiovascular disease', 'Digestive system disorder', 'Immune system disorder', 'Metabolic disorder', 'Neurological disorder', 'Other disease']);
const OR_EDGES = [1, 1.1, 1.25, 1.5, 2, 3, 5, 10];
const RAF_EDGES = linearEdges(0, 1, 0.1);
const SIZE_EDGES = [0, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6];
const SIZE_LABELS = ['<1K', '1K–3K', '3K–10K', '10K–30K', '30K–100K', '100K–300K', '300K–1M', '≥1M'];
const ASSOC_EDGES = [0, 1, 2, 4, 11, 31, 101, 301, 1001];
const ASSOC_LABELS = ['0', '1', '2–3', '4–10', '11–30', '31–100', '101–300', '301–1,000', '>1,000'];
const VARIANT_EDGES = logEdges(0, 8, 2);
const AUROC_EDGES = linearEdges(0.4, 1, 0.025);
const H2_EDGES = linearEdges(0, 0.8, 0.025);

const binOf = (edges, v) => {
  let i = 0;
  while (i + 1 < edges.length && v >= edges[i + 1]) i += 1;
  return i;
};
const ranked = (map, key) => [...map].map(([k, n]) => ({ [key]: k, n })).sort((a, b) => b.n - a.n || String(a[key]).localeCompare(String(b[key])));
const quartiles = (values, places) => {
  const s = [...values].sort((a, b) => a - b);
  return { q1: round(quantile(s, 0.25), places), median: round(quantile(s, 0.5), places), q3: round(quantile(s, 0.75), places) };
};
const middle = (values, places) => (values.length ? round(summary(values).median, places) : null);

export const traitParents = (text) => {
  const uris = new Map();
  const labels = new Map();
  for (const r of tsvRecords(text)) {
    const parent = r['Parent term'];
    if (!parent) continue;
    if (r['EFO URI'] && !uris.has(r['EFO URI'])) uris.set(own(r['EFO URI']), own(parent));
    const label = (r['EFO term'] ?? '').toLowerCase();
    if (label && !labels.has(label)) labels.set(own(label), own(parent));
  }
  return { uris, labels };
};

export const gwasAssociations = ({ parents }) => {
  let col = null;
  let associations = 0;
  let significant = 0;
  let placed = 0;
  let latestAdded = '';
  const perStudy = new Map();
  const years = new Map();
  const bins = new Map();
  const totals = {};
  const byYear = {};
  const traits = new Map();
  const genes = new Map();
  const context = {};
  const counts = OR_EDGES.map(() => RAF_EDGES.slice(1).map(() => 0));
  const ors = [];
  const parentOf = (uris, labels) => {
    for (const uri of uris.split(',')) {
      const p = parents.uris.get(uri.trim());
      if (p) return p;
    }
    for (const label of labels.split(',')) {
      const p = parents.labels.get(label.trim().toLowerCase());
      if (p) return p;
    }
    return 'Not mapped';
  };
  return {
    perStudy,
    add(line) {
      const c = line.split('\t');
      if (!col) {
        col = Object.fromEntries(c.map((h, i) => [h.trim(), i]));
        return;
      }
      associations += 1;
      const study = c[col['STUDY ACCESSION']];
      perStudy.set(perStudy.has(study) ? study : own(study), (perStudy.get(study) ?? 0) + 1);
      const added = c[col['DATE ADDED TO CATALOG']] ?? '';
      if (added > latestAdded) latestAdded = own(added);
      const year = /^\d{4}/.exec(c[col.DATE] ?? '')?.[0];
      const trait = c[col.MAPPED_TRAIT] ?? '';
      const parent = parentOf(c[col.MAPPED_TRAIT_URI] ?? '', trait);
      bump(totals, parent);
      if (year) {
        if (!years.has(year)) years.set(year, { associations: 0, studies: new Set() });
        const y = years.get(year);
        y.associations += 1;
        if (!y.studies.has(study)) y.studies.add(own(study));
        byYear[year] ??= {};
        bump(byYear[year], parent);
      }
      const t = traits.get(trait);
      if (t) t.n += 1;
      else traits.set(own(trait), { n: 1, parent });
      const gene = c[col.MAPPED_GENE] ?? '';
      for (const g of new Set(gene.split(/\s+-\s+|,\s*|;\s*/).map((x) => x.trim()).filter(Boolean))) genes.set(genes.has(g) ? g : own(g), (genes.get(g) ?? 0) + 1);
      const ctx = c[col.CONTEXT] ?? '';
      bump(context, !ctx ? 'Not given' : / x /.test(ctx) ? 'Interaction' : ctx.includes(';') ? 'Several variants' : own(ctx));
      const mlog = Number(c[col.PVALUE_MLOG]);
      if (Number.isFinite(mlog) && mlog >= MLOG_5E8) significant += 1;
      const chrom = c[col.CHR_ID];
      const pos = Number(c[col.CHR_POS]);
      if (CHROMS.has(chrom) && c[col.CHR_POS] && Number.isFinite(pos) && Number.isFinite(mlog)) {
        placed += 1;
        const bin = Math.floor(pos / 1e6);
        const key = `${chrom}:${bin}`;
        const e = bins.get(key);
        if (!e) bins.set(key, [chrom, bin, mlog, 1, own(trait), own(gene), own(c[col.SNPS] ?? '')]);
        else {
          e[3] += 1;
          if (mlog > e[2]) {
            e[2] = mlog;
            e[4] = own(trait);
            e[5] = own(gene);
            e[6] = own(c[col.SNPS] ?? '');
          }
        }
      }
      if (DISEASES.has(parent)) {
        const or = Number(c[col['OR or BETA']]);
        const raf = Number(c[col['RISK ALLELE FREQUENCY']]);
        if (or > 1 && raf > 0 && raf < 1 && !/increase|decrease|unit|s\.d\./i.test(c[col['95% CI (TEXT)']] ?? '')) {
          counts[binOf(OR_EDGES, or)][Math.min(RAF_EDGES.length - 2, Math.floor(raf * 10))] += 1;
          ors.push(or);
        }
      }
    },
    result() {
      const manhattan = {};
      for (const [chrom, bin, mlog, n, trait, gene, snp] of bins.values()) {
        manhattan[chrom] ??= [];
        manhattan[chrom].push([bin, round(mlog, 1), n, trait, gene, snp]);
      }
      for (const list of Object.values(manhattan)) list.sort((a, b) => a[0] - b[0]);
      return {
        associations,
        studies: perStudy.size,
        significant,
        placed,
        latestAdded,
        byYear: Object.fromEntries([...years].sort((a, b) => a[0].localeCompare(b[0])).map(([y, e]) => [y, { associations: e.associations, studies: e.studies.size }])),
        manhattan,
        categories: { totals, byYear },
        topTraits: [...traits].map(([trait, e]) => ({ trait, n: e.n, parent: e.parent })).sort((a, b) => b.n - a.n || a.trait.localeCompare(b.trait)).slice(0, 20),
        topGenes: ranked(genes, 'gene').slice(0, 20),
        context,
        effects: {
          orEdges: OR_EDGES, rafEdges: RAF_EDGES, counts, n: ors.length, median: middle(ors, 3),
        },
      };
    },
  };
};

export const studySizes = (perStudy, ancestryText) => {
  const lines = ancestryText.split('\n');
  const head = lines[0].split('\t');
  const at = (re) => head.findIndex((h) => re.test(h));
  const col = { acc: at(/^STUDY ACCESSION/), stage: at(/^STAGE/), n: at(/^NUMBER OF INDIV/) };
  const sizes = new Map();
  for (const line of lines.slice(1)) {
    if (!line) continue;
    const c = line.split('\t');
    const n = Number(c[col.n]);
    if (c[col.stage] !== 'initial' || !(n > 0)) continue;
    sizes.set(c[col.acc], (sizes.get(c[col.acc]) ?? 0) + n);
  }
  const counts = ASSOC_LABELS.map(() => SIZE_LABELS.map(() => 0));
  for (const [study, n] of sizes) counts[binOf(ASSOC_EDGES, perStudy.get(study) ?? 0)][binOf(SIZE_EDGES, n)] += 1;
  return {
    sizeLabels: SIZE_LABELS, assocLabels: ASSOC_LABELS, counts, studies: sizes.size,
  };
};

const FAMILIES = [
  [/sbayes/i, 'SBayes'],
  [/ldpred/i, 'LDpred'],
  [/prs-?cs/i, 'PRS-CS'],
  [/megaprs|ldak|bayesr/i, 'MegaPRS and LDAK'],
  [/snpnet|lasso|bigstatsr|penali[sz]ed|elastic.?net|ridge/i, 'Penalized regression'],
  [/prsice|p\s*\+\s*t|c\s*\+\s*t|prun|clump|threshold|plink/i, 'Pruning and thresholding'],
  [/genome-wide significant|gwas hits|significant (snps|variants)|snps passing|lead (snps|variants)|known (snps|variants)/i, 'Genome-wide significant variants'],
];

export const methodFamily = (name) => FAMILIES.find(([re]) => re.test(name ?? ''))?.[1] ?? 'Other';

const csvRecords = (text) => {
  const lines = text.split('\n').map((l) => l.replace(/\r$/, '')).filter(Boolean);
  if (!lines.length) return [];
  const head = csvLine(lines[0]).map((h) => h.trim());
  return lines.slice(1).map((line) => {
    const c = csvLine(line);
    return Object.fromEntries(head.map((h, i) => [h, (c[i] ?? '').trim()]));
  });
};

export const pgsScores = (text) => {
  const rows = csvRecords(text);
  const families = {};
  const byYear = {};
  const traits = new Map();
  const methods = new Set();
  const variants = [];
  for (const r of rows) {
    const family = methodFamily(r['PGS Development Method']);
    methods.add(r['PGS Development Method']);
    bump(families, family);
    const year = r['Release Date'].slice(0, 4);
    if (year) {
      byYear[year] ??= {};
      bump(byYear[year], family);
    }
    const v = Number(r['Number of Variants']);
    if (v > 0) variants.push(v);
    for (const t of new Set(r['Mapped Trait(s) (EFO label)'].split('|').map((x) => x.trim()).filter(Boolean))) traits.set(t, (traits.get(t) ?? 0) + 1);
  }
  return {
    scores: rows.length,
    methods: methods.size,
    families,
    byYear,
    variants: { edges: VARIANT_EDGES, counts: histogram(variants, VARIANT_EDGES), median: middle(variants, 1) },
    traits: ranked(traits, 'trait').slice(0, 20),
  };
};

const ancestryOf = (raw) => (raw.includes(',') ? 'Multiple categories' : raw || 'Not reported');

export const pgsAncestry = (devText, evalText) => {
  const out = { gwas: {}, development: {}, evaluation: {} };
  for (const r of csvRecords(devText)) {
    const n = Number(r['Number of Individuals']);
    if (n > 0) bump(out[/GWAS/.test(r['Stage of PGS Development']) ? 'gwas' : 'development'], ancestryOf(r['Broad Ancestry Category']), n);
  }
  for (const r of csvRecords(evalText)) {
    const n = Number(r['Number of Individuals']);
    if (n > 0) bump(out.evaluation, ancestryOf(r['Broad Ancestry Category']), n);
  }
  return out;
};

export const pgsCategories = (json) => {
  const map = new Map();
  for (const c of json?.results ?? []) {
    for (const t of c.efotraits ?? []) if (!map.has(t.id)) map.set(t.id, c.label);
  }
  return map;
};

const firstNumber = (s) => {
  const m = /^\s*(-?\d+(?:\.\d+)?)/.exec(s ?? '');
  return m ? Number(m[1]) : null;
};

export const pgsPerformance = (metricsText, scoresText, categories) => {
  const efo = new Map(csvRecords(scoresText).map((r) => [r['Polygenic Score (PGS) ID'], r['Mapped Trait(s) (EFO ID)'].split('|').map((x) => x.trim())]));
  const categoryOf = (pgs) => (efo.get(pgs) ?? []).map((id) => categories.get(id)).find(Boolean) ?? 'Not categorized';
  const rows = csvRecords(metricsText);
  const auroc = [];
  const byCategory = new Map();
  const cindex = [];
  const r2 = [];
  for (const r of rows) {
    const a = firstNumber(r['Area Under the Receiver-Operating Characteristic Curve (AUROC)']);
    if (a !== null && a > 0 && a <= 1) {
      auroc.push(a);
      const category = categoryOf(r['Evaluated Score']);
      if (!byCategory.has(category)) byCategory.set(category, []);
      byCategory.get(category).push(a);
    }
    const ci = firstNumber(r['Concordance Statistic (C-index)']);
    if (ci !== null && ci > 0 && ci <= 1) cindex.push(ci);
    const m = /^\s*R(?:²|2)\s*[:=]\s*(-?\d+(?:\.\d+)?)/i.exec(r['Other Metric(s)']);
    if (m) r2.push(Number(m[1]));
  }
  return {
    metrics: rows.length,
    auroc: {
      n: auroc.length,
      median: middle(auroc, 3),
      edges: AUROC_EDGES,
      counts: histogram(auroc, AUROC_EDGES),
      byCategory: [...byCategory].map(([category, v]) => ({ category, n: v.length, ...quartiles(v, 3) })).sort((a, b) => b.n - a.n || a.category.localeCompare(b.category)),
    },
    cindex: { n: cindex.length, median: middle(cindex, 3) },
    r2: { n: r2.length, median: middle(r2, 3) },
  };
};

const phenoKey = (r) => [r.trait_type, r.phenocode, r.pheno_sex, r.coding, r.modifier].join('|');

export const heritability = (h2Text, phenoText) => {
  const phenos = tsvRecords(phenoText);
  const names = new Map(phenos.map((p) => [phenoKey(p), p.description]));
  const byPop = {};
  const eur = [];
  for (const r of tsvRecords(h2Text)) {
    const raw = r['estimates.final.h2_liability'];
    const h2 = Number(raw);
    if (!raw || !Number.isFinite(h2)) continue;
    byPop[r.pop] ??= { estimates: 0, passing: 0 };
    byPop[r.pop].estimates += 1;
    if (r['qcflags.pass_all'] !== 'true') continue;
    byPop[r.pop].passing += 1;
    if (r.pop === 'EUR') {
      eur.push({
        description: names.get(phenoKey(r)) ?? r.phenocode, type: r.trait_type, h2, se: Number(r['estimates.final.h2_liability_se']),
      });
    }
  }
  const types = new Map();
  for (const e of eur) {
    if (!types.has(e.type)) types.set(e.type, []);
    types.get(e.type).push(e.h2);
  }
  return {
    phenotypes: phenos.length,
    byPop,
    eur: {
      passing: eur.length,
      median: middle(eur.map((e) => e.h2), 4),
      edges: H2_EDGES,
      counts: histogram(eur.map((e) => e.h2), H2_EDGES),
      byType: [...types].map(([type, v]) => ({ type, n: v.length, ...quartiles(v, 4) })).sort((a, b) => b.median - a.median || a.type.localeCompare(b.type)),
      top: [...eur].sort((a, b) => b.h2 - a.h2).slice(0, 20).map((e) => ({
        description: e.description, type: e.type, h2: round(e.h2, 4), se: round(e.se, 4),
      })),
    },
  };
};
