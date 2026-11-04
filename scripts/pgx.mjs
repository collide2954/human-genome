import { tsvRecords, bump, round } from './tables.mjs';

export const CPIC_GROUPS = ['African American/Afro-Caribbean', 'American', 'Central/South Asian', 'East Asian', 'European', 'Latino', 'Near Eastern', 'Oceanian', 'Sub-Saharan African'];
export const SCORE_BINS = ['0', '0.25', '0.5', '0.75', '1', '1.25', '1.5', '1.75', '2', '2.25', '2.5', '2.75', '≥3', 'n/a'];
export const FUNCTIONS = ['Normal function', 'Decreased function', 'No function', 'Increased function', 'Uncertain function', 'Unknown function'];
export const PHENOTYPES = ['Poor Metabolizer', 'Intermediate Metabolizer', 'Normal Metabolizer', 'Ultrarapid Metabolizer', 'Indeterminate'];

const HIGH_RISK = /priority\/high.?risk/i;
const present = (v) => v !== null && v !== undefined;
const list = (s) => String(s ?? '').split(';').map((x) => x.trim()).filter(Boolean);
const add6 = (a, b) => round((a ?? 0) + b, 6);

export const cpicPairs = (pairs) => {
  const current = pairs.filter((p) => !p.removed);
  const levels = {};
  const clinpgx = {};
  const testing = {};
  const genes = new Map();
  for (const p of current) {
    const level = p.cpiclevel ?? 'none';
    bump(levels, level);
    clinpgx[level] ??= {};
    bump(clinpgx[level], p.clinpgxlevel ?? 'none');
    bump(testing, p.pgxtesting ?? 'none');
    const g = genes.get(p.genesymbol) ?? { gene: p.genesymbol, total: 0, levels: {} };
    g.total += 1;
    bump(g.levels, level);
    genes.set(p.genesymbol, g);
  }
  return {
    pairs: current.length,
    removed: pairs.length - current.length,
    genes: genes.size,
    drugs: new Set(current.map((p) => p.drugid)).size,
    levels,
    byGene: [...genes.values()].sort((a, b) => (b.levels.A ?? 0) - (a.levels.A ?? 0) || b.total - a.total || a.gene.localeCompare(b.gene)),
    clinpgx,
    testing,
  };
};

export const cpicAlleles = (alleles) => {
  const genes = new Map();
  for (const a of alleles) {
    const g = genes.get(a.genesymbol) ?? {
      gene: a.genesymbol, total: 0, functions: {}, classes: [],
    };
    const status = a.clinicalfunctionalstatus;
    const known = status ? FUNCTIONS.find((f) => f.toLowerCase() === status.toLowerCase()) : null;
    g.total += 1;
    if (!status) bump(g.functions, 'Not assigned');
    else if (known) bump(g.functions, known);
    else {
      bump(g.functions, 'Gene-specific class');
      if (!g.classes.includes(status)) g.classes.push(status);
    }
    genes.set(a.genesymbol, g);
  }
  return { alleles: alleles.length, genes: [...genes.values()].sort((a, b) => b.total - a.total || a.gene.localeCompare(b.gene)) };
};

export const cpicPhenotypes = (results) => {
  const genes = new Map();
  for (const r of results) {
    const g = genes.get(r.genesymbol) ?? {
      gene: r.genesymbol, results: 0, highRisk: [], sums: {}, seen: new Set(),
    };
    const risky = HIGH_RISK.test(r.ehrpriority ?? '');
    g.results += 1;
    if (risky && !g.highRisk.includes(r.result)) g.highRisk.push(r.result);
    for (const [group, f] of Object.entries(r.frequency ?? {})) {
      if (!present(f)) continue;
      g.seen.add(group);
      if (risky) g.sums[group] = add6(g.sums[group], f);
    }
    genes.set(r.genesymbol, g);
  }
  return {
    groups: CPIC_GROUPS,
    genes: [...genes.values()].filter((g) => g.seen.size).sort((a, b) => a.gene.localeCompare(b.gene)).map((g) => ({
      gene: g.gene, results: g.results, highRisk: g.highRisk, shares: Object.fromEntries([...g.seen].map((group) => [group, g.sums[group] ?? 0])),
    })),
  };
};

const scoreBin = (raw) => {
  const s = String(raw ?? '').trim();
  const v = Number(s.replace('≥', ''));
  if (!s || !Number.isFinite(v)) return 'n/a';
  if (v >= 3) return '≥3';
  return SCORE_BINS.includes(String(v)) ? String(v) : 'n/a';
};

export const cyp2d6 = (alleles, results) => {
  const functions = {};
  const mix = {};
  const ranked = [];
  for (const a of alleles) {
    const fn = a.clinicalfunctionalstatus ?? 'Not assigned';
    const freqs = a.frequency ?? {};
    bump(functions, fn);
    for (const [group, f] of Object.entries(freqs)) {
      if (!present(f)) continue;
      mix[group] ??= {};
      mix[group][fn] = add6(mix[group][fn], f);
    }
    const values = Object.values(freqs).filter(present);
    const max = values.length ? Math.max(...values) : 0;
    if (max > 0) {
      ranked.push({
        name: a.name, function: fn, activity: a.activityvalue ?? null, freqs, max,
      });
    }
  }
  const phenotypes = {};
  const groups = {};
  for (const r of results) {
    const phenotype = PHENOTYPES.includes(r.result) ? r.result : 'Indeterminate';
    const bin = SCORE_BINS.indexOf(scoreBin(r.activityscore));
    for (const [group, f] of Object.entries(r.frequency ?? {})) {
      if (!present(f)) continue;
      phenotypes[group] ??= {};
      phenotypes[group][phenotype] = add6(phenotypes[group][phenotype], f);
      groups[group] ??= SCORE_BINS.map(() => 0);
      groups[group][bin] = add6(groups[group][bin], f);
    }
  }
  return {
    alleles: alleles.length,
    functions,
    top: ranked.sort((a, b) => b.max - a.max || a.name.localeCompare(b.name)).slice(0, 16).map((a) => ({
      name: a.name, function: a.function, activity: a.activity, freqs: a.freqs,
    })),
    phenotypes,
    scores: { bins: SCORE_BINS, groups },
    mix,
  };
};

const ENTITIES = {
  '&amp;': '&', '&lt;': '<', '&gt;': '>', '&quot;': '"', '&#39;': '\'', '&nbsp;': ' ',
};
const cellText = (html) => html.replace(/<[^>]+>/g, ' ').replace(/&(?:amp|lt|gt|quot|#39|nbsp);/g, (m) => ENTITIES[m]).replace(/\s+/g, ' ').trim();
const AREA_FIXES = { Onoclogy: 'Oncology' };
const rankRows = (map, key) => [...map].map(([k, n]) => ({ [key]: k, rows: n })).sort((a, b) => b.rows - a.rows || a[key].localeCompare(b[key]));
const tally = (map, key) => map.set(key, (map.get(key) ?? 0) + 1);

export const fdaBiomarkers = (html) => {
  const table = /<table[^>]*>([\s\S]*?)<\/table>/.exec(html)?.[1];
  if (!table) return null;
  const rows = [...table.matchAll(/<tr[^>]*>([\s\S]*?)<\/tr>/g)]
    .map((m) => [...m[1].matchAll(/<td[^>]*>([\s\S]*?)<\/td>/g)].map((c) => cellText(c[1])))
    .filter((cells) => cells.length >= 4);
  const areas = new Map();
  const markers = new Map();
  const sections = new Map();
  const drugs = new Set();
  let folded = 0;
  for (const [drug, rawArea, marker, labelSections] of rows) {
    const area = AREA_FIXES[rawArea] ?? rawArea;
    if (area !== rawArea) folded += 1;
    drugs.add(drug.replace(/\s*\(\d+\)$/, ''));
    tally(areas, area);
    tally(markers, marker);
    for (const s of new Set(labelSections.split(',').map((x) => x.trim()).filter(Boolean))) tally(sections, s);
  }
  return {
    rows: rows.length,
    drugs: drugs.size,
    biomarkers: markers.size,
    updated: /Content current as of[\s\S]{0,400}?datetime="(\d{4}-\d{2}-\d{2})/.exec(html)?.[1] ?? null,
    areas: rankRows(areas, 'area'),
    folded,
    top: rankRows(markers, 'biomarker').slice(0, 20),
    sections: rankRows(sections, 'section'),
  };
};

export const clinpgxAnnotations = (text) => {
  const rows = tsvRecords(text);
  const levels = {};
  const categories = {};
  const genes = new Map();
  let pediatric = 0;
  for (const r of rows) {
    const level = r['Level of Evidence'];
    bump(levels, level);
    for (const c of new Set(list(r['Phenotype Category']))) {
      categories[c] ??= {};
      bump(categories[c], level);
    }
    for (const gene of new Set(list(r.Gene))) {
      const g = genes.get(gene) ?? { gene, total: 0, levels: {} };
      g.total += 1;
      bump(g.levels, level);
      genes.set(gene, g);
    }
    if (/pediatric/i.test(r['Specialty Population'] ?? '')) pediatric += 1;
  }
  return {
    annotations: rows.length,
    levels,
    categories,
    genes: [...genes.values()].sort((a, b) => b.total - a.total || a.gene.localeCompare(b.gene)).slice(0, 20),
    pediatric,
  };
};

export const clinpgxLabels = (text) => {
  const rows = tsvRecords(text);
  const sources = {};
  const genes = new Map();
  let cancer = 0;
  let onFdaList = 0;
  for (const r of rows) {
    const source = r.Source || 'Unknown';
    sources[source] ??= {};
    bump(sources[source], r['Testing Level'] || 'No level given');
    for (const gene of new Set(list(r.Genes))) tally(genes, gene);
    if (r['Cancer Genome']) cancer += 1;
    if (r['Biomarker Flag'] === 'On FDA Biomarker List') onFdaList += 1;
  }
  return {
    labels: rows.length,
    sources,
    genes: genes.size,
    top: [...genes].map(([gene, n]) => ({ gene, labels: n })).sort((a, b) => b.labels - a.labels || a.gene.localeCompare(b.gene)).slice(0, 15),
    cancer,
    onFdaList,
  };
};

export const createdDate = (names) => names.map((n) => /CREATED_(\d{4}-\d{2}-\d{2})/.exec(n)?.[1]).find(Boolean) ?? null;
