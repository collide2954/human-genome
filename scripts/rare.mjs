import { histogram, summary } from '../src/stats.js';

const own = (s) => Buffer.from(s, 'utf8').toString('utf8');
const bump = (obj, key, by = 1) => {
  obj[key] = (obj[key] ?? 0) + by;
};
const tally = (map, key, by = 1) => map.set(map.has(key) ? key : own(key), (map.get(key) ?? 0) + by);
const ENTITIES = {
  '&lt;': '<', '&gt;': '>', '&amp;': '&', '&apos;': "'", '&quot;': '"', '&nbsp;': ' ',
};
const decode = (s) => s.replace(/&(?:lt|gt|amp|apos|quot|nbsp);/g, (m) => ENTITIES[m]).replace(/&#(\d+);/g, (_, n) => String.fromCharCode(Number(n)));
const byCount = (key) => (a, b) => b.n - a.n || String(a[key]).localeCompare(String(b[key]));

export const hpoTerms = (text) => {
  const terms = new Map();
  for (const stanza of text.split(/\n(?=\[)/)) {
    if (!stanza.startsWith('[Term]') || /^is_obsolete: true/m.test(stanza)) continue;
    const id = /^id: (\S+)/m.exec(stanza)?.[1];
    if (!id) continue;
    terms.set(own(id), { name: own(/^name: (.+)$/m.exec(stanza)?.[1] ?? ''), parents: [...stanza.matchAll(/^is_a: (\S+)/gm)].map((m) => own(m[1])) });
  }
  return {
    version: /^data-version:\s*(?:hp\/releases\/)?(\S+)/m.exec(text)?.[1] ?? null,
    terms,
    systems: [...terms].filter(([, t]) => t.parents.includes('HP:0000118')).map(([id]) => id),
  };
};

const PER_DISEASE = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];

export const hpoAnnotations = (text, { terms, systems }) => {
  const roots = new Set(systems);
  const memo = new Map();
  const systemsOf = (id) => {
    if (memo.has(id)) return memo.get(id);
    const out = new Set();
    memo.set(id, out);
    if (roots.has(id)) out.add(id);
    for (const p of terms.get(id)?.parents ?? []) for (const s of systemsOf(p)) out.add(s);
    return out;
  };
  let version = null;
  let col = null;
  let annotations = 0;
  const sources = {};
  const evidence = {};
  const diseases = new Map();
  const byTerm = new Map();
  const byInheritance = new Map();
  const add = (map, term, disease) => {
    if (!map.has(term)) map.set(own(term), new Set());
    map.get(term).add(disease);
  };
  for (const line of text.split('\n')) {
    if (!line) continue;
    if (line.startsWith('#')) {
      version = /^#version:\s*(\S+)/.exec(line)?.[1] ?? version;
      continue;
    }
    const c = line.split('\t');
    if (!col) {
      col = Object.fromEntries(c.map((h, i) => [h, i]));
      continue;
    }
    if (c[col.qualifier] === 'NOT') continue;
    const id = c[col.database_id];
    if (!diseases.has(id)) {
      diseases.set(own(id), { terms: new Set(), systems: new Set() });
      bump(sources, id.split(':')[0]);
    }
    const d = diseases.get(id);
    const term = c[col.hpo_id];
    if (c[col.aspect] === 'P') {
      annotations += 1;
      bump(evidence, c[col.evidence]);
      d.terms.add(term);
      for (const s of systemsOf(term)) d.systems.add(s);
      add(byTerm, term, id);
    } else if (c[col.aspect] === 'I') add(byInheritance, term, id);
  }
  const counts = [...diseases.values()].map((d) => d.terms.size).filter((n) => n > 0);
  const inSystem = new Map();
  for (const d of diseases.values()) for (const s of d.systems) inSystem.set(s, (inSystem.get(s) ?? 0) + 1);
  const ranked = (entries) => entries.map(([id, n]) => ({ id, name: terms.get(id)?.name ?? id, diseases: n }))
    .sort((a, b) => b.diseases - a.diseases || a.name.localeCompare(b.name));
  return {
    version,
    diseases: diseases.size,
    sources,
    annotations,
    evidence,
    perDisease: { edges: PER_DISEASE, counts: histogram(counts, PER_DISEASE), median: summary(counts).median },
    systems: ranked([...inSystem]),
    topTerms: ranked([...byTerm].map(([id, set]) => [id, set.size])).slice(0, 20),
    inheritance: ranked([...byInheritance].map(([id, set]) => [id, set.size])),
  };
};

const capIndex = (caps, n) => {
  const i = caps.findIndex((cap) => n <= cap);
  return i < 0 ? caps.length : i;
};
const BINS = ['1', '2', '3-5', '6-10', '>10'];
const binOf = (n) => capIndex([1, 2, 5, 10], n);

export const geneDiseases = (text) => {
  const lines = text.split('\n').filter(Boolean);
  const col = Object.fromEntries(lines[0].split('\t').map((h, i) => [h, i]));
  const types = {};
  const perGene = new Map();
  const diseases = new Set();
  for (const line of lines.slice(1)) {
    const c = line.split('\t');
    bump(types, c[col.association_type]);
    const gene = c[col.gene_symbol];
    if (c[col.association_type] !== 'MENDELIAN' || !gene || gene === '-') continue;
    if (!perGene.has(gene)) perGene.set(own(gene), new Set());
    perGene.get(gene).add(c[col.disease_id]);
    diseases.add(c[col.disease_id]);
  }
  const counts = new Array(BINS.length).fill(0);
  for (const set of perGene.values()) counts[binOf(set.size)] += 1;
  return {
    types,
    genes: perGene.size,
    diseases: diseases.size,
    perGene: { bins: BINS, counts },
    top: [...perGene].map(([gene, set]) => ({ gene, n: set.size })).sort(byCount('gene')).slice(0, 15).map(({ gene, n }) => ({ gene, diseases: n })),
  };
};

const blocks = (xml) => xml.split('<Disorder id="').slice(1);
const NAMED = (tag) => new RegExp(`<${tag} id="\\d+">\\s*<Name lang="en">([^<]*)</Name>`, 'g');
const namesIn = (block, tag) => [...block.matchAll(NAMED(tag))].map((m) => decode(m[1]));
const nameIn = (block, tag) => namesIn(block, tag)[0] ?? null;
const dateOf = (xml) => /<JDBOR date="(\d{4}-\d{2}-\d{2})/.exec(xml)?.[1] ?? null;
const isDisorder = (block) => nameIn(block, 'DisorderGroup') === 'Disorder';

export const orphaGenes = (xml) => {
  const genes = new Set();
  const associations = {};
  const counts = new Array(BINS.length).fill(0);
  let disorders = 0;
  for (const block of blocks(xml)) {
    const symbols = new Set();
    for (const a of block.split('<DisorderGeneAssociation>').slice(1)) {
      const symbol = /<Symbol>([^<]+)<\/Symbol>/.exec(a)?.[1];
      if (symbol) {
        symbols.add(symbol);
        genes.add(symbol);
      }
      const type = nameIn(a, 'DisorderGeneAssociationType');
      if (type) bump(associations, type);
    }
    if (!isDisorder(block) || !symbols.size) continue;
    disorders += 1;
    counts[binOf(symbols.size)] += 1;
  }
  return {
    date: dateOf(xml), disorders, genes: genes.size, associations, perDisorder: { bins: BINS, counts },
  };
};

export const orphaPrevalence = (xml) => {
  const classes = {};
  let disorders = 0;
  for (const block of blocks(xml)) {
    if (!isDisorder(block)) continue;
    disorders += 1;
    const list = block.split('<Prevalence id="').slice(1).map((p) => ({
      type: nameIn(p, 'PrevalenceType'), cls: nameIn(p, 'PrevalenceClass'), world: nameIn(p, 'PrevalenceGeographic') === 'Worldwide', valid: nameIn(p, 'PrevalenceValidationStatus') === 'Validated',
    }));
    const point = list.filter((p) => p.type === 'Point prevalence' && p.cls).sort((a, b) => b.world - a.world || b.valid - a.valid);
    if (point.length) bump(classes, point[0].cls);
    else if (list.length && list.every((p) => p.type === 'Cases/families')) bump(classes, 'Case reports only');
    else bump(classes, 'No point prevalence');
  }
  return { date: dateOf(xml), disorders, classes };
};

export const orphaAges = (xml) => {
  const onset = {};
  const inheritance = {};
  let disorders = 0;
  let noOnset = 0;
  let noInheritance = 0;
  for (const block of blocks(xml)) {
    if (!isDisorder(block)) continue;
    disorders += 1;
    const ages = namesIn(block, 'AverageAgeOfOnset');
    const types = namesIn(block, 'TypeOfInheritance');
    ages.forEach((a) => bump(onset, a));
    types.forEach((t) => bump(inheritance, t));
    if (!ages.length) noOnset += 1;
    if (!types.length) noInheritance += 1;
  }
  return {
    date: dateOf(xml), disorders, onset, inheritance, noOnset, noInheritance,
  };
};

const unquote = (s) => (s ?? '').replace(/^"|"$/g, '');
const GENE_BINS = ['0', '1', '2-10', '11-100', '101-1000', '>1000'];
const geneBin = (n) => capIndex([0, 1, 10, 100, 1000], n);

export const gtrTests = (text) => {
  const lines = text.split('\n');
  const col = Object.fromEntries(lines[0].split('\t').map((h, i) => [h, i]));
  let tests = 0;
  const labs = new Set();
  const countries = new Map();
  const methods = new Map();
  const genes = new Map();
  const development = {};
  const created = {};
  const counts = new Array(GENE_BINS.length).fill(0);
  for (const line of lines.slice(1)) {
    if (!line) continue;
    const c = line.split('\t');
    if (c[col.now_current] !== '1') continue;
    tests += 1;
    const lab = c[col.lab_id];
    labs.add(lab);
    const country = unquote(c[col.facility_country]) || 'not stated';
    if (!countries.has(country)) countries.set(own(country), { tests: 0, labs: new Set() });
    const e = countries.get(country);
    e.tests += 1;
    e.labs.add(lab);
    for (const m of unquote(c[col.method_categories]).split('|').filter(Boolean)) tally(methods, m);
    const list = unquote(c[col.genes]).split('|').filter(Boolean);
    counts[geneBin(list.length)] += 1;
    for (const g of list) tally(genes, g);
    bump(development, unquote(c[col.test_development]) || 'not stated');
    const year = /^(\d{4})/.exec(c[col.test_create_date] ?? '')?.[1];
    if (year) bump(created, year);
  }
  return {
    tests,
    labs: labs.size,
    countries: [...countries].map(([country, e]) => ({ country, tests: e.tests, labs: e.labs.size })).sort((a, b) => b.tests - a.tests || a.country.localeCompare(b.country)).slice(0, 15),
    methods: [...methods].map(([method, n]) => ({ method, n })).sort(byCount('method')).map(({ method, n }) => ({ method, tests: n })),
    genesPerTest: { bins: GENE_BINS, counts },
    topGenes: [...genes].map(([gene, n]) => ({ gene, n })).sort(byCount('gene')).slice(0, 20).map(({ gene, n }) => ({ gene, tests: n })),
    development,
    created,
  };
};

const SKIP = new Set(['cellular-gene-therapy-products', 'approved-cellular-and-gene-therapy-products']);

export const fdaProducts = (html) => {
  const start = html.indexOf('<main');
  const end = html.indexOf('</main>');
  const main = html.slice(Math.max(0, start), end > start ? end : html.length);
  const out = [];
  for (const m of main.matchAll(/<a href="([^"]+)"[^>]*>([^<]+)<\/a>/g)) {
    const path = m[1].replace(/^https:\/\/www\.fda\.gov/, '');
    if (!path.startsWith('/vaccines-blood-biologics/') || SKIP.has(path.replace(/\/$/, '').split('/').at(-1))) continue;
    const name = decode(m[2]).replace(/\s+/g, ' ').trim();
    out.push({ trade: /^([A-Z][A-Z0-9-]+)[™®]?(?=\s*\(|$)/.exec(name)?.[1] ?? name, name, url: `https://www.fda.gov${path}` });
  }
  return out;
};

const PHASES = ['NA', 'EARLY_PHASE1', 'PHASE1', 'PHASE2', 'PHASE3', 'PHASE4'];
const PHASE_LABEL = {
  NA: 'Not applicable', EARLY_PHASE1: 'Early phase 1', PHASE1: 'Phase 1', PHASE2: 'Phase 2', PHASE3: 'Phase 3', PHASE4: 'Phase 4',
};

export const trialYears = (studies) => {
  const byYear = {};
  const phases = {};
  let undated = 0;
  let recruiting = 0;
  for (const s of studies) {
    const p = s.protocolSection ?? {};
    const year = /^(\d{4})/.exec(p.statusModule?.startDateStruct?.date ?? '')?.[1];
    if (year) bump(byYear, year);
    else undated += 1;
    const top = (p.designModule?.phases ?? []).reduce((best, x) => (PHASES.indexOf(x) > PHASES.indexOf(best) ? x : best), null);
    bump(phases, top ? PHASE_LABEL[top] ?? top : 'Not stated');
    if (p.statusModule?.overallStatus === 'RECRUITING') recruiting += 1;
  }
  return {
    total: studies.length, byYear, undated, phases, recruiting,
  };
};
