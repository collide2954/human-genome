import { summary } from '../src/stats.js';
import { tsvRecords, bump, round } from './tables.mjs';

const ROLES = new Set(['Act', 'LoF']);

export const latestIntogen = (html) => {
  const files = [...html.matchAll(/file=(IntOGen-(Drivers|Cohorts)-(\d{8})\.zip)/g)];
  const newest = (kind) => files.filter((m) => m[2] === kind).sort((a, b) => b[3].localeCompare(a[3]))[0];
  const drivers = newest('Drivers');
  const cohorts = newest('Cohorts');
  if (!drivers || !cohorts) return null;
  const d = drivers[3];
  return { release: `${d.slice(0, 4)}-${d.slice(4, 6)}-${d.slice(6, 8)}`, drivers: drivers[1], cohorts: cohorts[1] };
};

const majority = (counts) => {
  const ranked = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  return ranked.length > 1 && ranked[0][1] === ranked[1][1] ? 'ambiguous' : ranked[0][0];
};

export const intogenDrivers = (text) => {
  const rows = tsvRecords(text);
  const geneRoles = new Map();
  const pairRoles = new Map();
  const typeGenes = new Map();
  const cohorts = new Set();
  for (const r of rows) {
    const role = ROLES.has(r.ROLE) ? r.ROLE : 'ambiguous';
    const pair = `${r.SYMBOL}\t${r.CANCER_TYPE}`;
    cohorts.add(r.COHORT);
    if (!geneRoles.has(r.SYMBOL)) geneRoles.set(r.SYMBOL, {});
    bump(geneRoles.get(r.SYMBOL), role);
    if (!pairRoles.has(pair)) pairRoles.set(pair, {});
    bump(pairRoles.get(pair), role);
    if (!typeGenes.has(r.CANCER_TYPE)) typeGenes.set(r.CANCER_TYPE, new Set());
    typeGenes.get(r.CANCER_TYPE).add(r.SYMBOL);
  }
  const roles = {};
  for (const counts of geneRoles.values()) bump(roles, majority(counts));
  const perGene = new Map();
  for (const [pair, counts] of pairRoles) {
    const [gene] = pair.split('\t');
    if (!perGene.has(gene)) {
      perGene.set(gene, {
        gene, types: 0, Act: 0, LoF: 0, ambiguous: 0,
      });
    }
    const g = perGene.get(gene);
    g.types += 1;
    g[majority(counts)] += 1;
  }
  return {
    rows: rows.length,
    genes: geneRoles.size,
    cancerTypes: typeGenes.size,
    cohorts: cohorts.size,
    roles,
    byType: [...typeGenes].map(([type, set]) => ({ type, genes: set.size })).sort((a, b) => b.genes - a.genes || a.type.localeCompare(b.type)),
    top: [...perGene.values()].sort((a, b) => b.types - a.types || a.gene.localeCompare(b.gene)).slice(0, 25),
  };
};

export const intogenCohorts = (text) => {
  const rows = tsvRecords(text);
  const platforms = {};
  const kinds = {};
  const ages = {};
  const names = {};
  const byCancer = new Map();
  let samples = 0;
  let mutations = 0;
  for (const r of rows) {
    const n = Number(r.SAMPLES);
    const m = Number(r.MUTATIONS);
    samples += n;
    mutations += m;
    bump(platforms, r.PLATFORM);
    bump(kinds, r.TYPE);
    bump(ages, r.AGE);
    names[r.CANCER] ??= r.CANCER_NAME;
    if (!byCancer.has(r.CANCER)) byCancer.set(r.CANCER, []);
    byCancer.get(r.CANCER).push({
      cohort: r.COHORT, platform: r.PLATFORM, age: r.AGE, samples: n, rate: n ? m / n : 0,
    });
  }
  const burden = [...byCancer].map(([cancer, list]) => ({
    cancer,
    name: names[cancer],
    median: round(summary(list.map((c) => c.rate)).median, 1),
    cohorts: [...list].sort((a, b) => a.rate - b.rate).map((c) => ({
      cohort: c.cohort, platform: c.platform, age: c.age, samples: c.samples, perSample: round(c.rate, 1),
    })),
  })).sort((a, b) => a.median - b.median || a.cancer.localeCompare(b.cancer));
  return {
    cohorts: rows.length, samples, mutations, platforms, kinds, ages, names, burden,
  };
};

export const gdcRelease = (status) => {
  const v = status?.data_release_version;
  return { release: v ? `${v.major}.${v.minor}` : null, date: v?.release_date ?? null };
};

export const gdcPrograms = (json) => {
  const hits = json?.data?.hits ?? [];
  const programs = new Map();
  let cases = 0;
  for (const p of hits) {
    const n = p.summary?.case_count ?? 0;
    const name = p.program?.name ?? 'Other';
    const entry = programs.get(name) ?? { program: name, cases: 0, projects: 0 };
    entry.cases += n;
    entry.projects += 1;
    programs.set(name, entry);
    cases += n;
  }
  return {
    projects: hits.length,
    cases,
    programs: [...programs.values()].sort((a, b) => b.cases - a.cases || a.program.localeCompare(b.program)),
    tcga: hits.filter((p) => p.program?.name === 'TCGA')
      .map((p) => ({ project: p.project_id, name: p.name, cases: p.summary?.case_count ?? 0 }))
      .sort((a, b) => a.project.localeCompare(b.project)),
  };
};

export const facetCounts = (json, field) => (json?.data?.aggregations?.[field]?.buckets ?? [])
  .map((b) => ({ key: b.key, n: b.doc_count }))
  .sort((a, b) => b.n - a.n || String(a.key).localeCompare(String(b.key)));

export const casesByProject = (hits) => {
  const sets = new Map();
  for (const h of hits) {
    const project = h.case?.project?.project_id;
    const id = h.case?.case_id;
    if (!project || !id) continue;
    if (!sets.has(project)) sets.set(project, new Set());
    sets.get(project).add(id);
  }
  return Object.fromEntries([...sets].map(([project, set]) => [project, set.size]));
};
