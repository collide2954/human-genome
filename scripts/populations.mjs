import { logEdges, histogram, summary } from '../src/stats.js';

const own = (s) => Buffer.from(s, 'utf8').toString('utf8');
const bump = (obj, key, by = 1) => {
  obj[key] = (obj[key] ?? 0) + by;
};
const ranked = (map, key) => [...map].map(([k, n]) => ({ [key]: k, n })).sort((a, b) => b.n - a.n || String(a[key]).localeCompare(String(b[key])));

export const gwasAncestry = (text) => {
  const lines = text.split('\n');
  const head = lines[0].split('\t');
  const at = (re) => head.findIndex((h) => re.test(h));
  const col = {
    acc: at(/^STUDY ACCESSION/), date: at(/^DATE/), stage: at(/^STAGE/), n: at(/^NUMBER OF INDIV/), cat: at(/^BROAD ANCESTRAL/), country: at(/^COUNTRY OF RECRUITMENT/),
  };
  let rows = 0;
  let participants = 0;
  const studies = new Map();
  const byYear = {};
  const stages = { initial: {}, replication: {} };
  const countries = new Map();
  for (const line of lines.slice(1)) {
    if (!line) continue;
    const c = line.split('\t');
    rows += 1;
    const year = (c[col.date] ?? '').slice(0, 4);
    if (!studies.has(c[col.acc])) studies.set(own(c[col.acc]), year);
    const n = Number(c[col.n]);
    if (!(n > 0)) continue;
    participants += n;
    const raw = (c[col.cat] ?? '').trim();
    const category = raw.includes(',') ? 'Multiple categories' : raw || 'NR';
    const stage = c[col.stage] === 'replication' ? 'replication' : 'initial';
    bump(stages[stage], category, n);
    if (stage === 'initial' && year) {
      byYear[year] ??= {};
      bump(byYear[year], category, n);
    }
    for (const country of new Set((c[col.country] ?? '').split(',').map((x) => x.trim()).filter((x) => x && x !== 'NR'))) {
      countries.set(countries.has(country) ? country : own(country), (countries.get(country) ?? 0) + n);
    }
  }
  const studiesByYear = {};
  for (const year of studies.values()) if (year) bump(studiesByYear, year);
  return {
    rows,
    studies: studies.size,
    participants,
    byYear,
    studiesByYear,
    stages,
    countries: ranked(countries, 'country').slice(0, 20).map(({ country, n }) => ({ country, participants: n })),
  };
};

export const majorHaplogroup = (call, kind) => {
  const s = String(call ?? '').trim();
  if (!/^[A-Z]/.test(s)) return null;
  if (kind === 'mt') return /^(L\d|HV|JT)/.exec(s)?.[1] ?? s[0];
  return s[0];
};

const AGE_EDGES = logEdges(2, 6, 4);
const SNP_EDGES = logEdges(3, 7, 4);
const column = (head, prefix) => head.findIndex((h) => h.replace(/^"/, '').startsWith(prefix));

export const aadrSummary = (text) => {
  const lines = text.split('\n');
  const head = lines[0].split('\t');
  const col = {
    id: column(head, 'Individual ID'),
    pub: column(head, 'First publication'),
    date: column(head, 'Date mean in BP'),
    country: column(head, 'Political Entity'),
    lat: column(head, 'Latitude'),
    lon: column(head, 'Longitude'),
    type: column(head, 'Data type'),
    snps: column(head, 'SNPs hit on autosomal targets (Computed using easystats on 1240k'),
    sex: column(head, 'Molecular Sex'),
    y: column(head, 'Y haplogroup  in ISOGG'),
    mt: column(head, 'mtDNA haplogroup if'),
    assessment: head.findIndex((h) => h.trim() === 'ASSESSMENT'),
  };
  const seen = new Set();
  let present = 0;
  const ages = [];
  const snps = [];
  let oldest = null;
  const published = {};
  const countries = new Map();
  const dataTypes = {};
  const sex = {};
  const y = new Map();
  const mt = new Map();
  const assessment = {};
  const cells = new Map();
  for (const line of lines.slice(1)) {
    if (!line) continue;
    const c = line.split('\t');
    const id = c[col.id];
    if (!id || seen.has(id)) continue;
    seen.add(own(id));
    const age = Number(c[col.date]);
    if (!(age > 0)) {
      present += 1;
      continue;
    }
    ages.push(age);
    const place = c[col.country] || 'not stated';
    if (!oldest || age > oldest.age) oldest = { id: own(id), age, place: own(place) };
    const year = /((?:19|20)\d{2})\D*$/.exec(c[col.pub] ?? '')?.[1];
    if (year) bump(published, year);
    countries.set(countries.has(place) ? place : own(place), (countries.get(place) ?? 0) + 1);
    bump(dataTypes, c[col.type] || 'not stated');
    bump(sex, (c[col.sex] || 'U')[0]);
    const yg = majorHaplogroup(c[col.y], 'y');
    if (yg) y.set(yg, (y.get(yg) ?? 0) + 1);
    const mg = majorHaplogroup(c[col.mt], 'mt');
    if (mg) mt.set(mg, (mt.get(mg) ?? 0) + 1);
    if (c[col.assessment]) bump(assessment, c[col.assessment]);
    const hit = Number(c[col.snps]);
    if (hit > 0) snps.push(hit);
    const lat = Number(c[col.lat]);
    const lon = Number(c[col.lon]);
    if (c[col.lat] && c[col.lon] && Number.isFinite(lat) && Number.isFinite(lon)) {
      const key = `${Math.floor(lat / 2) * 2},${Math.floor(lon / 2) * 2}`;
      if (!cells.has(key)) cells.set(key, { lat: Math.floor(lat / 2) * 2, lon: Math.floor(lon / 2) * 2, ages: [] });
      cells.get(key).ages.push(age);
    }
  }
  return {
    individuals: ages.length,
    present,
    ages: { edges: AGE_EDGES, counts: histogram(ages, AGE_EDGES), median: summary(ages)?.median ?? null },
    oldest,
    published,
    countries: ranked(countries, 'country').slice(0, 20),
    dataTypes,
    sex,
    y: ranked(y, 'group').slice(0, 15),
    mt: ranked(mt, 'group').slice(0, 15),
    assessment,
    snps: { edges: SNP_EDGES, counts: histogram(snps, SNP_EDGES), median: summary(snps)?.median ?? null },
    grid: [...cells.values()].sort((a, b) => a.lat - b.lat || a.lon - b.lon).map((g) => ({
      lat: g.lat, lon: g.lon, n: g.ages.length, age: summary(g.ages).median,
    })),
  };
};

export const aadrAnnoFile = (dataset) => {
  const latest = dataset?.data?.latestVersion ?? {};
  const file = (latest.files ?? []).map((f) => f.dataFile).find((f) => /_1240K\.aadr\.PUB\.anno$/.test(f?.filename ?? ''));
  if (!file) return null;
  return { id: file.id, version: /^(v[\d.]+(?:\.p\d+)?)_/.exec(file.filename)?.[1] ?? null, released: latest.releaseTime?.slice(0, 10) ?? null };
};

export const igsrSamples = (json, collection) => {
  const hits = json?.hits?.hits ?? [];
  const collections = new Map();
  const populations = new Map();
  const sex = {};
  let samples = 0;
  for (const { _source: s } of hits) {
    const titles = new Set((s.dataCollections ?? []).map((d) => d.title));
    for (const t of titles) collections.set(t, (collections.get(t) ?? 0) + 1);
    if (!titles.has(collection)) continue;
    samples += 1;
    bump(sex, s.sex ?? 'unknown');
    const p = s.populations?.[0];
    if (!p) continue;
    const e = populations.get(p.code) ?? {
      code: p.code, name: p.name, superpopulation: p.superpopulationCode, superpopulationName: p.superpopulationName, n: 0,
    };
    e.n += 1;
    populations.set(p.code, e);
  }
  const total = json?.hits?.total;
  return {
    total: typeof total === 'object' ? total.value : total ?? hits.length,
    collections: ranked(collections, 'title'),
    panel: { samples, populations: [...populations.values()].sort((a, b) => b.n - a.n || a.code.localeCompare(b.code)), sex },
  };
};

export const assemblySummary = (reports) => {
  const byYear = {};
  const levels = {};
  const submitters = new Map();
  const n50 = [];
  let best = null;
  for (const r of reports) {
    const info = r.assembly_info ?? {};
    const level = info.assembly_level ?? 'unknown';
    const date = info.release_date ?? '';
    const year = date.slice(0, 4);
    if (year) {
      byYear[year] ??= {};
      bump(byYear[year], level);
    }
    bump(levels, level);
    if (info.submitter) submitters.set(info.submitter, (submitters.get(info.submitter) ?? 0) + 1);
    const contig = Number(r.assembly_stats?.contig_n50);
    if (contig > 0 && date) {
      n50.push([date, contig, level, r.accession]);
      if (!best || contig > best.n50) {
        best = {
          accession: r.accession, name: info.assembly_name, n50: contig, submitter: info.submitter,
        };
      }
    }
  }
  return {
    total: reports.length, byYear, levels, submitters: ranked(submitters, 'submitter').slice(0, 12), n50, best,
  };
};
