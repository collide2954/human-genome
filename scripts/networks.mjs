import {
  logEdges, linearEdges, histogram, binIndex, summary, clusterOrder,
} from '../src/stats.js';

const own = (s) => Buffer.from(s, 'utf8').toString('utf8');
const bump = (map, key, by = 1) => {
  const v = map.get(key);
  map.set(v === undefined ? own(key) : key, (v ?? 0) + by);
};
const interner = () => {
  const seen = new Map();
  return (s) => {
    let v = seen.get(s);
    if (v === undefined) {
      v = own(s);
      seen.set(v, v);
    }
    return v;
  };
};
const ranked = (map, name) => [...map].map(([k, n]) => ({ [name]: k, n })).sort((a, b) => b.n - a.n || String(a[name]).localeCompare(String(b[name])));

const confidence = () => ({
  all: 0, medium: 0, high: 0, highest: 0,
});
const tally = (edges, score) => {
  edges.all += 1;
  if (score >= 400) edges.medium += 1;
  if (score >= 700) edges.high += 1;
  if (score >= 900) edges.highest += 1;
};

const CHANNELS = ['neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database', 'textmining'];
const SCORE_EDGES = linearEdges(150, 1000, 50);
const DEGREE_EDGES = logEdges(0, 4, 5);

export const stringNames = (text) => {
  const out = {};
  for (const line of text.split('\n')) {
    if (!line || line.startsWith('#')) continue;
    const [id, name] = line.split('\t');
    out[own(id)] = own(name);
  }
  return out;
};

export const stringSummary = ({ names = {} } = {}) => {
  const edges = confidence();
  const proteins = new Set();
  const degree = new Map();
  const scores = new Array(SCORE_EDGES.length - 1).fill(0);
  const channels = Object.fromEntries(CHANNELS.map((c) => [c, 0]));
  const only = Object.fromEntries(CHANNELS.map((c) => [c, 0]));
  const note = (p) => {
    if (!proteins.has(p)) proteins.add(own(p));
  };
  return {
    add(line) {
      if (!line || line.startsWith('protein1')) return;
      const a = line.indexOf(' ');
      const b = line.indexOf(' ', a + 1);
      const p1 = line.slice(0, a);
      const p2 = line.slice(a + 1, b);
      if (!(p1 < p2)) return;
      const f = line.slice(b + 1).split(' ').map(Number);
      const score = f[7];
      tally(edges, score);
      scores[binIndex(SCORE_EDGES, score)] += 1;
      note(p1);
      note(p2);
      if (score < 700) return;
      bump(degree, p1);
      bump(degree, p2);
      const present = CHANNELS.filter((_, i) => f[i] > 0);
      present.forEach((c) => { channels[c] += 1; });
      if (present.length === 1) only[present[0]] += 1;
    },
    result() {
      const degrees = [...degree.values()];
      return {
        edges,
        proteins: { all: proteins.size, high: degree.size },
        scores: { edges: SCORE_EDGES, counts: scores },
        degree: { edges: DEGREE_EDGES, counts: histogram(degrees, DEGREE_EDGES), summary: summary(degrees) },
        hubs: [...degree].sort((x, y) => y[1] - x[1] || x[0].localeCompare(y[0])).slice(0, 20)
          .map(([id, d]) => ({ id, name: names[id] ?? id.replace(/^9606\./, ''), degree: d })),
        channels,
        only,
      };
    },
  };
};

export const physicalSummary = () => {
  const edges = confidence();
  const all = new Set();
  const high = new Set();
  const note = (set, p) => {
    if (!set.has(p)) set.add(own(p));
  };
  return {
    add(line) {
      if (!line || line.startsWith('protein1')) return;
      const [p1, p2, s] = line.split(' ');
      if (!(p1 < p2)) return;
      const score = Number(s);
      tally(edges, score);
      note(all, p1);
      note(all, p2);
      if (score >= 700) {
        note(high, p1);
        note(high, p2);
      }
    },
    result: () => ({ edges, proteins: { all: all.size, high: high.size } }),
  };
};

const tidy = (v) => Number(v.toFixed(3));

export const reactomeSummary = ({ pathways, relations }) => {
  const names = new Map();
  for (const line of pathways.split('\n')) {
    const [id, name, species] = line.split('\t');
    if (species?.trim() === 'Homo sapiens') names.set(own(id), own(name));
  }
  const children = new Map();
  const hasParent = new Set();
  for (const line of relations.split('\n')) {
    const [parent, child] = line.split('\t').map((s) => s?.trim());
    if (!names.has(parent) || !names.has(child)) continue;
    if (!children.has(parent)) children.set(parent, []);
    children.get(parent).push(child);
    hasParent.add(child);
  }
  const tops = [...names.keys()].filter((id) => !hasParent.has(id));
  const depthOf = new Map(tops.map((t) => [t, 0]));
  const queue = [...tops];
  for (let i = 0; i < queue.length; i += 1) {
    for (const child of children.get(queue[i]) ?? []) {
      if (!depthOf.has(child)) {
        depthOf.set(child, depthOf.get(queue[i]) + 1);
        queue.push(child);
      }
    }
  }
  const topsOf = new Map();
  const below = new Map();
  for (const top of tops) {
    const seen = new Set([top]);
    const stack = [top];
    while (stack.length) {
      const id = stack.pop();
      if (!topsOf.has(id)) topsOf.set(id, new Set());
      topsOf.get(id).add(top);
      for (const child of children.get(id) ?? []) {
        if (!seen.has(child)) {
          seen.add(child);
          stack.push(child);
        }
      }
    }
    below.set(top, seen.size);
  }
  const members = new Map();
  const evidence = {};
  const canon = interner();
  return {
    add(line) {
      if (!line) return;
      const f = line.split('\t');
      if (f[5]?.trim() !== 'Homo sapiens' || !names.has(f[1])) return;
      evidence[f[4]] = (evidence[f[4]] ?? 0) + 1;
      if (!members.has(f[1])) members.set(f[1], new Set());
      members.get(f[1]).add(canon(f[0].split('-')[0]));
    },
    result() {
      const proteins = new Set();
      const topMembers = new Map(tops.map((t) => [t, new Set()]));
      for (const [id, set] of members) {
        for (const a of set) proteins.add(a);
        for (const t of topsOf.get(id) ?? []) for (const a of set) topMembers.get(t).add(a);
      }
      const topLevel = tops.map((id) => ({
        id, name: names.get(id), proteins: topMembers.get(id).size, pathways: below.get(id),
      })).sort((a, b) => b.proteins - a.proteins || a.name.localeCompare(b.name));
      const breadthOf = new Map();
      for (const set of topMembers.values()) for (const a of set) breadthOf.set(a, (breadthOf.get(a) ?? 0) + 1);
      const breadth = new Array(tops.length + 1).fill(0);
      for (const n of breadthOf.values()) breadth[n] += 1;
      const depth = [];
      for (const d of depthOf.values()) depth[d] = (depth[d] ?? 0) + 1;
      const sized = [...members].map(([id, set]) => ({ id, name: names.get(id), proteins: set.size }));
      const topIds = new Set(tops);
      const sets = topLevel.map((t) => topMembers.get(t.id));
      const matrix = sets.map((a) => sets.map((b) => {
        let both = 0;
        for (const x of a) if (b.has(x)) both += 1;
        const either = a.size + b.size - both;
        return either ? tidy(both / either) : 0;
      }));
      return {
        pathways: names.size,
        proteins: proteins.size,
        topLevel,
        depth: Array.from(depth, (n) => n ?? 0),
        breadth,
        sizes: { edges: DEGREE_EDGES, counts: histogram(sized.map((s) => s.proteins), DEGREE_EDGES) },
        largest: sized.filter((s) => !topIds.has(s.id)).sort((a, b) => b.proteins - a.proteins || a.name.localeCompare(b.name)).slice(0, 15),
        overlap: { ids: topLevel.map((t) => t.id), matrix, order: clusterOrder(matrix) },
        evidence,
      };
    },
  };
};

const NAMESPACES = ['biological_process', 'molecular_function', 'cellular_component'];

export const oboTerms = (text) => {
  const terms = new Map();
  const counts = Object.fromEntries(NAMESPACES.map((n) => [n, 0]));
  let obsolete = 0;
  for (const stanza of text.split(/\n(?=\[)/)) {
    if (!stanza.startsWith('[Term]')) continue;
    if (/^is_obsolete: true/m.test(stanza)) {
      obsolete += 1;
      continue;
    }
    const id = /^id: (\S+)/m.exec(stanza)?.[1];
    if (!id) continue;
    const namespace = NAMESPACES.find((n) => stanza.includes(`\nnamespace: ${n}`)) ?? null;
    terms.set(own(id), { name: own(/^name: (.+)$/m.exec(stanza)?.[1] ?? ''), namespace });
    if (namespace) counts[namespace] += 1;
  }
  return {
    version: /^data-version:\s*(\S+)/m.exec(text)?.[1] ?? null, counts, obsolete, terms,
  };
};

const EVIDENCE = {
  experimental: ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'],
  phylogenetic: ['IBA', 'IBD', 'IKR', 'IRD'],
  computational: ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA'],
  author: ['TAS', 'NAS', 'IC', 'ND'],
  electronic: ['IEA'],
};
const GROUP_OF = Object.fromEntries(Object.entries(EVIDENCE).flatMap(([g, codes]) => codes.map((c) => [c, g])));
const blank = () => Object.fromEntries(Object.keys(EVIDENCE).map((g) => [g, 0]));
const ASPECTS = ['P', 'F', 'C'];

export const goSummary = (terms = new Map()) => {
  let annotations = 0;
  let not = 0;
  const byAspect = Object.fromEntries(ASPECTS.map((a) => [a, blank()]));
  const perProtein = new Map();
  const perTerm = new Map();
  const assigned = new Map();
  const codes = new Map();
  const objects = new Map();
  const canon = interner();
  return {
    add(line) {
      if (!line || line.startsWith('!')) return;
      const f = line.split('\t');
      annotations += 1;
      bump(assigned, f[14]);
      if (f[3].startsWith('NOT')) {
        not += 1;
        return;
      }
      if (GROUP_OF[f[6]] && byAspect[f[8]]) byAspect[f[8]][GROUP_OF[f[6]]] += 1;
      bump(codes, f[6]);
      const type = canon(f[11] || 'unknown');
      if (!objects.has(type)) objects.set(type, new Set());
      objects.get(type).add(canon(`${f[0]}:${f[1]}`));
      if (f[0] !== 'UniProtKB' || type !== 'protein') return;
      const acc = canon(f[1]);
      const term = canon(f[4]);
      if (!perProtein.has(acc)) perProtein.set(acc, { symbol: own(f[2]), terms: new Set() });
      perProtein.get(acc).terms.add(term);
      if (!perTerm.has(term)) perTerm.set(term, { aspect: f[8] === 'P' ? 'P' : f[8] === 'F' ? 'F' : 'C', proteins: new Set() });
      perTerm.get(term).proteins.add(acc);
    },
    result() {
      const sizes = [...perProtein.values()].map((p) => p.terms.size);
      const topTerms = Object.fromEntries(ASPECTS.map((a) => [a, [...perTerm]
        .filter(([, t]) => t.aspect === a)
        .map(([id, t]) => ({ id, name: terms.get(id)?.name ?? null, n: t.proteins.size }))
        .sort((x, y) => y.n - x.n || x.id.localeCompare(y.id))
        .slice(0, 12)]));
      return {
        annotations,
        not,
        proteins: perProtein.size,
        types: Object.fromEntries([...objects].map(([t, set]) => [t, set.size]).sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))),
        terms: perTerm.size,
        byAspect,
        codes: ranked(codes, 'code'),
        perProtein: { edges: DEGREE_EDGES, counts: histogram(sizes, DEGREE_EDGES), summary: summary(sizes) },
        topProteins: [...perProtein].map(([acc, p]) => ({ acc, symbol: p.symbol, terms: p.terms.size }))
          .sort((a, b) => b.terms - a.terms || a.symbol.localeCompare(b.symbol)).slice(0, 15),
        topTerms,
        assignedBy: ranked(assigned, 'group').slice(0, 15),
      };
    },
  };
};
