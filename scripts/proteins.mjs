import {
  logEdges, histogram, summary, coverage, linearEdges, clusterOrder,
} from '../src/stats.js';

const own = (s) => (typeof s === 'string' ? Buffer.from(s, 'utf8').toString('utf8') : s);
const capped = (cap) => new Array(cap + 1).fill(0);
const bump = (map, key, by = 1) => map.set(key, (map.get(key) ?? 0) + by);
const ranked = (map, name) => [...map].map(([k, n]) => ({ [name]: k, n })).sort((a, b) => b.n - a.n || String(a[name]).localeCompare(String(b[name])));

const TOPOLOGY = /membrane protein|anchor|side$/i;

export const locationsOf = (text) => {
  if (!text) return [];
  const body = text.split('Note=')[0]
    .replace(/SUBCELLULAR LOCATION:\s*/g, '. ')
    .replace(/\{[^}]*\}/g, '')
    .replace(/\[[^\]]*\]:/g, '');
  const out = [];
  for (const seg of body.split(/[.;]/).map((s) => s.trim()).filter(Boolean)) {
    if (TOPOLOGY.test(seg)) continue;
    const main = seg.split(',')[0].trim();
    if (main && !out.includes(main)) out.push(main);
  }
  return out;
};

const notes = (text, feature) => [...String(text ?? '').matchAll(new RegExp(`${feature} [^;]+; /note="([^"]+)"`, 'g'))].map((m) => m[1]);

const disordered = (text) => [...String(text ?? '').matchAll(/REGION (\d+)\.\.(\d+); \/note="Disordered"/g)].map((m) => [Number(m[1]) - 1, Number(m[2])]);

export const uniprotSummary = () => {
  let col = null;
  const proteins = [];
  const existence = new Map();
  const pdbCounts = capped(50);
  const pdbTop = [];
  const tmCounts = capped(25);
  const locations = new Map();
  const pfam = new Map();
  const perProtein = capped(15);
  const share = new Array(10).fill(0);
  const ptm = {
    phosphoS: 0, phosphoT: 0, phosphoY: 0, acetyl: 0, methyl: 0, nGlyco: 0, oGlyco: 0, ubiquitin: 0, sumo: 0,
  };
  let signal = 0;
  let tmProteins = 0;
  let pdbAny = 0;
  let pfamAny = 0;
  let disorderResidues = 0;
  let disorderProteins = 0;
  return {
    add(line) {
      if (!line) return;
      const c = line.split('\t');
      if (!col) {
        col = Object.fromEntries(c.map((name, i) => [name, i]));
        return;
      }
      const acc = own(c[col.Entry]);
      const gene = own(c[col['Gene Names (primary)']] ?? '');
      const length = Number(c[col.Length]);
      proteins.push({ acc, gene, length });
      bump(existence, own(c[col['Protein existence']]));
      const pdb = (c[col.PDB] ?? '').split(';').filter(Boolean).length;
      pdbCounts[Math.min(pdb, 50)] += 1;
      if (pdb) {
        pdbAny += 1;
        pdbTop.push({ acc, gene, n: pdb });
      }
      const tm = (c[col.Transmembrane]?.match(/TRANSMEM /g) ?? []).length;
      tmCounts[Math.min(tm, 25)] += 1;
      if (tm) tmProteins += 1;
      if (/SIGNAL/.test(c[col['Signal peptide']] ?? '')) signal += 1;
      for (const loc of locationsOf(c[col['Subcellular location [CC]']])) bump(locations, own(loc));
      const domains = [...new Set((c[col.Pfam] ?? '').split(';').filter(Boolean))];
      domains.forEach((d) => bump(pfam, own(d)));
      perProtein[Math.min(domains.length, 15)] += 1;
      if (domains.length) pfamAny += 1;
      const dis = coverage(disordered(c[col.Region]));
      disorderResidues += dis;
      if (dis) disorderProteins += 1;
      share[Math.min(9, Math.floor((dis / length) * 10))] += 1;
      for (const note of notes(c[col['Modified residue']], 'MOD_RES')) {
        if (/^Phosphoserine/.test(note)) ptm.phosphoS += 1;
        else if (/^Phosphothreonine/.test(note)) ptm.phosphoT += 1;
        else if (/^Phosphotyrosine/.test(note)) ptm.phosphoY += 1;
        else if (/acetyl/i.test(note)) ptm.acetyl += 1;
        else if (/methyl/i.test(note)) ptm.methyl += 1;
      }
      for (const note of notes(c[col.Glycosylation], 'CARBOHYD')) {
        if (/^N-linked/.test(note)) ptm.nGlyco += 1;
        else if (/^O-linked/.test(note)) ptm.oGlyco += 1;
      }
      for (const note of notes(c[col['Cross-link']], 'CROSSLNK')) {
        if (/ubiquitin/i.test(note)) ptm.ubiquitin += 1;
        else if (/SUMO/.test(note)) ptm.sumo += 1;
      }
    },
    result() {
      const byLength = [...proteins].sort((a, b) => b.length - a.length);
      return {
        count: proteins.length,
        lengths: { edges: logEdges(1, 5, 10), counts: histogram(proteins.map((p) => p.length), logEdges(1, 5, 10)), summary: summary(proteins.map((p) => p.length)) },
        longest: byLength.slice(0, 10),
        shortest: byLength.slice(-10).reverse(),
        existence: Object.fromEntries(existence),
        pdb: { withAny: pdbAny, counts: pdbCounts, top: pdbTop.sort((a, b) => b.n - a.n).slice(0, 15) },
        transmem: { counts: tmCounts, proteins: tmProteins },
        signal,
        locations: ranked(locations, 'name').slice(0, 20),
        pfam: { top: ranked(pfam, 'id').slice(0, 25), withAny: pfamAny, perProtein },
        disorder: { residues: disorderResidues, proteins: disorderProteins, share, edges: linearEdges(0, 1, 0.1) },
        ptm,
      };
    },
  };
};

const unquote = (s) => s.replace(/^"|"$/g, '');
const listOf = (s) => (s ? s.split(/,\s*/).map((x) => x.trim()).filter(Boolean) : []);

const CONTEXTS = {
  tissue: 'RNA tissue specificity',
  singleCell: 'RNA single cell type specificity',
  blood: 'RNA blood cell specificity',
  brain: 'RNA brain regional specificity',
  cellLine: 'RNA cell line specificity',
  cancer: 'RNA cancer specificity',
};

const kindOf = (label) => {
  if (/^group enriched$/i.test(label)) return 'group';
  if (/enriched$/i.test(label)) return 'enriched';
  if (/enhanced$/i.test(label)) return 'enhanced';
  if (/^low .*specificity$/i.test(label)) return 'low';
  if (/^not detected$/i.test(label)) return 'none';
  return null;
};

const CONCENTRATION_EDGES = logEdges(0, 15, 1);
const concentration = (v) => {
  const n = Number(v);
  return v && Number.isFinite(n) && n > 0 ? n : null;
};

export const hpaSummary = (text) => {
  const lines = text.split('\n').filter((l) => l.length);
  const head = lines[0].split('\t').map(unquote);
  const col = Object.fromEntries(head.map((h, i) => [h, i]));
  const specificity = new Map();
  const enriched = new Map();
  const cells = new Map();
  const secretome = new Map();
  const locations = new Map();
  const classes = new Map();
  const contexts = Object.fromEntries(Object.keys(CONTEXTS).map((k) => [k, {
    enriched: 0, group: 0, enhanced: 0, low: 0, none: 0,
  }]));
  const plasma = [];
  for (const line of lines.slice(1)) {
    const c = line.replace(/\r$/, '').split('\t').map(unquote);
    for (const [key, name] of Object.entries(CONTEXTS)) {
      const kind = kindOf(c[col[name]] ?? '');
      if (kind) contexts[key][kind] += 1;
    }
    const im = concentration(c[col['Blood concentration - Conc. blood IM [pg/L]']]);
    const ms = concentration(c[col['Blood concentration - Conc. blood MS [pg/L]']]);
    if (im || ms) plasma.push({ gene: c[col.Gene], im, ms });
    const spec = c[col['RNA tissue specificity']];
    if (spec) bump(specificity, spec);
    if (spec === 'Tissue enriched') for (const t of (c[col['RNA tissue specific nTPM']] ?? '').split(';')) if (t.includes(':')) bump(enriched, t.split(':')[0].trim());
    if (c[col['RNA single cell type specificity']] === 'Cell type enriched') {
      for (const t of (c[col['RNA single cell type specific nCPM']] ?? '').split(';')) if (t.includes(':')) bump(cells, t.split(':')[0].trim());
    }
    if (c[col['Secretome location']]) bump(secretome, c[col['Secretome location']]);
    listOf(c[col['Subcellular main location']]).forEach((l) => bump(locations, l));
    listOf(c[col['Protein class']]).forEach((l) => bump(classes, l));
  }
  const levels = (key) => plasma.map((p) => p[key]).filter(Boolean);
  return {
    genes: lines.length - 1,
    contexts,
    blood: {
      edges: CONCENTRATION_EDGES,
      im: { count: levels('im').length, counts: histogram(levels('im'), CONCENTRATION_EDGES) },
      ms: { count: levels('ms').length, counts: histogram(levels('ms'), CONCENTRATION_EDGES) },
      top: plasma.sort((a, b) => Math.max(b.im ?? 0, b.ms ?? 0) - Math.max(a.im ?? 0, a.ms ?? 0)).slice(0, 15),
    },
    specificity: Object.fromEntries(specificity),
    enriched: ranked(enriched, 'tissue'),
    cellTypes: ranked(cells, 'cell').slice(0, 25),
    secretome: Object.fromEntries(secretome),
    locations: ranked(locations, 'name').slice(0, 20),
    classes: ranked(classes, 'name'),
  };
};

export const gtexSummary = (lines) => {
  let tissues = null;
  let genes = 0;
  let n = 0;
  let sums;
  let products;
  let expressed;
  let specific;
  for (const line of lines) {
    if (!line || line.startsWith('#') || /^\d+\t\d+$/.test(line)) continue;
    const c = line.split('\t');
    if (!tissues) {
      tissues = c.slice(2).map(own);
      const k = tissues.length;
      sums = new Float64Array(k);
      products = new Float64Array(k * k);
      expressed = new Array(k).fill(0);
      specific = new Array(k).fill(0);
      continue;
    }
    const tpm = c.slice(2).map(Number);
    genes += 1;
    let top = -1;
    let second = -1;
    let at = -1;
    for (let i = 0; i < tpm.length; i += 1) {
      if (tpm[i] >= 1) expressed[i] += 1;
      if (tpm[i] > top) {
        second = top;
        top = tpm[i];
        at = i;
      } else if (tpm[i] > second) second = tpm[i];
    }
    if (top < 1) continue;
    if (top >= 5 * Math.max(second, 0)) specific[at] += 1;
    n += 1;
    const logs = tpm.map((v) => Math.log10(v + 1));
    for (let i = 0; i < tpm.length; i += 1) {
      sums[i] += logs[i];
      for (let j = i; j < tpm.length; j += 1) products[i * tpm.length + j] += logs[i] * logs[j];
    }
  }
  const k = tissues.length;
  const mean = Array.from(sums, (s) => s / n);
  const cov = (i, j) => products[Math.min(i, j) * k + Math.max(i, j)] - n * mean[i] * mean[j];
  const correlation = tissues.map((_, i) => tissues.map((__, j) => Number((cov(i, j) / Math.sqrt(cov(i, i) * cov(j, j))).toFixed(3))));
  return {
    tissues, genes, correlated: n, expressed, specific, correlation, order: clusterOrder(correlation),
  };
};

const LEVELS = ['Tclin', 'Tchem', 'Tbio', 'Tdark'];
const facetOf = (json, name) => json?.data.targets.facets.find((f) => f.facet === name)?.values ?? [];

export const pharosSummary = (json, byLevel = {}) => {
  const { count } = json.data.targets;
  const levels = Object.fromEntries(facetOf(json, 'Target Development Level').map((v) => [v.name, v.value]));
  const pairs = (name) => facetOf(json, name).map((v) => [Number(v.name), v.value]).sort((a, b) => a[0] - b[0]);
  const families = new Map();
  for (const level of LEVELS) {
    for (const v of facetOf(byLevel[level], 'IDG Family')) {
      if (!families.has(v.name)) families.set(v.name, { family: v.name, ...Object.fromEntries(LEVELS.map((l) => [l, 0])) });
      families.get(v.name)[level] = v.value;
    }
  }
  const total = (f) => LEVELS.reduce((s, l) => s + f[l], 0);
  return {
    count,
    levels: Object.fromEntries(LEVELS.map((l) => [l, levels[l] ?? 0])),
    novelty: pairs('Log Novelty'),
    pubmed: pairs('Log PubMed Score'),
    families: [...families.values()].sort((a, b) => total(b) - total(a) || a.family.localeCompare(b.family)),
  };
};
