import { csvLine } from './variation.mjs';
import { CHANNELS, SUBSTITUTIONS } from '../src/lib.js';

const own = (s) => Buffer.from(s, 'utf8').toString('utf8');
const zeros = (n) => new Array(n).fill(0);
const spaced = (s) => String(s).replace(/_/g, ' ');

export const CLASSES = ['Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign/Likely benign', 'Benign', 'Conflicting', 'Other'];

const PATHOGENIC = new Set(['pathogenic', 'likely pathogenic']);
const BENIGN = new Set(['benign', 'likely benign']);
const UNCERTAIN = new Set(['uncertain significance', 'vus-high', 'vus-mid', 'vus-low']);

export const germlineClass = (raw) => {
  const text = spaced(raw ?? '').trim();
  if (!text || text === '-' || /^no classification/i.test(text)) return null;
  const primary = text.split(/[|;]|, (?!low penetrance)/)[0].trim();
  if (/^conflicting/i.test(primary)) return 'Conflicting';
  const terms = new Set(primary.split('/').map((t) => t.replace(/, low penetrance$/i, '').trim().toLowerCase()));
  const within = (set) => [...terms].every((t) => set.has(t));
  if (within(PATHOGENIC)) {
    if (terms.size > 1) return 'Pathogenic/Likely pathogenic';
    return terms.has('pathogenic') ? 'Pathogenic' : 'Likely pathogenic';
  }
  if (within(BENIGN)) {
    if (terms.size > 1) return 'Benign/Likely benign';
    return terms.has('benign') ? 'Benign' : 'Likely benign';
  }
  if (within(UNCERTAIN)) return 'Uncertain significance';
  return 'Other';
};

export const groupOf = (cls) => {
  if (cls === 'Pathogenic' || cls === 'Pathogenic/Likely pathogenic' || cls === 'Likely pathogenic') return 'P/LP';
  if (cls === 'Benign' || cls === 'Benign/Likely benign' || cls === 'Likely benign') return 'B/LB';
  if (cls === 'Uncertain significance') return 'VUS';
  if (cls === 'Conflicting') return 'Conflicting';
  return null;
};

export const CONSEQUENCES = ['splice site', 'nonsense', 'frameshift', 'stop lost', 'start lost', 'in-frame indel', 'missense', 'synonymous', 'UTR', 'intron', 'non-coding RNA', 'near gene', 'none'];

const CONSEQUENCE_OF = {
  splice_donor_variant: 'splice site',
  splice_acceptor_variant: 'splice site',
  nonsense: 'nonsense',
  stop_gained: 'nonsense',
  frameshift_variant: 'frameshift',
  stop_lost: 'stop lost',
  initiator_codon_variant: 'start lost',
  start_lost: 'start lost',
  inframe_deletion: 'in-frame indel',
  inframe_insertion: 'in-frame indel',
  inframe_indel: 'in-frame indel',
  missense_variant: 'missense',
  synonymous_variant: 'synonymous',
  '5_prime_UTR_variant': 'UTR',
  '3_prime_UTR_variant': 'UTR',
  intron_variant: 'intron',
  'non-coding_transcript_variant': 'non-coding RNA',
  genic_upstream_transcript_variant: 'near gene',
  genic_downstream_transcript_variant: 'near gene',
  '2KB_upstream_variant': 'near gene',
  '500B_downstream_variant': 'near gene',
};

export const consequenceOf = (mc) => {
  let best = CONSEQUENCES.length - 1;
  for (const part of String(mc ?? '').split(',')) {
    const cls = CONSEQUENCE_OF[part.split('|')[1]];
    if (cls) best = Math.min(best, CONSEQUENCES.indexOf(cls));
  }
  return CONSEQUENCES[best];
};

export const starsOf = (raw) => {
  const s = spaced(raw ?? '');
  if (!s || s === '-') return null;
  if (/practice guideline/.test(s)) return 4;
  if (/expert panel/.test(s)) return 3;
  if (/no assertion/.test(s)) return 0;
  if (/multiple submitters, no conflicts/.test(s)) return 2;
  if (/criteria provided/.test(s)) return 1;
  return 0;
};

const SUBS = SUBSTITUTIONS.map((s) => s.split('>'));
const NT = 'ACGT';
const COMP = {
  A: 'T', C: 'G', G: 'C', T: 'A',
};

export { CHANNELS };

export const channelOf = (b5, ref, b3, alt) => {
  if (![b5, ref, b3, alt].every((b) => COMP[b]) || ref === alt) return -1;
  const flip = ref === 'G' || ref === 'A';
  const [f5, r, f3, a] = flip ? [COMP[b3], COMP[ref], COMP[b5], COMP[alt]] : [b5, ref, b3, alt];
  return SUBS.findIndex(([x, y]) => x === r && y === a) * 16 + NT.indexOf(f5) * 4 + NT.indexOf(f3);
};

const infoOf = (info, keys) => {
  const out = {};
  for (const kv of info.split(';')) {
    const i = kv.indexOf('=');
    if (i > 0) {
      const k = kv.slice(0, i);
      if (keys.has(k)) out[k] = kv.slice(i + 1);
    }
  }
  return out;
};

const isoDate = (s) => (/^\d{8}$/.test(s) ? `${s.slice(0, 4)}-${s.slice(4, 6)}-${s.slice(6)}` : s);
const fileDate = (line) => {
  const m = /^##fileDate=(\S+)/.exec(line);
  return m ? isoDate(m[1]) : null;
};

const VCF_KEYS = new Set(['CLNSIG', 'CLNREVSTAT', 'CLNVC', 'MC', 'ORIGIN', 'ONC', 'SCI']);
const GERMLINE_BITS = 1 | 4 | 8 | 16 | 32 | 64 | 128;
const bump = (obj, key) => {
  obj[key] = (obj[key] ?? 0) + 1;
};

export const clinvarVcf = ({ base }) => {
  let release = null;
  let records = 0;
  const classes = Object.fromEntries([...CLASSES, 'none'].map((c) => [c, 0]));
  const consequences = Object.fromEntries(CONSEQUENCES.map((k) => [k, Object.fromEntries(CLASSES.map((c) => [c, 0]))]));
  const stars = Object.fromEntries(CLASSES.map((c) => [c, zeros(5)]));
  const types = {};
  const origin = {
    germline: 0, somatic: 0, both: 0, unknown: 0,
  };
  const oncogenicity = {};
  const clinicalImpact = {};
  const spectrum = {
    plp: zeros(96), blb: zeros(96), mismatches: 0, skipped: 0,
  };
  const snv = (c, group) => {
    const pos = Number(c[1]);
    const ref = c[3].toUpperCase();
    const seen = base(c[0], pos);
    if (seen === null || seen === 'N') {
      spectrum.skipped += 1;
      return;
    }
    if (seen !== ref) {
      spectrum.mismatches += 1;
      return;
    }
    const k = channelOf(base(c[0], pos - 1), ref, base(c[0], pos + 1), c[4].toUpperCase());
    if (k < 0) spectrum.skipped += 1;
    else spectrum[group === 'P/LP' ? 'plp' : 'blb'][k] += 1;
  };
  return {
    add(line) {
      if (!line) return;
      if (line.startsWith('#')) {
        release = fileDate(line) ?? release;
        return;
      }
      const c = line.split('\t');
      const info = infoOf(c[7] ?? '', VCF_KEYS);
      records += 1;
      bump(types, spaced(info.CLNVC ?? 'unknown'));
      const o = Number(info.ORIGIN ?? 0);
      const g = (o & GERMLINE_BITS) !== 0;
      const s = (o & 2) !== 0;
      origin[(g && s && 'both') || (g && 'germline') || (s && 'somatic') || 'unknown'] += 1;
      if (info.ONC && !/^no_classification/.test(info.ONC)) bump(oncogenicity, spaced(info.ONC));
      if (info.SCI && !/^no_classification/.test(info.SCI)) bump(clinicalImpact, spaced(info.SCI));
      const cls = germlineClass(info.CLNSIG);
      if (!cls) {
        classes.none += 1;
        return;
      }
      classes[cls] += 1;
      consequences[consequenceOf(info.MC)][cls] += 1;
      const st = starsOf(info.CLNREVSTAT);
      if (st !== null) stars[cls][st] += 1;
      const group = groupOf(cls);
      if ((group === 'P/LP' || group === 'B/LB') && info.CLNVC === 'single_nucleotide_variant' && c[3].length === 1 && c[4].length === 1) snv(c, group);
    },
    result: () => ({
      release, records, classes, consequences, stars, types, origin, oncogenicity, clinicalImpact, spectrum,
    }),
  };
};

export const vcfClassCounts = () => {
  let release = null;
  let records = 0;
  const classes = Object.fromEntries([...CLASSES, 'none'].map((c) => [c, 0]));
  return {
    add(line) {
      if (!line) return;
      if (line.startsWith('#')) {
        release = fileDate(line) ?? release;
        return;
      }
      records += 1;
      const m = /(?:^|;)CLNSIG=([^;]+)/.exec(line.slice(line.lastIndexOf('\t') + 1));
      classes[germlineClass(m?.[1]) ?? 'none'] += 1;
    },
    result: () => ({ release, records, classes }),
  };
};

export const earliestFile = (listing, re) => {
  const hits = [...listing.matchAll(new RegExp(re.source, 'g'))].map((m) => ({ file: m[0], version: m[1] }));
  return hits.length ? hits.reduce((best, h) => (h.version < best.version ? h : best)) : null;
};

const count = (v) => (/^\d+$/.test(v ?? '') ? Number(v) : 0);

export const geneSummary = (text, coding) => {
  const lines = text.split('\n');
  const date = /dated (.+?)\s*$/.exec(lines[0] ?? '')?.[1] ?? null;
  const col = Object.fromEntries(lines.find((l) => l.startsWith('#Symbol')).slice(1).split('\t').map((h, i) => [h, i]));
  const rows = [];
  for (const line of lines) {
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    if (!coding.has(c[col.Symbol])) continue;
    rows.push({
      gene: own(c[col.Symbol]),
      plp: count(c[col.Alleles_reported_Pathogenic_Likely_pathogenic]),
      vus: count(c[col.Number_uncertain]),
      conflicts: count(c[col.Number_with_conflicts]),
      alleles: count(c[col.Total_alleles]),
    });
  }
  const ranked = rows.filter((r) => r.plp > 0).sort((a, b) => b.plp - a.plp || a.gene.localeCompare(b.gene));
  const plp = ranked.reduce((s, r) => s + r.plp, 0);
  const cum = [0];
  ranked.forEach((r, i) => cum.push(cum[i] + r.plp));
  return {
    date,
    genes: rows.length,
    withPlp: ranked.length,
    plp,
    vus: rows.reduce((s, r) => s + r.vus, 0),
    top: ranked.slice(0, 20),
    lorenz: Array.from({ length: 101 }, (_, p) => Number((cum[Math.round((p / 100) * ranked.length)] / (plp || 1)).toFixed(4))),
  };
};

const METHODS = ['clinical testing', 'research', 'literature only', 'curation', 'phenotyping'];

export const organizations = (text) => {
  const lines = text.split('\n');
  const col = Object.fromEntries(lines.find((l) => l.startsWith('#organization')).slice(1).split('\t').map((h, i) => [h, i]));
  const types = new Map();
  const countries = new Map();
  const orgs = [];
  const methods = Object.fromEntries(METHODS.map((m) => [m, 0]));
  let submissions = 0;
  for (const line of lines) {
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    const n = Number(c[col['number of ClinVar submissions']]) || 0;
    const type = own(c[col['institution type']] || 'not stated');
    const country = own(c[col.country] || 'not stated');
    submissions += n;
    for (const [map, key] of [[types, type], [countries, country]]) {
      const v = map.get(key) ?? { orgs: 0, submissions: 0 };
      v.orgs += 1;
      v.submissions += n;
      map.set(key, v);
    }
    for (const m of METHODS) methods[m] += Number(c[col[`number of submissions from ${m}`]]) || 0;
    orgs.push({
      name: own(c[col.organization]), type, country, submissions: n,
    });
  }
  const ranked = (map, key) => [...map].map(([k, v]) => ({ [key]: k, ...v })).sort((a, b) => b.submissions - a.submissions || a[key].localeCompare(b[key]));
  return {
    organizations: orgs.length,
    submissions,
    types: ranked(types, 'type'),
    countries: ranked(countries, 'country').slice(0, 15),
    top: orgs.sort((a, b) => b.submissions - a.submissions || a.name.localeCompare(b.name)).slice(0, 15),
    methods,
  };
};

const kindOf = (aa, rest) => {
  if (aa === 'Ter' || rest === '?') return 'other';
  if (rest === '=') return 'synonymous';
  if (/fs/.test(rest) || /Ter|\*/.test(rest)) return 'truncating';
  if (/del|ins|dup/.test(rest)) return 'inframe';
  if (/^[A-Z][a-z]{2}$/.test(rest)) return rest === aa ? 'synonymous' : 'missense';
  return 'other';
};

export const proteinChange = (name) => {
  const m = /^(N[MR]_\d+\.\d+)\([^)]*\):c\.(\S+)(?: \(p\.(.+)\))?$/.exec(name ?? '');
  if (!m) return null;
  const [, transcript, cdna, protein] = m;
  if (protein) {
    const r = /^([A-Z][a-z]{2})(\d+)(.*)$/.exec(protein.replace(/[()]/g, ''));
    if (!r) return null;
    return {
      transcript, position: Number(r[2]), kind: kindOf(r[1], r[3]), aa: r[1],
    };
  }
  const splice = /^(\d+)[+-][12](?!\d)/.exec(cdna);
  return splice ? {
    transcript, position: Math.ceil(Number(splice[1]) / 3), kind: 'splice', aa: '',
  } : null;
};

const SUBMITTER_BINS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11+'];
const GROUPS = ['P/LP', 'VUS', 'B/LB', 'Conflicting'];

export const variantTable = ({ genes }) => {
  const wanted = new Set(genes);
  let col = null;
  let records = 0;
  const types = {};
  const submitters = { bins: SUBMITTER_BINS, ...Object.fromEntries(GROUPS.map((g) => [g, zeros(SUBMITTER_BINS.length)])) };
  const evaluated = Object.fromEntries(GROUPS.map((g) => [g, {}]));
  const perGene = new Map();
  const place = (gene, side, name) => {
    const g = perGene.get(gene) ?? { totals: { up: 0, down: 0 }, transcripts: new Map() };
    perGene.set(gene, g);
    g.totals[side] += 1;
    const tx = /^(N[MR]_\d+\.\d+)\(/.exec(name)?.[1];
    if (!tx) return;
    const t = g.transcripts.get(tx) ?? { rows: 0, up: new Map(), down: new Map() };
    g.transcripts.set(tx, t);
    t.rows += 1;
    const p = proteinChange(name);
    if (!p) return;
    const key = `${p.position}|${p.kind}|${p.aa}`;
    t[side].set(key, (t[side].get(key) ?? 0) + 1);
  };
  return {
    add(line) {
      if (!line) return;
      if (line.startsWith('#')) {
        col = Object.fromEntries(line.slice(1).split('\t').map((h, i) => [h, i]));
        return;
      }
      const c = line.split('\t');
      if (c[col.Assembly] !== 'GRCh38') return;
      records += 1;
      bump(types, c[col.Type]);
      const group = groupOf(germlineClass(c[col.ClinicalSignificance]));
      if (!group) return;
      const subs = Number(c[col.NumberSubmitters]);
      if (subs >= 1) submitters[group][Math.min(subs, SUBMITTER_BINS.length) - 1] += 1;
      bump(evaluated[group], /(\d{4})$/.exec(c[col.LastEvaluated])?.[1] ?? 'none');
      if ((group === 'P/LP' || group === 'B/LB') && wanted.has(c[col.GeneSymbol])) place(c[col.GeneSymbol], group === 'P/LP' ? 'up' : 'down', c[col.Name]);
    },
    result() {
      const lollipops = {};
      const listOf = (map) => [...map].map(([k, n]) => {
        const [pos, kind, aa] = k.split('|');
        return [Number(pos), kind, n, aa];
      }).sort((a, b) => a[0] - b[0] || a[1].localeCompare(b[1]));
      const placed = (list) => list.reduce((s, r) => s + r[2], 0);
      for (const [gene, g] of perGene) {
        const best = [...g.transcripts].sort((a, b) => b[1].rows - a[1].rows)[0];
        if (!best) continue;
        const up = listOf(best[1].up);
        const down = listOf(best[1].down);
        lollipops[gene] = {
          transcript: best[0], up, down, unplaced: { up: g.totals.up - placed(up), down: g.totals.down - placed(down) },
        };
      }
      return {
        records, types, submitters, evaluated, lollipops,
      };
    },
  };
};

export const pfamDomains = (json) => (json?.results ?? []).flatMap((r) => (r.proteins ?? []).flatMap((p) => (p.entry_protein_locations ?? []).flatMap((loc) => (loc.fragments ?? []).map((f) => ({
  id: r.metadata.accession, name: r.metadata.name, start: f.start, end: f.end,
}))))).sort((a, b) => a.start - b.start || a.end - b.end);

export const VALIDITY = ['Definitive', 'Strong', 'Moderate', 'Limited', 'Disputed', 'Refuted', 'No Known Disease Relationship'];
const SUPPORTED = new Set(['Definitive', 'Strong', 'Moderate']);

const bumpIn = (obj, key, sub) => {
  obj[key] ??= {};
  bump(obj[key], sub);
};

const moiLabel = (set) => {
  if (set.has('AD') && set.has('AR')) return 'both';
  return set.size === 1 ? [...set][0] : 'several';
};

const validityRows = (text) => {
  const lines = text.split('\n').map((l) => l.replace(/\r$/, ''));
  const at = lines.findIndex((l) => l.startsWith('"GENE SYMBOL"'));
  return {
    col: Object.fromEntries(csvLine(lines[at]).map((h, i) => [h, i])),
    rows: lines.slice(at + 1).filter((line) => line && !line.startsWith('"+++')).map(csvLine),
  };
};

export const supportedGenes = (text) => {
  const { col, rows } = validityRows(text);
  return new Set(rows.filter((c) => SUPPORTED.has(c[col.CLASSIFICATION])).map((c) => c[col['GENE ID (HGNC)']]));
};

export const validitySummary = (text) => {
  const { col, rows } = validityRows(text);
  const classes = Object.fromEntries(VALIDITY.map((c) => [c, 0]));
  const byMoi = {};
  const byYear = {};
  const gceps = new Map();
  const genes = new Set();
  const diseases = new Set();
  const moi = new Map();
  let curations = 0;
  for (const c of rows) {
    const gene = c[col['GENE SYMBOL']];
    const mode = c[col.MOI];
    const cls = c[col.CLASSIFICATION];
    curations += 1;
    genes.add(gene);
    diseases.add(c[col['DISEASE ID (MONDO)']]);
    if (cls in classes) classes[cls] += 1;
    bumpIn(byMoi, mode, cls);
    bumpIn(byYear, c[col['CLASSIFICATION DATE']].slice(0, 4), cls);
    gceps.set(c[col.GCEP], (gceps.get(c[col.GCEP]) ?? 0) + 1);
    if (SUPPORTED.has(cls)) moi.set(gene, (moi.get(gene) ?? new Set()).add(mode));
  }
  return {
    date: /FILE CREATED: (\d{4}-\d{2}-\d{2})/.exec(text)?.[1] ?? null,
    curations,
    genes: genes.size,
    diseases: diseases.size,
    classes,
    byMoi,
    byYear,
    gceps: [...gceps].map(([name, n]) => ({ name, n })).sort((a, b) => b.n - a.n || a.name.localeCompare(b.name)),
    moiGenes: new Map([...moi].map(([g, set]) => [g, moiLabel(set)])),
  };
};

const span = (lo, hi) => (Number.isFinite(lo) && Number.isFinite(hi) ? [lo, hi] : null);

const locationOf = (text) => {
  const t = String(text ?? '');
  if (/^5' (UTR|promoter)/i.test(t)) return "5' UTR";
  if (/^3' UTR/i.test(t)) return "3' UTR";
  if (/^Coding/i.test(t)) return 'coding exon';
  if (/^Intron/i.test(t)) return 'intron';
  if (/lncRNA|^Non-coding/i.test(t)) return 'non-coding';
  if (/^Exon/i.test(t)) return 'exon';
  return 'unknown';
};

export const repeatLoci = (json) => {
  const loci = (Array.isArray(json) ? json : []).map((l) => ({
    id: l.id,
    gene: l.gene,
    disease: l.disease_id,
    name: l.disease,
    evidence: l.evidence?.[0] ?? null,
    inheritance: (l.inheritance ?? []).join('/'),
    motif: l.pathogenic_motif_reference_orientation?.[0] ?? l.reference_motif_reference_orientation?.[0] ?? null,
    motifLength: l.motif_len ?? null,
    location: locationOf(l.location_in_gene),
    mechanism: l.mechanism ?? null,
    year: Number(/\b((?:19|20)\d{2})\b/.exec(l.year ?? '')?.[1]) || null,
    benign: span(l.benign_min, l.benign_max),
    intermediate: span(l.intermediate_min, l.intermediate_max),
    pathogenic: span(l.pathogenic_min, l.pathogenic_max),
    onset: span(l.typ_age_onset_min, l.typ_age_onset_max),
    genereviews: l.genereviews?.[0] ?? null,
  }));
  const evidence = {};
  for (const l of loci) if (l.evidence) bump(evidence, l.evidence);
  return { loci, evidence };
};

const CHROM_ORDER = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];
const score = (v) => (/^\d+$/.test(v ?? '') ? Number(v) : null);

export const regionSummary = (text) => {
  const lines = text.split('\n');
  const col = Object.fromEntries(lines.find((l) => l.startsWith('#ISCA ID')).slice(1).split('\t').map((h, i) => [h, i]));
  const hi = {};
  const ts = {};
  const sensitive = [];
  let regions = 0;
  for (const line of lines) {
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    regions += 1;
    const h = score(c[col['Haploinsufficiency Score']]);
    const t = score(c[col['Triplosensitivity Score']]);
    if (h !== null) bump(hi, h);
    if (t !== null) bump(ts, t);
    const loc = /^chr(\w+):(\d+)-(\d+)$/.exec((c[col['Genomic Location']] ?? '').trim());
    if (!loc || (h !== 3 && t !== 3)) continue;
    const name = own(c[col['ISCA Region Name']].trim());
    sensitive.push({
      id: own(c[col['ISCA ID']]), name, band: own(c[col.cytoBand]), chrom: loc[1], start: Number(loc[2]), end: Number(loc[3]), hi: h, ts: t, recurrent: /recurrent/i.test(name),
    });
  }
  sensitive.sort((a, b) => CHROM_ORDER.indexOf(a.chrom) - CHROM_ORDER.indexOf(b.chrom) || a.start - b.start);
  return {
    date: lines[1]?.startsWith('#') ? lines[1].slice(1).trim() : null, regions, hi, ts, sensitive,
  };
};
