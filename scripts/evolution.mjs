import { inflateSync } from 'node:zlib';
import {
  quantile, histogram, linearEdges, logEdges, binIndex, mergeIntervals,
} from '../src/stats.js';
import { round } from './tables.mjs';
import { CHROMS, bare, isStandard } from './parsers.mjs';
import { gtfRecord } from './sequence.mjs';

export const KINDS = [
  { key: 'primates', label: 'Primates', id: 9443 },
  { key: 'mammals', label: 'Other mammals', id: 40674 },
  { key: 'birds', label: 'Birds', id: 8782 },
  { key: 'reptiles', label: 'Reptiles', id: 8457 },
  { key: 'amphibians', label: 'Amphibians', id: 8292 },
  { key: 'rayfinned', label: 'Ray-finned fishes', id: 7898 },
  { key: 'cartilaginous', label: 'Cartilaginous fishes', id: 7777 },
  { key: 'other', label: 'Other', id: null },
];

export const kindOf = (parents, kinds = KINDS) => (kinds.find((k) => k.id === null || parents.includes(k.id)) ?? kinds.at(-1)).key;

export const taxonomyReports = (json) => new Map((json?.reports ?? []).filter((r) => r.taxonomy).map(({ taxonomy: t }) => [String(t.tax_id), {
  id: t.tax_id, name: t.current_scientific_name?.name ?? String(t.tax_id), group: t.group_name ?? '', parents: t.parents ?? [],
}]));

export const lineageOf = (parents, names) => parents.map((id) => ({ id, name: names.get(String(id))?.name ?? String(id) }));

export const sharedDepths = (taxa, humanParents) => {
  const out = new Map();
  for (const [id, t] of taxa) {
    const own = new Set(t.parents);
    let depth = -1;
    for (let i = 0; i < humanParents.length; i += 1) if (own.has(humanParents[i])) depth = i;
    if (depth >= 0) out.set(id, depth);
  }
  return out;
};

const humanRow = (line) => (line.startsWith('9606\t') ? line.split('\t') : null);

export const orthologPartners = ({ genes }) => {
  const partners = new Map();
  return {
    add(line) {
      const c = humanRow(line);
      if (c && genes.has(c[1])) partners.set(c[3], (partners.get(c[3]) ?? 0) + 1);
    },
    result: () => partners,
  };
};

export const orthologDepth = ({ genes, depthOf }) => {
  const deepest = new Map();
  return {
    add(line) {
      const c = humanRow(line);
      if (!c || !genes.has(c[1]) || !depthOf.has(c[3])) return;
      const d = depthOf.get(c[3]);
      if (!deepest.has(c[1]) || d < deepest.get(c[1])) deepest.set(c[1], d);
    },
    result: () => deepest,
  };
};

const nodeIndex = (depths) => new Map(depths.map((d, i) => [d, i]));
const closestFirst = (values) => [...new Set(values)].sort((a, b) => b - a);
const ascending = (values) => [...values].sort((a, b) => a - b);

export const vertebrateOrthologs = ({
  partners, taxa, lineage, depthOf, genes,
}) => {
  const placed = [...partners].filter(([id]) => taxa.has(id) && depthOf.has(id));
  const depths = closestFirst(placed.map(([id]) => depthOf.get(id)));
  const index = nodeIndex(depths);
  const species = placed.map(([id, n]) => [id, taxa.get(id).name, kindOf(taxa.get(id).parents), index.get(depthOf.get(id)), n])
    .sort((a, b) => a[3] - b[3] || b[4] - a[4] || a[1].localeCompare(b[1]));
  return {
    genes,
    nodes: depths.map((d, i) => {
      const counts = ascending(species.filter((s) => s[3] === i).map((s) => s[4]));
      return { name: lineage[d].name, species: counts.length, median: quantile(counts, 0.5) };
    }),
    species,
  };
};

const MODEL_ORDER = ['10090', '10116', '7955', '8364', '8355', '7227', '6239', '559292'];
const modelRank = (taxon) => (MODEL_ORDER.includes(taxon) ? MODEL_ORDER.indexOf(taxon) : MODEL_ORDER.length);

export const allianceOrthologs = (text, { genes, disease }) => {
  const species = new Map();
  const byGene = new Map();
  let head = null;
  for (const raw of text.split('\n')) {
    const line = raw.replace(/\r$/, '');
    if (!line || line.startsWith('#')) continue;
    const c = line.split('\t');
    if (!head) {
      head = Object.fromEntries(c.map((h, i) => [h, i]));
      continue;
    }
    const gene = c[head.Gene1ID];
    if (c[head.Gene1SpeciesTaxonID] !== 'NCBITaxon:9606' || !genes.has(gene)) continue;
    const taxon = c[head.Gene2SpeciesTaxonID].replace('NCBITaxon:', '');
    const of = Number(c[head.OutOfAlgorithms]);
    if (!species.has(taxon)) species.set(taxon, { taxon, name: c[head.Gene2SpeciesName], methods: 0, best: new Map() });
    const s = species.get(taxon);
    s.methods = Math.max(s.methods, of);
    s.best.set(gene, Math.max(s.best.get(gene) ?? 0, Number(c[head.AlgorithmsMatch]) / of));
    if (!byGene.has(gene)) byGene.set(gene, new Set());
    byGene.get(gene).add(taxon);
  }
  const diseaseGenes = [...disease].filter((g) => genes.has(g));
  return {
    version: /Alliance Database Version: (\S+)/.exec(text)?.[1] ?? null,
    generated: /Date file generated \(UTC\): (\d{4}-\d{2}-\d{2})/.exec(text)?.[1] ?? null,
    genes: genes.size,
    disease: diseaseGenes.length,
    species: [...species.values()].sort((a, b) => modelRank(a.taxon) - modelRank(b.taxon) || a.name.localeCompare(b.name)).map((s) => {
      const best = [...s.best.values()];
      return {
        taxon: s.taxon,
        name: s.name,
        methods: s.methods,
        genes: best.length,
        all: best.filter((f) => f >= 1).length,
        most: best.filter((f) => f >= 0.5 && f < 1).length,
        fewer: best.filter((f) => f < 0.5).length,
        disease: diseaseGenes.filter((g) => s.best.has(g)).length,
      };
    }),
    byGene,
  };
};

const EXAMPLES = 6;
const examplesOf = (symbols) => [...symbols].sort((a, b) => a.localeCompare(b)).slice(0, EXAMPLES);

export const depthSummary = ({
  genes, ncbi, alliance, species, depthOf, lineage,
}) => {
  const groups = new Map(closestFirst(species.filter((t) => depthOf.has(t)).map((t) => depthOf.get(t))).map((d) => [d, []]));
  const none = [];
  for (const g of genes) {
    let d = g.entrez && ncbi.has(g.entrez) ? ncbi.get(g.entrez) : Infinity;
    for (const t of alliance.get(g.hgnc) ?? []) if (depthOf.has(t)) d = Math.min(d, depthOf.get(t));
    if (d === Infinity) {
      none.push(g.symbol);
      continue;
    }
    if (!groups.has(d)) groups.set(d, []);
    groups.get(d).push(g.symbol);
  }
  return {
    genes: genes.length,
    nodes: [...groups].sort((a, b) => b[0] - a[0]).map(([d, list]) => ({ name: lineage[d].name, genes: list.length, examples: examplesOf(list) })),
    none: { genes: none.length, examples: examplesOf(none) },
  };
};

export const ensemblRelease = (html) => {
  const found = [...html.matchAll(/href="release-(\d+)\/"/g)].map((m) => Number(m[1]));
  return found.length ? Math.max(...found) : null;
};

export const ensemblTaxa = (text) => new Map(text.split('\n').filter((l) => l && !l.startsWith('#')).map((l) => l.split('\t'))
  .filter((c) => c[1] && c[3]).map((c) => [c[1], c[3].trim()]));

const TYPES = { ortholog_one2one: 'one2one', ortholog_one2many: 'one2many', ortholog_many2many: 'many2many' };

export const comparaSummary = () => {
  const bySpecies = new Map();
  return {
    add(line) {
      const c = line.split('\t');
      const type = TYPES[c[4]];
      if (!type) return;
      const first = c[2] === 'homo_sapiens';
      const partner = first ? c[7] : c[2];
      if ((!first && c[7] !== 'homo_sapiens') || partner === 'homo_sapiens') return;
      if (!bySpecies.has(partner)) {
        bySpecies.set(partner, {
          one2one: 0, one2many: 0, many2many: 0, highConfidence: 0, identities: [],
        });
      }
      const s = bySpecies.get(partner);
      s[type] += 1;
      if (type !== 'one2one') return;
      s.identities.push(Number(first ? c[3] : c[8]));
      if (c[13] === '1') s.highConfidence += 1;
    },
    result: () => bySpecies,
  };
};

const IDENTITY_EDGES = linearEdges(0, 100, 2);
const tenth = (v) => (Number.isFinite(v) ? round(v, 1) : null);

export const comparaSpecies = (summary, {
  taxa, reports, lineage, depthOf, keys = [],
}) => {
  const rows = [];
  for (const [species, s] of summary) {
    const id = taxa.get(species);
    if (!id || !reports.has(id) || !depthOf.has(id)) continue;
    rows.push({
      species, t: reports.get(id), depth: depthOf.get(id), s, sorted: ascending(s.identities),
    });
  }
  const depths = closestFirst(rows.map((r) => r.depth));
  const index = nodeIndex(depths);
  const species = rows.map((r) => ({
    species: r.species,
    name: r.t.name,
    kind: kindOf(r.t.parents),
    node: index.get(r.depth),
    one2one: r.s.one2one,
    one2many: r.s.one2many,
    many2many: r.s.many2many,
    highConfidence: r.s.one2one ? round(r.s.highConfidence / r.s.one2one, 3) : 0,
    median: tenth(quantile(r.sorted, 0.5)),
    q1: tenth(quantile(r.sorted, 0.25)),
    q3: tenth(quantile(r.sorted, 0.75)),
  })).sort((a, b) => a.node - b.node || (b.median ?? -1) - (a.median ?? -1) || a.name.localeCompare(b.name));
  return {
    nodes: depths.map((d, i) => {
      const medians = ascending(rows.filter((r) => index.get(r.depth) === i && r.sorted.length).map((r) => quantile(r.sorted, 0.5)));
      return { name: lineage[d].name, species: medians.length, median: tenth(quantile(medians, 0.5)) };
    }),
    species,
    identity: {
      edges: IDENTITY_EDGES,
      bySpecies: Object.fromEntries(keys.filter((k) => summary.has(k)).map((k) => [k, histogram(summary.get(k).identities, IDENTITY_EDGES)])),
    },
  };
};

export const GENOME_GROUPS = [
  { key: 'mammals', label: 'Mammals', id: 40674 },
  { key: 'birds', label: 'Birds', id: 8782 },
  { key: 'fishes', label: 'Ray-finned fishes', id: 7898 },
  { key: 'vertebrates', label: 'Other vertebrates', id: 7742 },
  { key: 'insects', label: 'Insects', id: 50557 },
  { key: 'animals', label: 'Other animals', id: 33208 },
  { key: 'plants', label: 'Green plants', id: 33090 },
  { key: 'fungi', label: 'Fungi', id: 4751 },
  { key: 'other', label: 'Other eukaryotes', id: null },
];

const isRefSeq = (accession) => String(accession ?? '').startsWith('GCF_');

export const genomePoints = (reports, taxa) => {
  const bySpecies = new Map();
  for (const r of reports) {
    const coding = Number(r.annotation_info?.stats?.gene_counts?.protein_coding);
    const length = Number(r.assembly_stats?.total_sequence_length);
    const id = String(r.organism?.tax_id ?? '');
    if (!(coding > 0) || !(length > 0) || !id) continue;
    const was = bySpecies.get(id);
    if (was && (isRefSeq(was.accession) || !isRefSeq(r.accession))) continue;
    bySpecies.set(id, {
      accession: r.accession, id, name: r.organism.organism_name, length, coding,
    });
  }
  const keys = GENOME_GROUPS.map((g) => g.key);
  const points = [...bySpecies.values()].sort((a, b) => a.length - b.length || a.name.localeCompare(b.name))
    .map((g) => [round(g.length / 1e6, 2), g.coding, keys.indexOf(kindOf(taxa.get(g.id)?.parents ?? [], GENOME_GROUPS)), g.name]);
  return { total: points.length, groups: GENOME_GROUPS.map((g) => g.label), points };
};

export const ASSEMBLY_LEVELS = ['Contig', 'Scaffold', 'Chromosome', 'Complete Genome'];
const HOMINOIDEA = 314295;
const binomial = (name) => name.split(' ').slice(0, 2).join(' ');

export const primateSummary = (reports, taxa) => {
  const names = [];
  const nameAt = new Map();
  const species = new Set();
  const levels = Object.fromEntries(ASSEMBLY_LEVELS.map((l) => [l, 0]));
  const points = [];
  const apes = new Map();
  let total = 0;
  for (const r of reports) {
    const name = r.organism?.organism_name ?? '';
    if (r.organism?.tax_id === 9606 || name.startsWith('Homo sapiens')) continue;
    total += 1;
    species.add(binomial(name));
    const level = r.assembly_info?.assembly_level ?? '';
    if (level in levels) levels[level] += 1;
    const date = r.assembly_info?.release_date ?? '';
    const n50 = Number(r.assembly_stats?.contig_n50);
    if (!(n50 > 0) || !date) continue;
    if (!nameAt.has(name)) nameAt.set(name, names.push(name) - 1);
    points.push([date, n50, ASSEMBLY_LEVELS.indexOf(level), nameAt.get(name)]);
    const key = binomial(name);
    if (taxa.get(String(r.organism.tax_id))?.parents.includes(HOMINOIDEA) && !(apes.get(key)?.n50 >= n50)) {
      apes.set(key, {
        name: key, accession: r.accession, level, n50, date,
      });
    }
  }
  return {
    total, species: species.size, levels, names, points, apes: [...apes.values()].sort((a, b) => b.n50 - a.n50),
  };
};

const BIGBED = 0x8789F2EB;

export const readBigBed = (buf) => {
  if (buf.readUInt32LE(0) !== BIGBED) throw new Error('Not a bigBed file');
  const big = (at) => Number(buf.readBigUInt64LE(at));
  const treeAt = big(8);
  const dataAt = big(16);
  const indexAt = big(24);
  const compressed = buf.readUInt32LE(52) > 0;
  const keySize = buf.readUInt32LE(treeAt + 8);
  const chroms = [];
  const walkChroms = (at) => {
    const leaf = buf[at] === 1;
    for (let k = 0; k < buf.readUInt16LE(at + 2); k += 1) {
      const item = at + 4 + k * (keySize + 8);
      if (leaf) {
        chroms.push({
          id: buf.readUInt32LE(item + keySize), name: buf.toString('latin1', item, item + keySize).replace(/\0+$/, ''), size: buf.readUInt32LE(item + keySize + 4),
        });
      } else walkChroms(big(item + keySize));
    }
  };
  const blocks = [];
  const walkIndex = (at) => {
    const leaf = buf[at] === 1;
    for (let k = 0; k < buf.readUInt16LE(at + 2); k += 1) {
      const item = at + 4 + k * (leaf ? 32 : 24);
      if (leaf) blocks.push([big(item + 16), big(item + 24)]);
      else walkIndex(big(item + 16));
    }
  };
  walkChroms(treeAt + 32);
  walkIndex(indexAt + 48);
  blocks.sort((a, b) => a[0] - b[0]);
  const names = new Map(chroms.map((c) => [c.id, c.name]));
  return {
    fieldCount: buf.readUInt16LE(32),
    count: big(dataAt),
    chroms,
    * records() {
      for (const [offset, size] of blocks) {
        const raw = buf.subarray(offset, offset + size);
        const block = compressed ? inflateSync(raw) : raw;
        let at = 0;
        while (at < block.length) {
          const end = block.indexOf(0, at + 12);
          yield [names.get(block.readUInt32LE(at)), block.readUInt32LE(at + 4), block.readUInt32LE(at + 8), block.toString('latin1', at + 12, end)];
          at = end + 1;
        }
      }
    },
  };
};

export const gapIntervals = (text) => {
  const out = {};
  for (const line of text.split('\n')) {
    const c = line.replace(/\r$/, '').split('\t');
    if (c.length >= 4 && isStandard(c[1])) (out[bare(c[1])] ??= []).push([Number(c[2]), Number(c[3])]);
  }
  for (const chrom of Object.keys(out)) out[chrom] = mergeIntervals(out[chrom]);
  return out;
};

const STANDARD = new Set(CHROMS);

export const regionLayers = () => {
  const out = {};
  let current = null;
  let raw = null;
  const flush = () => {
    if (!current) return;
    const cds = mergeIntervals(raw.cds);
    const coding = mergeIntervals([...cds, ...raw.pcExon]);
    const exons = mergeIntervals([...coding, ...raw.exon]);
    out[current] = {
      cds, coding, exons, genic: mergeIntervals([...exons, ...raw.gene]),
    };
  };
  return {
    add(line) {
      if (!line || line.startsWith('#')) return;
      const r = gtfRecord(line);
      if (!STANDARD.has(r.chrom)) return;
      if (r.chrom !== current) {
        if (out[r.chrom]) throw new Error(`GTF lines for chr${r.chrom} are not grouped`);
        flush();
        current = r.chrom;
        raw = {
          cds: [], pcExon: [], exon: [], gene: [],
        };
      }
      const iv = [r.start - 1, r.end];
      if (r.feature === 'CDS') raw.cds.push(iv);
      else if (r.feature === 'exon') (r.attrs.gene_type === 'protein_coding' ? raw.pcExon : raw.exon).push(iv);
      else if (r.feature === 'gene') raw.gene.push(iv);
    },
    result() {
      flush();
      current = null;
      return out;
    },
  };
};

const overlapBases = (a, b) => {
  let i = 0;
  let j = 0;
  let n = 0;
  while (i < a.length && j < b.length) {
    const lo = Math.max(a[i][0], b[j][0]);
    const hi = Math.min(a[i][1], b[j][1]);
    if (hi > lo) n += hi - lo;
    if (a[i][1] < b[j][1]) i += 1;
    else j += 1;
  }
  return n;
};

const covered = (list) => list.reduce((n, [s, e]) => n + e - s, 0);
const LENGTH_EDGES = logEdges(0, 4, 10);
const LENGTH_CAP = 10000;
const NO_LAYERS = {
  cds: [], coding: [], exons: [], genic: [],
};

export const conservedSummary = ({ layers, gaps, lengths }) => {
  const classes = Object.fromEntries(['cds', 'utr', 'ncExon', 'intron', 'intergenic'].map((k) => [k, { bases: 0, conserved: 0 }]));
  const byChrom = {};
  const exact = new Uint32Array(LENGTH_CAP + 1);
  const counts = new Array(LENGTH_EDGES.length - 1).fill(0);
  const done = new Set();
  let elements = 0;
  let longest = null;
  let current = null;
  let merged = [];
  const tally = (key, bases, conserved) => {
    classes[key].bases += bases;
    classes[key].conserved += conserved;
  };
  const flush = (chrom) => {
    const l = layers[chrom] ?? NO_LAYERS;
    const g = gaps[chrom] ?? [];
    const outside = mergeIntervals([...l.genic, ...g]);
    const total = covered(merged);
    const [o1, o2, o3, o4, og] = [l.cds, l.coding, l.exons, l.genic, outside].map((layer) => overlapBases(merged, layer));
    const [c1, c2, c3, c4] = [l.cds, l.coding, l.exons, l.genic].map(covered);
    const size = lengths[chrom] ?? 0;
    tally('cds', c1, o1);
    tally('utr', c2 - c1, o2 - o1);
    tally('ncExon', c3 - c2, o3 - o2);
    tally('intron', c4 - c3, o4 - o3);
    tally('intergenic', size - covered(outside), total - og);
    byChrom[chrom] = { length: size - covered(g), conserved: total };
    done.add(chrom);
    merged = [];
  };
  return {
    add(chrom, start, end) {
      if (chrom !== current) {
        if (current !== null) flush(current);
        if (done.has(chrom)) throw new Error(`Conserved elements on chr${chrom} are not grouped`);
        current = chrom;
      }
      const last = merged.at(-1);
      if (last && start < last[0]) throw new Error(`Conserved elements on chr${chrom} are not sorted`);
      elements += 1;
      const length = end - start;
      exact[Math.min(length, LENGTH_CAP)] += 1;
      counts[binIndex(LENGTH_EDGES, length)] += 1;
      if (!longest || length > longest.length) {
        longest = {
          chrom, start: start + 1, end, length,
        };
      }
      if (last && start <= last[1]) last[1] = Math.max(last[1], end);
      else merged.push([start, end]);
    },
    result() {
      if (current !== null) flush(current);
      current = null;
      for (const chrom of Object.keys(lengths)) if (!done.has(chrom)) flush(chrom);
      let median = null;
      for (let len = 0, seen = 0; len <= LENGTH_CAP && median === null; len += 1) {
        seen += exact[len];
        if (elements && seen >= elements / 2) median = len;
      }
      return {
        elements,
        bases: Object.values(byChrom).reduce((n, c) => n + c.conserved, 0),
        byChrom,
        classes,
        median,
        longest,
        lengths: { edges: LENGTH_EDGES, counts },
      };
    },
  };
};
