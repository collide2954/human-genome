import { quantile, histogram, linearEdges } from '../src/stats.js';
import { round } from './tables.mjs';

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
