import {
  exact, pct, codonOf, perLitre,
} from './lib.js';

const nuclear = (list) => list.filter((c) => c.chrom !== 'MT');
const sum = (list, key) => list.reduce((n, x) => n + (x[key] ?? 0), 0);
const role = (asm, name) => asm.roles.find((r) => r.role === name);

export const assemblyNumbers = (snap) => {
  const { grch38, t2t } = snap.assembly;
  const g = nuclear(grch38.chromosomes);
  const t = nuclear(t2t.chromosomes);
  const grch38Length = sum(g, 'length');
  const grch38Ungapped = sum(g, 'ungapped');
  const t2tLength = sum(t, 'length');
  const perChrom = t.map((c) => ({ chrom: c.chrom, bases: c.length - g.find((x) => x.chrom === c.chrom).ungapped }));
  return {
    grch38Length,
    grch38Ungapped,
    grch38N: grch38Length - grch38Ungapped,
    t2tLength,
    added: t2tLength - grch38Ungapped,
    patches: role(grch38, 'fix-patch').count + role(grch38, 'novel-patch').count,
    mostAdded: perChrom.reduce((a, b) => (b.bases > a.bases ? b : a)),
    perChrom,
  };
};

const chromosomes = (snap) => {
  const list = nuclear(snap.assembly.grch38.chromosomes);
  const longest = list.reduce((a, b) => (b.length > a.length ? b : a));
  const shortest = list.reduce((a, b) => (b.length < a.length ? b : a));
  return `Chromosome ${longest.chrom} spans ${exact(longest.length)} bp in GRCh38; chromosome ${shortest.chrom}, the shortest, ${exact(shortest.length)}`;
};

const gaps = (snap) => {
  const n = assemblyNumbers(snap);
  return `GRCh38 leaves ${exact(n.grch38N)} bases of its chromosomes as N, ${pct(n.grch38N / n.grch38Length)} of their length`;
};

const scaffolds = (snap) => {
  const { grch38 } = snap.assembly;
  return `Besides its chromosomes, ${grch38.name} holds ${exact(role(grch38, 'alt-scaffold').count)} alternate loci in ${grch38.altLoci} regions and ${exact(assemblyNumbers(snap).patches)} patches`;
};

const t2t = (snap) => {
  const n = assemblyNumbers(snap);
  return `T2T-CHM13 sequences ${exact(n.added)} more bases than GRCh38 does; ${exact(n.mostAdded.bases)} of them are on chromosome ${n.mostAdded.chrom}`;
};

const gcOf = (c) => (c.C + c.G) / (c.A + c.C + c.G + c.T);

export const sequenceNumbers = (snap) => {
  const list = nuclear(snap.sequence.chromosomes);
  const acgt = list.reduce((n, c) => n + c.A + c.C + c.G + c.T, 0);
  const C = sum(list, 'C');
  const G = sum(list, 'G');
  const ranked = list.map((c) => ({ chrom: c.chrom, gc: gcOf(c) })).sort((a, b) => b.gc - a.gc);
  return {
    acgt,
    gc: (C + G) / acgt,
    cpgOe: sum(list, 'cpg') / ((C * G) / acgt),
    g4: list.reduce((n, c) => n + c.g4.plus + c.g4.minus, 0),
    richest: ranked[0],
    poorest: ranked.at(-1),
  };
};

const composition = (snap) => {
  const n = sequenceNumbers(snap);
  return `GRCh38 is ${pct(n.gc)} G and C; chromosome ${n.richest.chrom} is the richest at ${pct(n.richest.gc)}, chromosome ${n.poorest.chrom} the poorest at ${pct(n.poorest.gc)}`;
};

const cpg = (snap) => `CpG turns up at ${pct(sequenceNumbers(snap).cpgOe)} of the rate the C and G content predicts`;

const motifs = (snap) => {
  const runs = Object.entries(snap.sequence.homopolymers).map(([base, r]) => ({ base, ...r }));
  const longest = runs.reduce((a, b) => (b.length > a.length ? b : a));
  return `${exact(sequenceNumbers(snap).g4)} G-quadruplex motifs; the longest run of one letter is ${longest.length} ${longest.base}s on chromosome ${longest.chrom}`;
};

const transposons = (snap) => {
  const { acgt } = sequenceNumbers(snap);
  const { rmsk } = snap.repeats;
  const l1 = rmsk.families.find((f) => f.family === 'L1');
  return `RepeatMasker marks ${pct(rmsk.total.bp / acgt)} of GRCh38’s sequenced bases as repeats; L1 elements alone cover ${pct(l1.bp / acgt)}`;
};

const alu = (snap) => {
  const { alu: a, l1hs, families } = snap.repeats.rmsk;
  const total = families.find((f) => f.family === 'Alu').count;
  return `${exact(total)} Alu copies, ${exact(a.AluY.count)} in the youngest subfamily, AluY; ${exact(l1hs.count)} L1HS copies, ${exact(l1hs.fullLength)} of them at least 6 kb long`;
};

const segdups = (snap) => {
  const { superDups } = snap.repeats;
  const covered = Object.values(superDups.covered).reduce((a, b) => a + b, 0);
  return `Segmental duplications cover ${exact(covered)} bases, ${pct(covered / sequenceNumbers(snap).acgt)} of GRCh38; ${exact(superDups.kinds.inter)} of ${exact(superDups.pairs)} pairs join different chromosomes`;
};

const tandem = (snap) => {
  const { units, motifs: top } = snap.repeats.tandem;
  return `Tandem Repeats Finder calls ${exact(units.reduce((n, u) => n + u.count, 0))} tandem repeats; the most frequent units are ${top[0].motif}, then ${top[1].motif}`;
};

export const codingCounts = (snap) => [
  { source: 'GENCODE', n: snap.genes.genes.byType.protein_coding },
  { source: 'NCBI Gene', n: snap.catalogs.geneInfo.types['protein-coding'] },
  { source: 'HGNC', n: snap.catalogs.hgnc.locusGroups['protein-coding gene'] },
  { source: 'MANE Select', n: snap.catalogs.mane.select },
];

const genecount = (snap) => {
  const list = codingCounts(snap);
  const lo = list.reduce((a, b) => (b.n < a.n ? b : a));
  const hi = list.reduce((a, b) => (b.n > a.n ? b : a));
  return `Protein-coding gene counts run from ${exact(lo.n)} (${lo.source}) to ${exact(hi.n)} (${hi.source})`;
};

const perMb = (v) => (v >= 100 ? exact(v) : v.toFixed(1));

export const geneDensity = (snap) => nuclear(snap.sequence.chromosomes)
  .filter((c) => snap.genes.genes.byChrom[c.chrom])
  .map((c) => ({ chrom: c.chrom, genes: snap.genes.genes.byChrom[c.chrom].protein_coding ?? 0, perMb: ((snap.genes.genes.byChrom[c.chrom].protein_coding ?? 0) / (c.A + c.C + c.G + c.T)) * 1e6 }));

const density = (snap) => {
  const list = geneDensity(snap);
  const dense = list.reduce((a, b) => (b.perMb > a.perMb ? b : a));
  const sparse = list.reduce((a, b) => (b.perMb < a.perMb ? b : a));
  return `Chromosome ${dense.chrom} carries ${perMb(dense.perMb)} protein-coding genes per sequenced megabase; chromosome ${sparse.chrom}, ${perMb(sparse.perMb)}`;
};

const genesize = (snap) => {
  const { spans } = snap.genes;
  return `The median protein-coding gene spans ${exact(spans.pc.summary.median)} bp; ${spans.longest[0].name}, the longest, ${exact(spans.longest[0].span)}`;
};

const exons = (snap) => {
  const { mane } = snap.genes;
  return `A MANE Select transcript has a median of ${mane.exons.summary.median} exons; ${mane.exons.top[0].name} has ${mane.exons.top[0].exons}; ${exact(mane.singleExon)} have one`;
};

const territory = (snap) => {
  const t = snap.genes.territory;
  const { acgt } = sequenceNumbers(snap);
  return `Coding sequence covers ${pct(t.cds / acgt)} of GRCh38, exons of all genes ${pct((t.cds + t.utr + t.ncExon) / acgt)} and introns ${pct(t.intron / acgt)}`;
};

const isoforms = (snap) => {
  const iso = snap.genes.isoforms;
  return `Protein-coding genes have a median of ${iso.summary.median} GENCODE transcripts; ${iso.top[0].name} has ${exact(iso.top[0].n)}`;
};

const names = (snap) => {
  const { hgnc } = snap.catalogs;
  return `HGNC has approved ${exact(hgnc.approved)} gene symbols; ${hgnc.spreadsheet.count} were changed so spreadsheets would stop reading them as dates`;
};

const pseudogenes = (snap) => {
  const { hgnc } = snap.catalogs;
  const [a, b] = hgnc.pseudoParents;
  return `${exact(hgnc.locusGroups.pseudogene)} named pseudogenes; ${a.parent} has ${exact(a.n)} copies, ${b.parent} ${exact(b.n)}`;
};

const immune = (snap) => {
  const { V, D, J } = snap.genes.immune.IGH;
  return `${V} IGHV, ${D} IGHD and ${J} IGHJ segments make ${exact(V * D * J)} heavy-chain V-D-J combinations before junctional diversity`;
};

const splicing = (snap) => {
  const c = snap.sequence.mane.classes;
  const total = Object.values(c).reduce((a, b) => a + b, 0);
  return `${pct(c['GT-AG'] / total)} of MANE Select introns start GT and end AG; ${exact(c['GC-AG'])} start GC and ${exact(c['AT-AC'])} ${c['AT-AC'] === 1 ? 'runs' : 'run'} from AT to AC`;
};

const translation = (snap) => {
  const { kozak, stops, coding } = snap.sequence.mane;
  const [a, c, g, t] = kozak[3];
  return `A purine sits three bases before the start codon in ${pct((a + g) / (a + c + g + t))} of MANE coding sequences; ${pct(stops.TGA / coding)} end in TGA`;
};

const STOP_CODONS = new Set([48, 50, 56]);

export const senseCodons = (counts) => counts.map((n, i) => ({ codon: codonOf(i), n })).filter((_, i) => !STOP_CODONS.has(i));

const codons = (snap) => {
  const sense = senseCodons(snap.sequence.mane.codons);
  const total = sense.reduce((s, c) => s + c.n, 0);
  const top = sense.reduce((a, b) => (b.n > a.n ? b : a));
  const rare = sense.reduce((a, b) => (b.n < a.n ? b : a));
  const perK = (n) => ((n / total) * 1000).toFixed(1);
  return `${top.codon} is the most used codon, ${perK(top.n)} per thousand; ${rare.codon} the least used sense codon, ${perK(rare.n)}`;
};

export const proteinLabel = (p) => {
  const genes = String(p.gene ?? '').split(/;\s*/).filter(Boolean);
  if (!genes.length) return p.acc;
  return genes.length === 1 ? genes[0] : `${genes[0]} and ${genes.length - 1} more`;
};

const proteome = (snap) => {
  const P = snap.proteins;
  const [top] = P.longest;
  return `UniProt ${P.release} lists ${exact(P.count)} reviewed human proteins; the median is ${exact(P.lengths.summary.median)} residues long and ${proteinLabel(top)}, the longest, ${exact(top.length)}`;
};

const structures = (snap) => {
  const { count, pdb } = snap.proteins;
  return `${exact(pdb.withAny)} reviewed human proteins, ${pct(pdb.withAny / count)}, have at least one PDB structure; ${proteinLabel(pdb.top[0])} appears in ${exact(pdb.top[0].n)}`;
};

const membrane = (snap) => {
  const { transmem, signal } = snap.proteins;
  return `${exact(transmem.proteins)} proteins cross a membrane: ${exact(transmem.counts[1])} once and ${exact(transmem.counts[7])} seven times; ${exact(signal)} carry a signal peptide`;
};

const domains = (snap) => {
  const { count, pfam } = snap.proteins;
  const [top] = pfam.top;
  return `${pct(pfam.withAny / count)} of proteins carry a Pfam domain; ${exact(top.n)} carry ${top.short} (${top.id}), the most common`;
};

export const proteomeResidues = (snap) => snap.proteins.lengths.summary.mean * snap.proteins.lengths.summary.n;

const modifications = (snap) => {
  const { disorder, ptm } = snap.proteins;
  return `${pct(disorder.residues / proteomeResidues(snap))} of residues sit in disordered regions; UniProt lists ${exact(ptm.phosphoS)} phosphoserines, ${exact(ptm.phosphoT)} phosphothreonines and ${exact(ptm.phosphoY)} phosphotyrosines`;
};

const targets = (snap) => {
  const { count, levels } = snap.targets;
  return `Pharos counts ${exact(levels.Tclin)} targets of approved drugs (Tclin) among ${exact(count)} human proteins; ${exact(levels.Tdark)} are Tdark`;
};

const specificity = (snap) => {
  const { hpa } = snap.expression;
  const [a, b] = hpa.enriched;
  return `HPA calls ${exact(hpa.specificity['Tissue enriched'])} genes tissue enriched; ${a.tissue} has ${exact(a.n)} of them, ${b.tissue} ${exact(b.n)}`;
};

export const gtexTissues = (snap) => {
  const { gtex } = snap.expression;
  return gtex.tissues.map((id, i) => ({
    ...(gtex.sites.find((s) => s.id === id) ?? { id, name: id }),
    expressed: gtex.expressed[i],
    specific: gtex.specific?.[i] ?? null,
    index: i,
  }));
};

const gtex = (snap) => {
  const list = gtexTissues(snap);
  const most = list.reduce((a, b) => (b.expressed > a.expressed ? b : a));
  const least = list.reduce((a, b) => (b.expressed < a.expressed ? b : a));
  return `In GTEx ${snap.expression.gtex.release}, ${exact(most.expressed)} genes reach a median of 1 TPM in ${most.name}; ${exact(least.expressed)} do in ${least.name}`;
};

export const plasmaLevel = (p, pick = Math.max) => pick(...[p.im, p.ms].filter((v) => v !== null && v !== undefined));

const plasma = (snap) => {
  const { blood } = snap.expression.hpa;
  const [top] = blood.top;
  return `Plasma levels in HPA run from ${perLitre(plasmaLevel(top))} for ${top.gene} to ${perLitre(plasmaLevel(blood.lowest, Math.min))} for ${blood.lowest.gene}, across ${exact(blood.proteins)} proteins`;
};

const interactions = (snap) => {
  const { version, edges, proteins, hubs } = snap.networks;
  return `STRING ${version} links ${exact(proteins.high)} human proteins by ${exact(edges.high)} interactions scored 0.7 or more; ${hubs[0].name} has the most partners, ${exact(hubs[0].degree)}`;
};

const pathways = (snap) => {
  const R = snap.pathways;
  const [top] = R.topLevel;
  return `Reactome ${R.release} places ${exact(R.proteins)} human proteins in ${exact(R.pathways)} pathways; ${top.name}, the largest of ${R.topLevel.length} top-level pathways, holds ${exact(top.proteins)}`;
};

export const goTotals = (snap) => {
  const groups = Object.values(snap.ontology.byAspect);
  const total = groups.reduce((n, g) => n + Object.values(g).reduce((a, b) => a + b, 0), 0);
  return { total, electronic: groups.reduce((n, g) => n + g.electronic, 0) };
};

const ontology = (snap) => {
  const { total, electronic } = goTotals(snap);
  return `${exact(snap.ontology.proteins)} human proteins carry ${exact(total)} GO annotations; ${pct(electronic / total)} are electronic inferences (IEA)`;
};

export const findings = (snap) => ({
  interactions: interactions(snap),
  pathways: pathways(snap),
  ontology: ontology(snap),
  proteome: proteome(snap),
  structures: structures(snap),
  membrane: membrane(snap),
  domains: domains(snap),
  modifications: modifications(snap),
  targets: targets(snap),
  specificity: specificity(snap),
  gtex: gtex(snap),
  plasma: plasma(snap),
  splicing: splicing(snap),
  translation: translation(snap),
  codons: codons(snap),
  genecount: genecount(snap),
  density: density(snap),
  genesize: genesize(snap),
  exons: exons(snap),
  territory: territory(snap),
  isoforms: isoforms(snap),
  names: names(snap),
  pseudogenes: pseudogenes(snap),
  immune: immune(snap),
  transposons: transposons(snap),
  alu: alu(snap),
  segdups: segdups(snap),
  tandem: tandem(snap),
  chromosomes: chromosomes(snap),
  gaps: gaps(snap),
  scaffolds: scaffolds(snap),
  t2t: t2t(snap),
  composition: composition(snap),
  cpg: cpg(snap),
  motifs: motifs(snap),
});
