import { exact, pct } from './lib.js';

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

export const findings = (snap) => ({
  chromosomes: chromosomes(snap),
  gaps: gaps(snap),
  scaffolds: scaffolds(snap),
  t2t: t2t(snap),
  composition: composition(snap),
  cpg: cpg(snap),
  motifs: motifs(snap),
});
