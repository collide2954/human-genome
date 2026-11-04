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

export const findings = (snap) => ({
  chromosomes: chromosomes(snap),
  gaps: gaps(snap),
  scaffolds: scaffolds(snap),
  t2t: t2t(snap),
});
