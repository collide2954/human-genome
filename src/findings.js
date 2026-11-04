import {
  exact, pct, codonOf, perLitre, GENETIC_CODE, AA3, bases, longDate,
} from './lib.js';
import { reverseComplement } from './stats.js';
import {
  GENE_THERAPIES, DIAGNOSTIC_YIELD, CFDNA_SCREENING, GNOMAD_V4,
} from './facts.js';

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

export const trnaPairs = (counts, anticodons) => {
  const sense = senseCodons(counts);
  const total = sense.reduce((s, c) => s + c.n, 0);
  return sense.map(({ codon, n }) => {
    const aa = AA3[GENETIC_CODE[codon]];
    const anticodon = reverseComplement(codon);
    return {
      codon, aa, anticodon, perK: (n / total) * 1000, copies: anticodons[`${aa}-${anticodon}`] ?? 0,
    };
  });
};

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
  const { proteins, types } = snap.ontology;
  const others = Object.values(types).reduce((a, b) => a + b, 0) - proteins;
  return `GO holds ${exact(total)} annotations of ${exact(proteins)} human proteins and ${exact(others)} other gene products; ${pct(electronic / total)} are electronic (IEA)`;
};

const ccres = (snap) => {
  const { count, classes } = snap.regulation.ccre;
  const bp = Object.values(classes).reduce((n, c) => n + c.bp, 0);
  return `ENCODE’s registry lists ${exact(count)} candidate cis-regulatory elements covering ${pct(bp / sequenceNumbers(snap).acgt)} of GRCh38; ${exact(classes.dELS.count)} are distal enhancer-like`;
};

const cpgislands = (snap) => {
  const { count, bp, promoters } = snap.regulation.cpg;
  return `${exact(count)} CpG islands cover ${pct(bp / sequenceNumbers(snap).acgt)} of GRCh38; ${pct(promoters.inIsland / promoters.tss)} of MANE Select transcripts start inside one`;
};

const trna = (snap) => {
  const pairs = trnaPairs(snap.sequence.mane.codons, snap.trnas.anticodons);
  return `GtRNAdb lists ${exact(snap.trnas.genes)} high-confidence tRNA genes; ${pairs.filter((p) => !p.copies).length} of the ${pairs.length} sense codons have none with an exactly pairing anticodon`;
};

const encode = (snap) => {
  const { total, assays } = snap.encode;
  return `ENCODE has released ${exact(total)} human experiments; ${assays[0].name} is the most common assay, with ${exact(assays[0].n)}`;
};

const mito = (snap) => {
  const { genes, variants } = snap.structure.mito;
  const kinds = (k) => genes.filter((g) => g.kind === k).length;
  const length = snap.sequence.chromosomes.find((c) => c.chrom === 'MT').length;
  return `The mitochondrial genome’s ${exact(length)} bp hold ${genes.length} genes, ${kinds('protein')} for proteins, ${kinds('tRNA')} tRNAs and ${kinds('rRNA')} rRNAs; gnomAD finds variants at ${exact(variants.positions)} positions`;
};

const censat = (snap) => {
  const { total, classes, arrays } = snap.structure.censat;
  const top = arrays.hor.reduce((a, b) => (b.bp > a.bp ? b : a));
  return `T2T-CHM13 maps ${bases(total)} of satellite arrays; active α-satellite HORs make up ${bases(classes.hor)}, the most on chromosome ${top.chrom} (${bases(top.bp)})`;
};

const sexchrom = (snap) => {
  const by = snap.genes.genes.byChrom;
  const y = snap.assembly.t2t.chromosomes.find((c) => c.chrom === 'Y').length;
  const sat = Object.values(snap.structure.censat.byChrom.Y ?? {}).reduce((a, b) => a + b, 0);
  return `Chromosome X carries ${exact(by.X.protein_coding)} protein-coding genes and Y ${exact(by.Y.protein_coding)}; in T2T-CHM13 satellite arrays make up ${pct(sat / y)} of Y`;
};

const dbsnp = (snap) => {
  const d = snap.variation.dbsnp;
  return `dbSNP build ${d.build} holds ${exact(d.total)} reference SNP records from ${exact(d.ss.total)} submitted variants; ${exact(d.status['Live RS'])} are live`;
};

const sv = (snap) => {
  const c = snap.variation.sv.classes;
  return `dbVar’s non-redundant sets hold ${exact(c.deletion)} deletions, ${exact(c.insertion)} insertions and ${exact(c.duplication)} duplications on GRCh38`;
};

const constraint = (snap) => {
  const { gnomad, byDosage } = snap.constraint;
  return `${exact(gnomad.below06)} of ${exact(gnomad.genes)} genes in gnomAD v4.1 have a LOEUF below 0.6; ClinGen haploinsufficient genes have a median of ${byDosage.sufficient.median.toFixed(2)}`;
};

const dosage = (snap) => {
  const { hi, ts } = snap.constraint.dosage;
  return `ClinGen finds sufficient evidence of haploinsufficiency for ${exact(hi[3] ?? 0)} genes and of triplosensitivity for ${exact(ts[3] ?? 0)}`;
};

const PLP = ['Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic'];
const total = (list) => list.reduce((a, b) => a + b, 0);

export const titv = (channels) => {
  const ti = total(channels.slice(32, 48)) + total(channels.slice(64, 80));
  return ti / (total(channels) - ti);
};

const clinvar = (snap) => {
  const { records, classes } = snap.clinvar;
  const plp = total(PLP.map((k) => classes[k] ?? 0));
  return `ClinVar’s GRCh38 file lists ${exact(records)} variants; ${pct(classes['Uncertain significance'] / records)} are of uncertain significance, ${pct(plp / records)} pathogenic or likely pathogenic`;
};

const review = (snap) => {
  const panel = total(Object.values(snap.clinvar.stars).map((s) => s[3]));
  const groups = Object.entries(snap.clinvarTable.submitters).filter(([k]) => k !== 'bins').map(([, v]) => v);
  return `Expert panels have reviewed ${exact(panel)} ClinVar variants; ${pct(total(groups.map((g) => g[0])) / total(groups.map(total)))} of classified variants rest on one submitter`;
};

const growth = (snap) => {
  const { points } = snap.clinvarGrowth;
  const [first, last] = [points[0], points.at(-1)];
  const fold = last.classes['Uncertain significance'] / first.classes['Uncertain significance'];
  return `ClinVar grew from ${exact(first.records)} variants in ${longDate(first.release.slice(0, 7))} to ${exact(last.records)}; those of uncertain significance rose ${fold.toFixed(1)}-fold`;
};

const kindShare = (lollipop, kind) => total(lollipop.up.filter((r) => r[1] === kind).map((r) => r[2])) / total(lollipop.up.map((r) => r[2]));

const lollipops = (snap) => {
  const l = snap.clinvarTable.lollipops;
  return `${pct(kindShare(l.BRCA1, 'truncating'))} of BRCA1’s P/LP variants placed on the protein truncate it; in MYH7 ${pct(kindShare(l.MYH7, 'missense'))} are missense`;
};

const spectrum = (snap) => {
  const s = snap.clinvar.spectrum;
  return `P/LP SNVs in ClinVar carry ${titv(s.plp).toFixed(2)} transitions per transversion; B/LB SNVs carry ${titv(s.blb).toFixed(2)}`;
};

const validity = (snap) => {
  const v = snap.validity;
  return `ClinGen has curated ${exact(v.curations)} gene–disease pairs in ${exact(v.genes)} genes; ${exact(v.classes.Definitive)} are definitive, ${exact(v.classes.Disputed + v.classes.Refuted)} disputed or refuted`;
};

const cnvregions = (snap) => {
  const r = snap.regions;
  return `Of ${exact(r.regions)} genomic regions ClinGen has curated, ${exact(r.hi[3] ?? 0)} are haploinsufficient and ${exact(r.ts[3] ?? 0)} triplosensitive`;
};

const expansions = (snap) => {
  const e = snap.expansions;
  return `STRchive lists ${exact(e.loci.length)} disease-linked tandem repeat loci, ${exact(e.evidence.Definitive ?? 0)} with definitive evidence; ${exact(e.loci.filter((l) => l.year >= 2015).length)} were found since 2015`;
};

const hpo = (snap) => {
  const h = snap.hpo;
  const top = h.systems[0];
  return `HPO annotates ${exact(h.diseases)} rare diseases, a median of ${h.perDisease.median} phenotype terms each; ${pct(top.diseases / h.diseases)} involve ${/^Abnormality of /.test(top.name) ? top.name.replace(/^Abnormality of /, '') : top.name.toLowerCase()}`;
};

const orphanet = (snap) => {
  const { prevalence, ages } = snap.orphanet;
  return `${exact(prevalence.classes['<1 / 1 000 000'] ?? 0)} of ${exact(prevalence.disorders)} Orphanet disorders affect fewer than 1 in a million people; ${exact(ages.onset.Neonatal ?? 0)} can start in the newborn period`;
};

const gtr = (snap) => {
  const g = snap.gtr;
  return `The Genetic Testing Registry lists ${exact(g.tests)} current tests from ${exact(g.labs)} laboratories; ${pct(g.genesPerTest.counts[1] / g.tests)} look at a single gene`;
};

const therapy = (snap) => {
  const { fda } = snap;
  const known = new Map(GENE_THERAPIES.map((g) => [g.trade, g]));
  const present = fda.filter((p) => known.has(p.trade));
  return `FDA lists ${exact(fda.length)} approved cell and gene therapy products; ${exact(present.length)} deliver or modify genes, ${exact(present.filter((p) => known.get(p.trade).group === 'AAV').length)} of them with AAV vectors`;
};

const diagnosticYield = () => {
  const of = (key) => DIAGNOSTIC_YIELD.find((y) => y.key === key).value;
  return `Genome sequencing diagnosed ${of('wgs')}% of children with a suspected genetic disease, exome ${of('wes')}% and microarray ${of('cma')}%`;
};

const prenatal = () => {
  const t21 = CFDNA_SCREENING.find((s) => s.condition === 'Trisomy 21');
  return `Cell-free DNA screening detects ${t21.detection}% of trisomy 21 pregnancies with a false-positive rate of ${t21.fpr}%`;
};

const ancestry = (snap) => {
  const initial = snap.ancestry.stages.initial;
  const all = total(Object.values(initial));
  return `${pct(initial.European / all)} of GWAS discovery participants in the GWAS Catalog are listed as European; ${pct((initial['African American or Afro-Caribbean'] ?? 0) / all)} as African American or Afro-Caribbean`;
};

const ancient = (snap) => {
  const a = snap.ancient;
  return `The AADR ${a.version} holds ${exact(a.individuals)} ancient individuals; the oldest, ${a.oldest.id} from ${a.oldest.place}, dates to ${exact(a.oldest.age)} years before present`;
};

const panels = (snap) => {
  const { panel } = snap.igsr;
  return `The 1000 Genomes 30x set holds ${exact(panel.samples)} genomes from ${exact(panel.populations.length)} populations; gnomAD v4 pools ${exact(GNOMAD_V4.people)} people`;
};

const assemblies = (snap) => {
  const a = snap.assemblies;
  const year = Object.keys(a.byYear).sort().at(-1);
  return `NCBI lists ${exact(a.total)} human genome assemblies, ${exact(total(Object.values(a.byYear[year])))} released in ${year}; the longest contig N50 is ${bases(a.best.n50)}`;
};

const perTumour = (v) => (v < 10 ? v.toFixed(1) : exact(v));
const extremes = (shares) => {
  const entries = Object.entries(shares);
  const [minGroup, min] = entries.reduce((a, b) => (b[1] < a[1] ? b : a));
  const [maxGroup, max] = entries.reduce((a, b) => (b[1] > a[1] ? b : a));
  return {
    minGroup, min, maxGroup, max,
  };
};

export const tcgaShares = (snap) => {
  const { genes, projects, mutated } = snap.gdc.tcga;
  return genes.flatMap((gene) => projects.filter((p) => p.tested > 0).map((p) => {
    const n = mutated[gene]?.[p.project] ?? 0;
    return {
      gene, project: p.project, name: p.name, n, tested: p.tested, share: n / p.tested,
    };
  }));
};

export const exomeBurden = (snap) => snap.intogen.cohorts.burden.map((row) => {
  const exomes = row.cohorts.filter((c) => c.platform === 'WXS');
  const samples = exomes.reduce((n, c) => n + c.samples, 0);
  return {
    cancer: row.cancer, name: row.name, cohorts: exomes.length, samples, perTumour: samples ? exomes.reduce((n, c) => n + c.perSample * c.samples, 0) / samples : null,
  };
}).filter((r) => r.perTumour !== null);

const drivers = (snap) => {
  const d = snap.intogen.drivers;
  return `IntOGen calls ${exact(d.genes)} genes cancer drivers across ${exact(d.cancerTypes)} cancer types; ${d.top[0].gene} is a driver in ${exact(d.top[0].types)} of them`;
};

const tcga = (snap) => {
  const cells = tcgaShares(snap);
  const best = cells.reduce((a, b) => (b.share > a.share ? b : a));
  return `${best.gene} has a somatic mutation in ${pct(best.share)} of TCGA ${best.name} cases, the highest of ${snap.gdc.tcga.genes.length} drivers across ${new Set(cells.map((c) => c.project)).size} projects`;
};

const burden = (snap) => {
  const rows = exomeBurden(snap);
  const low = rows.reduce((a, b) => (b.perTumour < a.perTumour ? b : a));
  const high = rows.reduce((a, b) => (b.perTumour > a.perTumour ? b : a));
  return `Exome cohorts in IntOGen average ${perTumour(low.perTumour)} mutations per tumour in ${low.name} and ${perTumour(high.perTumour)} in ${high.name}`;
};

const somatic = (snap) => `The GDC holds ${exact(snap.gdc.cases)} cases from ${exact(snap.gdc.projects)} projects; ClinVar lists ${exact(total(Object.values(snap.clinvar.oncogenicity)))} somatic oncogenicity classifications`;

const cpic = (snap) => {
  const { levels, byGene } = snap.cpic.pairs;
  const withA = byGene.filter((g) => g.levels.A);
  const top = withA.reduce((a, b) => (b.levels.A > a.levels.A ? b : a));
  return `CPIC rates ${exact(levels.A)} gene–drug pairs level A across ${exact(withA.length)} genes; ${top.gene} has the most, ${exact(top.levels.A)}`;
};

const phenotypes = (snap) => {
  const widest = snap.cpic.phenotypes.genes.map((g) => ({ gene: g.gene, ...extremes(g.shares) })).reduce((a, b) => (b.max - b.min > a.max - a.min ? b : a));
  return `Predicted high-risk ${widest.gene} phenotypes run from ${pct(widest.min)} of ${widest.minGroup} to ${pct(widest.max)} of ${widest.maxGroup} people in CPIC’s tables`;
};

const cyp2d6 = (snap) => {
  const c = snap.cpic.cyp2d6;
  const poor = Object.fromEntries(Object.entries(c.phenotypes).filter(([, p]) => p['Poor Metabolizer'] !== undefined).map(([group, p]) => [group, p['Poor Metabolizer']]));
  const e = extremes(poor);
  return `CPIC lists ${exact(c.alleles)} CYP2D6 alleles; predicted poor metabolizers run from ${pct(e.min)} of ${e.minGroup} to ${pct(e.max)} of ${e.maxGroup} people`;
};

const fdalabels = (snap) => {
  const f = snap.fdaPgx;
  const oncology = f.areas.find((a) => a.area === 'Oncology')?.rows ?? 0;
  return `FDA’s biomarker table pairs ${exact(f.drugs)} drugs with ${exact(f.biomarkers)} biomarkers in ${exact(f.rows)} rows; ${pct(oncology / f.rows)} are in oncology`;
};

const clinpgx = (snap) => {
  const a = snap.clinpgx.annotations;
  return `ClinPGx holds ${exact(a.annotations)} clinical annotations; ${exact(a.levels['1A'] ?? 0)} reach level 1A and ${exact(a.levels[3] ?? 0)} sit at level 3`;
};

const sentence = (s) => `${s[0].toUpperCase()}${s.slice(1)}`;
const inline = (s) => (/^[A-Z][a-z]/.test(s) ? `${s[0].toLowerCase()}${s.slice(1)}` : s);
const shareOf = (counts, key) => counts[key] / total(Object.values(counts));

const manhattan = (snap) => {
  const g = snap.gwas;
  return `The GWAS Catalog holds ${exact(g.associations)} associations from ${exact(g.studies)} studies; ${pct(g.significant / g.associations)} reach p < 5×10⁻⁸`;
};

const gwasgrowth = (snap) => {
  const years = Object.entries(snap.gwas.byYear);
  const [top, peak] = years.reduce((a, b) => (b[1].associations > a[1].associations ? b : a));
  const [base, then] = years.find(([y]) => Number(y) === Number(top) - 10) ?? years[0];
  return `Studies published in ${top} report ${exact(peak.associations)} catalog associations, ${Math.round(peak.associations / then.associations)} times as many as those published in ${base}`;
};

const gwastraits = (snap) => {
  const [trait] = snap.gwas.topTraits;
  const [gene] = snap.gwas.topGenes;
  return `${sentence(trait.trait)} has the most catalog associations, ${exact(trait.n)}; ${gene.gene} is the gene mapped to the most, ${exact(gene.n)}`;
};

const effects = (snap) => {
  const e = snap.gwas.effects;
  const above = total(e.counts.slice(4).flat());
  return `Among ${exact(e.n)} disease associations with odds ratios, the median OR is ${e.median.toFixed(2)}; ${pct(above / e.n)} exceed 2`;
};

const pgs = (snap) => {
  const p = snap.pgs;
  return `The PGS Catalog holds ${exact(p.catalog.scores)} scores for ${exact(p.catalog.traits)} traits; the median score uses ${exact(p.variants.median)} variants`;
};

const pgsancestry = (snap) => {
  const { gwas, evaluation } = snap.pgs.ancestry;
  return `${pct(shareOf(gwas, 'European'))} of individuals in the GWAS behind polygenic scores are listed as European; ${pct(shareOf(evaluation, 'European'))} of those in evaluations`;
};

const pgsperformance = (snap) => {
  const { auroc, r2 } = snap.pgs.performance;
  return `Polygenic scores report a median AUROC of ${auroc.median.toFixed(3)} across ${exact(auroc.n)} evaluations; the median R² is ${r2.median.toFixed(2)}`;
};

const heritability = (snap) => {
  const { eur } = snap.panukb;
  return `Median SNP heritability of ${exact(eur.passing)} UK Biobank traits in European-ancestry people is ${eur.median.toFixed(2)}; ${inline(eur.top[0].description)} reaches ${eur.top[0].h2.toFixed(2)}`;
};

export const findings = (snap) => ({
  manhattan: manhattan(snap),
  gwasgrowth: gwasgrowth(snap),
  gwastraits: gwastraits(snap),
  effects: effects(snap),
  pgs: pgs(snap),
  pgsancestry: pgsancestry(snap),
  pgsperformance: pgsperformance(snap),
  heritability: heritability(snap),
  drivers: drivers(snap),
  tcga: tcga(snap),
  burden: burden(snap),
  somatic: somatic(snap),
  cpic: cpic(snap),
  phenotypes: phenotypes(snap),
  cyp2d6: cyp2d6(snap),
  fdalabels: fdalabels(snap),
  clinpgx: clinpgx(snap),
  ancestry: ancestry(snap),
  ancient: ancient(snap),
  panels: panels(snap),
  assemblies: assemblies(snap),
  hpo: hpo(snap),
  orphanet: orphanet(snap),
  gtr: gtr(snap),
  therapy: therapy(snap),
  yield: diagnosticYield(),
  prenatal: prenatal(),
  clinvar: clinvar(snap),
  review: review(snap),
  growth: growth(snap),
  lollipops: lollipops(snap),
  spectrum: spectrum(snap),
  validity: validity(snap),
  cnvregions: cnvregions(snap),
  expansions: expansions(snap),
  dbsnp: dbsnp(snap),
  sv: sv(snap),
  constraint: constraint(snap),
  dosage: dosage(snap),
  mito: mito(snap),
  censat: censat(snap),
  sexchrom: sexchrom(snap),
  ccres: ccres(snap),
  cpgislands: cpgislands(snap),
  trna: trna(snap),
  encode: encode(snap),
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
