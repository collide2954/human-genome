import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  findings, assemblyNumbers, sequenceNumbers, proteinLabel, trnaPairs,
} from '../../src/findings.js';

test('proteinLabel names a protein by its first gene and counts the rest', () => {
  assert.equal(proteinLabel({ acc: 'P61769', gene: 'B2M' }), 'B2M');
  assert.equal(proteinLabel({ acc: 'P62805', gene: 'H4C1; H4C2; H4C3' }), 'H4C1 and 2 more');
  assert.equal(proteinLabel({ acc: 'P01858', gene: '' }), 'P01858');
});
import snapshot from '../../src/snapshot.js';

const chrom = (c, length, ungapped) => ({ chrom: c, length, ungapped, gc: 41, spanned: 1, unspanned: 1, refseq: null, genbank: null });

const fixture = {
  generated: '2026-11-08T06:17:00.000Z',
  assembly: {
    grch38: {
      name: 'GRCh38.p14',
      altLoci: 35,
      roles: [
        { role: 'assembled-molecule', count: 4, length: 400016569 },
        { role: 'unlocalized-scaffold', count: 42, length: 6978808 },
        { role: 'unplaced-scaffold', count: 127, length: 4485509 },
        { role: 'alt-scaffold', count: 261, length: 109535387 },
        { role: 'fix-patch', count: 164, length: 65946498 },
        { role: 'novel-patch', count: 90, length: 23679459 },
      ],
      chromosomes: [chrom('1', 248956422, 230481121), chrom('21', 46709983, 40088623), chrom('Y', 57227415, 26415048), chrom('MT', 16569, 16569)],
    },
    t2t: {
      name: 'T2T-CHM13v2.0',
      chromosomes: [chrom('1', 248387328, 248387328), chrom('21', 45090682, 45090682), chrom('Y', 62460029, 62460029)],
    },
    centromeres: { 1: 2905865, 21: 2050048, Y: 227095 },
    gaps: { 1: { contig: { count: 2, bp: 100 } } },
  },
  sequence: {
    chromosomes: [
      {
        chrom: '1', length: 1000, A: 300, C: 200, G: 200, T: 300, N: 0, cpg: 10, masked: 500, g4: { plus: 3, minus: 4 },
      },
      {
        chrom: '19', length: 1000, A: 250, C: 250, G: 250, T: 250, N: 0, cpg: 20, masked: 600, g4: { plus: 5, minus: 5 },
      },
      {
        chrom: 'MT', length: 100, A: 25, C: 25, G: 25, T: 25, N: 0, cpg: 4, masked: 0, g4: { plus: 0, minus: 0 },
      },
    ],
    homopolymers: {
      A: { length: 72, chrom: '4', start: 1 }, C: { length: 29, chrom: '7', start: 1 }, G: { length: 81, chrom: '2', start: 1 }, T: { length: 90, chrom: '12', start: 66057593 },
    },
  },
};

fixture.repeats = {
  rmsk: {
    total: { count: 50, bp: 1000 },
    families: [{ cls: 'LINE', family: 'L1', count: 10, bp: 400 }, { cls: 'SINE', family: 'Alu', count: 30, bp: 300 }],
    alu: { AluY: { count: 5, bp: 50 }, AluS: { count: 20, bp: 200 }, AluJ: { count: 5, bp: 50 } },
    l1hs: { count: 8, fullLength: 3 },
  },
  superDups: { pairs: 4, kinds: { intra: 1, inter: 3 }, covered: { 1: 60, 19: 40 } },
  tandem: {
    units: [{ unit: '1', count: 7, bp: 70 }, { unit: '2', count: 5, bp: 50 }],
    motifs: [{ motif: 'A', count: 7, bp: 70 }, { motif: 'AC', count: 4, bp: 40 }],
  },
};

fixture.genes = {
  release: { gencode: '50', ensembl: '116', date: '2026-04-08' },
  genes: {
    total: 100, byType: { protein_coding: 40, lncRNA: 30, processed_pseudogene: 20, miRNA: 10 }, byChrom: { 1: { protein_coding: 30 }, 19: { protein_coding: 10 } },
  },
  transcripts: { total: 500, byType: { protein_coding: 200, nonsense_mediated_decay: 50 } },
  spans: { pc: { summary: { median: 29403 } }, longest: [{ name: 'RBFOX1', span: 2473539 }] },
  isoforms: { summary: { median: 15 }, top: [{ name: 'ZBTB20', n: 361 }] },
  mane: { transcripts: 38, singleExon: 3, exons: { summary: { median: 8 }, top: [{ name: 'TTN', exons: 363 }] } },
  territory: {
    cds: 20, utr: 40, ncExon: 60, intron: 1000,
  },
  immune: {
    IGH: {
      V: 61, D: 37, J: 6, C: 9, pseudo: 95,
    },
  },
};
fixture.catalogs = {
  hgnc: {
    approved: 45289, locusGroups: { 'protein-coding gene': 19299, pseudogene: 14671 }, spreadsheet: { count: 26 }, pseudoParents: [{ parent: 'RN7SK', n: 293 }, { parent: 'RPL21', n: 137 }],
  },
  geneInfo: { types: { 'protein-coding': 20600 } },
  mane: { select: 19363 },
};

const codonFixture = Array.from({ length: 64 }, (_, i) => ([48, 50, 56].includes(i) ? 0 : 10));
codonFixture[34] = 400;
codonFixture[25] = 2;
fixture.sequence.mane = {
  classes: {
    'GT-AG': 90, 'GC-AG': 8, 'AT-AC': 1, other: 1,
  },
  kozak: Array.from({ length: 10 }, (_, i) => (i === 3 ? [6, 1, 2, 1] : [1, 1, 1, 1])),
  stops: { TAA: 3, TAG: 2, TGA: 5 },
  coding: 10,
  codons: codonFixture,
};

fixture.proteins = {
  release: '2026_03',
  count: 20431,
  lengths: { summary: { n: 20431, median: 415, mean: 558.8 } },
  longest: [{ acc: 'Q8WZ42', gene: 'TTN', length: 34350 }],
  pdb: { withAny: 8856, top: [{ acc: 'P61769', gene: 'B2M', n: 1345 }] },
  transmem: { counts: [15200, 2391, 325, 159, 413, 77, 277, 945], proteins: 5231 },
  signal: 3616,
  pfam: { withAny: 19457, top: [{ id: 'PF00096', n: 681, name: 'Zinc finger, C2H2 type', short: 'zf-C2H2' }] },
  disorder: { residues: 1141682, proteins: 11304 },
  ptm: {
    phosphoS: 33295, phosphoT: 6065, phosphoY: 2153, nGlyco: 16263,
  },
};
fixture.targets = {
  count: 20412,
  levels: {
    Tclin: 704, Tchem: 1904, Tbio: 12303, Tdark: 5501,
  },
};
fixture.expression = {
  hpa: {
    specificity: { 'Tissue enriched': 3132 },
    enriched: [{ tissue: 'testis', n: 932 }, { tissue: 'brain', n: 475 }],
    blood: { proteins: 4520, top: [{ gene: 'ALB', im: 4e13, ms: null }], lowest: { gene: 'IL6', im: 1200, ms: 5000 } },
  },
  gtex: {
    release: 'v10',
    tissues: ['Testis', 'Whole_Blood', 'Liver'],
    expressed: [25533, 11947, 15000],
    sites: [{ id: 'Testis', name: 'Testis' }, { id: 'Whole_Blood', name: 'Whole Blood' }, { id: 'Liver', name: 'Liver' }],
  },
};

fixture.networks = {
  version: '12.5',
  edges: {
    all: 6000000, medium: 900000, high: 250000, highest: 120000,
  },
  proteins: { all: 19500, high: 16000 },
  hubs: [{ id: '9606.ENSP1', name: 'UBC', degree: 2500 }],
};
fixture.pathways = {
  release: '98',
  pathways: 2700,
  proteins: 11000,
  topLevel: [
    {
      id: 'R-HSA-162582', name: 'Signal Transduction', proteins: 2700, pathways: 500,
    },
    {
      id: 'R-HSA-1430728', name: 'Metabolism', proteins: 2100, pathways: 600,
    },
    {
      id: 'R-HSA-168256', name: 'Immune System', proteins: 2000, pathways: 200,
    },
  ],
};
fixture.ontology = {
  annotations: 1005,
  not: 5,
  proteins: 19800,
  types: { protein: 19800, miRNA: 3000, protein_complex: 2000 },
  byAspect: {
    P: {
      experimental: 300, phylogenetic: 100, computational: 50, author: 38, electronic: 112,
    },
    F: {
      experimental: 150, phylogenetic: 0, computational: 0, author: 0, electronic: 100,
    },
    C: {
      experimental: 50, phylogenetic: 0, computational: 0, author: 0, electronic: 100,
    },
  },
};

test('findings state the network, pathway and ontology sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.interactions, 'STRING 12.5 links 16,000 human proteins by 250,000 interactions scored 0.7 or more; UBC has the most partners, 2,500');
  assert.equal(f.pathways, 'Reactome 98 places 11,000 human proteins in 2,700 pathways; Signal Transduction, the largest of 3 top-level pathways, holds 2,700');
  assert.equal(f.ontology, 'GO holds 1,000 annotations of 19,800 human proteins and 5,000 other gene products; 31.2% are electronic (IEA)');
});

test('findings state the protein sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.proteome, 'UniProt 2026_03 lists 20,431 reviewed human proteins; the median is 415 residues long and TTN, the longest, 34,350');
  assert.equal(f.structures, '8,856 reviewed human proteins, 43.3%, have at least one PDB structure; B2M appears in 1,345');
  assert.equal(f.membrane, '5,231 proteins cross a membrane: 2,391 once and 945 seven times; 3,616 carry a signal peptide');
  assert.equal(f.domains, '95.2% of proteins carry a Pfam domain; 681 carry zf-C2H2 (PF00096), the most common');
  assert.equal(f.modifications, '10% of residues sit in disordered regions; UniProt lists 33,295 phosphoserines, 6,065 phosphothreonines and 2,153 phosphotyrosines');
  assert.equal(f.targets, 'Pharos counts 704 targets of approved drugs (Tclin) among 20,412 human proteins; 5,501 are Tdark');
});

test('findings state the expression sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.specificity, 'HPA calls 3,132 genes tissue enriched; testis has 932 of them, brain 475');
  assert.equal(f.gtex, 'In GTEx v10, 25,533 genes reach a median of 1 TPM in Testis; 11,947 do in Whole Blood');
  assert.equal(f.plasma, 'Plasma levels in HPA run from 40 g/L for ALB to 1.2 ng/L for IL6, across 4,520 proteins');
});

fixture.regulation = {
  ccre: {
    count: 1000,
    classes: {
      dELS: { count: 600, bp: 120 }, pELS: { count: 200, bp: 50 }, PLS: { count: 100, bp: 20 }, 'CA-CTCF': { count: 100, bp: 10 },
    },
  },
  cpg: { count: 30, bp: 20, promoters: { tss: 200, inIsland: 120 } },
};
fixture.trnas = { genes: 432, anticodons: { 'Ala-AGC': 26, 'Glu-CTC': 8, 'Met-CAT': 11 } };
fixture.encode = { total: 20450, assays: [{ name: 'TF ChIP-seq', n: 3502 }] };

const mtGene = (kind, i) => ({
  name: `${kind}${i}`, start: i, end: i + 1, strand: '+', kind, product: null,
});
fixture.structure = {
  mito: {
    genes: [...Array.from({ length: 13 }, (_, i) => mtGene('protein', i)), ...Array.from({ length: 22 }, (_, i) => mtGene('tRNA', i)), ...Array.from({ length: 2 }, (_, i) => mtGene('rRNA', i))],
    variants: { positions: 8793 },
  },
  censat: {
    total: 240521368,
    classes: { hor: 70325827, hsat3: 69334314 },
    arrays: { hor: [{ chrom: '18', bp: 5479654 }, { chrom: '21', bp: 432213 }] },
    byChrom: { Y: { hsat3: 21685400, hsat1B: 14206138 } },
  },
};
fixture.genes.genes.byChrom.X = { protein_coding: 867 };
fixture.genes.genes.byChrom.Y = { protein_coding: 61 };

fixture.variation = {
  dbsnp: {
    build: 157, total: 1206053617, ss: { total: 4849775973 }, status: { 'Live RS': 1172689405 },
  },
  sv: {
    classes: {
      deletion: 2178603, duplication: 634313, cnv: 153161, insertion: 1651014,
    },
  },
};
fixture.constraint = {
  gnomad: { genes: 17331, below06: 3784 },
  dosage: { hi: { 3: 422 }, ts: { 3: 3 } },
  byDosage: { sufficient: { n: 314, median: 0.279 }, other: { n: 17017, median: 0.972 } },
};

test('findings state the variation and constraint sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.dbsnp, 'dbSNP build 157 holds 1,206,053,617 reference SNP records from 4,849,775,973 submitted variants; 1,172,689,405 are live');
  assert.equal(f.sv, 'dbVar’s non-redundant sets hold 2,178,603 deletions, 1,651,014 insertions and 634,313 duplications on GRCh38');
  assert.equal(f.constraint, '3,784 of 17,331 genes in gnomAD v4.1 have a LOEUF below 0.6; ClinGen haploinsufficient genes have a median of 0.28');
  assert.equal(f.dosage, 'ClinGen finds sufficient evidence of haploinsufficiency for 422 genes and of triplosensitivity for 3');
});

test('findings state the chromosome structure sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.mito, 'The mitochondrial genome’s 100 bp hold 37 genes, 13 for proteins, 22 tRNAs and 2 rRNAs; gnomAD finds variants at 8,793 positions');
  assert.equal(f.censat, 'T2T-CHM13 maps 240.5 Mb of satellite arrays; active α-satellite HORs make up 70.3 Mb, the most on chromosome 18 (5.5 Mb)');
  assert.equal(f.sexchrom, 'Chromosome X carries 867 protein-coding genes and Y 61; in T2T-CHM13 satellite arrays make up 57.5% of Y');
});

test('findings state the regulation sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.ccres, 'ENCODE’s registry lists 1,000 candidate cis-regulatory elements covering 10% of GRCh38; 600 are distal enhancer-like');
  assert.equal(f.cpgislands, '30 CpG islands cover 1% of GRCh38; 60% of MANE Select transcripts start inside one');
  assert.equal(f.trna, 'GtRNAdb lists 432 high-confidence tRNA genes; 58 of the 61 sense codons have none with an exactly pairing anticodon');
  assert.equal(f.encode, 'ENCODE has released 20,450 human experiments; TF ChIP-seq is the most common assay, with 3,502');
});

test('trnaPairs sets each sense codon against the tRNA genes whose anticodon pairs with it', () => {
  const pairs = trnaPairs(codonFixture, {
    'Ala-AGC': 26, 'Glu-CTC': 8, 'Met-CAT': 11, 'iMet-CAT': 9, 'SeC-TCA': 1,
  });
  assert.equal(pairs.length, 61);
  const of = (codon) => pairs.find((p) => p.codon === codon);
  assert.equal(of('GCT').copies, 26);
  assert.equal(of('GCT').anticodon, 'AGC');
  assert.equal(of('GAG').copies, 8);
  assert.equal(of('GAG').aa, 'Glu');
  assert.equal(of('GAG').perK.toFixed(1), '403.2');
  assert.equal(of('ATG').copies, 11);
  assert.equal(of('GCC').copies, 0);
});

test('findings state the splicing, translation and codon sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.splicing, '90% of MANE Select introns start GT and end AG; 8 start GC and 1 runs from AT to AC');
  assert.equal(f.translation, 'A purine sits three bases before the start codon in 80% of MANE coding sequences; 50% end in TGA');
  assert.equal(f.codons, 'GAG is the most used codon, 403.2 per thousand; CGC the least used sense codon, 2.0');
});

test('findings state the gene sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.genecount, 'Protein-coding gene counts run from 40 (GENCODE) to 20,600 (NCBI Gene)');
  assert.equal(f.density, 'Chromosome 1 carries 30,000 protein-coding genes per sequenced megabase; chromosome 19, 10,000');
  assert.equal(f.genesize, 'The median protein-coding gene spans 29,403 bp; RBFOX1, the longest, 2,473,539');
  assert.equal(f.exons, 'A MANE Select transcript has a median of 8 exons; TTN has 363; 3 have one');
  assert.equal(f.territory, 'Coding sequence covers 1% of GRCh38, exons of all genes 6% and introns 50%');
  assert.equal(f.isoforms, 'Protein-coding genes have a median of 15 GENCODE transcripts; ZBTB20 has 361');
  assert.equal(f.names, 'HGNC has approved 45,289 gene symbols; 26 were changed so spreadsheets would stop reading them as dates');
  assert.equal(f.pseudogenes, '14,671 named pseudogenes; RN7SK has 293 copies, RPL21 137');
  assert.equal(f.immune, '61 IGHV, 37 IGHD and 6 IGHJ segments make 13,542 heavy-chain V-D-J combinations before junctional diversity');
});

test('findings state the repeat sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.transposons, 'RepeatMasker marks 50% of GRCh38’s sequenced bases as repeats; L1 elements alone cover 20%');
  assert.equal(f.alu, '30 Alu copies, 5 in the youngest subfamily, AluY; 8 L1HS copies, 3 of them at least 6 kb long');
  assert.equal(f.segdups, 'Segmental duplications cover 100 bases, 5% of GRCh38; 3 of 4 pairs join different chromosomes');
  assert.equal(f.tandem, 'Tandem Repeats Finder calls 12 tandem repeats; the most frequent units are A, then AC');
});

test('sequenceNumbers totals the nuclear chromosomes and the CpG ratio', () => {
  const n = sequenceNumbers(fixture);
  assert.equal(n.acgt, 2000);
  assert.equal(n.gc, 0.45);
  assert.equal(n.cpgOe, 30 / ((450 * 450) / 2000));
  assert.equal(n.g4, 17);
  assert.deepEqual(n.richest, { chrom: '19', gc: 0.5 });
  assert.deepEqual(n.poorest, { chrom: '1', gc: 0.4 });
});

test('findings state the composition and motif sections as numbers', () => {
  const f = findings(fixture);
  assert.equal(f.composition, 'GRCh38 is 45% G and C; chromosome 19 is the richest at 50%, chromosome 1 the poorest at 40%');
  assert.equal(f.cpg, 'CpG turns up at 29.6% of the rate the C and G content predicts');
  assert.equal(f.motifs, '17 G-quadruplex motifs; the longest run of one letter is 90 Ts on chromosome 12');
});

test('assemblyNumbers totals the nuclear chromosomes of both assemblies', () => {
  const n = assemblyNumbers(fixture);
  assert.equal(n.grch38Length, 352893820);
  assert.equal(n.grch38Ungapped, 296984792);
  assert.equal(n.grch38N, 55909028);
  assert.equal(n.t2tLength, 355938039);
  assert.equal(n.added, 58953247);
  assert.equal(n.patches, 254);
  assert.deepEqual(n.mostAdded, { chrom: 'Y', bases: 36044981 });
});

test('findings state each assembly section as a number', () => {
  const f = findings(fixture);
  assert.equal(f.chromosomes, 'Chromosome 1 spans 248,956,422 bp in GRCh38; chromosome 21, the shortest, 46,709,983');
  assert.equal(f.gaps, 'GRCh38 leaves 55,909,028 bases of its chromosomes as N, 15.8% of their length');
  assert.equal(f.scaffolds, 'Besides its chromosomes, GRCh38.p14 holds 261 alternate loci in 35 regions and 254 patches');
  assert.equal(f.t2t, 'T2T-CHM13 sequences 58,953,247 more bases than GRCh38 does; 36,044,981 of them are on chromosome Y');
});

test('the real snapshot gives every finding as plain text without filler words', () => {
  const f = findings(snapshot);
  assert.ok(Object.keys(f).length >= 4);
  for (const [key, text] of Object.entries(f)) {
    assert.equal(typeof text, 'string', key);
    assert.ok(text.length > 10 && text.length < 140, `${key}: ${text}`);
    assert.doesNotMatch(text, /NaN|undefined|null|Infinity/, `${key}: ${text}`);
    assert.doesNotMatch(text, /\b(only|just|mere|merely|soar\w*|plung\w*|stagger\w*|shock\w*|remarkabl\w*|astonish\w*)\b/i, `${key}: ${text}`);
  }
});
