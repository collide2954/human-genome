import { test } from 'node:test';
import assert from 'node:assert/strict';
import { deflateSync } from 'node:zlib';
import {
  KINDS, kindOf, taxonomyReports, lineageOf, sharedDepths, orthologPartners, orthologDepth, vertebrateOrthologs, allianceOrthologs, depthSummary,
  ensemblRelease, ensemblTaxa, comparaSummary, comparaSpecies, GENOME_GROUPS, genomePoints, ASSEMBLY_LEVELS, primateSummary,
  readBigBed, gapIntervals, regionLayers, conservedSummary,
} from '../../scripts/evolution.mjs';
import { linearEdges, logEdges } from '../../src/stats.js';

const HUMAN = [1, 2759, 33154, 33213, 7742, 32524, 40674, 9443, 9604, 207598, 9605];
const NAMES = {
  1: 'root', 2759: 'Eukaryota', 33154: 'Opisthokonta', 33213: 'Bilateria', 7742: 'Vertebrata', 32524: 'Amniota', 40674: 'Mammalia', 9443: 'Primates', 9604: 'Hominidae', 207598: 'Homininae', 9605: 'Homo',
};
const report = (id, name, parents, group = '') => ({ taxonomy: { tax_id: id, current_scientific_name: { name }, group_name: group, parents } });
const TAXA = {
  reports: [
    report(9606, 'Homo sapiens', HUMAN, 'primates'),
    report(9598, 'Pan troglodytes', [1, 2759, 33154, 33213, 7742, 32524, 40674, 9443, 9604, 207598, 9596], 'primates'),
    report(9597, 'Pan paniscus', [1, 2759, 33154, 33213, 7742, 32524, 40674, 9443, 9604, 207598, 9596], 'primates'),
    report(10090, 'Mus musculus', [1, 2759, 33154, 33213, 7742, 32524, 40674, 10066], 'rodents'),
    report(9031, 'Gallus gallus', [1, 2759, 33154, 33213, 7742, 32524, 8457, 8782, 9030], 'birds'),
    report(7955, 'Danio rerio', [1, 2759, 33154, 33213, 7742, 7898, 7953], 'ray-finned fishes'),
    report(7227, 'Drosophila melanogaster', [1, 2759, 33154, 33208, 33213, 6656, 50557], 'flies'),
    report(559292, 'Saccharomyces cerevisiae S288C', [1, 2759, 33154, 4751, 4932], 'ascomycete fungi'),
  ],
};
const NAME_REPORTS = { reports: Object.entries(NAMES).map(([id, name]) => report(Number(id), name, [])) };
const taxa = taxonomyReports(TAXA);
const lineage = lineageOf(HUMAN, taxonomyReports(NAME_REPORTS));
const depthOf = sharedDepths(taxa, HUMAN);

test('taxonomyReports keys NCBI Datasets taxonomy reports by tax ID with name, group and lineage', () => {
  assert.equal(taxa.size, 8);
  assert.deepEqual(taxa.get('9598'), {
    id: 9598, name: 'Pan troglodytes', group: 'primates', parents: [1, 2759, 33154, 33213, 7742, 32524, 40674, 9443, 9604, 207598, 9596],
  });
});

test('lineageOf names the human lineage from the root down', () => {
  assert.equal(lineage.length, 11);
  assert.deepEqual(lineage[0], { id: 1, name: 'root' });
  assert.deepEqual(lineage.at(-2), { id: 207598, name: 'Homininae' });
});

test('sharedDepths places each taxon at the deepest node of the human lineage it shares', () => {
  assert.equal(depthOf.get('9598'), 9);
  assert.equal(depthOf.get('10090'), 6);
  assert.equal(depthOf.get('9031'), 5);
  assert.equal(depthOf.get('7955'), 4);
  assert.equal(depthOf.get('7227'), 3);
  assert.equal(depthOf.get('559292'), 2);
});

test('kindOf sorts species into primates, other mammals, birds, reptiles, amphibians, fishes and other', () => {
  assert.deepEqual(KINDS.map((k) => k.key), ['primates', 'mammals', 'birds', 'reptiles', 'amphibians', 'rayfinned', 'cartilaginous', 'other']);
  assert.equal(kindOf(taxa.get('9598').parents), 'primates');
  assert.equal(kindOf(taxa.get('10090').parents), 'mammals');
  assert.equal(kindOf(taxa.get('9031').parents), 'birds');
  assert.equal(kindOf(taxa.get('7955').parents), 'rayfinned');
  assert.equal(kindOf(taxa.get('7227').parents), 'other');
});

const ORTHO = [
  '#tax_id\tGeneID\trelationship\tOther_tax_id\tOther_GeneID',
  '9606\t1\tOrtholog\t9598\t101',
  '9606\t1\tOrtholog\t10090\t201',
  '9606\t1\tOrtholog\t7955\t301',
  '9606\t2\tOrtholog\t9598\t102',
  '9606\t2\tOrtholog\t10090\t202',
  '9606\t3\tOrtholog\t9598\t103',
  '9606\t9\tOrtholog\t9598\t109',
  '10090\t201\tOrtholog\t10116\t401',
  '9606\t4\tOrtholog\t9031\t501',
  '9606\t1\tOrtholog\t9597\t601',
  '9606\t2\tOrtholog\t9597\t602',
];
const CODING = new Set(['1', '2', '3', '4', '5']);

const feed = (acc, lines) => {
  for (const line of lines) acc.add(line);
  return acc.result();
};

test('orthologPartners counts human protein-coding genes with an NCBI ortholog in each species', () => {
  assert.deepEqual([...feed(orthologPartners({ genes: CODING }), ORTHO)], [['9598', 3], ['10090', 2], ['7955', 1], ['9031', 1], ['9597', 2]]);
});

test('orthologDepth keeps the deepest shared node at which each human gene has an ortholog', () => {
  assert.deepEqual(Object.fromEntries(feed(orthologDepth({ genes: CODING, depthOf }), ORTHO)), {
    1: 4, 2: 6, 3: 9, 4: 5,
  });
});

test('vertebrateOrthologs groups species by the deepest taxon they share with human', () => {
  const v = vertebrateOrthologs({
    partners: feed(orthologPartners({ genes: CODING }), ORTHO), taxa, lineage, depthOf, genes: 5,
  });
  assert.equal(v.genes, 5);
  assert.deepEqual(v.nodes, [
    { name: 'Homininae', species: 2, median: 2.5 },
    { name: 'Mammalia', species: 1, median: 2 },
    { name: 'Amniota', species: 1, median: 1 },
    { name: 'Vertebrata', species: 1, median: 1 },
  ]);
  assert.deepEqual(v.species, [
    ['9598', 'Pan troglodytes', 'primates', 0, 3],
    ['9597', 'Pan paniscus', 'primates', 0, 2],
    ['10090', 'Mus musculus', 'mammals', 1, 2],
    ['9031', 'Gallus gallus', 'birds', 2, 1],
    ['7955', 'Danio rerio', 'rayfinned', 3, 1],
  ]);
});

const ALLIANCE_HEAD = 'Gene1ID\tGene1Symbol\tGene1SpeciesTaxonID\tGene1SpeciesName\tGene2ID\tGene2Symbol\tGene2SpeciesTaxonID\tGene2SpeciesName\tAlgorithms\tAlgorithmsMatch\tOutOfAlgorithms\tIsBestScore\tIsBestRevScore';
const pair = (g1, s1, t1, n1, g2, s2, t2, n2, match, of) => [g1, s1, `NCBITaxon:${t1}`, n1, g2, s2, `NCBITaxon:${t2}`, n2, 'A|B', match, of, 'Yes', 'Yes'].join('\t');
const ALLIANCE = [
  '##########',
  '# Alliance Database Version: 9.0.0',
  '# Date file generated (UTC): 2026-04-05 03:33',
  '##########',
  ALLIANCE_HEAD,
  pair('HGNC:1', 'G1', 9606, 'Homo sapiens', 'FB:1', 'g1', 7227, 'Drosophila melanogaster', 2, 9),
  pair('HGNC:1', 'G1', 9606, 'Homo sapiens', 'FB:2', 'g1b', 7227, 'Drosophila melanogaster', 9, 9),
  pair('HGNC:2', 'G2', 9606, 'Homo sapiens', 'SGD:1', 'Y2', 559292, 'Saccharomyces cerevisiae', 5, 9),
  pair('HGNC:6', 'G6', 9606, 'Homo sapiens', 'MGI:6', 'g6', 10090, 'Mus musculus', 1, 10),
  pair('HGNC:7', 'G7', 9606, 'Homo sapiens', 'MGI:7', 'g7', 10090, 'Mus musculus', 10, 10),
  pair('MGI:6', 'g6', 10090, 'Mus musculus', 'HGNC:6', 'G6', 9606, 'Homo sapiens', 1, 10),
  pair('HGNC:2', 'G2', 9606, 'Homo sapiens', 'MGI:2', 'g2', 10090, 'Mus musculus', 10, 10),
].join('\n');
const HGNC_CODING = new Set(['HGNC:1', 'HGNC:2', 'HGNC:3', 'HGNC:4', 'HGNC:5', 'HGNC:6']);

test('allianceOrthologs counts human genes with model organism orthologs by how many methods agree', () => {
  const a = allianceOrthologs(ALLIANCE, { genes: HGNC_CODING, disease: new Set(['HGNC:2', 'HGNC:3', 'HGNC:99']) });
  assert.equal(a.version, '9.0.0');
  assert.equal(a.generated, '2026-04-05');
  assert.equal(a.genes, 6);
  assert.equal(a.disease, 2);
  assert.deepEqual(a.species, [
    {
      taxon: '10090', name: 'Mus musculus', methods: 10, genes: 2, all: 1, most: 0, fewer: 1, disease: 1,
    },
    {
      taxon: '7227', name: 'Drosophila melanogaster', methods: 9, genes: 1, all: 1, most: 0, fewer: 0, disease: 0,
    },
    {
      taxon: '559292', name: 'Saccharomyces cerevisiae', methods: 9, genes: 1, all: 0, most: 1, fewer: 0, disease: 1,
    },
  ]);
  assert.deepEqual([...a.byGene].map(([g, s]) => [g, [...s].sort()]), [['HGNC:1', ['7227']], ['HGNC:2', ['10090', '559292']], ['HGNC:6', ['10090']]]);
});

test('depthSummary counts human genes by the most distant taxon in which they have an ortholog', () => {
  const genes = [
    { hgnc: 'HGNC:1', entrez: '1', symbol: 'G1' },
    { hgnc: 'HGNC:2', entrez: '2', symbol: 'G2' },
    { hgnc: 'HGNC:3', entrez: '3', symbol: 'G3' },
    { hgnc: 'HGNC:4', entrez: '4', symbol: 'G4' },
    { hgnc: 'HGNC:5', entrez: '5', symbol: 'G5' },
    { hgnc: 'HGNC:6', entrez: null, symbol: 'G6' },
  ];
  const d = depthSummary({
    genes,
    ncbi: feed(orthologDepth({ genes: CODING, depthOf }), ORTHO),
    alliance: allianceOrthologs(ALLIANCE, { genes: HGNC_CODING, disease: new Set() }).byGene,
    species: ['9598', '9597', '10090', '7955', '9031', '7227', '559292'],
    depthOf,
    lineage,
  });
  assert.equal(d.genes, 6);
  assert.deepEqual(d.nodes, [
    { name: 'Homininae', genes: 1, examples: ['G3'] },
    { name: 'Mammalia', genes: 1, examples: ['G6'] },
    { name: 'Amniota', genes: 1, examples: ['G4'] },
    { name: 'Vertebrata', genes: 0, examples: [] },
    { name: 'Bilateria', genes: 1, examples: ['G1'] },
    { name: 'Opisthokonta', genes: 1, examples: ['G2'] },
  ]);
  assert.deepEqual(d.none, { genes: 1, examples: ['G5'] });
});

test('ensemblRelease finds the newest numbered release in the FTP listing', () => {
  assert.equal(ensemblRelease('<a href="release-99/">release-99/</a><a href="release-115/">x</a><a href="release-116/">y</a><a href="current_fasta/">z</a>'), 116);
  assert.equal(ensemblRelease('<a href="current_fasta/">z</a>'), null);
});

const SPECIES = [
  '#name\tspecies\tdivision\ttaxonomy_id\tassembly\tassembly_accession\tgenebuild\tvariation\tmicroarray\tpan_compara\tpeptide_compara\tgenome_alignments\tother_alignments\tcore_db\tspecies_id',
  'Chimpanzee\tpan_troglodytes\tEnsemblVertebrates\t9598\tPan_tro_3.0\tGCA_000001515.5\t2018\tN\tN\tN\tY\tY\tY\tpan_troglodytes_core_116_3\t1\t',
  'Mouse\tmus_musculus\tEnsemblVertebrates\t10090\tGRCm39\tGCA_000001635.9\t2020\tY\tY\tY\tY\tY\tY\tmus_musculus_core_116_39\t1\t',
].join('\n');

test('ensemblTaxa maps Ensembl species names to NCBI tax IDs', () => {
  assert.deepEqual([...ensemblTaxa(SPECIES)], [['pan_troglodytes', '9598'], ['mus_musculus', '10090']]);
});

const COMPARA = [
  'gene_stable_id\tprotein_stable_id\tspecies\tidentity\thomology_type\thomology_gene_stable_id\thomology_protein_stable_id\thomology_species\thomology_identity\tdn\tds\tgoc_score\twga_coverage\tis_high_confidence\thomology_id',
  'ENSG1\tENSP1\thomo_sapiens\t99.5\tortholog_one2one\tENSPTRG1\tENSPTRP1\tpan_troglodytes\t99.1\tNULL\tNULL\t100\t100\t1\t1',
  'ENSG2\tENSP2\thomo_sapiens\t97\tortholog_one2one\tENSPTRG2\tENSPTRP2\tpan_troglodytes\t96\tNULL\tNULL\t75\t90\t0\t2',
  'ENSG3\tENSP3\thomo_sapiens\t80\tortholog_one2many\tENSPTRG3\tENSPTRP3\tpan_troglodytes\t70\tNULL\tNULL\t50\tNULL\t1\t3',
  'ENSG4\tENSP4\thomo_sapiens\t40\tother_paralog\tENSG5\tENSP5\thomo_sapiens\t41\tNULL\tNULL\tNULL\tNULL\tNULL\t4',
  'ENSMUSG1\tENSMUSP1\tmus_musculus\t85\tortholog_one2one\tENSG1\tENSP1\thomo_sapiens\t84\tNULL\tNULL\t100\t100\t1\t5',
  'ENSMUSG2\tENSMUSP2\tmus_musculus\t60\tortholog_many2many\tENSG6\tENSP6\thomo_sapiens\t62\tNULL\tNULL\t0\t0\t0\t6',
  'ENSMUSG3\tENSMUSP3\tmus_musculus\t70\tortholog_one2one\tENSRNOG1\tENSRNOP1\trattus_norvegicus\t71\tNULL\tNULL\t100\t100\t1\t7',
];

test('comparaSummary counts human orthologs per Ensembl species by type with the identity of each one-to-one pair', () => {
  const c = feed(comparaSummary(), COMPARA);
  assert.deepEqual([...c.keys()], ['pan_troglodytes', 'mus_musculus']);
  assert.deepEqual(c.get('pan_troglodytes'), {
    one2one: 2, one2many: 1, many2many: 0, highConfidence: 1, identities: [99.5, 97],
  });
  assert.deepEqual(c.get('mus_musculus'), {
    one2one: 1, one2many: 0, many2many: 1, highConfidence: 1, identities: [84],
  });
});

test('comparaSpecies places Ensembl species on the human lineage with identity quartiles and histograms', () => {
  const s = comparaSpecies(feed(comparaSummary(), COMPARA), {
    taxa: ensemblTaxa(SPECIES), reports: taxa, lineage, depthOf, keys: ['pan_troglodytes', 'mus_musculus'],
  });
  assert.deepEqual(s.nodes, [{ name: 'Homininae', species: 1, median: 98.3 }, { name: 'Mammalia', species: 1, median: 84 }]);
  assert.deepEqual(s.species, [
    {
      species: 'pan_troglodytes', name: 'Pan troglodytes', kind: 'primates', node: 0, one2one: 2, one2many: 1, many2many: 0, highConfidence: 0.5, median: 98.3, q1: 97.6, q3: 98.9,
    },
    {
      species: 'mus_musculus', name: 'Mus musculus', kind: 'mammals', node: 1, one2one: 1, one2many: 0, many2many: 1, highConfidence: 1, median: 84, q1: 84, q3: 84,
    },
  ]);
  assert.deepEqual(s.identity.edges, linearEdges(0, 100, 2));
  assert.equal(s.identity.bySpecies.pan_troglodytes[48], 1);
  assert.equal(s.identity.bySpecies.pan_troglodytes[49], 1);
  assert.equal(s.identity.bySpecies.mus_musculus[42], 1);
  assert.equal(s.identity.bySpecies.mus_musculus.reduce((a, b) => a + b, 0), 1);
});

const genome = (accession, taxId, name, length, coding) => ({
  accession, organism: { tax_id: taxId, organism_name: name }, assembly_stats: { total_sequence_length: String(length) }, annotation_info: { stats: { gene_counts: { protein_coding: coding } } },
});

test('genomePoints keeps one annotated reference genome per species with its length, protein-coding genes and group', () => {
  const reports = [
    genome('GCF_000001405.40', 9606, 'Homo sapiens', 3099734149, 20080),
    genome('GCA_000001405.29', 9606, 'Homo sapiens', 3099734149, 20000),
    genome('GCF_000001635.27', 10090, 'Mus musculus', 2728222451, 22000),
    genome('GCF_000146045.2', 559292, 'Saccharomyces cerevisiae S288C', 12157105, 6000),
    genome('GCA_000001215.4', 7227, 'Drosophila melanogaster', 143726002, 0),
    { accession: 'GCA_2', organism: { tax_id: 6239, organism_name: 'Caenorhabditis elegans' }, assembly_stats: { total_sequence_length: '100286401' } },
  ];
  assert.deepEqual(GENOME_GROUPS.map((g) => g.key), ['mammals', 'birds', 'fishes', 'vertebrates', 'insects', 'animals', 'plants', 'fungi', 'other']);
  const g = genomePoints(reports, taxa);
  assert.equal(g.total, 3);
  assert.deepEqual(g.groups, GENOME_GROUPS.map((x) => x.label));
  assert.deepEqual(g.points, [
    [12.16, 6000, 7, 'Saccharomyces cerevisiae S288C'],
    [2728.22, 22000, 0, 'Mus musculus'],
    [3099.73, 20080, 0, 'Homo sapiens'],
  ]);
});

const assembly = (accession, taxId, name, level, date, n50) => ({
  accession, organism: { tax_id: taxId, organism_name: name }, assembly_info: { assembly_level: level, release_date: date }, assembly_stats: { contig_n50: n50 },
});

test('primateSummary charts non-human primate assemblies and the most contiguous assembly of each ape', () => {
  const apes = taxonomyReports({
    reports: [
      report(9598, 'Pan troglodytes', [9443, 314295, 9604]),
      report(9595, 'Gorilla gorilla gorilla', [9443, 314295, 9604, 9592, 9593]),
      report(9544, 'Macaca mulatta', [9443, 9526]),
    ],
  });
  const p = primateSummary([
    assembly('GCA_1', 9606, 'Homo sapiens', 'Complete Genome', '2022-01-24', 150617247),
    assembly('GCA_2', 9598, 'Pan troglodytes', 'Complete Genome', '2025-04-09', 147700000),
    assembly('GCA_3', 9598, 'Pan troglodytes', 'Scaffold', '2006-03-16', 39000),
    assembly('GCA_4', 9595, 'Gorilla gorilla gorilla', 'Chromosome', '2018-01-01', 23000000),
    assembly('GCA_5', 9544, 'Macaca mulatta', 'Contig', '2010-05-05', 0),
    assembly('GCA_6', 9544, 'Macaca mulatta', 'Complete Genome', '2025-07-24', 140000000),
  ], apes);
  assert.deepEqual(ASSEMBLY_LEVELS, ['Contig', 'Scaffold', 'Chromosome', 'Complete Genome']);
  assert.equal(p.total, 5);
  assert.equal(p.species, 3);
  assert.deepEqual(p.levels, {
    Contig: 1, Scaffold: 1, Chromosome: 1, 'Complete Genome': 2,
  });
  assert.deepEqual(p.names, ['Pan troglodytes', 'Gorilla gorilla gorilla', 'Macaca mulatta']);
  assert.deepEqual(p.points, [['2025-04-09', 147700000, 3, 0], ['2006-03-16', 39000, 1, 0], ['2018-01-01', 23000000, 2, 1], ['2025-07-24', 140000000, 3, 2]]);
  assert.deepEqual(p.apes, [
    {
      name: 'Pan troglodytes', accession: 'GCA_2', level: 'Complete Genome', n50: 147700000, date: '2025-04-09',
    },
    {
      name: 'Gorilla gorilla', accession: 'GCA_4', level: 'Chromosome', n50: 23000000, date: '2018-01-01',
    },
  ]);
});

const u8 = (v) => Buffer.from([v]);
const u16 = (v) => {
  const b = Buffer.alloc(2);
  b.writeUInt16LE(v);
  return b;
};
const u32 = (v) => {
  const b = Buffer.alloc(4);
  b.writeUInt32LE(v);
  return b;
};
const u64 = (v) => {
  const b = Buffer.alloc(8);
  b.writeBigUInt64LE(BigInt(v));
  return b;
};

const bigBed = (chroms, blocks) => {
  const keySize = Math.max(...chroms.map(([n]) => n.length));
  const key = (name) => Buffer.concat([Buffer.from(name), Buffer.alloc(keySize - name.length)]);
  const leafA = chroms.slice(0, 2);
  const leafB = chroms.slice(2);
  const leafAAt = 64 + 32 + 4 + 2 * (keySize + 8);
  const leafBAt = leafAAt + 4 + leafA.length * (keySize + 8);
  const tree = Buffer.concat([
    u32(0x78CA8C91), u32(2), u32(keySize), u32(8), u64(chroms.length), u64(0),
    u8(0), u8(0), u16(2), key(leafA[0][0]), u64(leafAAt), key(leafB[0][0]), u64(leafBAt),
    u8(1), u8(0), u16(leafA.length), ...leafA.flatMap(([n, size], i) => [key(n), u32(i), u32(size)]),
    u8(1), u8(0), u16(leafB.length), ...leafB.flatMap(([n, size], i) => [key(n), u32(i + 2), u32(size)]),
  ]);
  const dataAt = 64 + tree.length;
  const packed = blocks.map((recs) => deflateSync(Buffer.concat(recs.flatMap(([c, s, e, rest]) => [u32(c), u32(s), u32(e), Buffer.from(`${rest}\0`)]))));
  const offsets = [];
  let at = dataAt + 8;
  for (const p of packed) {
    offsets.push([at, p.length]);
    at += p.length;
  }
  const data = Buffer.concat([u64(blocks.reduce((n, b) => n + b.length, 0)), ...packed]);
  const indexAt = dataAt + data.length;
  const span = (recs) => [recs[0][0], recs[0][1], recs.at(-1)[0], recs.at(-1)[2]];
  const index = Buffer.concat([
    u32(0x2468ACE0), u32(2), u64(blocks.length), u32(0), u32(0), u32(1), u32(500), u64(indexAt), u32(512), u32(0),
    u8(0), u8(0), u16(1), u32(0), u32(0), u32(1), u32(500), u64(indexAt + 48 + 28),
    u8(1), u8(0), u16(blocks.length), ...blocks.flatMap((recs, i) => [...span(recs).map(u32), u64(offsets[i][0]), u64(offsets[i][1])]),
  ]);
  const header = Buffer.concat([u32(0x8789F2EB), u16(4), u16(0), u64(64), u64(dataAt), u64(indexAt), u16(5), u16(5), u64(0), u64(0), u32(16384), u64(0)]);
  return Buffer.concat([header, tree, data, index]);
};

test('readBigBed walks the chromosome tree and the R-tree index and inflates each data block', () => {
  const file = readBigBed(bigBed([['chr1', 1000], ['chr2', 500], ['chrUn_x', 300]], [
    [[0, 130, 160, 'lod=10\t500'], [0, 310, 320, 'lod=5\t400']],
    [[1, 10, 20, 'lod=7\t450'], [1, 15, 30, 'lod=8\t460']],
  ]));
  assert.equal(file.fieldCount, 5);
  assert.equal(file.count, 4);
  assert.deepEqual(file.chroms, [{ id: 0, name: 'chr1', size: 1000 }, { id: 1, name: 'chr2', size: 500 }, { id: 2, name: 'chrUn_x', size: 300 }]);
  assert.deepEqual([...file.records()], [
    ['chr1', 130, 160, 'lod=10\t500'], ['chr1', 310, 320, 'lod=5\t400'], ['chr2', 10, 20, 'lod=7\t450'], ['chr2', 15, 30, 'lod=8\t460'],
  ]);
});

test('gapIntervals reads UCSC gap rows on the standard chromosomes', () => {
  assert.deepEqual(gapIntervals(['585\tchr1\t900\t1000\t1\tN\t100\ttelomere\tno', '1\tchrUn_x\t0\t10\t1\tN\t10\tcontig\tno'].join('\n')), { 1: [[900, 1000]] });
});

const GTF = [
  ['chr1', 'gene', 101, 400, 'gene_id "G1"; gene_type "protein_coding"; gene_name "G1";'],
  ['chr1', 'transcript', 101, 400, 'gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";'],
  ['chr1', 'exon', 101, 150, 'gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";'],
  ['chr1', 'exon', 301, 400, 'gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";'],
  ['chr1', 'CDS', 121, 150, 'gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";'],
  ['chr1', 'CDS', 301, 350, 'gene_id "G1"; transcript_id "T1"; gene_type "protein_coding";'],
  ['chr1', 'gene', 601, 700, 'gene_id "G2"; gene_type "lncRNA"; gene_name "G2";'],
  ['chr1', 'exon', 601, 650, 'gene_id "G2"; transcript_id "T2"; gene_type "lncRNA";'],
  ['chrM', 'gene', 1, 70, 'gene_id "G3"; gene_type "Mt_tRNA"; gene_name "G3";'],
].map(([chrom, feature, start, end, attrs]) => [chrom, 'HAVANA', feature, start, end, '.', '+', '.', attrs].join('\t'));

const layers = feed(regionLayers(), ['##description: test', ...GTF]);

test('regionLayers merges coding, exon and gene territory per chromosome as nested layers', () => {
  assert.deepEqual(layers, {
    1: {
      cds: [[120, 150], [300, 350]], coding: [[100, 150], [300, 400]], exons: [[100, 150], [300, 400], [600, 650]], genic: [[100, 400], [600, 700]],
    },
  });
});

test('conservedSummary splits conserved bases by chromosome and by coding, UTR, non-coding exon, intron and intergenic territory', () => {
  const acc = conservedSummary({ layers, gaps: { 1: [[900, 1000]] }, lengths: { 1: 1000, 2: 500 } });
  for (const [chrom, start, end] of [['1', 130, 160], ['1', 310, 320], ['1', 360, 380], ['1', 620, 630], ['1', 680, 690], ['1', 800, 850], ['2', 10, 20], ['2', 15, 30]]) acc.add(chrom, start, end);
  const c = acc.result();
  assert.equal(c.elements, 8);
  assert.equal(c.bases, 150);
  assert.deepEqual(c.byChrom, { 1: { length: 900, conserved: 130 }, 2: { length: 500, conserved: 20 } });
  assert.deepEqual(c.classes, {
    cds: { bases: 80, conserved: 30 },
    utr: { bases: 70, conserved: 20 },
    ncExon: { bases: 50, conserved: 10 },
    intron: { bases: 200, conserved: 20 },
    intergenic: { bases: 1000, conserved: 70 },
  });
  assert.equal(c.median, 10);
  assert.deepEqual(c.longest, {
    chrom: '1', start: 801, end: 850, length: 50,
  });
  assert.deepEqual(c.lengths.edges, logEdges(0, 4, 10));
  assert.equal(c.lengths.counts.reduce((a, b) => a + b, 0), 8);
});
