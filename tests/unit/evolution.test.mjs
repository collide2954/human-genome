import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  KINDS, kindOf, taxonomyReports, lineageOf, sharedDepths, orthologPartners, orthologDepth, vertebrateOrthologs, allianceOrthologs, depthSummary,
  ensemblRelease, ensemblTaxa, comparaSummary, comparaSpecies,
} from '../../scripts/evolution.mjs';
import { linearEdges } from '../../src/stats.js';

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
