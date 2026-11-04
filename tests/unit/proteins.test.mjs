import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  uniprotSummary, locationsOf, hpaSummary, gtexSummary, pharosSummary,
} from '../../scripts/proteins.mjs';

const HEAD = ['Entry', 'Gene Names (primary)', 'Length', 'Protein existence', 'PDB', 'Transmembrane', 'Signal peptide', 'Subcellular location [CC]', 'Pfam', 'Region', 'Modified residue', 'Glycosylation', 'Cross-link'].join('\t');
const entry = (o) => [o.acc, o.gene ?? '', o.length, o.pe ?? 'Evidence at protein level', o.pdb ?? '', o.tm ?? '', o.signal ?? '', o.loc ?? '', o.pfam ?? '', o.region ?? '', o.mod ?? '', o.glyco ?? '', o.cross ?? ''].join('\t');

const LINES = [
  HEAD,
  entry({
    acc: 'P04637', gene: 'TP53', length: 393, pdb: '1A1U;1AIE;1TUP;', pfam: 'PF00870;PF08563;', loc: 'SUBCELLULAR LOCATION: Cytoplasm {ECO:1}. Nucleus {ECO:2}. Note=Moves to the nucleus.',
    region: 'REGION 50..96; /note="Disordered"; /evidence="x"; REGION 66..110; /note="Interaction with WWOX"; REGION 90..100; /note="Disordered"',
    mod: 'MOD_RES 9; /note="Phosphoserine; by HIPK4"; MOD_RES 18; /note="Phosphothreonine"; MOD_RES 120; /note="N6-acetyllysine"; MOD_RES 333; /note="Asymmetric dimethylarginine"',
    cross: 'CROSSLNK 291; /note="Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)"; CROSSLNK 386; /note="Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO)"',
  }),
  entry({
    acc: 'P08588', gene: 'ADRB1', length: 477, pfam: 'PF00001;', signal: '',
    tm: Array.from({ length: 7 }, (_, i) => `TRANSMEM ${i * 40 + 1}..${i * 40 + 20}; /note="Helical"`).join('; '),
    loc: 'SUBCELLULAR LOCATION: Cell membrane {ECO:1}; Multi-pass membrane protein {ECO:2}. Early endosome {ECO:3}.',
    glyco: 'CARBOHYD 15; /note="N-linked (GlcNAc...) asparagine"; CARBOHYD 20; /note="O-linked (GalNAc...) threonine"',
  }),
  entry({
    acc: 'Q8WZ42', gene: 'TTN', length: 34350, pe: 'Evidence at protein level', pdb: '1BPV;', pfam: 'PF07679;PF07679;PF00041;', loc: 'SUBCELLULAR LOCATION: Cytoplasm, myofibril, sarcomere {ECO:1}. Nucleus {ECO:2}.',
  }),
  entry({
    acc: 'P01308', gene: 'INS', length: 110, pe: 'Evidence at transcript level', signal: 'SIGNAL 1..24', loc: 'SUBCELLULAR LOCATION: Secreted {ECO:1}.', pfam: 'PF00001;',
  }),
];

const run = () => {
  const acc = uniprotSummary();
  LINES.forEach((l) => acc.add(l));
  return acc.result();
};

test('locationsOf keeps main compartments and drops topology, notes and evidence', () => {
  assert.deepEqual(locationsOf('SUBCELLULAR LOCATION: Cell membrane {ECO:1}; Multi-pass membrane protein {ECO:2}. Early endosome {ECO:3}. Note=Found in the Golgi.'), ['Cell membrane', 'Early endosome']);
  assert.deepEqual(locationsOf('SUBCELLULAR LOCATION: Cytoplasm, myofibril, sarcomere {ECO:1}. Nucleus {ECO:2}. [Isoform 2]: Cytoplasm {ECO:3}.'), ['Cytoplasm', 'Nucleus']);
  assert.deepEqual(locationsOf(''), []);
});

test('uniprotSummary counts proteins, lengths, evidence and structures', () => {
  const u = run();
  assert.equal(u.count, 4);
  assert.deepEqual(u.longest[0], { acc: 'Q8WZ42', gene: 'TTN', length: 34350 });
  assert.equal(u.lengths.summary.median, 435);
  assert.deepEqual(u.existence, { 'Evidence at protein level': 3, 'Evidence at transcript level': 1 });
  assert.equal(u.pdb.withAny, 2);
  assert.deepEqual(u.pdb.top[0], { acc: 'P04637', gene: 'TP53', n: 3 });
  assert.deepEqual(u.pdb.counts.slice(0, 4), [2, 1, 0, 1]);
});

test('uniprotSummary counts membrane spans, signal peptides, locations and domains', () => {
  const u = run();
  assert.deepEqual(u.transmem.counts.slice(0, 8), [3, 0, 0, 0, 0, 0, 0, 1]);
  assert.equal(u.transmem.proteins, 1);
  assert.equal(u.signal, 1);
  assert.deepEqual(u.locations.slice(0, 2), [{ name: 'Cytoplasm', n: 2 }, { name: 'Nucleus', n: 2 }]);
  assert.deepEqual(u.pfam.top[0], { id: 'PF00001', n: 2 });
  assert.equal(u.pfam.withAny, 4);
  assert.deepEqual(u.pfam.perProtein.slice(0, 4), [0, 2, 2, 0]);
});

test('uniprotSummary measures disorder once per residue and counts modifications by kind', () => {
  const u = run();
  assert.equal(u.disorder.residues, 51);
  assert.equal(u.disorder.proteins, 1);
  assert.equal(u.disorder.share.reduce((a, b) => a + b, 0), 4);
  assert.deepEqual(u.ptm, {
    phosphoS: 1, phosphoT: 1, phosphoY: 0, acetyl: 1, methyl: 1, nGlyco: 1, oGlyco: 1, ubiquitin: 1, sumo: 1,
  });
});

const HPA_HEAD = ['Gene', 'Ensembl', 'Protein class', 'Evidence', 'RNA tissue specificity', 'RNA tissue specific nTPM', 'RNA single cell type specificity', 'RNA single cell type specific nCPM', 'Secretome location', 'Subcellular main location'].map((h) => `"${h}"`).join('\t');
const hpa = (o) => [o.gene, 'ENSG1', o.cls ?? '', o.ev ?? 'Evidence at protein level', o.spec, o.ntpm ?? '', o.sc ?? '', o.ncpm ?? '', o.secretome ?? '', o.loc ?? ''].map((v) => `"${v}"`).join('\t');

test('hpaSummary counts specificity categories, enriched tissues and cell types, and locations', () => {
  const h = hpaSummary([
    HPA_HEAD,
    hpa({ gene: 'A', spec: 'Tissue enriched', ntpm: 'testis: 120.5', sc: 'Cell type enriched', ncpm: 'Late spermatids: 300', cls: 'Enzymes, Predicted membrane proteins', loc: 'Nucleoplasm, Cytosol' }),
    hpa({ gene: 'B', spec: 'Tissue enriched', ntpm: 'brain: 50', loc: 'Cytosol' }),
    hpa({ gene: 'C', spec: 'Group enriched', ntpm: 'brain: 10;retina: 20', secretome: 'Secreted to blood' }),
    hpa({ gene: 'D', spec: 'Low tissue specificity' }),
  ].join('\n'));
  assert.deepEqual(h.specificity, { 'Tissue enriched': 2, 'Group enriched': 1, 'Low tissue specificity': 1 });
  assert.deepEqual(h.enriched, [{ tissue: 'brain', n: 1 }, { tissue: 'testis', n: 1 }]);
  assert.deepEqual(h.cellTypes, [{ cell: 'Late spermatids', n: 1 }]);
  assert.deepEqual(h.secretome, { 'Secreted to blood': 1 });
  assert.deepEqual(h.locations[0], { name: 'Cytosol', n: 2 });
  assert.deepEqual(h.classes.find((c) => c.name === 'Enzymes'), { name: 'Enzymes', n: 1 });
});

test('gtexSummary counts expressed and tissue-specific genes per tissue', () => {
  const lines = [
    '#1.2', '3\t3',
    'Name\tDescription\tLiver\tTestis\tBrain_Cortex',
    'ENSG1.1\tALB\t5000\t0.5\t0.1',
    'ENSG2.1\tPRM1\t0\t900\t0',
    'ENSG3.1\tACTB\t500\t400\t300',
  ];
  const g = gtexSummary(lines);
  assert.deepEqual(g.tissues, ['Liver', 'Testis', 'Brain_Cortex']);
  assert.deepEqual(g.expressed, [2, 2, 1]);
  assert.deepEqual(g.specific, [1, 1, 0]);
  assert.equal(g.genes, 3);
  assert.equal(g.correlation.length, 3);
  assert.equal(g.correlation[0][0], 1);
});

test('pharosSummary reads target development levels and novelty', () => {
  const p = pharosSummary({
    data: {
      targets: {
        count: 10,
        facets: [
          { facet: 'Target Development Level', values: [{ name: 'Tbio', value: 5 }, { name: 'Tdark', value: 3 }, { name: 'Tclin', value: 1 }, { name: 'Tchem', value: 1 }] },
          { facet: 'Log Novelty', values: [{ name: '-2', value: 4 }, { name: '0', value: 6 }] },
        ],
      },
    },
  });
  assert.deepEqual(p.levels, {
    Tclin: 1, Tchem: 1, Tbio: 5, Tdark: 3,
  });
  assert.equal(p.count, 10);
  assert.deepEqual(p.novelty, [[-2, 4], [0, 6]]);
});
