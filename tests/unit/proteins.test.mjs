import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  uniprotSummary, locationsOf, hpaSummary, gtexSummary, pharosSummary, tissueSites, pfamNames, alignSites,
} from '../../scripts/proteins.mjs';

test('alignSites names every GTEx column, trusting an API label only when it spells the column name', () => {
  const site = (id, name, organ, samples) => ({
    id, name, organ, samples, female: 1, male: 1, eGenes: 5, sGenes: 2,
  });
  const out = alignSites(['Liver', 'Brain_Anterior_cingulate_cortex_BA24', 'Liver_Hepatocyte', 'Small_Intestine_Terminal_Ileum_Mixed_Cell', 'Adipose_Visceral_Omentum'], [
    site('Liver', 'Liver', 'Liver', 262),
    site('Brain_Anterior_cingulate_cortex_BA24', 'Brain - Anterior qcingulate cortex (BA24)', 'Brain', 200),
    site('Small_Intestine_Terminal_Ileum', 'Small Intestine - Terminal Ileum', 'Small Intestine', 190),
    site('Adipose_Visceral_Omentum', 'Adipose - Visceral (Omentum)', 'Adipose Tissue', 587),
  ]);
  assert.deepEqual(out.map((s) => s.name), ['Liver', 'Brain Anterior cingulate cortex BA24', 'Liver Hepatocyte', 'Small Intestine Terminal Ileum Mixed Cell', 'Adipose - Visceral (Omentum)']);
  assert.deepEqual(out.map((s) => s.organ), ['Liver', 'Brain', 'Liver', 'Small Intestine', 'Adipose Tissue']);
  assert.deepEqual(out.map((s) => s.samples), [262, 200, null, null, 587]);
  assert.deepEqual(out.map((s) => s.id), ['Liver', 'Brain_Anterior_cingulate_cortex_BA24', 'Liver_Hepatocyte', 'Small_Intestine_Terminal_Ileum_Mixed_Cell', 'Adipose_Visceral_Omentum']);
  assert.equal(out[2].eGenes, null);
});

test('tissueSites keeps the name, organ, sample counts by sex and eGene counts of each GTEx tissue', () => {
  const rows = [
    {
      tissueSiteDetailId: 'Liver', tissueSiteDetail: 'Liver', tissueSite: 'Liver', eGeneCount: 9000, sGeneCount: 2000, rnaSeqSampleSummary: { totalCount: 262, female: { count: 80 }, male: { count: 182 } },
    },
    {
      tissueSiteDetailId: 'Liver_Hepatocyte', tissueSiteDetail: 'Liver - Hepatocyte', tissueSite: 'Liver', eGeneCount: null, rnaSeqSampleSummary: { totalCount: 30, female: { count: 10 }, male: { count: 20 } },
    },
  ];
  assert.deepEqual(tissueSites(rows), [
    {
      id: 'Liver', name: 'Liver', organ: 'Liver', samples: 262, female: 80, male: 182, eGenes: 9000, sGenes: 2000,
    },
    {
      id: 'Liver_Hepatocyte', name: 'Liver - Hepatocyte', organ: 'Liver', samples: 30, female: 10, male: 20, eGenes: null, sGenes: null,
    },
  ]);
});

test('pfamNames attaches InterPro names to the top domains', () => {
  const top = [{ id: 'PF00069', n: 3 }, { id: 'PF99999', n: 1 }];
  const named = pfamNames(top, { PF00069: { metadata: { accession: 'PF00069', name: { name: 'Protein kinase domain', short: 'Pkinase' } } } });
  assert.deepEqual(named, [
    {
      id: 'PF00069', n: 3, name: 'Protein kinase domain', short: 'Pkinase',
    },
    {
      id: 'PF99999', n: 1, name: null, short: null,
    },
  ]);
});

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

const HPA_HEAD = ['Gene', 'Ensembl', 'Protein class', 'Evidence', 'RNA tissue specificity', 'RNA tissue specific nTPM', 'RNA single cell type specificity', 'RNA single cell type specific nCPM', 'RNA blood cell specificity', 'RNA brain regional specificity', 'RNA cell line specificity', 'RNA cancer specificity', 'Secretome location', 'Subcellular main location', 'Blood concentration - Conc. blood IM [pg/L]', 'Blood concentration - Conc. blood MS [pg/L]'].map((h) => (h.includes(' ') ? `"${h}"` : h)).join('\t');
const hpa = (o) => [o.gene, 'ENSG1', o.cls ?? '', o.ev ?? 'Evidence at protein level', o.spec, o.ntpm ?? '', o.sc ?? '', o.ncpm ?? '', o.blood ?? '', o.brain ?? '', o.line ?? '', o.cancer ?? '', o.secretome ?? '', o.loc ?? '', o.im ?? 'NA', o.ms ?? 'NA'].map((v) => (String(v).includes(' ') ? `"${v}"` : v)).join('\t');

test('hpaSummary sorts each atlas into enriched, group enriched, enhanced, low and not detected', () => {
  const h = hpaSummary([
    HPA_HEAD,
    hpa({ gene: 'A', spec: 'Tissue enriched', sc: 'Cell type enhanced', blood: 'Immune cell enriched', brain: 'Region enriched', line: 'Cell line enhanced', cancer: 'Low cancer specificity' }),
    hpa({ gene: 'B', spec: 'Group enriched', sc: 'Group enriched', blood: 'Not detected in immune cells', brain: 'Low region specificity', line: 'Not detected', cancer: 'Cancer enhanced' }),
    hpa({ gene: 'C', spec: 'Not detected', sc: 'Low cell type specificity', blood: 'Low immune cell specificity', brain: 'Not detected', cancer: 'NA' }),
  ].join('\n'));
  assert.deepEqual(h.contexts.tissue, {
    enriched: 1, group: 1, enhanced: 0, low: 0, none: 1,
  });
  assert.deepEqual(h.contexts.singleCell, {
    enriched: 0, group: 1, enhanced: 1, low: 1, none: 0,
  });
  assert.deepEqual(h.contexts.blood, {
    enriched: 1, group: 0, enhanced: 0, low: 1, none: 1,
  });
  assert.deepEqual(h.contexts.brain, {
    enriched: 1, group: 0, enhanced: 0, low: 1, none: 1,
  });
  assert.deepEqual(h.contexts.cellLine, {
    enriched: 0, group: 0, enhanced: 1, low: 0, none: 1,
  });
  assert.deepEqual(h.contexts.cancer, {
    enriched: 0, group: 0, enhanced: 1, low: 1, none: 0,
  });
  assert.equal(h.genes, 3);
});

test('hpaSummary bins plasma concentrations by decade for immunoassay and mass spectrometry', () => {
  const h = hpaSummary([
    HPA_HEAD,
    hpa({ gene: 'ALB', spec: 'Tissue enriched', im: '40000000000000' }),
    hpa({ gene: 'IL6', spec: 'Low tissue specificity', im: '100000' }),
    hpa({ gene: 'APOA1', spec: 'Tissue enriched', ms: '2.5e10' }),
  ].join('\n'));
  assert.equal(h.blood.edges[0], 1);
  assert.equal(h.blood.edges.at(-1), 1e15);
  assert.equal(h.blood.im.count, 2);
  assert.equal(h.blood.im.counts[13], 1);
  assert.equal(h.blood.im.counts[5], 1);
  assert.equal(h.blood.ms.count, 1);
  assert.equal(h.blood.ms.counts[10], 1);
  assert.deepEqual(h.blood.top[0], { gene: 'ALB', im: 4e13, ms: null });
  assert.deepEqual(h.blood.top.map((t) => t.gene), ['ALB', 'APOA1', 'IL6']);
  assert.equal(h.blood.proteins, 3);
  assert.deepEqual(h.blood.lowest, { gene: 'IL6', im: 100000, ms: null });
});

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

test('gtexSummary correlates tissues over genes reaching 1 TPM somewhere and orders similar tissues together', () => {
  const head = ['#1.2', '6\t4', 'Name\tDescription\tLiver\tTestis\tBrain_Cortex\tBrain_Cerebellum'];
  const rows = [
    'ENSG1.1\tALB\t5000\t0.5\t0.1\t0.2',
    'ENSG2.1\tPRM1\t0\t900\t0\t0',
    'ENSG3.1\tACTB\t500\t400\t300\t350',
    'ENSG4.1\tGFAP\t1\t2\t800\t600',
  ];
  const a = gtexSummary([...head, ...rows]);
  const b = gtexSummary([...head, ...rows, 'ENSG5.1\tOFF\t0\t0\t0\t0', 'ENSG6.1\tLOW\t0.2\t0.1\t0\t0.3']);
  assert.deepEqual(b.correlation, a.correlation);
  assert.equal(b.genes, 6);
  assert.equal(b.correlated, 4);
  const at = (t) => a.order.indexOf(a.tissues.indexOf(t));
  assert.equal(Math.abs(at('Brain_Cortex') - at('Brain_Cerebellum')), 1);
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

test('pharosSummary splits each development level by IDG family, largest family first', () => {
  const family = (values) => ({ data: { targets: { count: values.reduce((n, v) => n + v.value, 0), facets: [{ facet: 'IDG Family', values }] } } });
  const p = pharosSummary({ data: { targets: { count: 4, facets: [] } } }, {
    Tclin: family([{ name: 'GPCR', value: 1 }]),
    Tdark: family([{ name: 'Other', value: 1 }, { name: 'GPCR', value: 2 }]),
  });
  assert.deepEqual(p.families, [
    {
      family: 'GPCR', Tclin: 1, Tchem: 0, Tbio: 0, Tdark: 2,
    },
    {
      family: 'Other', Tclin: 0, Tchem: 0, Tbio: 0, Tdark: 1,
    },
  ]);
});
