import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  SOURCES, SOURCE_TITLES, SOURCE_GROUPS, EXAMPLES, PARS, GENE_THERAPIES, GENE_THERAPIES_CHECKED, DIAGNOSTIC_YIELD, CFDNA_SCREENING, CARRIER_TIER3, BIOBANKS, GNOMAD_V4,
  ZOONOMIA, DE_NOVO, ARCHAIC_GENOMES, ARCHAIC_ANCESTRY, APE_GENOMES,
} from '../../src/facts.js';

test('evolution sources are cited in their own group after complex traits, with the license of the Alliance data', () => {
  for (const key of ['ncbiOrthologs', 'sayers2025', 'ncbiTaxonomy', 'schoch2020', 'alliance', 'alliance2024', 'allianceLicense', 'ensemblCompara', 'dyer2025', 'vilella2009',
    'ucscCons470', 'siepel2005', 'zoonomia2020', 'christmas2023', 'sullivan2023', 'keough2023', 'xue2023', 'jonsson2017', 'green2010', 'meyer2012', 'prufer2014', 'prufer2017',
    'mafessoni2020', 'sankararaman2016', 'yoo2025']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Evolution', key);
  }
  assert.equal(SOURCE_GROUPS.indexOf('Evolution'), SOURCE_GROUPS.indexOf('Complex traits') + 1);
  assert.equal(SOURCES.allianceLicense, 'https://creativecommons.org/licenses/by/4.0/');
});

test('Zoonomia, de novo mutation, archaic and ape genome numbers are fixed values from their papers', () => {
  assert.deepEqual(ZOONOMIA, {
    species: 240, constrainedBases: 332000000, constrainedShare: 10.7, significantBases: 101000000, significantShare: 3.3, outsideCoding: 80, unannotated: 50, ultraconserved: 4552, zooHARs: 312, hcondels: 10032, hcondelLength: 2.56, hcondelActive: 800,
  });
  assert.deepEqual(DE_NOVO.slopes, [
    {
      key: 'father', label: "Father's age", value: 1.51, lo: 1.45, hi: 1.57,
    },
    {
      key: 'mother', label: "Mother's age", value: 0.37, lo: 0.32, hi: 0.43,
    },
  ]);
  assert.deepEqual(DE_NOVO.maternal, [
    {
      key: 'cpg', label: 'CpG>TpG', value: -0.26, lo: -0.33, hi: -0.19,
    },
    {
      key: 'cg', label: 'C>G', value: 0.33, lo: 0.28, hi: 0.38,
    },
  ]);
  assert.equal(DE_NOVO.trios, 1548);
  assert.equal(DE_NOVO.mutations, 108778);
  assert.equal(DE_NOVO.phased, 42961);
  assert.equal(ARCHAIC_GENOMES.length, 6);
  for (const g of ARCHAIC_GENOMES) {
    assert.ok(SOURCES[g.source], g.name);
    assert.ok(['Neandertal', 'Denisovan'].includes(g.group), g.name);
    assert.ok(g.coverage > 0 && g.year >= 2010 && g.year <= 2020, g.name);
  }
  assert.deepEqual(ARCHAIC_GENOMES.find((g) => g.year === 2014 && g.coverage === 52), {
    name: 'Altai Neandertal, Denisova Cave toe phalanx', group: 'Neandertal', year: 2014, coverage: 52, source: 'prufer2014', note: null,
  });
  for (const a of ARCHAIC_ANCESTRY) {
    assert.ok(SOURCES[a.source], a.label);
    assert.ok(a.value !== null || (a.lo !== null && a.hi !== null), a.label);
  }
  assert.deepEqual(ARCHAIC_ANCESTRY.map((a) => [a.label, a.kind, a.value, a.lo, a.hi]), [
    ['East Asians', 'Neandertal', null, 2.3, 2.6],
    ['People in Western Eurasia', 'Neandertal', null, 1.8, 2.4],
    ['Non-Africans outside Oceania', 'Neandertal', null, 1.8, 2.6],
    ['Papuans', 'Denisovan', 1.12, null, null],
    ['Sherpas', 'Denisovan', 0.1, null, null],
  ]);
  assert.deepEqual(APE_GENOMES, {
    species: 6, gapless: 215, errorEvery: 2700000, unaligned: [12.5, 27.3], gapFold: [5, 15], inaccessible: [10, 15], mhc: { deletions: 328, insertions: 422 }, spacer: { bonobo: 95.2, siamang: 80.6 },
  });
});

test('chromosome structure sources are cited', () => {
  for (const key of ['anderson1981', 'andrews1999', 'laricchia2022', 'gnomadMito', 'altemose2022', 'rhie2023', 't2tAnnotation']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Assemblies and sequence', key);
  }
});

test('variation and constraint sources are cited', () => {
  for (const key of ['dbsnp', 'sherry2001', 'dbvar', 'gnomad', 'karczewski2020']) assert.equal(SOURCE_TITLES[key]?.[0], 'Variation and populations', key);
  for (const key of ['clingen', 'rehm2015', 'riggs2020']) assert.equal(SOURCE_TITLES[key]?.[0], 'Clinical genetics', key);
  for (const key of ['impc', 'dickinson2016']) assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
});

test('the GRCh38 pseudoautosomal regions match on X and Y', () => {
  assert.deepEqual(PARS.map((p) => p.name), ['PAR1', 'PAR2']);
  for (const p of PARS) assert.equal(p.x[1] - p.x[0], p.y[1] - p.y[0], p.name);
  assert.deepEqual(PARS[0].x, [10001, 2781479]);
  assert.deepEqual(PARS[1].y, [56887903, 57217415]);
});
import { identify } from '../../src/lib.js';

test('clinical variant sources are cited', () => {
  for (const key of ['clinvar', 'landrum2018', 'landrum2025', 'richards2015', 'strande2017', 'clingenValidity', 'dosageMap', 'riggs2012', 'strchive', 'hiatt2025']) assert.equal(SOURCE_TITLES[key]?.[0], 'Clinical genetics', key);
  for (const key of ['mclaren2016', 'alexandrov2020']) assert.equal(SOURCE_TITLES[key]?.[0], 'Variation and populations', key);
});

test('rare disease and testing sources are cited', () => {
  for (const key of ['hpo', 'gargano2024', 'orphanet', 'nguengang2020', 'gtr', 'rubinstein2013', 'clinicaltrials', 'fdaCgt', 'miller2010', 'clark2018', 'srivastava2019', 'gil2017', 'gregg2021']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Clinical genetics', key);
  }
});

test('gene therapies list each FDA product with its first approval year, vector and group', () => {
  assert.equal(GENE_THERAPIES.length, 34);
  assert.equal(new Set(GENE_THERAPIES.map((g) => g.trade)).size, 34);
  const groups = new Set(['AAV', 'Other in vivo vectors', 'CAR-T and TCR-T cells', 'Blood stem cells', 'Other gene-modified cells']);
  for (const g of GENE_THERAPIES) {
    assert.ok(groups.has(g.group), g.trade);
    assert.ok(g.year >= 2015 && g.year <= 2026, g.trade);
    assert.ok(g.vector && g.condition && g.proper, g.trade);
  }
  assert.deepEqual(GENE_THERAPIES.find((g) => g.trade === 'LUXTURNA'), {
    trade: 'LUXTURNA', proper: 'voretigene neparvovec-rzyl', year: 2017, vector: 'AAV2', group: 'AAV', condition: 'biallelic RPE65 retinal dystrophy',
  });
  assert.equal(GENE_THERAPIES.find((g) => g.trade === 'IMLYGIC').year, 2015);
  assert.equal(GENE_THERAPIES.filter((g) => g.group === 'AAV').length, 11);
  assert.equal(GENE_THERAPIES_CHECKED, '2026-11-04');
});

test('diagnostic yields, cfDNA screening and carrier screening come from cited meta-analyses and guidance', () => {
  for (const y of DIAGNOSTIC_YIELD) {
    assert.ok(['miller2010', 'clark2018', 'srivastava2019'].includes(y.source), y.key);
    assert.ok(y.value !== null || (y.lo !== null && y.hi !== null), y.key);
    if (y.lo !== null) assert.ok(y.lo <= (y.value ?? y.lo) && (y.value ?? y.hi) <= y.hi, y.key);
  }
  assert.deepEqual(DIAGNOSTIC_YIELD.find((y) => y.key === 'wgs'), {
    key: 'wgs', test: 'Genome sequencing', measure: 'diagnostic', value: 41, lo: 34, hi: 48, source: 'clark2018',
  });
  assert.equal(CFDNA_SCREENING.length, 5);
  assert.deepEqual(CFDNA_SCREENING[0], {
    condition: 'Trisomy 21', detection: 99.7, lo: 99.1, hi: 99.9, fpr: 0.04, fprLo: 0.02, fprHi: 0.07, cases: 1963, unaffected: 223932,
  });
  assert.deepEqual(CARRIER_TIER3, { genes: 97, definitive: 39, frequency: '1 in 200' });
});

test('population and ancient genome sources are cited', () => {
  for (const key of ['gwasCatalog', 'cerezo2025', 'popejoy2016', 'aadr', 'mallick2024', 'igsr', 'fairley2020', 'thousandGenomes2015', 'byrskaBishop2022', 'ncbiDatasets', 'oleary2024', 'liao2023', 'gnomadV4', 'bycroft2018', 'allOfUs2024', 'kurki2023', 'nagai2017', 'chen2011', 'verma2024', 'ziyatdinov2023', 'feng2022']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Variation and populations', key);
  }
});

test('cancer genomics sources are cited in their own group after clinical genetics', () => {
  for (const key of ['intogen', 'martinezJimenez2020', 'gdc', 'heath2021', 'tcga2013', 'bailey2018', 'lawrence2013', 'alexandrov2013', 'horak2022', 'li2017']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Cancer genomics', key);
  }
  assert.equal(SOURCE_GROUPS.indexOf('Cancer genomics'), SOURCE_GROUPS.indexOf('Clinical genetics') + 1);
});

test('pharmacogenomics sources are cited, with the license of the ClinPGx data', () => {
  for (const key of ['cpic', 'relling2011', 'caudle2014', 'caudle2020', 'gaedigk2017', 'huddart2019', 'fdaPgx', 'clinpgx', 'whirlCarrillo2021', 'clinpgxLicense']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Pharmacogenomics', key);
  }
  assert.equal(SOURCES.clinpgxLicense, 'https://creativecommons.org/licenses/by-sa/4.0/');
  assert.match(SOURCE_TITLES.clinpgx[1], /PharmGKB/);
  assert.match(SOURCE_TITLES.cpic[1], /CPIC®/);
});

test('complex trait sources are cited, with the license of the Pan-UKB results', () => {
  for (const key of ['gwasRelease', 'manolio2009', 'visscher2017', 'pgsCatalog', 'lambert2021', 'lambert2024', 'martin2019', 'panukb', 'karczewski2024', 'bulikSullivan2015', 'panukbLicense']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Complex traits', key);
  }
  assert.equal(SOURCES.panukbLicense, 'https://creativecommons.org/licenses/by/4.0/');
  assert.match(SOURCE_TITLES.karczewski2024[1], /medRxiv/);
});

test('biobank sizes and gnomAD v4 totals are fixed values from their papers and release post', () => {
  assert.equal(BIOBANKS.length, 8);
  for (const b of BIOBANKS) {
    assert.ok(SOURCES[b.source], b.name);
    assert.ok(b.people > 100000 && b.measure, b.name);
  }
  assert.deepEqual(BIOBANKS[0], {
    name: 'Million Veteran Program', people: 635969, qualifier: null, measure: 'participants in its GWAS', source: 'verma2024',
  });
  assert.deepEqual(BIOBANKS.find((b) => b.name === 'UK Biobank'), {
    name: 'UK Biobank', people: 500000, qualifier: 'about', measure: 'participants', source: 'bycroft2018',
  });
  assert.deepEqual(GNOMAD_V4, {
    people: 807162, exomes: 730947, genomes: 76215, ukBiobank: 416555, nonEuropean: 138000,
  });
});

test('every source is an https URL with a title in a known group', () => {
  assert.ok(Object.keys(SOURCES).length >= 6);
  for (const [key, url] of Object.entries(SOURCES)) {
    assert.match(url, /^https:\/\//, key);
    const [group, title] = SOURCE_TITLES[key] ?? [];
    assert.ok(SOURCE_GROUPS.includes(group), `${key}: ${group}`);
    assert.ok(title?.length > 5 && !/^https?:|^www\./.test(title), key);
  }
  for (const key of Object.keys(SOURCE_TITLES)) assert.ok(SOURCES[key], `title without source: ${key}`);
});

test('protein, target and expression sources are cited in their groups', () => {
  const proteins = ['uniprot', 'uniprot2025', 'interpro', 'pfam2021', 'interpro2025', 'pdb', 'burley2019', 'pharos', 'pharos2023', 'oprea2018'];
  const expression = ['hpa', 'uhlen2015', 'karlsson2021', 'uhlen2019', 'gtex', 'gtex2020'];
  for (const key of proteins) assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
  for (const key of expression) assert.equal(SOURCE_TITLES[key]?.[0], 'Expression and regulation', key);
  for (const key of [...proteins, ...expression]) assert.ok(SOURCES[key], key);
});

test('network, pathway and ontology sources are cited', () => {
  for (const key of ['string', 'szklarczyk2023', 'reactome', 'milacic2024', 'geneOntology', 'go2023', 'ashburner2000']) {
    assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
    assert.ok(SOURCES[key], key);
  }
});

test('regulation and RNA gene sources are cited', () => {
  for (const key of ['encode', 'encode2020', 'encode2012', 'gardinerGarden1987']) assert.equal(SOURCE_TITLES[key]?.[0], 'Expression and regulation', key);
  for (const key of ['gtrnadb', 'chan2016']) assert.equal(SOURCE_TITLES[key]?.[0], 'Genes and proteins', key);
});

test('every decoder example is recognised, unique and described', () => {
  assert.ok(EXAMPLES.length >= 18);
  assert.equal(new Set(EXAMPLES.map((e) => e.id)).size, EXAMPLES.length);
  for (const e of EXAMPLES) {
    assert.ok(identify(e.id), e.id);
    assert.ok(e.is.length > 3, e.id);
  }
  const types = new Set(EXAMPLES.map((e) => identify(e.id).type));
  assert.ok(types.size >= 15, [...types].join(', '));
});
