import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  gwasAncestry, aadrSummary, majorHaplogroup, igsrSamples, assemblySummary,
} from '../../scripts/populations.mjs';

const GWAS_HEAD = ['STUDY ACCESSION', 'PUBMEDID', 'FIRST AUTHOR', 'DATE', 'INITIAL SAMPLE DESCRIPTION', 'REPLICATION SAMPLE DESCRIPTION', 'STAGE', 'NUMBER OF INDIVDUALS', 'BROAD ANCESTRAL CATEGORY', 'COUNTRY OF ORIGIN', 'COUNTRY OF RECRUITMENT', 'ADDITONAL ANCESTRY DESCRIPTION'];
const gwas = (acc, date, stage, n, category, recruited) => [acc, '1', 'A', date, 'x', 'NA', stage, n, category, 'NR', recruited, ''].join('\t');

test('gwasAncestry sums participants by ancestral category, year, stage and country of recruitment', () => {
  const g = gwasAncestry([
    GWAS_HEAD.join('\t'),
    gwas('GCST1', '2016-12-13', 'initial', '15124', 'Hispanic or Latin American', 'U.S., NR'),
    gwas('GCST1', '2016-12-13', 'initial', '11731', 'African American or Afro-Caribbean', 'U.S.'),
    gwas('GCST2', '2017-01-05', 'initial', '400000', 'European', 'U.K.'),
    gwas('GCST2', '2017-01-05', 'replication', '5000', 'East Asian', 'Japan'),
    gwas('GCST3', '2017-03-01', 'initial', '1000', 'European, East Asian', 'U.K., China'),
    gwas('GCST4', '2018-02-02', 'initial', '', 'NR', 'NR'),
  ].join('\n'));
  assert.equal(g.rows, 6);
  assert.equal(g.studies, 4);
  assert.equal(g.participants, 432855);
  assert.deepEqual(g.byYear, {
    2016: { 'Hispanic or Latin American': 15124, 'African American or Afro-Caribbean': 11731 },
    2017: { European: 400000, 'Multiple categories': 1000 },
  });
  assert.deepEqual(g.studiesByYear, { 2016: 1, 2017: 2, 2018: 1 });
  assert.deepEqual(g.stages.initial, {
    'Hispanic or Latin American': 15124, 'African American or Afro-Caribbean': 11731, European: 400000, 'Multiple categories': 1000,
  });
  assert.deepEqual(g.stages.replication, { 'East Asian': 5000 });
  assert.deepEqual(g.countries, [
    { country: 'U.K.', participants: 401000 }, { country: 'U.S.', participants: 26855 }, { country: 'Japan', participants: 5000 }, { country: 'China', participants: 1000 },
  ]);
});

const AADR_HEAD = [
  'Genetic ID (suffices: ".DG" is a high coverage shotgun genome)', 'Persistent Genetic ID', 'Individual ID', 'First publication: Abbreviation for earliest paper',
  '"Date mean in BP in years before 1950 CE [OxCal mu]"', 'Political Entity', 'Latitude', 'Longitude', 'Data type',
  'SNPs hit on autosomal targets (Computed using easystats on 1240k snpset)', 'Molecular Sex',
  'Y haplogroup  in ISOGG notation automatically called based on Yfull 12.03', 'mtDNA haplogroup if >2x or published', 'ASSESSMENT',
];
const aadr = (...cells) => cells.join('\t');

test('majorHaplogroup keeps the top clade of a Y or mtDNA call', () => {
  assert.equal(majorHaplogroup('R1b1a1b', 'y'), 'R');
  assert.equal(majorHaplogroup('n/a (female)', 'y'), null);
  assert.equal(majorHaplogroup('..', 'mt'), null);
  assert.equal(majorHaplogroup('L3e2', 'mt'), 'L3');
  assert.equal(majorHaplogroup('HV0a', 'mt'), 'HV');
  assert.equal(majorHaplogroup('U5a1a1', 'mt'), 'U');
});

test('aadrSummary counts each ancient individual once by age, first publication, place, data, sex and haplogroup', () => {
  const a = aadrSummary([
    AADR_HEAD.join('\t'),
    aadr('I1.SG', 'I1.SG', 'I1', 'RasmussenNature2010', '4000', 'Greenland', '70.0', '-50.0', 'Shotgun', '900000', 'M', 'Q1a', 'D2a1', 'PASS'),
    aadr('I2.AG', 'I2.AG', 'I2', 'MathiesonNature2015', '7500', 'Hungary', '47.5', '19.1', '1240K', '500000', 'F', 'n/a (female)', 'U5a1', 'PASS'),
    aadr('I2.SG', 'I2.SG', 'I2', 'MathiesonNature2015', '7500', 'Hungary', '47.5', '19.1', 'Shotgun', '700000', 'F', 'n/a (female)', 'U5a1', 'PASS'),
    aadr('I3.DG', 'I3.DG', 'I3', 'FuNature2014', '45000', 'Russia', '57.0', '73.0', 'Shotgun', '1000000', 'M', 'K2a', 'R*', 'PASS'),
    aadr('I4.HO', 'I4.HO', 'I4', 'LazaridisNature2014', '0', 'France', '48.8', '2.3', 'Genotyping', '..', 'F', 'n/a (female)', '..', 'PASS'),
    aadr('I5.AG', 'I5.AG', 'I5', 'PattersonNature2022', '3000', 'United Kingdom', '51.5', '-0.1', '1240K', '30000', 'U', '..', '..', 'QUESTIONABLE'),
  ].join('\n'));
  assert.equal(a.individuals, 4);
  assert.equal(a.present, 1);
  assert.equal(a.ages.counts.reduce((x, y) => x + y, 0), 4);
  assert.equal(a.ages.edges[0], 100);
  assert.deepEqual(a.oldest, { id: 'I3', age: 45000, place: 'Russia' });
  assert.deepEqual(a.published, {
    2010: 1, 2014: 1, 2015: 1, 2022: 1,
  });
  assert.deepEqual(a.countries.map((c) => c.country), ['Greenland', 'Hungary', 'Russia', 'United Kingdom']);
  assert.deepEqual(a.dataTypes, { Shotgun: 2, '1240K': 2 });
  assert.deepEqual(a.sex, { M: 2, F: 1, U: 1 });
  assert.deepEqual(a.y, [{ group: 'K', n: 1 }, { group: 'Q', n: 1 }]);
  assert.deepEqual(a.mt, [{ group: 'D', n: 1 }, { group: 'R', n: 1 }, { group: 'U', n: 1 }]);
  assert.deepEqual(a.assessment, { PASS: 3, QUESTIONABLE: 1 });
  assert.equal(a.snps.counts.reduce((x, y) => x + y, 0), 4);
  assert.deepEqual(a.grid, [
    { lat: 46, lon: 18, n: 1, age: 7500 }, { lat: 50, lon: -2, n: 1, age: 3000 }, { lat: 56, lon: 72, n: 1, age: 45000 }, { lat: 70, lon: -50, n: 1, age: 4000 },
  ]);
});

const sample = (name, sex, pop, superpop, collections) => ({
  _source: {
    name, sex, populations: [{ code: pop, name: `${pop} name`, superpopulationCode: superpop, superpopulationName: `${superpop} ancestry` }], dataCollections: collections.map((title) => ({ title })),
  },
});

test('igsrSamples counts samples per population in one data collection and samples per collection', () => {
  const s = igsrSamples({
    hits: {
      total: 4,
      hits: [
        sample('HG1', 'female', 'GBR', 'EUR', ['1000 Genomes on GRCh38', '1000 Genomes 30x on GRCh38']),
        sample('HG2', 'male', 'GBR', 'EUR', ['1000 Genomes 30x on GRCh38']),
        sample('NA1', 'female', 'YRI', 'AFR', ['1000 Genomes 30x on GRCh38', 'Human Genome Structural Variation Consortium, Phase 3']),
        sample('HGDP1', 'male', 'Han', 'EAS', ['Human Genome Diversity Project']),
      ],
    },
  }, '1000 Genomes 30x on GRCh38');
  assert.equal(s.total, 4);
  assert.deepEqual(s.collections, [
    { title: '1000 Genomes 30x on GRCh38', n: 3 }, { title: '1000 Genomes on GRCh38', n: 1 }, { title: 'Human Genome Diversity Project', n: 1 }, { title: 'Human Genome Structural Variation Consortium, Phase 3', n: 1 },
  ]);
  assert.equal(s.panel.samples, 3);
  assert.deepEqual(s.panel.populations, [
    {
      code: 'GBR', name: 'GBR name', superpopulation: 'EUR', superpopulationName: 'EUR ancestry', n: 2,
    },
    {
      code: 'YRI', name: 'YRI name', superpopulation: 'AFR', superpopulationName: 'AFR ancestry', n: 1,
    },
  ]);
  assert.deepEqual(s.panel.sex, { female: 2, male: 1 });
});

const assembly = (accession, level, date, submitter, n50, project = 'PRJNA1') => ({
  accession,
  assembly_info: {
    assembly_level: level, release_date: date, submitter, assembly_name: `${accession} name`, bioproject_accession: project,
  },
  assembly_stats: { contig_n50: n50, total_sequence_length: '3100000000' },
});

test('assemblySummary counts human assemblies by year, level and submitter and keeps each contig N50', () => {
  const a = assemblySummary([
    assembly('GCA_1', 'Chromosome', '2009-02-27', 'Genome Reference Consortium', 38440852),
    assembly('GCA_2', 'Complete Genome', '2022-01-24', 'T2T Consortium', 154259566),
    assembly('GCA_3', 'Chromosome', '2023-05-10', 'Human Pangenome Reference Consortium', 40000000, 'PRJNA730823'),
    assembly('GCA_4', 'Scaffold', '2023-07-01', 'Human Pangenome Reference Consortium', 35000000, 'PRJNA730823'),
    assembly('GCA_5', 'Contig', '2015-03-03', 'Some Lab', 1200000),
  ]);
  assert.equal(a.total, 5);
  assert.deepEqual(a.byYear, {
    2009: { Chromosome: 1 }, 2015: { Contig: 1 }, 2022: { 'Complete Genome': 1 }, 2023: { Chromosome: 1, Scaffold: 1 },
  });
  assert.deepEqual(a.levels, {
    Chromosome: 2, 'Complete Genome': 1, Scaffold: 1, Contig: 1,
  });
  assert.deepEqual(a.submitters[0], { submitter: 'Human Pangenome Reference Consortium', n: 2 });
  assert.equal(a.n50.length, 5);
  assert.deepEqual(a.n50[1], ['2022-01-24', 154259566, 'Complete Genome', 'GCA_2']);
  assert.deepEqual(a.best, { accession: 'GCA_2', name: 'GCA_2 name', n50: 154259566, submitter: 'T2T Consortium' });
});
