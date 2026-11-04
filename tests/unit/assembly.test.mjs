import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  parseAssemblyReport, summarizeRoles, parseAssemblyStats, parseCytobands, parseGaps, parseCentromeres, buildAssembly,
} from '../../scripts/parsers.mjs';

const REPORT = [
  '# Assembly name:  GRCh38.p14',
  '# Date:           2022-02-03',
  '# GenBank assembly accession: GCA_000001405.29',
  '# RefSeq assembly accession: GCF_000001405.40',
  '## Assembly-Unit names:',
  '## Primary Assembly',
  '## PATCHES',
  '## ALT_REF_LOCI_1',
  '## ALT_REF_LOCI_2',
  '## non-nuclear',
  '# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name',
  '1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\tNC_000001.11\tPrimary Assembly\t248956422\tchr1',
  'HSCHR1_CTG1_UNLOCALIZED\tunlocalized-scaffold\t1\tChromosome\tKI270706.1\t=\tNT_187361.1\tPrimary Assembly\t175055\tchr1_KI270706v1_random',
  'HSCHRUN_RANDOM_CTG1\tunplaced-scaffold\tna\tna\tKI270302.1\t=\tNT_187395.1\tPrimary Assembly\t2274\tchrUn_KI270302v1',
  'HSCHR6_MHC_APD_CTG1\talt-scaffold\t6\tChromosome\tGL000250.2\t=\tNT_167244.2\tALT_REF_LOCI_2\t4672374\tchr6_GL000250v2_alt',
  'HG2266_PATCH\tfix-patch\t1\tChromosome\tKN538360.1\t=\tNW_009646194.1\tPATCHES\t460100\tchr1_KN538360v1_fix',
  'HSCHR1_5_CTG3\tnovel-patch\t1\tChromosome\tKQ983255.1\t=\tNW_015495298.1\tPATCHES\t165050\tchr1_KQ983255v1_alt',
  'MT\tassembled-molecule\tMT\tMitochondrion\tJ01415.2\t=\tNC_012920.1\tnon-nuclear\t16569\tchrM',
].join('\n');

test('assembly reports give the header and one row per sequence', () => {
  const r = parseAssemblyReport(REPORT);
  assert.equal(r.name, 'GRCh38.p14');
  assert.equal(r.date, '2022-02-03');
  assert.equal(r.genbank, 'GCA_000001405.29');
  assert.equal(r.refseq, 'GCF_000001405.40');
  assert.equal(r.altLoci, 2);
  assert.equal(r.sequences.length, 7);
  assert.deepEqual(r.sequences[0], {
    name: '1', role: 'assembled-molecule', molecule: '1', genbank: 'CM000663.2', refseq: 'NC_000001.11', unit: 'Primary Assembly', length: 248956422, ucsc: 'chr1',
  });
  assert.equal(r.sequences.at(-1).molecule, 'MT');
});

test('summarizeRoles counts sequences and bases by role in a fixed order', () => {
  const roles = summarizeRoles(parseAssemblyReport(REPORT).sequences);
  assert.deepEqual(roles.map((x) => x.role), ['assembled-molecule', 'unlocalized-scaffold', 'unplaced-scaffold', 'alt-scaffold', 'fix-patch', 'novel-patch']);
  assert.deepEqual(roles[0], { role: 'assembled-molecule', count: 2, length: 248972991 });
  assert.deepEqual(roles[3], { role: 'alt-scaffold', count: 1, length: 4672374 });
});

const STATS = [
  '# Assembly name:  GRCh38.p14',
  '# unit-name\tmolecule-name\tmolecule-type/loc\tsequence-type\tstatistic\tvalue',
  'all\tall\tall\tall\ttotal-length\t3298430636',
  'all\tall\tall\tall\ttotal-gap-length\t161611139',
  'all\tall\tall\tall\tcontig-N50\t54806562',
  'all\tall\tall\tall\tgc-perc\t41',
  'Primary Assembly\t1\tChromosome\tall\ttotal-length\t249698942',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\ttotal-length\t248956422',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\tungapped-length\t230481121',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\tspanned-gaps\t64',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\tunspanned-gaps\t13',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\tgc-perc\t41.5',
  'Primary Assembly\t1\tChromosome\tunlocalized-scaffold\ttotal-length\t742520',
  'PATCHES\t1\tChromosome\tall\ttotal-length\t7866676',
  'non-nuclear\tMT\tMitochondrion\tassembled-molecule\ttotal-length\t16569',
  'non-nuclear\tMT\tMitochondrion\tassembled-molecule\tungapped-length\t16569',
].join('\n');

test('assembly stats give genome totals and assembled chromosomes', () => {
  const s = parseAssemblyStats(STATS);
  assert.deepEqual(s.total, { totalLength: 3298430636, totalGapLength: 161611139, contigN50: 54806562, gcPerc: 41 });
  assert.deepEqual(s.chromosomes, [
    { chrom: '1', totalLength: 248956422, ungappedLength: 230481121, spannedGaps: 64, unspannedGaps: 13, gcPerc: 41.5 },
    { chrom: 'MT', totalLength: 16569, ungappedLength: 16569 },
  ]);
});

const CYTO = [
  'chr1\t0\t2300000\tp36.33\tgneg',
  'chr1\t121700000\t123400000\tp11.1\tacen',
  'chr1_KI270706v1_random\t0\t175055\t\tgneg',
  'chrM\t0\t16569\t\tgneg',
  'chr21\t0\t3100000\tp13\tgvar',
  'chr21\t3100000\t7000000\tp12\tstalk',
].join('\n');

test('cytobands keep named bands on the 24 chromosomes', () => {
  assert.deepEqual(parseCytobands(CYTO), [
    { chrom: '1', start: 0, end: 2300000, name: 'p36.33', stain: 'gneg' },
    { chrom: '1', start: 121700000, end: 123400000, name: 'p11.1', stain: 'acen' },
    { chrom: '21', start: 0, end: 3100000, name: 'p13', stain: 'gvar' },
    { chrom: '21', start: 3100000, end: 7000000, name: 'p12', stain: 'stalk' },
  ]);
});

const GAPS = [
  '585\tchr1\t0\t10000\t1\tN\t10000\ttelomere\tno',
  '586\tchr1\t207666\t257666\t5\tN\t50000\tcontig\tno',
  '587\tchr1\t297968\t347968\t7\tN\t50000\tcontig\tno',
  '588\tchr1\t125184587\t143184587\t9\tN\t18000000\theterochromatin\tno',
  '10\tchr1_KI270706v1_random\t0\t100\t1\tN\t100\tcontig\tno',
  '11\tchrY\t10316945\t10544039\t3\tN\t227094\tcontig\tno',
].join('\n');

test('gaps sum by chromosome and type', () => {
  assert.deepEqual(parseGaps(GAPS), {
    1: { telomere: { count: 1, bp: 10000 }, contig: { count: 2, bp: 100000 }, heterochromatin: { count: 1, bp: 18000000 } },
    Y: { contig: { count: 1, bp: 227094 } },
  });
});

const CEN = [
  '189\tchr1\t122026459\t122224535\tGJ211836.1',
  '189\tchr1\t122224635\t122503147\tGJ211837.1',
  '23\tchr2\t92188145\t94090557\tGJ211855.1',
].join('\n');

test('centromere models sum to bases per chromosome', () => {
  assert.deepEqual(parseCentromeres(CEN), { 1: 476588, 2: 1902412 });
});

const T2T_REPORT = [
  '# Assembly name:  T2T-CHM13v2.0',
  '# Date:           2022-01-24',
  '# GenBank assembly accession: GCA_009914755.4',
  '# RefSeq assembly accession: GCF_009914755.1',
  '1\tassembled-molecule\t1\tChromosome\tCP068277.2\t=\tNC_060925.1\tPrimary Assembly\t248387328\tchr1',
  'MT\tassembled-molecule\tMT\tMitochondrion\tCP068254.1\t<>\tna\tnon-nuclear\t16569\tchrM',
].join('\n');

const T2T_STATS = [
  'all\tall\tall\tall\ttotal-length\t3117292070',
  'all\tall\tall\tall\ttotal-gap-length\t0',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\ttotal-length\t248387328',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\tungapped-length\t248387328',
  'Primary Assembly\t1\tChromosome\tassembled-molecule\tgc-perc\t41.5',
].join('\n');

test('buildAssembly joins reports and stats into one record per chromosome', () => {
  const a = buildAssembly({
    report38: REPORT, stats38: STATS, reportT2t: T2T_REPORT, statsT2t: T2T_STATS, cyto: CYTO, gaps: GAPS, cen: CEN,
  });
  assert.equal(a.grch38.name, 'GRCh38.p14');
  assert.equal(a.grch38.altLoci, 2);
  assert.equal(a.grch38.roles.length, 6);
  assert.equal(a.grch38.total.totalGapLength, 161611139);
  assert.deepEqual(a.grch38.chromosomes, [
    { chrom: '1', length: 248956422, ungapped: 230481121, gc: 41.5, spanned: 64, unspanned: 13, refseq: 'NC_000001.11', genbank: 'CM000663.2' },
    { chrom: 'MT', length: 16569, ungapped: 16569, gc: null, spanned: null, unspanned: null, refseq: 'NC_012920.1', genbank: 'J01415.2' },
  ]);
  assert.equal(a.t2t.name, 'T2T-CHM13v2.0');
  assert.deepEqual(a.t2t.chromosomes, [
    { chrom: '1', length: 248387328, ungapped: 248387328, gc: 41.5, spanned: null, unspanned: null, refseq: 'NC_060925.1', genbank: 'CP068277.2' },
  ]);
  assert.equal(a.bands.length, 4);
  assert.deepEqual(a.gaps.Y, { contig: { count: 1, bp: 227094 } });
  assert.equal(a.centromeres[1], 476588);
});
