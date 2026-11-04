import { test } from 'node:test';
import assert from 'node:assert/strict';
import { hgncSummary, geneInfoSummary, maneSummary } from '../../scripts/catalogs.mjs';

const HEAD = ['hgnc_id', 'symbol', 'name', 'locus_group', 'locus_type', 'status', 'location', 'alias_symbol', 'alias_name', 'prev_symbol', 'prev_name',
  'gene_group', 'gene_group_id', 'date_approved_reserved', 'date_symbol_changed', 'date_name_changed', 'date_modified', 'entrez_id', 'mane_select'];
const hg = (o) => HEAD.map((k) => o[k] ?? '').join('\t');

const HGNC = [
  HEAD.join('\t'),
  hg({ hgnc_id: 'HGNC:1', symbol: 'SEPTIN7', locus_group: 'protein-coding gene', locus_type: 'gene with protein product', status: 'Approved', prev_symbol: 'SEPT7', gene_group: 'Septins', date_approved_reserved: '1999-01-21', date_symbol_changed: '2020-03-04', alias_symbol: 'CDC10' }),
  hg({ hgnc_id: 'HGNC:2', symbol: 'MARCHF1', locus_group: 'protein-coding gene', locus_type: 'gene with protein product', status: 'Approved', prev_symbol: 'MARCH1', gene_group: 'Ring finger proteins|MARCH family', date_approved_reserved: '2004-05-01', date_symbol_changed: '2020-06-10' }),
  hg({ hgnc_id: 'HGNC:3', symbol: 'C1orf112', locus_group: 'protein-coding gene', locus_type: 'gene with protein product', status: 'Approved', date_approved_reserved: '2004-01-01' }),
  hg({ hgnc_id: 'HGNC:4', symbol: 'FAM20A', locus_group: 'protein-coding gene', locus_type: 'gene with protein product', status: 'Approved', date_approved_reserved: '1999-07-01', alias_symbol: 'A|B|C' }),
  hg({ hgnc_id: 'HGNC:5', symbol: 'GAPDHP1', locus_group: 'pseudogene', locus_type: 'pseudogene', status: 'Approved', date_approved_reserved: '1999-02-02' }),
  hg({ hgnc_id: 'HGNC:6', symbol: 'GAPDHP2', locus_group: 'pseudogene', locus_type: 'pseudogene', status: 'Approved', date_approved_reserved: '2001-02-02' }),
  hg({ hgnc_id: 'HGNC:7', symbol: 'GAPDH', locus_group: 'protein-coding gene', locus_type: 'gene with protein product', status: 'Approved', date_approved_reserved: '1986-01-01', gene_group: 'Ring finger proteins' }),
  hg({ hgnc_id: 'HGNC:8', symbol: 'MIR21', locus_group: 'non-coding RNA', locus_type: 'RNA, micro', status: 'Approved', date_approved_reserved: '2003-01-01', prev_symbol: 'MIRN21|MIR-21' }),
].join('\n');

test('hgncSummary counts approved symbols by group and year', () => {
  const h = hgncSummary(HGNC);
  assert.equal(h.approved, 8);
  assert.deepEqual(h.locusGroups, { 'protein-coding gene': 5, pseudogene: 2, 'non-coding RNA': 1 });
  assert.deepEqual(h.perYear.approved, [[1986, 1], [1999, 3], [2001, 1], [2003, 1], [2004, 2]]);
  assert.deepEqual(h.perYear.symbolChanged, [[2020, 2]]);
});

test('hgncSummary finds placeholders, groups, renames, pseudogene parents and spreadsheet renames', () => {
  const h = hgncSummary(HGNC);
  assert.deepEqual(h.placeholders, { orf: 1, fam: 1, kiaa: 0 });
  assert.deepEqual(h.groups.slice(0, 1), [{ name: 'Ring finger proteins', n: 2 }]);
  assert.deepEqual(h.renamed[0], { symbol: 'MIR21', prev: ['MIRN21', 'MIR-21'] });
  assert.deepEqual(h.aliases[0], { symbol: 'FAM20A', n: 3 });
  assert.deepEqual(h.pseudoParents, [{ parent: 'GAPDH', n: 2 }]);
  assert.deepEqual(h.spreadsheet, { count: 2, examples: ['SEPT7 → SEPTIN7', 'MARCH1 → MARCHF1'] });
});

const INFO = [
  '#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\tchromosome\tmap_location\tdescription\ttype_of_gene\tSymbol_from_nomenclature_authority\tFull_name_from_nomenclature_authority\tNomenclature_status\tOther_designations\tModification_date\tFeature_type',
  '9606\t1\tA1BG\t-\tA1B|ABG|GAB|HYST2477\tMIM:138670\t19\t19q13.43\talpha-1-B glycoprotein\tprotein-coding\tA1BG\talpha-1-B glycoprotein\tO\t-\t20261011\t-',
  '9606\t2\tA2M\t-\tA2MD\t-\t12\t12p13.31\talpha-2-macroglobulin\tprotein-coding\tA2M\t-\tO\t-\t20261011\t-',
  '9606\t3\tLOC1\t-\t-\t-\t1\t-\tuncharacterized LOC1\tncRNA\t-\t-\t-\t-\t20261011\t-',
  '9606\t4\tP1\t-\t-\t-\t1\t-\tpseudo\tpseudo\t-\t-\t-\t-\t20261011\t-',
].join('\n');

test('geneInfoSummary counts NCBI gene records by type and finds the most synonyms', () => {
  const g = geneInfoSummary(INFO);
  assert.equal(g.total, 4);
  assert.deepEqual(g.types, { 'protein-coding': 2, ncRNA: 1, pseudo: 1 });
  assert.deepEqual(g.synonyms[0], { symbol: 'A1BG', n: 4 });
  assert.equal(g.uncharacterized, 1);
});

test('maneSummary counts MANE Select and Plus Clinical transcripts', () => {
  const text = ['#NCBI_GeneID\tEnsembl_Gene\tHGNC_ID\tsymbol\tname\tRefSeq_nuc\tRefSeq_prot\tEnsembl_nuc\tEnsembl_prot\tMANE_status\tGRCh38_chr\tchr_start\tchr_end\tchr_strand',
    'GeneID:1\tE\tHGNC:5\tA1BG\tx\tNM_1\tNP_1\tENST1\tENSP1\tMANE Select\tNC_000019.10\t1\t2\t-',
    'GeneID:2\tE\tHGNC:6\tA2M\tx\tNM_2\tNP_2\tENST2\tENSP2\tMANE Select\tNC_000012.12\t1\t2\t-',
    'GeneID:2\tE\tHGNC:6\tA2M\tx\tNM_3\tNP_3\tENST3\tENSP3\tMANE Plus Clinical\tNC_000012.12\t1\t2\t-'].join('\n');
  assert.deepEqual(maneSummary(text), { select: 2, plusClinical: 1, genes: 2 });
});
