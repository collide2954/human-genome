const table = (text) => {
  const lines = text.split('\n').filter((l) => l.length);
  const head = lines[0].replace(/^#/, '').split('\t');
  return {
    col: Object.fromEntries(head.map((name, i) => [name, i])),
    rows: lines.slice(1).map((l) => l.replace(/\r$/, '').split('\t')),
  };
};

const list = (v) => (v && v !== '-' ? v.replace(/^"|"$/g, '').split('|').map((s) => s.trim()).filter(Boolean) : []);

const tallyBy = (rows, key) => {
  const out = new Map();
  for (const r of rows) {
    const k = key(r);
    if (k !== null && k !== undefined && k !== '') out.set(k, (out.get(k) ?? 0) + 1);
  }
  return out;
};

const byYear = (rows, i) => [...tallyBy(rows, (r) => (r[i] ? Number(r[i].slice(0, 4)) : null))].sort((a, b) => a[0] - b[0]);

const PLACEHOLDERS = { orf: /^C(\d+|X|Y)orf\d+$/, fam: /^FAM\d+/, kiaa: /^KIAA\d+/ };
const SPREADSHEET = /^(SEPT|MARCH|DEC)\d+$/;

export const codingSymbols = (text) => {
  const { col, rows } = table(text);
  return new Set(rows.filter((r) => r[col.status] === 'Approved' && r[col.locus_group] === 'protein-coding gene').map((r) => r[col.symbol]));
};

export const codingGenes = (text) => {
  const { col, rows } = table(text);
  return rows.filter((r) => r[col.status] === 'Approved' && r[col.locus_group] === 'protein-coding gene')
    .map((r) => ({ hgnc: r[col.hgnc_id], entrez: r[col.entrez_id] || null, symbol: r[col.symbol] }));
};

export const hgncSummary = (text) => {
  const { col, rows: all } = table(text);
  const rows = all.filter((r) => r[col.status] === 'Approved');
  const symbols = new Set(rows.map((r) => r[col.symbol]));
  const groups = new Map();
  for (const r of rows) for (const g of list(r[col.gene_group])) groups.set(g, (groups.get(g) ?? 0) + 1);
  const parents = tallyBy(rows.filter((r) => r[col.locus_group] === 'pseudogene'), (r) => {
    const stem = r[col.symbol].match(/^(.+?)P\d+$/)?.[1];
    return stem && symbols.has(stem) ? stem : null;
  });
  const spreadsheet = rows.flatMap((r) => list(r[col.prev_symbol]).filter((p) => SPREADSHEET.test(p)).map((p) => `${p} → ${r[col.symbol]}`));
  const ranked = (fn, n) => rows.map(fn).filter((x) => x.n > 0).sort((a, b) => b.n - a.n).slice(0, n);
  return {
    approved: rows.length,
    locusGroups: Object.fromEntries(tallyBy(rows, (r) => r[col.locus_group])),
    locusTypes: Object.fromEntries([...tallyBy(rows, (r) => r[col.locus_type])].sort((a, b) => b[1] - a[1])),
    perYear: {
      approved: byYear(rows, col.date_approved_reserved),
      symbolChanged: byYear(rows, col.date_symbol_changed),
      nameChanged: byYear(rows, col.date_name_changed),
    },
    placeholders: Object.fromEntries(Object.entries(PLACEHOLDERS).map(([k, re]) => [k, rows.filter((r) => re.test(r[col.symbol])).length])),
    groups: [...groups].map(([name, n]) => ({ name, n })).sort((a, b) => b.n - a.n).slice(0, 30),
    renamed: ranked((r) => ({ symbol: r[col.symbol], prev: list(r[col.prev_symbol]), n: list(r[col.prev_symbol]).length }), 12).map(({ symbol, prev }) => ({ symbol, prev })),
    aliases: ranked((r) => ({ symbol: r[col.symbol], n: list(r[col.alias_symbol]).length }), 10),
    pseudoParents: [...parents].map(([parent, n]) => ({ parent, n })).sort((a, b) => b.n - a.n).slice(0, 15),
    spreadsheet: { count: spreadsheet.length, examples: spreadsheet.slice(0, 6) },
  };
};

export const geneInfoSummary = (text) => {
  const { col, rows } = table(text);
  return {
    total: rows.length,
    types: Object.fromEntries([...tallyBy(rows, (r) => r[col.type_of_gene])].sort((a, b) => b[1] - a[1])),
    synonyms: rows.map((r) => ({ symbol: r[col.Symbol], n: list(r[col.Synonyms]).length })).sort((a, b) => b.n - a.n).slice(0, 10),
    uncharacterized: rows.filter((r) => /^uncharacterized/i.test(r[col.description])).length,
  };
};

export const maneSummary = (text) => {
  const { col, rows } = table(text);
  const status = tallyBy(rows, (r) => r[col.MANE_status]);
  return { select: status.get('MANE Select') ?? 0, plusClinical: status.get('MANE Plus Clinical') ?? 0, genes: new Set(rows.map((r) => r[col.NCBI_GeneID])).size };
};
