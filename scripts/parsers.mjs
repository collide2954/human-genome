export const CHROMS = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];
const STANDARD = new Set(CHROMS);

export const bare = (name) => String(name).replace(/^chr/, '');
export const isStandard = (name) => STANDARD.has(bare(name));

const rows = (text) => text.split('\n').filter((l) => l && !l.startsWith('#')).map((l) => l.replace(/\r$/, '').split('\t'));

const header = (text, key) => text.match(new RegExp(`^# ${key}:\\s*(.+?)\\s*$`, 'm'))?.[1] ?? null;

export const toModule = (value) => `export default ${JSON.stringify(value, null, 2)};\n`;

export const parseAssemblyReport = (text) => ({
  name: header(text, 'Assembly name'),
  date: header(text, 'Date'),
  genbank: header(text, 'GenBank assembly accession'),
  refseq: header(text, 'RefSeq assembly accession'),
  altLoci: text.split('\n').filter((l) => /^## ALT_REF_LOCI_\d+/.test(l)).length,
  sequences: rows(text).filter((c) => c.length >= 10).map((c) => ({
    name: c[0], role: c[1], molecule: c[2], genbank: c[4], refseq: c[6], unit: c[7], length: Number(c[8]), ucsc: c[9],
  })),
});

const ROLES = ['assembled-molecule', 'unlocalized-scaffold', 'unplaced-scaffold', 'alt-scaffold', 'fix-patch', 'novel-patch'];

export const summarizeRoles = (sequences) => ROLES.map((role) => {
  const hits = sequences.filter((s) => s.role === role);
  return { role, count: hits.length, length: hits.reduce((n, s) => n + s.length, 0) };
});

const camel = (s) => s.replace(/-([a-z0-9])/gi, (_, c) => c.toUpperCase());

export const parseAssemblyStats = (text) => {
  const total = {};
  const chromosomes = new Map();
  for (const [unit, molecule, , kind, stat, value] of rows(text)) {
    if (unit === 'all' && molecule === 'all' && kind === 'all') total[camel(stat)] = Number(value);
    else if ((unit === 'Primary Assembly' || unit === 'non-nuclear') && kind === 'assembled-molecule') {
      if (!chromosomes.has(molecule)) chromosomes.set(molecule, { chrom: molecule });
      chromosomes.get(molecule)[camel(stat)] = Number(value);
    }
  }
  return { total, chromosomes: [...chromosomes.values()] };
};

export const parseCytobands = (text) => rows(text)
  .filter((c) => c.length >= 5 && isStandard(c[0]) && c[3])
  .map((c) => ({ chrom: bare(c[0]), start: Number(c[1]), end: Number(c[2]), name: c[3], stain: c[4] }));

export const parseGaps = (text) => {
  const out = {};
  for (const c of rows(text)) {
    const chrom = bare(c[1]);
    if (!STANDARD.has(chrom)) continue;
    out[chrom] ??= {};
    out[chrom][c[7]] ??= { count: 0, bp: 0 };
    out[chrom][c[7]].count += 1;
    out[chrom][c[7]].bp += Number(c[3]) - Number(c[2]);
  }
  return out;
};

export const parseCentromeres = (text) => {
  const out = {};
  for (const c of rows(text)) {
    const chrom = bare(c[1]);
    if (!STANDARD.has(chrom)) continue;
    out[chrom] = (out[chrom] ?? 0) + Number(c[3]) - Number(c[2]);
  }
  return out;
};
