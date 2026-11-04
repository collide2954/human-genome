export const CHROMS = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];
const STANDARD = new Set(CHROMS);

export const bare = (name) => String(name).replace(/^chr/, '');
export const isStandard = (name) => STANDARD.has(bare(name));

const rows = (text) => text.split('\n').filter((l) => l && !l.startsWith('#')).map((l) => l.replace(/\r$/, '').split('\t'));

const header = (text, key) => text.match(new RegExp(`^# ${key}:\\s*(.+?)\\s*$`, 'm'))?.[1] ?? null;

export const toModule = (value) => `export default ${JSON.stringify(value, null, 2)};\n`;

const newer = (a, b) => {
  const x = a.split('.').map(Number);
  const y = b.split('.').map(Number);
  for (let i = 0; i < Math.max(x.length, y.length); i += 1) {
    if ((x[i] ?? 0) !== (y[i] ?? 0)) return (x[i] ?? 0) > (y[i] ?? 0);
  }
  return false;
};

export const latestFile = (listing, re) => {
  const hits = [...listing.matchAll(new RegExp(re.source, 'g'))].map((m) => ({ file: m[0], version: m[1] }));
  return hits.length ? hits.reduce((best, h) => (newer(h.version, best.version) ? h : best)) : null;
};

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

const joinChromosomes = (report, stats, nuclearOnly) => {
  const byChrom = new Map(stats.chromosomes.map((c) => [c.chrom, c]));
  return report.sequences
    .filter((s) => s.role === 'assembled-molecule' && !(nuclearOnly && s.unit === 'non-nuclear'))
    .map((s) => {
      const st = byChrom.get(s.molecule) ?? {};
      return {
        chrom: s.molecule,
        length: s.length,
        ungapped: st.ungappedLength ?? null,
        gc: st.gcPerc ?? null,
        spanned: st.spannedGaps ?? null,
        unspanned: st.unspannedGaps ?? null,
        refseq: s.refseq === 'na' ? null : s.refseq,
        genbank: s.genbank,
      };
    });
};

export const buildAssembly = ({ report38, stats38, reportT2t, statsT2t, cyto, gaps, cen }) => {
  const r38 = parseAssemblyReport(report38);
  const s38 = parseAssemblyStats(stats38);
  const rT2t = parseAssemblyReport(reportT2t);
  const sT2t = parseAssemblyStats(statsT2t);
  return {
    grch38: {
      name: r38.name, date: r38.date, genbank: r38.genbank, refseq: r38.refseq, altLoci: r38.altLoci,
      roles: summarizeRoles(r38.sequences), total: s38.total, chromosomes: joinChromosomes(r38, s38, false),
    },
    t2t: {
      name: rT2t.name, date: rT2t.date, genbank: rT2t.genbank, refseq: rT2t.refseq, total: sT2t.total, chromosomes: joinChromosomes(rT2t, sT2t, true),
    },
    bands: parseCytobands(cyto),
    gaps: parseGaps(gaps),
    centromeres: parseCentromeres(cen),
  };
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
