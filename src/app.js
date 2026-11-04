import snapshot from './snapshot.js';
import {
  SOURCES, SOURCE_TITLES, SOURCE_GROUPS, EXAMPLES,
} from './facts.js';
import {
  exact, bases, coords, signed, identify, matches, LENSES, readLens, inLens, lensSearch, toCsv, pct, longDate, information, GENETIC_CODE, AMINO_ACIDS, codonOf, perLitre,
} from './lib.js';
import { pearson, linearFit } from './stats.js';
import { FILES } from './queries.js';
import {
  findings, assemblyNumbers, sequenceNumbers, codingCounts, geneDensity, senseCodons, gtexTissues, plasmaLevel, proteomeResidues, proteinLabel,
} from './findings.js';
import {
  ideogram, hbars, dotplot, strips, matrix, rampLegend, treemap, stackedArea, waffle, columns, histo, scatter, logo, codonTable, heatmap, observeFigures, countUp, tipLines,
} from './charts.js';

const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];

const h = (tag, props = {}, ...kids) => {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') node.className = v;
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const kid of kids.flat()) {
    if (kid !== null && kid !== undefined && kid !== false) node.append(kid.nodeType ? kid : document.createTextNode(String(kid)));
  }
  return node;
};

const ext = (href, text, props = {}) => h('a', { href, target: '_blank', rel: 'noopener', ...props }, text);
const fill = (sel, ...nodes) => $(sel).replaceChildren(...nodes);
const CAT = ['var(--s1)', 'var(--s2)', 'var(--s3)', 'var(--s4)', 'var(--s5)', 'var(--data)', 'var(--data-2)'];
const SITE = 'https://collide2954.github.io/human-genome/';
const REPO = 'https://github.com/collide2954/human-genome';
const ucsc = (chrom, start, end) => `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${encodeURIComponent(`chr${chrom}:${start}-${end}`)}`;

const F = findings(snapshot);
const A = snapshot.assembly;
const N = assemblyNumbers(snapshot);
const NUCLEAR = A.grch38.chromosomes.filter((c) => c.chrom !== 'MT');

const state = { query: '', lens: 'all' };

const attachData = (mountSel, name, heads, rows) => {
  const caption = $(mountSel).closest('.fig').querySelector('figcaption');
  caption.querySelectorAll('.figdata').forEach((n) => n.remove());
  const table = h('table', { class: 'data' },
    h('thead', {}, h('tr', {}, heads.map((c) => h('th', { scope: 'col' }, c)))),
    h('tbody', {}, rows.map((r) => h('tr', { 'data-search': r.join(' ') }, r.map((v) => h('td', { class: typeof v === 'number' ? 'n' : null }, v === null ? '—' : typeof v === 'number' ? exact(v) : v))))));
  caption.append(
    h('span', { class: 'figdata' }, ' ', h('details', { class: 'data' }, h('summary', {}, 'Table'), h('div', { class: 'scroll' }, table))),
    h('a', { class: 'csv figdata', download: `${name}.csv`, href: `data:text/csv;charset=utf-8,${encodeURIComponent(toCsv(heads, rows))}` }, 'CSV'),
  );
};

const day = (iso) => new Date(iso).toLocaleDateString('en-GB', { day: 'numeric', month: 'long', year: 'numeric', timeZone: 'UTC' });

const renderHead = () => {
  $('#stamp').textContent = `Data as of ${day(snapshot.generated)}. Rebuilt every Sunday from NCBI, UCSC, GENCODE and EBI files.`;
  for (const node of $$('[data-finding]')) node.textContent = F[node.dataset.finding];
};

const renderLede = () => {
  countUp($('[data-hero="bases"]'), N.t2tLength);
  const b = (text) => h('b', {}, text);
  fill('#lede-text',
    'That is T2T-CHM13v2.0, the first gapless assembly. GRCh38.p14, the reference the databases annotate, spans ', b(`${exact(N.grch38Length)} bases`),
    ' across chromosomes 1 to 22, X and Y, and ', b(`${exact(N.grch38N)} of them are N`), ': sequence it could not read.');
  const items = [
    ['GRCh38 contig N50', bases(A.grch38.total.contigN50)],
    ['GRCh38 scaffold N50', bases(A.grch38.total.scaffoldN50)],
    ['Alternate loci', exact(A.grch38.roles.find((r) => r.role === 'alt-scaffold').count)],
    ['Patches in GRCh38.p14', exact(N.patches)],
    ['Giemsa bands', exact(A.bands.length)],
    ['Gaps in T2T-CHM13', exact(A.t2t.total.totalGapLength ?? 0)],
  ];
  fill('#ticker', ...items.map(([label, value]) => h('div', {}, h('dt', {}, label), h('dd', {}, value))));
};

const STAINS = {
  gneg: 'Giemsa-negative, light band', gpos25: 'Giemsa-positive, 25% stain', gpos50: 'Giemsa-positive, 50% stain', gpos75: 'Giemsa-positive, 75% stain',
  gpos100: 'Giemsa-positive, full stain', acen: 'Centromere', gvar: 'Variable heterochromatin', stalk: 'Acrocentric stalk',
};

const renderChromosomes = () => {
  const acro = [...new Set(A.bands.filter((b) => b.stain === 'stalk').map((b) => b.chrom))];
  $('#chromosomes-dek').textContent = `${exact(A.bands.length)} bands at the 850-band level, from 1p36.33 to Yq12. Chromosomes ${acro.slice(0, -1).join(', ')} and ${acro.at(-1)} are acrocentric: their short arms are stalks and satellites holding the ribosomal RNA genes.`;
  ideogram({
    mount: $('#ideogram-chart'),
    chromosomes: NUCLEAR,
    bands: A.bands,
    aria: 'Ideogram of the 24 GRCh38 chromosomes with their Giemsa bands',
    describe: (band) => tipLines(`${band.chrom}${band.name}`, coords(band.chrom, band.start + 1, band.end),
      `${STAINS[band.stain] ?? band.stain}, ${bases(band.end - band.start)}`, { href: ucsc(band.chrom, band.start + 1, band.end), text: 'Open in the UCSC browser' }),
  });
  attachData('#ideogram-chart', 'cytobands', ['Chromosome', 'Band', 'Start', 'End', 'Stain'], A.bands.map((b) => [b.chrom, b.name, b.start + 1, b.end, b.stain]));
};

const GAP_LABELS = {
  heterochromatin: 'heterochromatin', short_arm: 'short arm', contig: 'between contigs', scaffold: 'between scaffolds', telomere: 'telomere', centromere: 'centromere',
};

const renderGaps = () => {
  const types = [...new Set(Object.values(A.gaps).flatMap((t) => Object.keys(t)))]
    .map((type) => ({ type, bp: Object.values(A.gaps).reduce((n, t) => n + (t[type]?.bp ?? 0), 0) }))
    .sort((a, b) => b.bp - a.bp);
  const segments = types.map((t, i) => ({ key: t.type, label: GAP_LABELS[t.type] ?? t.type, color: CAT[i] }));
  const rows = NUCLEAR.map((c) => ({
    key: c.chrom, label: c.chrom, values: Object.fromEntries(types.map((t) => [t.type, A.gaps[c.chrom]?.[t.type]?.bp ?? 0])),
  }));
  const top = rows.reduce((a, r) => (Object.values(r.values).reduce((x, y) => x + y, 0) > Object.values(a.values).reduce((x, y) => x + y, 0) ? r : a));
  $('#gaps-dek').textContent = `The ${exact(N.grch38N)} N bases fall in ${exact(Object.values(A.gaps).reduce((n, t) => n + Object.values(t).reduce((m, g) => m + g.count, 0), 0))} gaps. ${GAP_LABELS[types[0].type][0].toUpperCase()}${GAP_LABELS[types[0].type].slice(1)} accounts for ${exact(types[0].bp)} bases; chromosome ${top.key} holds the most Ns.`;
  hbars({
    mount: $('#gaps-chart'),
    rows,
    segments,
    aria: 'N bases per GRCh38 chromosome by gap type',
    describe: (r) => {
      const total = Object.values(r.values).reduce((a, b) => a + b, 0);
      return tipLines(`${exact(total)} N`, `chromosome ${r.key}`,
        ...types.filter((t) => r.values[t.type]).map((t) => `${GAP_LABELS[t.type] ?? t.type}: ${exact(r.values[t.type])} in ${A.gaps[r.key][t.type].count} gaps`));
    },
  });
  attachData('#gaps-chart', 'gaps', ['Chromosome', 'Gap type', 'Gaps', 'Bases'],
    NUCLEAR.flatMap((c) => Object.entries(A.gaps[c.chrom] ?? {}).map(([type, g]) => [c.chrom, type, g.count, g.bp])));
  const ROLE_LABELS = {
    'assembled-molecule': 'Chromosomes', 'unlocalized-scaffold': 'Unlocalized scaffolds', 'unplaced-scaffold': 'Unplaced scaffolds',
    'alt-scaffold': 'Alternate loci', 'fix-patch': 'Fix patches', 'novel-patch': 'Novel patches',
  };
  dotplot({
    mount: $('#roles-chart'),
    rows: A.grch38.roles.map((r) => ({ key: r.role, label: ROLE_LABELS[r.role], value: r.length, hi: r.role === 'alt-scaffold' })),
    format: bases,
    aria: 'GRCh38.p14 sequences by role, total length',
    describe: (r) => {
      const role = A.grch38.roles.find((x) => x.role === r.key);
      return tipLines(bases(role.length), r.label, `${exact(role.count)} sequences`, `${exact(role.length)} bp`);
    },
  });
  attachData('#roles-chart', 'sequence-roles', ['Role', 'Sequences', 'Bases'], A.grch38.roles.map((r) => [ROLE_LABELS[r.role], r.count, r.length]));
};

const renderT2t = () => {
  const rows = N.perChrom.map((p) => ({ key: p.chrom, label: p.chrom, values: { d: p.bases } }));
  const acro = ['13', '14', '15', '21', '22'].reduce((n, c) => n + N.perChrom.find((p) => p.chrom === c).bases, 0);
  $('#t2t-dek').textContent = `Chromosomes 13, 14, 15, 21 and 22, whose short arms GRCh38 left blank, gain ${exact(acro)} bases between them. Chromosome X is the one that shrinks, by ${exact(-Math.min(...N.perChrom.map((p) => p.bases)))} bases.`;
  hbars({
    mount: $('#t2t-chart'),
    rows,
    segments: [{ key: 'd', label: 'Bases gained', color: 'var(--accent)' }],
    showSign: true,
    aria: 'Bases each T2T-CHM13 chromosome adds over the sequenced bases of GRCh38',
    describe: (r) => {
      const t = A.t2t.chromosomes.find((c) => c.chrom === r.key);
      const g = A.grch38.chromosomes.find((c) => c.chrom === r.key);
      return tipLines(`${signed(r.values.d)} bases`, `chromosome ${r.key}`, `T2T-CHM13 ${exact(t.length)} bp, no gaps`, `GRCh38 ${exact(g.length)} bp, ${exact(g.ungapped)} sequenced`);
    },
  });
  attachData('#t2t-chart', 't2t-vs-grch38', ['Chromosome', 'T2T-CHM13 bases', 'GRCh38 bases', 'GRCh38 sequenced bases', 'Difference'],
    N.perChrom.map((p) => {
      const g = A.grch38.chromosomes.find((c) => c.chrom === p.chrom);
      return [p.chrom, p.bases + g.ungapped, g.length, g.ungapped, p.bases];
    }));
};

const S = snapshot.sequence;
const SN = sequenceNumbers(snapshot);
const ACGT = 'ACGT';
const acgtOf = (c) => c.A + c.C + c.G + c.T;
const gcPct = (v) => `${(v * 100).toFixed(1)}%`;

const renderComposition = () => {
  const nuclear = S.chromosomes.filter((c) => c.chrom !== 'MT');
  const cells = nuclear.flatMap((c) => S.bins[c.chrom].map((v, j) => ({ chrom: c.chrom, j, v }))).filter((b) => b.v !== null);
  const rich = cells.reduce((a, b) => (b.v > a.v ? b : a));
  const poor = cells.reduce((a, b) => (b.v < a.v ? b : a));
  const masked = nuclear.reduce((n, c) => n + c.masked, 0);
  const binOf = (b) => coords(b.chrom, b.j * 1e6 + 1, Math.min((b.j + 1) * 1e6, S.chromosomes.find((c) => c.chrom === b.chrom).length));
  $('#composition-dek').textContent = `Of ${exact(SN.acgt)} sequenced bases, ${pct(masked / SN.acgt)} are soft-masked as repeats. The richest megabase, ${binOf(rich)}, is ${gcPct(rich.v)} GC; the poorest, ${binOf(poor)}, is ${gcPct(poor.v)}.`;
  const lo = Math.min(...cells.map((b) => b.v));
  const hi = Math.max(...cells.map((b) => b.v));
  strips({
    mount: $('#gc-chart'),
    rows: nuclear.map((c) => ({ key: c.chrom, bins: S.bins[c.chrom], right: gcPct((c.C + c.G) / acgtOf(c)) })),
    lo,
    hi,
    legend: rampLegend(lo, hi, gcPct),
    aria: 'GC content of each megabase of each chromosome',
    describe: (b) => tipLines(`${gcPct(b.value)} GC`, binOf({ chrom: b.row.key, j: b.index }), `chromosome ${b.row.key}: ${b.row.right}`),
  });
  attachData('#gc-chart', 'gc-per-megabase', ['Chromosome', 'Start', 'End', 'GC'],
    nuclear.flatMap((c) => S.bins[c.chrom].map((v, j) => [c.chrom, j * 1e6 + 1, Math.min((j + 1) * 1e6, c.length), v])));
  const freq = ACGT.split('').map((b) => nuclear.reduce((n, c) => n + c[b], 0) / SN.acgt);
  const pairs = S.dinucs.reduce((a, b) => a + b, 0);
  const oe = (i, j) => S.dinucs[i * 4 + j] / (freq[i] * freq[j] * pairs);
  matrix({
    mount: $('#dinuc-chart'),
    rows: [...ACGT],
    cols: [...ACGT],
    values: [0, 1, 2, 3].map((i) => [0, 1, 2, 3].map((j) => oe(i, j))),
    format: (v) => v.toFixed(2),
    aria: 'Dinucleotide observed over expected ratios',
    describe: (i, j) => tipLines(`${oe(i, j).toFixed(2)} observed/expected`, `${ACGT[i]}p${ACGT[j]}: ${ACGT[i]} followed by ${ACGT[j]}`,
      `${exact(S.dinucs[i * 4 + j])} in GRCh38`, `${exact(freq[i] * freq[j] * pairs)} expected from base frequencies`),
  });
  attachData('#dinuc-chart', 'dinucleotides', ['Dinucleotide', 'Observed', 'Expected', 'Ratio'],
    [0, 1, 2, 3].flatMap((i) => [0, 1, 2, 3].map((j) => [`${ACGT[i]}${ACGT[j]}`, S.dinucs[i * 4 + j], Math.round(freq[i] * freq[j] * pairs), Number(oe(i, j).toFixed(4))])));
  const cpgOe = (c) => c.cpg / ((c.C * c.G) / acgtOf(c));
  hbars({
    mount: $('#cpg-chart'),
    rows: S.chromosomes.map((c) => ({ key: c.chrom, label: c.chrom, values: { oe: cpgOe(c) } })),
    segments: [{ key: 'oe', label: 'CpG observed/expected', color: 'var(--s2)' }],
    format: (v) => v.toFixed(3),
    tickFormat: (v) => String(v),
    aria: 'CpG observed over expected by chromosome',
    describe: (r) => {
      const c = S.chromosomes.find((x) => x.chrom === r.key);
      return tipLines(`${cpgOe(c).toFixed(3)} observed/expected`, r.key === 'MT' ? 'mitochondrial genome' : `chromosome ${r.key}`, `${exact(c.cpg)} CpG`, `${gcPct((c.C + c.G) / acgtOf(c))} GC`);
    },
  });
  attachData('#cpg-chart', 'cpg-by-chromosome', ['Chromosome', 'CpG', 'C', 'G', 'Sequenced bases', 'Observed over expected'],
    S.chromosomes.map((c) => [c.chrom, c.cpg, c.C, c.G, acgtOf(c), Number(cpgOe(c).toFixed(4))]));
  hbars({
    mount: $('#g4-chart'),
    rows: nuclear.map((c) => ({ key: c.chrom, label: c.chrom, values: { plus: c.g4.plus, minus: c.g4.minus } })),
    segments: [{ key: 'plus', label: 'G-rich strand', color: 'var(--s1)' }, { key: 'minus', label: 'C-rich strand', color: 'var(--s2)' }],
    aria: 'G-quadruplex motifs by chromosome and strand',
    describe: (r) => {
      const c = nuclear.find((x) => x.chrom === r.key);
      return tipLines(`${exact(c.g4.plus + c.g4.minus)} motifs`, `chromosome ${r.key}`, `${exact(c.g4.plus)} G-rich, ${exact(c.g4.minus)} C-rich`,
        `${((c.g4.plus + c.g4.minus) / (acgtOf(c) / 1e6)).toFixed(1)} per sequenced megabase`);
    },
  });
  attachData('#g4-chart', 'g-quadruplex-motifs', ['Chromosome', 'G-rich strand', 'C-rich strand', 'Sequenced bases'], nuclear.map((c) => [c.chrom, c.g4.plus, c.g4.minus, acgtOf(c)]));
};

const R = snapshot.repeats;
const CLASS_COLORS = {
  LINE: 'var(--s1)', SINE: 'var(--s2)', LTR: 'var(--s3)', DNA: 'var(--s4)', Satellite: 'var(--s5)',
};
const CLASS_NAMES = {
  LINE: 'LINE', SINE: 'SINE', LTR: 'LTR', DNA: 'DNA transposons', Satellite: 'satellites', Simple_repeat: 'simple repeats', Low_complexity: 'low complexity', Retroposon: 'SVA',
};

const renderTransposons = () => {
  const { rmsk } = R;
  const top = rmsk.families.slice(0, 18);
  const rest = rmsk.families.slice(18);
  const items = [
    ...top.map((f) => ({
      key: `${f.cls}/${f.family}`, label: f.family.replace(/_/g, ' '), sub: pct(f.bp / SN.acgt), value: f.bp, color: CLASS_COLORS[f.cls] ?? 'var(--data)', f,
    })),
    {
      key: 'other', label: 'Other', sub: pct(rest.reduce((n, f) => n + f.bp, 0) / SN.acgt), value: rest.reduce((n, f) => n + f.bp, 0), color: 'var(--data-2)',
      f: { cls: `${rest.length} families`, family: 'Other', count: rest.reduce((n, f) => n + f.count, 0), bp: rest.reduce((n, f) => n + f.bp, 0) },
    },
  ];
  $('#transposons-dek').textContent = `${rmsk.classes.slice(0, 5).map((c) => `${CLASS_NAMES[c.cls] ?? c.cls} ${pct(c.bp / SN.acgt)}`).join(', ')}. ${exact(rmsk.total.count)} annotations in ${exact(rmsk.families.length)} families.`;
  treemap({
    mount: $('#repeat-chart'),
    items,
    aria: 'Bases in each repeat family',
    describe: (item) => tipLines(bases(item.f.bp), `${item.f.family} (${item.f.cls})`, `${exact(item.f.count)} copies`, `${pct(item.f.bp / SN.acgt)} of sequenced bases`),
  });
  attachData('#repeat-chart', 'repeat-families', ['Class', 'Family', 'Copies', 'Bases'], rmsk.families.map((f) => [f.cls, f.family, f.count, f.bp]));
};

const renderLandscape = () => {
  const { landscape, alu: subs } = R.rmsk;
  const years = landscape.edges.slice(0, -1);
  const peakOf = (cls) => landscape.bp[cls].indexOf(Math.max(...landscape.bp[cls]));
  $('#landscape-dek').textContent = `SINE bases peak at ${peakOf('SINE')}–${peakOf('SINE') + 1}% divergence, the Alu expansion; LINE bases at ${peakOf('LINE')}–${peakOf('LINE') + 1}%, DNA transposons at ${peakOf('DNA')}–${peakOf('DNA') + 1}%. Copies under 2% divergence are the recent, sometimes still active, L1 and Alu.`;
  stackedArea({
    mount: $('#landscape-chart'),
    years,
    series: landscape.classes.map((cls) => ({ key: cls, label: CLASS_NAMES[cls], color: CLASS_COLORS[cls], values: landscape.bp[cls] })),
    format: bases,
    xFormat: (v) => `${v}–${v + 1}% divergence`,
    xTick: (v) => `${v}%`,
    aria: 'Transposon bases by divergence from consensus',
  });
  attachData('#landscape-chart', 'repeat-landscape', ['Divergence from', 'Divergence to', ...landscape.classes],
    years.map((v, i) => [v, v + 1, ...landscape.classes.map((cls) => landscape.bp[cls][i])]));
  const parts = [['AluJ', 'var(--s3)'], ['AluS', 'var(--s2)'], ['AluY', 'var(--s1)']].map(([key, color]) => ({ label: key, value: subs[key].count, color }));
  waffle({ mount: $('#alu-chart'), parts, aria: 'Alu copies by subfamily' });
  attachData('#alu-chart', 'alu-subfamilies', ['Subfamily', 'Copies', 'Bases'], ['AluJ', 'AluS', 'AluY'].map((k) => [k, subs[k].count, subs[k].bp]));
};

const renderSegdups = () => {
  const { superDups } = R;
  const { edges, counts } = superDups.identity;
  const rows = counts.map((count, i) => ({ lo: edges[i], hi: edges[i + 1], count }));
  const high = rows.filter((r) => r.lo >= 0.99).reduce((n, r) => n + r.count, 0);
  const covered = Object.values(superDups.covered).reduce((a, b) => a + b, 0);
  const share = (c) => (superDups.covered[c.chrom] ?? 0) / acgtOf(S.chromosomes.find((x) => x.chrom === c.chrom));
  const most = NUCLEAR.reduce((a, b) => (share(b) > share(a) ? b : a));
  $('#segdups-dek').textContent = `${exact(high)} of the ${exact(superDups.pairs)} pairs are at least 99% identical. ${exact(superDups.kinds.intra)} pairs sit on one chromosome. Chromosome ${most.chrom} has the largest share in duplications, ${pct(share(most))} of its sequenced bases; ${bases(covered)} are duplicated across the genome.`;
  const pctLabel = (v) => `${(v * 100).toFixed(1)}`;
  columns({
    mount: $('#segdup-chart'),
    rows,
    key: (r) => r.lo,
    value: (r) => r.count,
    attr: 'bin',
    label: (r) => pctLabel(r.lo),
    height: 220,
    aria: 'Segmental duplication pairs by identity',
    describe: (r) => tipLines(`${exact(r.count)} pairs`, `${pctLabel(r.lo)}–${pctLabel(r.hi)}% identity`),
  });
  attachData('#segdup-chart', 'segdup-identity', ['Identity from', 'Identity to', 'Pairs'], rows.map((r) => [r.lo, r.hi, r.count]));
  hbars({
    mount: $('#segdup-chrom-chart'),
    rows: NUCLEAR.map((c) => ({ key: c.chrom, label: c.chrom, values: { bp: superDups.covered[c.chrom] ?? 0 } })),
    segments: [{ key: 'bp', label: 'Duplicated bases', color: 'var(--s2)' }],
    aria: 'Bases in segmental duplications by chromosome',
    describe: (r) => tipLines(bases(r.values.bp), `chromosome ${r.key}`, `${pct(share({ chrom: r.key }))} of its sequenced bases`),
  });
  attachData('#segdup-chrom-chart', 'segdup-by-chromosome', ['Chromosome', 'Duplicated bases', 'Sequenced bases'],
    NUCLEAR.map((c) => [c.chrom, superDups.covered[c.chrom] ?? 0, acgtOf(S.chromosomes.find((x) => x.chrom === c.chrom))]));
};

const renderTandem = () => {
  const { units, motifs: top } = R.tandem;
  const unitLabel = (u) => `${u} bp`;
  const micro = units.slice(0, 6).reduce((n, u) => n + u.count, 0);
  $('#tandem-dek').textContent = `${exact(micro)} calls have units of 1 to 6 bases, the microsatellites; ${exact(units.at(-1).count)} have units over 50 bases. (AC)n alone is called ${exact(top.find((m) => m.motif === 'AC').count)} times.`;
  hbars({
    mount: $('#tandem-chart'),
    rows: units.map((u) => ({ key: u.unit, label: unitLabel(u.unit), values: { count: u.count } })),
    segments: [{ key: 'count', label: 'Calls', color: 'var(--data)' }],
    labelShare: 0.22,
    aria: 'Tandem repeats by unit length',
    describe: (r) => {
      const u = units.find((x) => x.unit === r.key);
      return tipLines(`${exact(u.count)} repeats`, `units of ${unitLabel(u.unit)}`, `${bases(u.bp)} called`);
    },
  });
  attachData('#tandem-chart', 'tandem-units', ['Unit length', 'Repeats', 'Bases'], units.map((u) => [u.unit, u.count, u.bp]));
  const shown = top.slice(0, 15);
  dotplot({
    mount: $('#motif-chart'),
    rows: shown.map((m, i) => ({ key: m.motif, label: `(${m.motif})n`, value: m.count, hi: i === 0 })),
    aria: 'The most frequent tandem repeat units',
    describe: (r) => {
      const m = shown.find((x) => x.motif === r.key);
      return tipLines(`${exact(m.count)} repeats`, `(${m.motif})n with its rotations and reverse complement`, `${bases(m.bp)} called`);
    },
  });
  attachData('#motif-chart', 'tandem-motifs', ['Unit', 'Repeats', 'Bases'], top.map((m) => [m.motif, m.count, m.bp]));
};

const G = snapshot.genes;
const CAT2 = snapshot.catalogs;
const typeLabel = (t) => (t === 'protein_coding' ? 'protein-coding' : t.replace(/_/g, ' '));
const lengthOfChrom = (chrom) => S.chromosomes.find((c) => c.chrom === chrom)?.length ?? A.grch38.chromosomes.find((c) => c.chrom === chrom).length;
const megabase = (chrom, j) => coords(chrom, j * 1e6 + 1, Math.min((j + 1) * 1e6, lengthOfChrom(chrom)));
const shares = (counts) => {
  const total = counts.reduce((a, b) => a + b, 0) || 1;
  return counts.map((c) => (c / total) * 100);
};
const binsOf = (edges) => edges.slice(0, -1).map((lo, i) => ({ lo, hi: edges[i + 1] }));
const shortBases = (v) => (v >= 1e6 ? `${Number((v / 1e6).toPrecision(3))}M` : v >= 1e3 ? `${Number((v / 1e3).toPrecision(3))}k` : String(Math.round(v)));

const renderGenecount = () => {
  const counts = codingCounts(snapshot);
  const notes = {
    GENCODE: `release ${G.release.gencode}, Ensembl ${G.release.ensembl}, ${longDate(G.release.date)}`,
    'NCBI Gene': 'records of type protein-coding in gene_info',
    HGNC: 'approved symbols in the protein-coding locus group',
    'MANE Select': `release ${CAT2.mane.release}, one transcript per gene`,
  };
  $('#genecount-dek').textContent = `GENCODE ${G.release.gencode} (Ensembl ${G.release.ensembl}, ${longDate(G.release.date)}) lists ${exact(G.genes.total)} genes and ${exact(G.transcripts.total)} transcripts on the chromosomes and the mitochondrial genome. NCBI Gene also keeps ${exact(CAT2.geneInfo.types['biological-region'] ?? 0)} records for regions that are not genes, and ${exact(CAT2.geneInfo.uncharacterized)} of its records are described as uncharacterized.`;
  dotplot({
    mount: $('#source-count-chart'),
    rows: counts.map((c) => ({ key: c.source, label: c.source, value: c.n, hi: c.source === 'GENCODE' })),
    scale: 'linear',
    labelShare: 0.3,
    aria: 'Protein-coding gene counts by catalog',
    describe: (r) => tipLines(exact(r.value), `protein-coding genes in ${r.label}`, notes[r.key]),
  });
  attachData('#source-count-chart', 'protein-coding-counts', ['Source', 'Protein-coding genes'], counts.map((c) => [c.source, c.n]));
  const types = Object.entries(G.genes.byType).slice(0, 16);
  dotplot({
    mount: $('#biotype-chart'),
    rows: types.map(([t, n], i) => ({ key: t, label: typeLabel(t), value: n, hi: i === 1 })),
    aria: 'GENCODE genes by biotype',
    describe: (r) => tipLines(exact(r.value), `${r.label} genes`, `${exact(G.transcripts.byType[r.key] ?? 0)} transcripts of this type`),
  });
  attachData('#biotype-chart', 'gene-biotypes', ['Biotype', 'Genes', 'Transcripts of the same type'], Object.entries(G.genes.byType).map(([t, n]) => [t, n, G.transcripts.byType[t] ?? 0]));
};

const renderDensity = () => {
  const { bins } = G.density;
  const dens = geneDensity(snapshot);
  const max = Math.max(...NUCLEAR.flatMap((c) => bins[c.chrom] ?? []));
  strips({
    mount: $('#density-chart'),
    rows: NUCLEAR.map((c) => ({ key: c.chrom, bins: bins[c.chrom] ?? [], right: (dens.find((d) => d.chrom === c.chrom)?.perMb ?? 0).toFixed(1) })),
    lo: 0,
    hi: max,
    legend: rampLegend(0, max, (v) => String(Math.round(v))),
    aria: 'Protein-coding genes per megabase',
    describe: (b) => tipLines(`${b.value} protein-coding genes`, megabase(b.row.key, b.index), `chromosome ${b.row.key}: ${b.row.right} per sequenced megabase`),
  });
  attachData('#density-chart', 'genes-per-megabase', ['Chromosome', 'Start', 'End', 'Protein-coding genes'],
    NUCLEAR.flatMap((c) => (bins[c.chrom] ?? []).map((n, j) => [c.chrom, j * 1e6 + 1, Math.min((j + 1) * 1e6, lengthOfChrom(c.chrom)), n])));
  const points = NUCLEAR.flatMap((c) => (bins[c.chrom] ?? []).map((n, j) => ({
    key: `${c.chrom}:${j}`, chrom: c.chrom, j, gc: S.bins[c.chrom]?.[j], y: n,
  }))).filter((p) => p.gc !== null && p.gc !== undefined).map((p) => ({ ...p, x: Number((p.gc * 100).toFixed(1)) }));
  const fit = linearFit(points.map((p) => p.x), points.map((p) => p.y));
  const r = pearson(points.map((p) => p.x), points.map((p) => p.y));
  const top = points.reduce((a, b) => (b.y > a.y ? b : a));
  top.hi = true;
  $('#density-dek').textContent = `Across ${exact(points.length)} megabases, gene count and GC content correlate with Pearson r = ${r.toFixed(2)}; each extra point of GC goes with ${fit.slope.toFixed(1)} more protein-coding genes per megabase. The densest megabase, ${megabase(top.chrom, top.j)}, holds ${top.y}.`;
  scatter({
    mount: $('#density-gc-chart'),
    points,
    xFormat: (v) => `${v}%`,
    yFormat: (v) => String(v),
    xLabel: 'GC content of the megabase',
    yLabel: 'Protein-coding genes',
    trend: fit,
    aria: 'Genes per megabase against GC content',
    describe: (p) => tipLines(`${p.y} protein-coding genes`, megabase(p.chrom, p.j), `${p.x.toFixed(1)}% GC`),
  });
  attachData('#density-gc-chart', 'genes-and-gc-per-megabase', ['Chromosome', 'Start', 'GC percent', 'Protein-coding genes'], points.map((p) => [p.chrom, p.j * 1e6 + 1, p.x, p.y]));
};

const renderGenesize = () => {
  const { pc, lnc, longest } = G.spans;
  const pcShare = shares(pc.counts);
  const lncShare = shares(lnc.counts);
  const bins = binsOf(pc.edges);
  const overMb = bins.reduce((n, b, i) => n + (b.lo >= 1e6 ? pc.counts[i] : 0), 0);
  $('#genesize-dek').textContent = `Quartiles: ${exact(pc.summary.q1)}, ${exact(pc.summary.median)} and ${exact(pc.summary.q3)} bp for protein-coding genes; ${exact(lnc.summary.q1)}, ${exact(lnc.summary.median)} and ${exact(lnc.summary.q3)} bp for lncRNA genes. ${exact(overMb)} protein-coding genes span a megabase or more.`;
  histo({
    mount: $('#span-chart'),
    bins,
    series: [{ key: 'pc', label: 'Protein-coding', color: 'var(--s1)', values: pcShare }, { key: 'lnc', label: 'lncRNA', color: 'var(--s2)', values: lncShare }],
    yFormat: (v) => `${v}%`,
    tickLabel: (b) => shortBases(b.lo),
    aria: 'Gene spans of protein-coding and lncRNA genes',
    describe: (i) => tipLines(`${bases(bins[i].lo)} to ${bases(bins[i].hi)}`, 'genomic span',
      `Protein-coding: ${exact(pc.counts[i])} genes, ${pcShare[i].toFixed(1)}%`, `lncRNA: ${exact(lnc.counts[i])} genes, ${lncShare[i].toFixed(1)}%`),
  });
  attachData('#span-chart', 'gene-spans', ['Span from', 'Span to', 'Protein-coding genes', 'lncRNA genes'], bins.map((b, i) => [Math.round(b.lo), Math.round(b.hi), pc.counts[i], lnc.counts[i]]));
  dotplot({
    mount: $('#longest-chart'),
    rows: longest.map((g, i) => ({ key: g.name, label: g.name, value: g.span, hi: i === 0 })),
    scale: 'linear',
    format: bases,
    aria: 'The longest protein-coding genes',
    describe: (r) => {
      const g = longest.find((x) => x.name === r.key);
      return tipLines(bases(g.span), g.name, coords(g.chrom, g.start, g.end));
    },
  });
  attachData('#longest-chart', 'longest-genes', ['Gene', 'Chromosome', 'Start', 'End', 'Span'], longest.map((g) => [g.name, g.chrom, g.start, g.end, g.span]));
};

const renderStructure = () => {
  const m = G.mane;
  const ex = m.exons;
  $('#structure-dek').textContent = `Introns have a median length of ${exact(m.introns.summary.median)} bp and exons ${exact(m.exonLengths.summary.median)} bp. The longest MANE intron, in ${m.introns.longest.name}, runs ${exact(m.introns.longest.length)} bp. Median untranslated regions: ${exact(m.utr5.summary.median)} bp before the coding sequence and ${exact(m.utr3.summary.median)} bp after it.`;
  const exonBins = ex.counts.map((n, i) => ({ n: i + 1 }));
  histo({
    mount: $('#exon-count-chart'),
    bins: exonBins,
    series: [{ key: 'tx', label: 'MANE Select transcripts', color: 'var(--data)', values: ex.counts }],
    tickEvery: 1,
    centered: true,
    tickLabel: (b, i) => (i === 0 || (i + 1) % 10 === 0 ? (i === 99 ? '100+' : String(i + 1)) : null),
    aria: 'MANE Select transcripts by number of exons',
    describe: (i) => tipLines(`${exact(ex.counts[i])} transcripts`, i === 99 ? '100 or more exons' : `${i + 1} exon${i ? 's' : ''}`),
  });
  attachData('#exon-count-chart', 'exons-per-transcript', ['Exons', 'MANE Select transcripts'], ex.counts.map((n, i) => [i === 99 ? '100+' : i + 1, n]));
  const grid = binsOf(Array.from({ length: 71 }, (_, i) => Number((10 ** (i / 10)).toPrecision(12))));
  const exonCounts = grid.map((_, i) => m.exonLengths.counts[i] ?? 0);
  const intronCounts = grid.map((_, i) => (i >= 10 ? m.introns.counts[i - 10] ?? 0 : 0));
  const exonShare = shares(exonCounts);
  const intronShare = shares(intronCounts);
  histo({
    mount: $('#intron-chart'),
    bins: grid,
    series: [{ key: 'exon', label: 'Exons', color: 'var(--s1)', values: exonShare }, { key: 'intron', label: 'Introns', color: 'var(--s2)', values: intronShare }],
    yFormat: (v) => `${v}%`,
    tickLabel: (b) => shortBases(b.lo),
    aria: 'Lengths of MANE Select exons and introns',
    describe: (i) => tipLines(`${bases(grid[i].lo)} to ${bases(grid[i].hi)}`, 'length', `Exons: ${exact(exonCounts[i])}, ${exonShare[i].toFixed(1)}%`, `Introns: ${exact(intronCounts[i])}, ${intronShare[i].toFixed(1)}%`),
  });
  attachData('#intron-chart', 'exon-and-intron-lengths', ['Length from', 'Length to', 'Exons', 'Introns'], grid.map((b, i) => [Math.round(b.lo), Math.round(b.hi), exonCounts[i], intronCounts[i]]));
  const t = G.territory;
  const intergenic = SN.acgt - t.cds - t.utr - t.ncExon - t.intron;
  const parts = [
    { label: 'Coding sequence', value: t.cds, color: 'var(--s1)' },
    { label: 'UTRs and other exons of protein-coding genes', value: t.utr, color: 'var(--s3)' },
    { label: 'Exons of other genes', value: t.ncExon, color: 'var(--s4)' },
    { label: 'Introns', value: t.intron, color: 'var(--data)' },
    { label: 'Between genes', value: intergenic, color: 'var(--data-2)' },
  ];
  waffle({ mount: $('#territory-chart'), parts, aria: 'Sequenced bases by annotation', format: bases });
  attachData('#territory-chart', 'genome-territory', ['Category', 'Bases'], parts.map((p) => [p.label, p.value]));
};

const renderIsoforms = () => {
  const iso = G.isoforms;
  const tx = G.transcripts.byType;
  $('#isoforms-dek').textContent = `${exact(tx.nonsense_mediated_decay ?? 0)} transcripts are annotated for nonsense-mediated decay and ${exact(tx.retained_intron ?? 0)} retain an intron. Most transcripts per gene: ${iso.top.slice(0, 3).map((g) => `${g.name} ${exact(g.n)}`).join(', ')}.`;
  histo({
    mount: $('#isoform-chart'),
    bins: iso.counts.map((n, i) => ({ n: i + 1 })),
    series: [{ key: 'genes', label: 'Protein-coding genes', color: 'var(--data)', values: iso.counts }],
    tickEvery: 1,
    centered: true,
    tickLabel: (b, i) => (i === 0 || (i + 1) % 10 === 0 ? (i === 49 ? '50+' : String(i + 1)) : null),
    aria: 'Protein-coding genes by number of transcripts',
    describe: (i) => tipLines(`${exact(iso.counts[i])} genes`, i === 49 ? '50 or more transcripts' : `${i + 1} transcript${i ? 's' : ''}`),
  });
  attachData('#isoform-chart', 'transcripts-per-gene', ['Transcripts', 'Protein-coding genes'], iso.counts.map((n, i) => [i === 49 ? '50+' : i + 1, n]));
  const types = Object.entries(tx).slice(0, 12);
  dotplot({
    mount: $('#tx-biotype-chart'),
    rows: types.map(([type, n]) => ({ key: type, label: typeLabel(type), value: n, hi: type === 'nonsense_mediated_decay' })),
    aria: 'GENCODE transcripts by biotype',
    describe: (r) => tipLines(exact(r.value), `${r.label} transcripts`),
  });
  attachData('#tx-biotype-chart', 'transcript-biotypes', ['Biotype', 'Transcripts'], Object.entries(tx));
};

const renderNames = () => {
  const H = CAT2.hgnc;
  const first = H.perYear.approved[0];
  $('#names-dek').textContent = `${exact(H.placeholders.orf)} genes still carry open-reading-frame placeholders such as C1orf112, ${exact(H.placeholders.fam)} are FAM family names and ${exact(H.placeholders.kiaa)} keep a KIAA clone number. The spreadsheet renames include ${H.spreadsheet.examples.slice(0, 3).join(', ')}. ${exact(first[1])} of today’s symbols date from ${first[0]}.`;
  const yearly = (sel, name, list, what) => {
    columns({
      mount: $(sel),
      rows: list.map(([year, n]) => ({ year, n })),
      key: (r) => r.year,
      value: (r) => r.n,
      attr: 'bin',
      label: (r) => r.year,
      height: 200,
      aria: `HGNC symbols ${what} each year`,
      describe: (r) => tipLines(`${exact(r.n)} symbols`, `${what} in ${r.year}`),
    });
    attachData(sel, name, ['Year', 'Symbols'], list);
  };
  yearly('#approvals-chart', 'hgnc-approvals', H.perYear.approved, 'approved or reserved');
  yearly('#renames-chart', 'hgnc-symbol-changes', H.perYear.symbolChanged, 'last renamed');
  const groups = H.groups.slice(0, 15);
  dotplot({
    mount: $('#groups-chart'),
    rows: groups.map((g, i) => ({ key: g.name, label: g.name, value: g.n, hi: i === 0 })),
    scale: 'linear',
    labelShare: 0.5,
    aria: 'The largest HGNC gene groups',
    describe: (r) => tipLines(`${exact(r.value)} genes`, r.label),
  });
  attachData('#groups-chart', 'hgnc-gene-groups', ['Group', 'Genes'], H.groups.map((g) => [g.name, g.n]));
};

const renderPseudogenes = () => {
  const types = Object.entries(G.genes.byType).filter(([t]) => /pseudogene$/.test(t));
  const total = types.reduce((n, [, v]) => n + v, 0);
  const processed = G.genes.byType.processed_pseudogene ?? 0;
  $('#pseudogenes-dek').textContent = `GENCODE annotates ${exact(total)} pseudogenes; ${pct(processed / total)} are processed copies of mRNA. HGNC names copies after their parent gene, and ribosomal protein genes dominate the list.`;
  dotplot({
    mount: $('#pseudotype-chart'),
    rows: types.map(([t, n], i) => ({ key: t, label: typeLabel(t).replace(/ pseudogene$/, ''), value: n, hi: i === 0 })),
    aria: 'GENCODE pseudogenes by type',
    describe: (r) => tipLines(exact(r.value), `${typeLabel(r.key)}s`),
  });
  attachData('#pseudotype-chart', 'pseudogene-types', ['Type', 'Pseudogenes'], types);
  const parents = CAT2.hgnc.pseudoParents;
  dotplot({
    mount: $('#parents-chart'),
    rows: parents.map((p, i) => ({ key: p.parent, label: p.parent, value: p.n, hi: i === 0 })),
    scale: 'linear',
    aria: 'Genes with the most named pseudogenes',
    describe: (r) => tipLines(`${exact(r.value)} pseudogenes`, `named ${r.key}P1, ${r.key}P2 and so on`),
  });
  attachData('#parents-chart', 'pseudogene-parents', ['Parent gene', 'Named pseudogenes'], parents.map((p) => [p.parent, p.n]));
};

const LOCI = ['IGH', 'IGK', 'IGL', 'TRA', 'TRB', 'TRG', 'TRD'];

const renderImmune = () => {
  const im = G.immune;
  $('#immune-dek').textContent = `Light chains pair one V and one J: ${im.IGK.V} × ${im.IGK.J} = ${exact(im.IGK.V * im.IGK.J)} kappa and ${im.IGL.V} × ${im.IGL.J} = ${exact(im.IGL.V * im.IGL.J)} lambda combinations. T-cell receptors: ${im.TRA.V} × ${im.TRA.J} = ${exact(im.TRA.V * im.TRA.J)} alpha and ${im.TRB.V} × ${im.TRB.D} × ${im.TRB.J} = ${exact(im.TRB.V * Math.max(1, im.TRB.D) * im.TRB.J)} beta. GENCODE’s counts include open reading frames that IMGT may class apart.`;
  hbars({
    mount: $('#immune-chart'),
    rows: LOCI.map((l) => ({ key: l, label: l, values: im[l] })),
    segments: [['V', 'var(--s1)'], ['D', 'var(--s2)'], ['J', 'var(--s3)'], ['C', 'var(--s4)'], ['pseudo', 'var(--data-2)']].map(([key, color]) => ({ key, label: key === 'pseudo' ? 'Pseudogenes' : `${key} genes`, color })),
    labelShare: 0.1,
    aria: 'Immune receptor gene segments by locus',
    describe: (r) => tipLines(`${exact(Object.values(r.values).reduce((a, b) => a + b, 0))} genes`, `${r.key} locus`, `V ${r.values.V}, D ${r.values.D}, J ${r.values.J}, C ${r.values.C}, pseudogenes ${r.values.pseudo}`),
  });
  attachData('#immune-chart', 'immune-receptor-segments', ['Locus', 'V', 'D', 'J', 'C', 'Pseudogenes'], LOCI.map((l) => [l, im[l].V, im[l].D, im[l].J, im[l].C, im[l].pseudo]));
};

const MN = S.mane;
const signedPos = (v) => (v < 0 ? `−${-v}` : `+${v}`);

const drawLogo = (sel, name, counts, labels, split, aria) => {
  const positions = information(counts);
  logo({
    mount: $(sel),
    positions,
    labels,
    split,
    aria,
    describe: (i) => tipLines(`${positions[i].bits.toFixed(2)} bits`, `position ${labels[i]}`,
      ...[...positions[i].letters].reverse().map((l) => `${l.base} ${(l.p * 100).toFixed(1)}%`), `${exact(positions[i].n)} sequences`),
  });
  attachData(sel, name, ['Position', 'A', 'C', 'G', 'T', 'Bits'], counts.map((row, i) => [labels[i], ...row, Number(positions[i].bits.toFixed(4))]));
};

const renderSplicing = () => {
  const c = MN.classes;
  const donor = information(MN.donor);
  $('#splicing-dek').textContent = `${exact(MN.introns)} introns from ${exact(MN.transcripts)} MANE Select transcripts. The intron’s first G and T are fixed in ${pct(MN.donor[3][2] / donor[3].n)} and ${pct(MN.donor[4][3] / donor[4].n)} of them; the acceptor’s A and G in ${pct(MN.acceptor[18][0] / information(MN.acceptor)[18].n)} and ${pct(MN.acceptor[19][2] / information(MN.acceptor)[19].n)}. ${exact(c.other)} introns use other pairs.`;
  drawLogo('#donor-chart', 'splice-donor', MN.donor, [-3, -2, -1, 1, 2, 3, 4, 5, 6].map(signedPos), 3, '5 prime splice site logo');
  drawLogo('#acceptor-chart', 'splice-acceptor', MN.acceptor, [...Array.from({ length: 20 }, (_, i) => i - 20), 1, 2, 3].map(signedPos), 20, '3 prime splice site logo');
  const classes = Object.entries(c);
  dotplot({
    mount: $('#intron-class-chart'),
    rows: classes.map(([k, n], i) => ({ key: k, label: k === 'other' ? 'Other' : k, value: n, hi: i === 0 })),
    aria: 'Introns by terminal dinucleotides',
    describe: (r) => tipLines(`${exact(r.value)} introns`, r.label, pct(r.value / MN.introns)),
  });
  attachData('#intron-class-chart', 'intron-classes', ['Ends', 'Introns'], classes);
};

const renderTranslation = () => {
  const k = information(MN.kozak);
  const [a, , g] = MN.kozak[3];
  $('#translation-dek').textContent = `At −3, A in ${pct(a / k[3].n)} and G in ${pct(g / k[3].n)}; at +4, G in ${pct(MN.kozak[9][2] / k[9].n)}. ${exact(MN.nonAtg)} MANE coding sequences start with a codon other than ATG, and ${exact(MN.broken)} could not be read as a whole number of codons ending in a stop.`;
  drawLogo('#kozak-chart', 'kozak-context', MN.kozak, [-6, -5, -4, -3, -2, -1, 1, 2, 3, 4].map(signedPos), 6, 'Start codon context logo');
  const parts = [['TGA', 'var(--s1)'], ['TAA', 'var(--s2)'], ['TAG', 'var(--s3)']].map(([codon, color]) => ({ label: codon, value: MN.stops[codon], color }));
  waffle({ mount: $('#stop-chart'), parts, aria: 'Stop codons of MANE coding sequences' });
  attachData('#stop-chart', 'stop-codons', ['Stop codon', 'Coding sequences'], parts.map((p) => [p.label, p.value]));
};

const AA3 = {
  A: 'Ala', R: 'Arg', N: 'Asn', D: 'Asp', C: 'Cys', Q: 'Gln', E: 'Glu', G: 'Gly', H: 'His', I: 'Ile', L: 'Leu', K: 'Lys', M: 'Met', F: 'Phe', P: 'Pro', S: 'Ser', T: 'Thr', W: 'Trp', Y: 'Tyr', V: 'Val', '*': 'Stop',
};

const renderCodons = () => {
  const sense = senseCodons(MN.codons);
  const total = sense.reduce((s, c) => s + c.n, 0);
  const byCodon = Object.fromEntries(MN.codons.map((n, i) => [codonOf(i), n]));
  const aaTotals = {};
  for (const c of sense) aaTotals[GENETIC_CODE[c.codon]] = (aaTotals[GENETIC_CODE[c.codon]] ?? 0) + c.n;
  const cells = Object.fromEntries(Object.entries(GENETIC_CODE).map(([codon, aa]) => [codon, {
    aa, aa3: AA3[aa], perK: aa === '*' ? 0 : (byCodon[codon] / total) * 1000, share: aa === '*' ? null : byCodon[codon] / aaTotals[aa],
  }]));
  const max = Math.max(...Object.values(cells).map((c) => c.perK));
  const leu = Object.entries(GENETIC_CODE).filter(([, aa]) => aa === 'L').map(([codon]) => [codon, cells[codon].share]).sort((x, y) => y[1] - x[1]);
  const arg = Object.entries(GENETIC_CODE).filter(([, aa]) => aa === 'R').map(([codon]) => [codon, cells[codon].share]).sort((x, y) => y[1] - x[1]);
  const rarest = [...sense].sort((x, y) => x.n - y.n);
  const cgRun = rarest.findIndex((c) => !c.codon.includes('CG'));
  const cgNote = cgRun >= 2 ? ` The ${cgRun} least used codons, ${rarest.slice(0, cgRun).map((c) => c.codon).join(', ')}, all contain CG.` : '';
  $('#codons-dek').textContent = `${exact(total)} sense codons from ${exact(MN.coding)} coding sequences. ${leu[0][0]} carries ${pct(leu[0][1])} of leucines and ${leu.at(-1)[0]} ${pct(leu.at(-1)[1])}; among arginine codons, ${arg[0][0]} ${pct(arg[0][1])} and ${arg.at(-1)[0]} ${pct(arg.at(-1)[1])}.${cgNote}`;
  codonTable({
    mount: $('#codon-chart'),
    cells,
    max,
    aria: 'Codon usage per thousand codons',
    describe: (codon) => {
      const c = cells[codon];
      return c.aa === '*'
        ? tipLines(`${exact(MN.stops[codon])} coding sequences`, `${codon}, a stop codon`)
        : tipLines(`${c.perK.toFixed(2)} per thousand`, `${codon}, ${AMINO_ACIDS[c.aa]}`, `${pct(c.share)} of ${AA3[c.aa]} codons`, `${exact(byCodon[codon])} codons`);
    },
  });
  attachData('#codon-chart', 'codon-usage', ['Codon', 'Amino acid', 'Count', 'Per thousand', 'Share of amino acid'],
    Object.entries(cells).map(([codon, c]) => [codon, AA3[c.aa], c.aa === '*' ? MN.stops[codon] : byCodon[codon], c.aa === '*' ? null : Number(c.perK.toFixed(3)), c.share === null ? null : Number(c.share.toFixed(4))]));
  const aas = Object.entries(aaTotals).sort((x, y) => y[1] - x[1]);
  dotplot({
    mount: $('#amino-chart'),
    rows: aas.map(([aa, n], i) => ({ key: aa, label: `${AMINO_ACIDS[aa]} (${AA3[aa]})`, value: (n / total) * 100, hi: i === 0 })),
    scale: 'linear',
    format: (v) => `${v.toFixed(1)}%`,
    labelShare: 0.36,
    aria: 'Amino acid composition of MANE coding sequences',
    describe: (r) => tipLines(`${r.value.toFixed(2)}%`, r.label, `${exact(aaTotals[r.key])} residues`),
  });
  attachData('#amino-chart', 'amino-acids', ['Amino acid', 'Residues', 'Percent'], aas.map(([aa, n]) => [AMINO_ACIDS[aa], n, Number(((n / total) * 100).toFixed(3))]));
};

const P = snapshot.proteins;
const TG = snapshot.targets;
const HPA = snapshot.expression.hpa;
const GT = snapshot.expression.gtex;
const uniprotEntry = (acc) => `https://www.uniprot.org/uniprotkb/${acc}/entry`;
const proteinName = proteinLabel;
const genesOf = (p) => (String(p.gene).includes(';') ? `Genes: ${p.gene}` : null);
const sumOf = (list) => list.reduce((a, b) => a + b, 0);
const capped = (i, last, unit) => (i === last ? `${i + 1} or more ${unit}` : `${i + 1} ${unit}`);

const EXISTENCE = ['Evidence at protein level', 'Evidence at transcript level', 'Inferred from homology', 'Predicted', 'Uncertain'];

const renderProteome = () => {
  const { summary: s, edges, counts } = P.lengths;
  const [short] = P.shortest;
  $('#proteome-dek').textContent = `Quartiles: ${exact(s.q1)}, ${exact(s.median)} and ${exact(s.q3)} amino acids; the mean is ${exact(s.mean)}. ${pct(P.existence['Evidence at protein level'] / P.count)} of entries have evidence at protein level and ${exact(P.existence.Uncertain ?? 0)} are flagged uncertain. The shortest, ${proteinName(short)}, has ${short.length}.`;
  const bins = binsOf(edges);
  const last = bins.length - 1;
  const span = (i) => (i === 0 ? `fewer than ${exact(bins[0].hi)} amino acids` : i === last ? `${exact(bins[i].lo)} or more amino acids` : `${exact(bins[i].lo)} to ${exact(bins[i].hi)} amino acids`);
  histo({
    mount: $('#protein-length-chart'),
    bins,
    series: [{ key: 'proteins', label: 'Reviewed proteins', color: 'var(--data)', values: counts }],
    tickLabel: (b) => shortBases(b.lo),
    aria: 'Reviewed human proteins by length',
    describe: (i) => tipLines(`${exact(counts[i])} proteins`, span(i)),
  });
  attachData('#protein-length-chart', 'protein-lengths', ['Length from', 'Length to', 'Proteins'], bins.map((b, i) => [Math.round(b.lo), Math.round(b.hi), counts[i]]));
  const levels = EXISTENCE.filter((l) => P.existence[l]);
  dotplot({
    mount: $('#existence-chart'),
    rows: levels.map((l, i) => ({ key: l, label: l, value: P.existence[l], hi: i === 0 })),
    labelShare: 0.46,
    aria: 'Reviewed human proteins by existence level',
    describe: (r) => tipLines(`${exact(r.value)} proteins`, r.label, pct(r.value / P.count)),
  });
  attachData('#existence-chart', 'protein-existence', ['Existence level', 'Proteins'], levels.map((l) => [l, P.existence[l]]));
  const longest = P.longest.slice(0, 10);
  dotplot({
    mount: $('#longest-protein-chart'),
    rows: longest.map((p, i) => ({ key: p.acc, label: proteinName(p), value: p.length, hi: i === 0 })),
    scale: 'linear',
    aria: 'The longest reviewed human proteins',
    href: (r) => uniprotEntry(r.key),
    linkText: () => 'Open in UniProt',
    describe: (r) => tipLines(`${exact(r.value)} amino acids`, `${r.label}, UniProt ${r.key}`, genesOf(longest.find((p) => p.acc === r.key))),
  });
  attachData('#longest-protein-chart', 'longest-proteins', ['Gene', 'UniProt', 'Amino acids'], P.longest.map((p) => [p.gene, p.acc, p.length]));
};

const renderStructures = () => {
  const { counts, top, withAny } = P.pdb;
  const tail = counts.slice(1);
  const last = tail.length - 1;
  $('#structures-dek').textContent = `${exact(counts[0])} reviewed proteins have no PDB entry. Of the ${exact(withAny)} that do, ${exact(counts[1])} appear in one entry and ${exact(counts.at(-1))} in 50 or more.`;
  histo({
    mount: $('#pdb-count-chart'),
    bins: tail.map((n, i) => ({ n: i + 1 })),
    series: [{ key: 'proteins', label: 'Proteins', color: 'var(--data)', values: tail }],
    tickEvery: 1,
    centered: true,
    tickLabel: (b, i) => (i === 0 || (i + 1) % 10 === 0 ? (i === last ? '50+' : String(i + 1)) : null),
    aria: 'Proteins by number of PDB entries',
    describe: (i) => tipLines(`${exact(tail[i])} proteins`, capped(i, last, 'PDB entries').replace(/^1 PDB entries$/, '1 PDB entry')),
  });
  attachData('#pdb-count-chart', 'pdb-entries-per-protein', ['PDB entries', 'Proteins'], counts.map((n, i) => [i === counts.length - 1 ? `${i}+` : String(i), n]));
  dotplot({
    mount: $('#pdb-top-chart'),
    rows: top.map((p, i) => ({ key: p.acc, label: proteinName(p), value: p.n, hi: i === 0 })),
    scale: 'linear',
    aria: 'Proteins in the most PDB entries',
    href: (r) => uniprotEntry(r.key),
    linkText: () => 'Open in UniProt',
    describe: (r) => tipLines(`${exact(r.value)} PDB entries`, `${r.label}, UniProt ${r.key}`, genesOf(top.find((p) => p.acc === r.key))),
  });
  attachData('#pdb-top-chart', 'most-solved-proteins', ['Gene', 'UniProt', 'PDB entries'], top.map((p) => [p.gene, p.acc, p.n]));
};

const renderMembrane = () => {
  const tm = P.transmem.counts;
  const spans = tm.slice(1);
  const last = spans.length - 1;
  const [nuc] = HPA.locations;
  $('#membrane-dek').textContent = `${exact(tm[1])} proteins have one transmembrane segment and ${exact(sumOf(tm.slice(2)))} have two or more; ${exact(tm[7])} have exactly seven. ${exact(P.signal)} carry an N-terminal signal peptide. By immunofluorescence, HPA finds ${exact(nuc.n)} proteins mainly in the ${nuc.name.toLowerCase()}.`;
  histo({
    mount: $('#tm-chart'),
    bins: spans.map((n, i) => ({ n: i + 1 })),
    series: [{ key: 'proteins', label: 'Proteins', color: 'var(--data)', values: spans }],
    tickEvery: 1,
    centered: true,
    tickLabel: (b, i) => (i === 0 || (i + 1) % 5 === 0 ? (i === last ? '25+' : String(i + 1)) : null),
    aria: 'Proteins by number of transmembrane segments',
    describe: (i) => tipLines(`${exact(spans[i])} proteins`, capped(i, last, 'transmembrane segments').replace(/^1 transmembrane segments$/, '1 transmembrane segment')),
  });
  attachData('#tm-chart', 'transmembrane-segments', ['Segments', 'Proteins'], tm.map((n, i) => [i === tm.length - 1 ? `${i}+` : String(i), n]));
  const where = (sel, name, list, aria, what) => {
    dotplot({
      mount: $(sel),
      rows: list.slice(0, 20).map((l, i) => ({ key: l.name, label: l.name, value: l.n, hi: i === 0 })),
      scale: 'linear',
      labelShare: 0.46,
      aria,
      describe: (r) => tipLines(`${exact(r.value)} proteins`, r.label, what),
    });
    attachData(sel, name, ['Location', 'Proteins'], list.map((l) => [l.name, l.n]));
  };
  where('#location-chart', 'uniprot-locations', P.locations, 'Subcellular locations in UniProt', 'curated UniProt location');
  where('#hpa-location-chart', 'hpa-locations', HPA.locations, 'Main subcellular locations in the Human Protein Atlas', 'main location by immunofluorescence');
};

const renderDomains = () => {
  const { top, perProtein } = P.pfam;
  const last = perProtein.length - 1;
  $('#domains-dek').textContent = `${exact(perProtein[0])} reviewed proteins carry no Pfam domain, ${exact(perProtein[1])} carry one family and ${exact(sumOf(perProtein.slice(5)))} carry five or more. The next most common after ${top[0].short} are ${top.slice(1, 4).map((p) => `${p.short} (${exact(p.n)})`).join(', ')}.`;
  dotplot({
    mount: $('#pfam-chart'),
    rows: top.map((p, i) => ({ key: p.id, label: p.name ?? p.id, value: p.n, hi: i === 0 })),
    scale: 'linear',
    labelShare: 0.36,
    labelMax: 330,
    aria: 'The most common Pfam families in human proteins',
    href: (r) => `https://www.ebi.ac.uk/interpro/entry/pfam/${r.key}/`,
    linkText: () => 'Open in InterPro',
    describe: (r) => {
      const p = top.find((x) => x.id === r.key);
      return tipLines(`${exact(p.n)} proteins`, `${p.id} · ${p.short ?? ''}`, p.name);
    },
  });
  attachData('#pfam-chart', 'pfam-families', ['Pfam', 'Short name', 'Name', 'Proteins'], top.map((p) => [p.id, p.short, p.name, p.n]));
  histo({
    mount: $('#pfam-per-protein-chart'),
    bins: perProtein.map((n, i) => ({ n: i })),
    series: [{ key: 'proteins', label: 'Proteins', color: 'var(--data)', values: perProtein }],
    tickEvery: 1,
    centered: true,
    tickLabel: (b, i) => (i % 5 === 0 ? (i === last ? `${i}+` : String(i)) : null),
    aria: 'Proteins by number of distinct Pfam families',
    describe: (i) => tipLines(`${exact(perProtein[i])} proteins`, i === last ? `${i} or more Pfam families` : `${i} Pfam ${i === 1 ? 'family' : 'families'}`),
  });
  attachData('#pfam-per-protein-chart', 'pfam-families-per-protein', ['Pfam families', 'Proteins'], perProtein.map((n, i) => [i === last ? `${i}+` : String(i), n]));
};

const PTMS = [
  ['phosphoS', 'Phosphoserine'], ['phosphoT', 'Phosphothreonine'], ['phosphoY', 'Phosphotyrosine'], ['acetyl', 'Acetylation'], ['methyl', 'Methylation'],
  ['nGlyco', 'N-linked glycosylation'], ['oGlyco', 'O-linked glycosylation'], ['ubiquitin', 'Ubiquitin cross-link'], ['sumo', 'SUMO cross-link'],
];

const renderModifications = () => {
  const { share, proteins } = P.disorder;
  $('#modifications-dek').textContent = `${exact(P.disorder.residues)} of ${exact(proteomeResidues(snapshot))} residues, in ${exact(proteins)} proteins, fall in disordered regions; ${exact(share.at(-1))} proteins are disordered over 90% or more of their length. N-linked glycans are annotated at ${exact(P.ptm.nGlyco)} sites and O-linked at ${exact(P.ptm.oGlyco)}; ${exact(P.ptm.sumo)} lysines carry a SUMO cross-link.`;
  histo({
    mount: $('#disorder-chart'),
    bins: share.map((n, i) => ({ lo: i / 10 })),
    series: [{ key: 'proteins', label: 'Proteins', color: 'var(--data)', values: share }],
    tickEvery: 1,
    tickLabel: (b, i) => (i % 2 === 0 ? `${i * 10}%` : null),
    aria: 'Proteins by share of residues in disordered regions',
    describe: (i) => tipLines(`${exact(share[i])} proteins`, `${i * 10}% to ${(i + 1) * 10}% of residues disordered`),
  });
  attachData('#disorder-chart', 'disorder-share', ['Disordered from', 'Disordered to', 'Proteins'], share.map((n, i) => [`${i * 10}%`, `${(i + 1) * 10}%`, n]));
  dotplot({
    mount: $('#ptm-chart'),
    rows: PTMS.map(([key, label], i) => ({ key, label, value: P.ptm[key], hi: i === 0 })),
    labelShare: 0.42,
    aria: 'Annotated modification sites by kind',
    describe: (r) => tipLines(`${exact(r.value)} sites`, r.label),
  });
  attachData('#ptm-chart', 'modification-sites', ['Modification', 'Sites'], PTMS.map(([key, label]) => [label, P.ptm[key]]));
};

const LEVELS = [['Tclin', 'var(--s1)'], ['Tchem', 'var(--s2)'], ['Tbio', 'var(--s3)'], ['Tdark', 'var(--data-2)']];

const halfBins = (sel, name, pairs, aria, what) => {
  histo({
    mount: $(sel),
    bins: pairs.map(([x]) => ({ x })),
    series: [{ key: 'targets', label: 'Targets', color: 'var(--data)', values: pairs.map(([, n]) => n) }],
    tickEvery: 1,
    tickLabel: (b) => (Number.isInteger(b.x) && b.x % 2 === 0 ? signed(b.x, String) : null),
    aria,
    describe: (i) => tipLines(`${exact(pairs[i][1])} targets`, `${what} ${signed(pairs[i][0], String)} to ${signed(pairs[i][0] + 0.5, String)}`),
  });
  attachData(sel, name, [`${what} from`, 'Targets'], pairs.map(([x, n]) => [x, n]));
};

const renderTargets = () => {
  const { levels, count, families } = TG;
  const fam = (name) => families.find((f) => f.family === name);
  const famTotal = (f) => LEVELS.reduce((n, [l]) => n + f[l], 0);
  const kin = fam('Kinase');
  const gpcr = fam('GPCR');
  const olf = fam('oGPCR');
  $('#targets-dek').textContent = `${pct(levels.Tclin / count)} of human proteins are approved-drug targets and ${pct(levels.Tchem / count)} have potent ligands without one. ${exact(kin.Tclin)} of ${exact(famTotal(kin))} kinases and ${exact(gpcr.Tclin)} of ${exact(famTotal(gpcr))} non-olfactory GPCRs are Tclin; ${exact(olf.Tdark)} of ${exact(famTotal(olf))} olfactory receptors are Tdark.`;
  const parts = LEVELS.map(([l, color]) => ({ label: l, value: levels[l], color }));
  waffle({ mount: $('#tdl-chart'), parts, aria: 'Human proteins by target development level' });
  attachData('#tdl-chart', 'target-development-levels', ['Level', 'Targets'], parts.map((p) => [p.label, p.value]));
  hbars({
    mount: $('#family-chart'),
    rows: families.map((f) => ({
      key: f.family, label: f.family, f, values: Object.fromEntries(LEVELS.map(([l]) => [l, (f[l] / famTotal(f)) * 100])),
    })),
    segments: LEVELS.map(([key, color]) => ({ key, label: key, color })),
    labelShare: 0.22,
    max: 100,
    R: 96,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => `${exact(famTotal(r.f))} targets`,
    aria: 'Target development level by protein family, as shares',
    describe: (r) => tipLines(`${exact(famTotal(r.f))} targets`, r.key, ...LEVELS.map(([l]) => `${l}: ${exact(r.f[l])}, ${pct(r.f[l] / famTotal(r.f))}`)),
  });
  attachData('#family-chart', 'target-levels-by-family', ['Family', 'Tclin', 'Tchem', 'Tbio', 'Tdark'], families.map((f) => [f.family, f.Tclin, f.Tchem, f.Tbio, f.Tdark]));
  halfBins('#pubmed-score-chart', 'pharos-pubmed-score', TG.pubmed, 'Targets by log PubMed score', 'log PubMed score');
  halfBins('#novelty-chart', 'pharos-novelty', TG.novelty, 'Targets by log novelty', 'log novelty');
};

const CONTEXTS = [['tissue', 'Tissues'], ['singleCell', 'Single cell types'], ['blood', 'Blood cells'], ['brain', 'Brain regions'], ['cellLine', 'Cell lines'], ['cancer', 'Cancers']];
const KINDS = [['enriched', 'Enriched', 'var(--s1)'], ['group', 'Group enriched', 'var(--s2)'], ['enhanced', 'Enhanced', 'var(--s3)'], ['low', 'Low specificity', 'var(--data-2)'], ['none', 'Not detected', 'var(--s5)']];

const renderSpecificity = () => {
  const c = HPA.contexts;
  $('#specificity-dek').textContent = `The single cell type section calls ${exact(c.singleCell.enriched)} genes cell type enriched and the brain section ${exact(c.brain.enriched)} region enriched. ${exact(c.tissue.none)} genes are not detected in any tissue, and ${exact(c.blood.none)} in no immune cell type. The most cell-type-enriched genes belong to ${HPA.cellTypes[0].cell.toLowerCase()} (${exact(HPA.cellTypes[0].n)}).`;
  hbars({
    mount: $('#contexts-chart'),
    rows: CONTEXTS.map(([key, label]) => ({ key, label, values: c[key] })),
    segments: KINDS.map(([key, label, color]) => ({ key, label, color })),
    aria: 'Genes by specificity call in six HPA sections',
    describe: (r) => {
      const total = sumOf(Object.values(r.values));
      return tipLines(`${exact(total)} genes`, r.label, ...KINDS.map(([k, label]) => `${label}: ${exact(r.values[k])}, ${pct(r.values[k] / total)}`));
    },
  });
  attachData('#contexts-chart', 'hpa-specificity', ['Section', ...KINDS.map(([, label]) => label)], CONTEXTS.map(([key, label]) => [label, ...KINDS.map(([k]) => c[key][k])]));
  const ranked = (sel, name, list, field, aria, what) => {
    dotplot({
      mount: $(sel),
      rows: list.slice(0, 20).map((t, i) => ({ key: t[field], label: t[field], value: t.n, hi: i === 0 })),
      scale: 'linear',
      labelShare: 0.44,
      aria,
      describe: (r) => tipLines(`${exact(r.value)} genes`, `${what} ${r.label}`),
    });
    attachData(sel, name, [field === 'tissue' ? 'Tissue' : 'Cell type', 'Enriched genes'], list.map((t) => [t[field], t.n]));
  };
  ranked('#enriched-chart', 'tissue-enriched-genes', HPA.enriched, 'tissue', 'Tissues with the most tissue-enriched genes', 'tissue enriched in');
  ranked('#celltype-chart', 'cell-type-enriched-genes', HPA.cellTypes, 'cell', 'Cell types with the most cell-type-enriched genes', 'cell type enriched in');
};

const shortTissue = (name) => name.replace(/\s*\([^)]*\)/g, '');

const renderGtex = () => {
  const tissues = gtexTissues(snapshot);
  const byExpressed = [...tissues].sort((a, b) => b.expressed - a.expressed);
  const bySpecific = [...tissues].sort((a, b) => b.specific - a.specific);
  const eq = tissues.filter((t) => t.eGenes && t.samples);
  const r = pearson(eq.map((t) => t.samples), eq.map((t) => t.eGenes));
  $('#gtex-dek').textContent = `${exact(GT.genes)} genes in the median expression table, ${exact(GT.correlated)} of them at 1 TPM in at least one tissue. ${exact(bySpecific[0].specific)} genes are at least five times higher in ${bySpecific[0].name.toLowerCase()} than in any other tissue; ${bySpecific[1].name.toLowerCase()} is next with ${exact(bySpecific[1].specific)}. Across ${eq.length} tissues with eQTL results, eGene count and sample size correlate with r = ${r.toFixed(2)}.`;
  const names = tissues.map((t) => shortTissue(t.name));
  const off = GT.correlation.flatMap((row, i) => row.filter((_, j) => j !== i)).sort((a, b) => a - b);
  const lo = off[Math.floor(off.length * 0.02)];
  heatmap({
    mount: $('#tissue-corr-chart'),
    labels: names,
    values: GT.correlation,
    order: GT.order,
    lo,
    hi: 1,
    color: 'var(--s2)',
    legend: rampLegend(lo, 1, (v) => (v <= lo ? `≤${v.toFixed(2)}` : v.toFixed(2)), 'var(--s2)'),
    aria: 'Correlation of median expression between GTEx tissues',
    describe: (i, j) => tipLines(`r = ${GT.correlation[i][j].toFixed(3)}`, `${tissues[i].name} and ${tissues[j].name}`, `over ${exact(GT.correlated)} genes`),
  });
  attachData('#tissue-corr-chart', 'gtex-tissue-correlation', ['Tissue', ...GT.order.map((i) => tissues[i].id)], GT.order.map((i) => [tissues[i].name, ...GT.order.map((j) => GT.correlation[i][j])]));
  const people = (t) => (t.samples ? `${exact(t.samples)} RNA-seq samples, ${exact(t.female)} female and ${exact(t.male)} male` : 'sample counts not listed by the GTEx API');
  dotplot({
    mount: $('#tissue-genes-chart'),
    rows: byExpressed.map((t, i) => ({ key: t.id, label: shortTissue(t.name), value: t.expressed, hi: i === 0 })),
    scale: 'linear',
    rowH: 18,
    labelShare: 0.3,
    labelMax: 330,
    aria: 'Genes at 1 TPM or more in each GTEx tissue',
    describe: (row) => {
      const t = tissues.find((x) => x.id === row.key);
      return tipLines(`${exact(t.expressed)} genes`, t.name, `${exact(t.specific)} at least five times higher here than in any other tissue`, people(t));
    },
  });
  attachData('#tissue-genes-chart', 'gtex-expressed-genes', ['Tissue', 'GTEx ID', 'Genes at 1 TPM', 'Tissue-specific genes', 'RNA-seq samples', 'Female', 'Male'],
    byExpressed.map((t) => [t.name, t.id, t.expressed, t.specific, t.samples, t.female, t.male]));
  const top = eq.reduce((a, b) => (b.eGenes > a.eGenes ? b : a));
  scatter({
    mount: $('#egene-chart'),
    points: eq.map((t) => ({
      key: t.id, x: t.samples, y: t.eGenes, hi: t === top, t,
    })),
    xFormat: exact,
    xLabel: 'RNA-seq samples',
    yLabel: 'eGenes',
    trend: linearFit(eq.map((t) => t.samples), eq.map((t) => t.eGenes)),
    aria: 'eGenes against sample size across GTEx tissues',
    describe: (p) => tipLines(`${exact(p.y)} eGenes`, p.t.name, `${exact(p.x)} samples`, p.t.sGenes ? `${exact(p.t.sGenes)} sGenes` : null),
  });
  attachData('#egene-chart', 'gtex-egenes', ['Tissue', 'RNA-seq samples', 'eGenes', 'sGenes'], eq.map((t) => [t.name, t.samples, t.eGenes, t.sGenes]));
};

const SECRETOME_ORDER = (a, b) => b[1] - a[1];

const renderPlasma = () => {
  const { blood } = HPA;
  const bins = binsOf(blood.edges);
  $('#plasma-dek').textContent = `${exact(blood.im.count)} proteins have an immunoassay level and ${exact(blood.ms.count)} a mass spectrometry level; ${exact(blood.im.count + blood.ms.count - blood.proteins)} have both. ${exact(HPA.secretome['Secreted to blood'] ?? 0)} genes are classed as secreted to blood.`;
  histo({
    mount: $('#plasma-chart'),
    bins,
    series: [
      { key: 'im', label: 'Immunoassay', color: 'var(--s1)', values: blood.im.counts },
      { key: 'ms', label: 'Mass spectrometry', color: 'var(--s2)', values: blood.ms.counts },
    ],
    tickEvery: 1,
    tickLabel: (b, i) => (i % 3 === 0 ? perLitre(b.lo) : null),
    aria: 'Proteins by plasma concentration and method',
    describe: (i) => tipLines(`${perLitre(bins[i].lo)} to ${perLitre(bins[i].hi)}`, 'plasma concentration', `Immunoassay: ${exact(blood.im.counts[i])} proteins`, `Mass spectrometry: ${exact(blood.ms.counts[i])} proteins`),
  });
  attachData('#plasma-chart', 'plasma-concentrations', ['From pg per L', 'To pg per L', 'Immunoassay', 'Mass spectrometry'], bins.map((b, i) => [b.lo, b.hi, blood.im.counts[i], blood.ms.counts[i]]));
  dotplot({
    mount: $('#plasma-top-chart'),
    rows: blood.top.map((p, i) => ({ key: p.gene, label: p.gene, value: plasmaLevel(p), hi: i === 0 })),
    format: perLitre,
    tickFormat: perLitre,
    aria: 'Proteins with the highest plasma concentration',
    describe: (r) => {
      const p = blood.top.find((x) => x.gene === r.key);
      return tipLines(perLitre(r.value), p.gene, p.im ? `Immunoassay: ${perLitre(p.im)}` : 'No immunoassay level', p.ms ? `Mass spectrometry: ${perLitre(p.ms)}` : 'No mass spectrometry level');
    },
  });
  attachData('#plasma-top-chart', 'plasma-top-proteins', ['Gene', 'Immunoassay pg per L', 'Mass spectrometry pg per L'], blood.top.map((p) => [p.gene, p.im, p.ms]));
  const sec = Object.entries(HPA.secretome).sort(SECRETOME_ORDER);
  dotplot({
    mount: $('#secretome-chart'),
    rows: sec.map(([k, n], i) => ({ key: k, label: k, value: n, hi: i === 0 })),
    scale: 'linear',
    labelShare: 0.5,
    labelMax: 330,
    aria: 'Secreted proteins by destination',
    describe: (r) => tipLines(`${exact(r.value)} genes`, r.label),
  });
  attachData('#secretome-chart', 'hpa-secretome', ['Secretome location', 'Genes'], sec);
};

const decode = () => {
  const value = $('#decoder-input').value;
  if (!value.trim()) {
    fill('#decoder-out', h('p', { class: 'none' }, 'Type an identifier, or pick an example below.'));
    return;
  }
  const hit = identify(value, { bands: A.bands });
  if (!hit) {
    fill('#decoder-out', h('p', { class: 'none' }, `No match for “${value.trim()}”.`));
    return;
  }
  const known = EXAMPLES.find((e) => identify(e.id)?.id === hit.id)?.is;
  fill('#decoder-out', h('dl', {},
    h('dt', {}, 'Type'), h('dd', {}, h('span', { class: 'type' }, hit.type)),
    h('dt', {}, 'ID'), h('dd', {}, h('code', {}, hit.id)),
    h('dt', {}, 'Source'), h('dd', {}, hit.db),
    h('dt', {}, 'Format'), h('dd', {}, hit.note),
    ...hit.parts.flatMap(([k, v]) => [h('dt', {}, k), h('dd', {}, v)]),
    known && h('dt', {}, 'Known as'), known && h('dd', {}, known),
    h('dt', {}, 'Open'), h('dd', {}, ext(hit.url, hit.url))));
};

const renderDecoder = () => {
  fill('#examples', ...EXAMPLES.map((e) => h('button', { type: 'button', class: 'example', 'data-id': e.id, title: e.is }, e.id)));
  $('#examples').addEventListener('click', (event) => {
    const btn = event.target.closest('[data-id]');
    if (!btn) return;
    $('#decoder-input').value = btn.dataset.id;
    decode();
  });
  $('#decoder-input').addEventListener('input', decode);
  decode();
};

const DATA_SOURCES = [
  ['NCBI GRCh38.p14 assembly report', FILES.grch38Report],
  ['NCBI GRCh38.p14 assembly statistics', FILES.grch38Stats],
  ['NCBI T2T-CHM13v2.0 assembly report', FILES.t2tReport],
  ['NCBI T2T-CHM13v2.0 assembly statistics', FILES.t2tStats],
  ['UCSC hg38 cytogenetic bands', FILES.cytobands],
  ['UCSC hg38 gaps', FILES.gaps],
  ['UCSC hg38 centromere models', FILES.centromeres],
  ['UCSC hg38 2bit sequence', FILES.hg38TwoBit],
  ['NCBI MANE release folder', FILES.maneDir],
  ['UCSC hg38 RepeatMasker table', FILES.rmsk],
  ['UCSC hg38 segmental duplications', FILES.superDups],
  ['UCSC hg38 Tandem Repeats Finder table', FILES.simpleRepeat],
  ['GENCODE latest release folder', FILES.gencodeDir],
  ['HGNC complete set', FILES.hgnc],
  ['NCBI Homo sapiens gene_info', FILES.geneInfo],
  ['UniProt reviewed human entries, streamed as TSV', FILES.uniprot],
  ['InterPro Pfam entry API', FILES.interproPfam],
  ['Pharos GraphQL API', FILES.pharos],
  ['Human Protein Atlas proteinatlas.tsv', FILES.hpa],
  ['GTEx v10 gene median TPM', FILES.gtexMedian],
  ['GTEx Portal tissue summary API', FILES.gtexTissues],
];

const renderSources = () => {
  const host = (url) => new URL(url).hostname.replace(/^www\./, '');
  const item = (title, url) => h('li', { 'data-search': `${title} ${url}` }, ext(url, title), ' ', h('small', {}, host(url)));
  const group = (title, items) => h('div', { class: 'srcgroup' }, h('h3', {}, title), h('ul', {}, items));
  const titled = (g) => Object.entries(SOURCES).filter(([key]) => SOURCE_TITLES[key]?.[0] === g).map(([key, url]) => item(SOURCE_TITLES[key][1], url));
  fill('#source-list',
    group('Data files', DATA_SOURCES.map(([name, url]) => item(name, url))),
    ...SOURCE_GROUPS.map((g) => [g, titled(g)]).filter(([, items]) => items.length).map(([g, items]) => group(g, items)));
  const when = day(snapshot.generated);
  fill('#data-use',
    'Every number on this page, as of ', when, ', is in ', h('a', { href: 'snapshot.json', download: '' }, 'snapshot.json'),
    '. The code and the weekly rebuild are on ', ext(REPO, 'GitHub'),
    `. Cite as: Human Genome in Numbers, snapshot of ${when}. ${SITE}`);
  fill('#foot',
    h('p', {}, `Built from public NCBI, UCSC, GENCODE and EBI data. Snapshot ${snapshot.generated.slice(0, 16).replace('T', ' ')} UTC. Not an official NIH, NCBI, EMBL-EBI or UCSC site.`),
    h('p', { class: 'made' }, 'Made by a person, with help from AI tools. Every line, chart and figure was reviewed before publishing; the numbers now refresh on their own each Sunday.'));
};

const SEARCHABLE = ['#source-list li', '.fig details.data tbody tr'];

const applySearch = () => {
  for (const sel of SEARCHABLE) {
    for (const node of $$(sel)) node.hidden = !matches(node.dataset.search ?? node.textContent, state.query);
  }
};

let searchTimer;
const onSearch = (event) => {
  clearTimeout(searchTimer);
  searchTimer = setTimeout(() => {
    state.query = event.target.value;
    applySearch();
  }, 120);
};

const THEME_KEY = 'human-genome-theme';
const darkQuery = matchMedia('(prefers-color-scheme: dark)');
const theme = () => document.documentElement.dataset.theme || (darkQuery.matches ? 'dark' : 'light');

const store = {
  get: (key) => {
    try {
      return localStorage.getItem(key);
    } catch {
      return null;
    }
  },
  set: (key, value) => {
    try {
      localStorage.setItem(key, value);
    } catch {
      return false;
    }
    return true;
  },
};

const drawCharts = () => {
  renderChromosomes();
  renderGaps();
  renderT2t();
  renderComposition();
  renderTransposons();
  renderLandscape();
  renderSegdups();
  renderTandem();
  renderGenecount();
  renderDensity();
  renderGenesize();
  renderStructure();
  renderIsoforms();
  renderNames();
  renderPseudogenes();
  renderImmune();
  renderSplicing();
  renderTranslation();
  renderCodons();
  renderProteome();
  renderStructures();
  renderMembrane();
  renderDomains();
  renderModifications();
  renderTargets();
  renderSpecificity();
  renderGtex();
  renderPlasma();
  applySearch();
};

const paintTheme = () => {
  $('#theme').textContent = theme() === 'dark' ? 'Light' : 'Dark';
  $('#theme').setAttribute('aria-pressed', String(theme() === 'dark'));
};

const initTheme = () => {
  const saved = store.get(THEME_KEY);
  if (saved) document.documentElement.dataset.theme = saved;
  paintTheme();
  darkQuery.addEventListener('change', () => {
    paintTheme();
    drawCharts();
  });
  $('#theme').addEventListener('click', () => {
    const next = theme() === 'dark' ? 'light' : 'dark';
    document.documentElement.dataset.theme = next;
    store.set(THEME_KEY, next);
    paintTheme();
    drawCharts();
  });
};

const applyLens = (lens, push) => {
  state.lens = lens;
  $$('#lenses button').forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.lens === lens)));
  $$('main > section.story').forEach((sec) => { sec.hidden = !inLens(sec.dataset.lens, lens); });
  for (const ch of $$('main > .chapter')) {
    const members = [];
    for (let n = ch.nextElementSibling; n && !n.classList.contains('chapter'); n = n.nextElementSibling) members.push(n);
    ch.hidden = members.every((n) => n.hidden);
  }
  $$('#jump a').forEach((a) => { a.hidden = $(`#${a.dataset.for}`).hidden; });
  $('#lens-hint').textContent = LENSES.find((l) => l.id === lens).hint;
  if (push) {
    history.replaceState(null, '', `${location.pathname}${lensSearch(location.search, lens)}${location.hash}`);
    drawCharts();
  }
};

const renderLenses = () => {
  fill('#lenses', ...LENSES.map((l) => h('button', { type: 'button', 'data-lens': l.id, 'aria-pressed': 'false', title: l.hint }, l.label)));
  fill('#jump', ...$$('main > .chapter, main > section[data-nav]').map((n) => (n.classList.contains('chapter')
    ? h('a', { href: `#${n.id}`, 'data-for': n.id, class: 'jump-ch' }, $('.chapter-title', n).textContent)
    : h('a', { href: `#${n.id}`, 'data-for': n.id }, n.dataset.nav))));
  $('#lenses').addEventListener('click', (event) => {
    const btn = event.target.closest('button[data-lens]');
    if (btn) applyLens(btn.dataset.lens, true);
  });
};

let resizeTimer;
let lastWidth = 0;
const onResize = () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    if (window.innerWidth === lastWidth) return;
    lastWidth = window.innerWidth;
    drawCharts();
  }, 160);
};

const init = () => {
  initTheme();
  renderLenses();
  renderHead();
  renderLede();
  renderDecoder();
  renderSources();
  drawCharts();
  applyLens(readLens(location.search), false);
  observeFigures();
  lastWidth = window.innerWidth;
  window.addEventListener('resize', onResize);
  $('#q').addEventListener('input', onSearch);
};

init();
