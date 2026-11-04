import snapshot from './snapshot.js';
import {
  SOURCES, SOURCE_TITLES, SOURCE_GROUPS, EXAMPLES, PARS, GENE_THERAPIES, GENE_THERAPIES_CHECKED, DIAGNOSTIC_YIELD, CFDNA_SCREENING, CARRIER_TIER3, BIOBANKS, GNOMAD_V4,
} from './facts.js';
import {
  exact, cell, compact, bases, coords, signed, identify, matches, LENSES, readLens, inLens, lensSearch, toCsv, pct, longDate, information, GENETIC_CODE, AMINO_ACIDS, AA3, codonOf, perLitre, CHANNELS, SUBSTITUTIONS, wordy,
} from './lib.js';
import {
  pearson, linearFit, linearEdges, summary,
} from './stats.js';
import {
  FILES, TRIAL_MODALITIES, assemblyPage, gdcSsmPage,
} from './queries.js';
import {
  findings, assemblyNumbers, sequenceNumbers, codingCounts, geneDensity, senseCodons, gtexTissues, plasmaLevel, proteomeResidues, proteinLabel, trnaPairs, titv,
} from './findings.js';
import {
  ideogram, hbars, dotplot, strips, matrix, rampLegend, treemap, stackedArea, waffle, columns, histo, scatter, dotmap, dotRows, manhattan, logo, codonTable, heatmap, genomeMap, dumbbell, grid as shareGrid, lollipop as lollipopChart, spectrum, ranges, intervals, lines, observeFigures, countUp, tipLines, measureWidths, forgetWidths,
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
  const scroll = h('div', { class: 'scroll' });
  const details = h('details', { class: 'data' }, h('summary', {}, 'Table'), scroll);
  details.addEventListener('toggle', () => {
    if (!details.open || scroll.firstChild) return;
    scroll.append(h('table', { class: 'data' },
      h('thead', {}, h('tr', {}, heads.map((c) => h('th', { scope: 'col' }, c)))),
      h('tbody', {}, rows.map((r) => h('tr', { 'data-search': r.join(' '), hidden: !matches(r.join(' '), state.query) }, r.map((v) => h('td', { class: typeof v === 'number' ? 'n' : null }, v === null ? '—' : typeof v === 'number' ? cell(v) : v)))))));
  });
  const csv = h('a', { class: 'csv figdata', download: `${name}.csv`, href: '#' }, 'CSV');
  let written = false;
  const write = () => {
    if (written) return;
    written = true;
    csv.setAttribute('href', `data:text/csv;charset=utf-8,${encodeURIComponent(toCsv(heads, rows))}`);
  };
  for (const type of ['pointerdown', 'focus', 'click']) csv.addEventListener(type, write);
  caption.append(h('span', { class: 'figdata' }, ' ', details), csv);
};

const day = (iso) => new Date(iso).toLocaleDateString('en-GB', { day: 'numeric', month: 'long', year: 'numeric', timeZone: 'UTC' });

const renderHead = () => {
  $('#stamp').textContent = `Data as of ${day(snapshot.generated)}. Rebuilt every Sunday from public genome and clinical files.`;
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

const NET = snapshot.networks;
const RX = snapshot.pathways;
const GO = snapshot.ontology;
const trimmed = (counts) => counts.slice(0, counts.findLastIndex((n) => n > 0) + 1);
const wholeRange = (b, unit) => {
  const lo = Math.ceil(b.lo);
  const hi = Math.ceil(b.hi) - 1;
  return lo >= hi ? `${exact(lo)} ${unit}` : `${exact(lo)} to ${exact(hi)} ${unit}`;
};

const CHANNEL_NAMES = {
  textmining: 'Text mining', coexpression: 'Co-expression', experimental: 'Experiments', database: 'Curated databases', neighborhood: 'Gene neighbourhood', cooccurence: 'Co-occurrence', fusion: 'Gene fusion',
};

const renderInteractions = () => {
  const { degree, edges, physical } = NET;
  const counts = trimmed(degree.counts);
  const bins = binsOf(degree.edges).slice(0, counts.length);
  $('#interactions-dek').textContent = `${exact(edges.all)} human protein pairs carry a STRING score; ${exact(edges.medium)} reach 0.4 and ${exact(edges.highest)} reach 0.9. At 0.7 the median protein has ${exact(degree.summary.median)} partners and ${exact(degree.counts[0])} have one. The physical subnetwork, limited to evidence that proteins bind, holds ${exact(physical.edges.high)} links at 0.7 among ${exact(physical.proteins.high)} proteins.`;
  histo({
    mount: $('#degree-chart'),
    bins,
    series: [{ key: 'proteins', label: 'Proteins', color: 'var(--data)', values: counts }],
    yLog: true,
    tickEvery: 5,
    tickLabel: (b) => shortBases(b.lo),
    aria: 'Proteins by number of high-confidence partners, log scales',
    describe: (i) => tipLines(`${exact(counts[i])} proteins`, wholeRange(bins[i], 'partners')),
  });
  attachData('#degree-chart', 'string-degree', ['Partners from', 'Partners to', 'Proteins'], bins.map((b, i) => [Math.ceil(b.lo), Math.ceil(b.hi) - 1, counts[i]]));
  dotplot({
    mount: $('#hub-chart'),
    rows: NET.hubs.map((p, i) => ({ key: p.id, label: p.name, value: p.degree, hi: i === 0 })),
    scale: 'linear',
    aria: 'Proteins with the most high-confidence partners',
    href: (r) => `https://string-db.org/network/${r.key}`,
    linkText: () => 'Open in STRING',
    describe: (r) => tipLines(`${exact(r.value)} partners`, r.label, `STRING ${r.key}`),
  });
  attachData('#hub-chart', 'string-hubs', ['Protein', 'STRING ID', 'Partners at 0.7'], NET.hubs.map((p) => [p.name, p.id, p.degree]));
  const channels = Object.keys(CHANNEL_NAMES).sort((a, b) => NET.channels[b] - NET.channels[a]);
  hbars({
    mount: $('#channel-chart'),
    rows: channels.map((c) => ({ key: c, label: CHANNEL_NAMES[c], values: { only: NET.only[c], shared: NET.channels[c] - NET.only[c] } })),
    segments: [{ key: 'only', label: 'This channel only', color: 'var(--s1)' }, { key: 'shared', label: 'With other channels', color: 'var(--data-2)' }],
    labelShare: 0.3,
    aria: 'High-confidence links by evidence channel',
    describe: (r) => tipLines(`${exact(NET.channels[r.key])} links`, r.label, `${exact(NET.only[r.key])} rest on this channel alone`, pct(NET.channels[r.key] / edges.high)),
  });
  attachData('#channel-chart', 'string-channels', ['Channel', 'Links at 0.7 with this channel', 'Links on this channel alone'], channels.map((c) => [CHANNEL_NAMES[c], NET.channels[c], NET.only[c]]));
  const scoreBins = binsOf(NET.scores.edges);
  const CUTS = [150, 400, 700, 900];
  histo({
    mount: $('#score-chart'),
    bins: scoreBins,
    series: [{ key: 'links', label: 'Links', color: 'var(--data)', values: NET.scores.counts }],
    yLog: true,
    tickEvery: 1,
    tickLabel: (b) => (CUTS.includes(b.lo) ? String(b.lo / 1000) : null),
    aria: 'Human links by combined score, log scale',
    describe: (i) => tipLines(`${exact(NET.scores.counts[i])} links`, `score ${(scoreBins[i].lo / 1000).toFixed(2)} to ${(scoreBins[i].hi / 1000).toFixed(2)}`),
  });
  attachData('#score-chart', 'string-scores', ['Score from', 'Score to', 'Links'], scoreBins.map((b, i) => [b.lo / 1000, b.hi / 1000, NET.scores.counts[i]]));
};

const reactomeLink = (id) => `https://reactome.org/PathwayBrowser/#/${id}`;

const renderPathways = () => {
  const tops = RX.topLevel;
  const ev = RX.evidence;
  $('#pathways-dek').textContent = `The hierarchy runs ${RX.depth.length} levels deep, from ${exact(RX.depth[0])} top-level pathways down. ${exact(RX.breadth[1])} proteins sit under a single top-level pathway and ${exact(sumOf(RX.breadth.slice(5)))} under five or more. ${pct((ev.TAS ?? 0) / sumOf(Object.values(ev)))} of the human mappings are curated (TAS); the rest are inferred from other species (IEA).`;
  dotplot({
    mount: $('#toplevel-chart'),
    rows: tops.map((t, i) => ({ key: t.id, label: t.name, value: t.proteins, hi: i === 0 })),
    scale: 'linear',
    rowH: 20,
    labelShare: 0.46,
    aria: 'Proteins in each top-level Reactome pathway',
    href: (r) => reactomeLink(r.key),
    linkText: () => 'Open in Reactome',
    describe: (r) => {
      const t = tops.find((x) => x.id === r.key);
      return tipLines(`${exact(t.proteins)} proteins`, t.name, `${exact(t.pathways)} pathways, itself included`);
    },
  });
  attachData('#toplevel-chart', 'reactome-top-level', ['Pathway', 'Reactome ID', 'Proteins', 'Pathways below and including'], tops.map((t) => [t.name, t.id, t.proteins, t.pathways]));
  const { matrix: jac, order } = RX.overlap;
  const off = jac.flatMap((row, i) => row.filter((_, j) => j !== i));
  const hi = Math.max(...off);
  heatmap({
    mount: $('#pathway-overlap-chart'),
    labels: tops.map((t) => t.name),
    values: jac,
    order,
    lo: 0,
    hi,
    color: 'var(--s4)',
    legend: rampLegend(0, hi, (v) => (v >= hi ? `≥${v.toFixed(2)}` : v.toFixed(2)), 'var(--s4)'),
    aria: 'Shared proteins between top-level Reactome pathways',
    describe: (i, j) => tipLines(`Jaccard index ${jac[i][j].toFixed(3)}`, `${tops[i].name} and ${tops[j].name}`, `${exact(tops[i].proteins)} and ${exact(tops[j].proteins)} proteins`),
  });
  attachData('#pathway-overlap-chart', 'reactome-overlap', ['Pathway', ...order.map((i) => tops[i].id)], order.map((i) => [tops[i].name, ...order.map((j) => jac[i][j])]));
  const breadth = RX.breadth.slice(1, RX.breadth.findLastIndex((v) => v > 0) + 1);
  histo({
    mount: $('#breadth-chart'),
    bins: breadth.map((n, i) => ({ n: i + 1 })),
    series: [{ key: 'proteins', label: 'Proteins', color: 'var(--data)', values: breadth }],
    yLog: true,
    centered: true,
    tickEvery: 1,
    tickLabel: (b, i) => (i === 0 || (i + 1) % 5 === 0 ? String(i + 1) : null),
    aria: 'Proteins by number of top-level pathways, log scale',
    describe: (i) => tipLines(`${exact(breadth[i])} proteins`, `in ${i + 1} top-level pathway${i ? 's' : ''}`),
  });
  attachData('#breadth-chart', 'reactome-breadth', ['Top-level pathways', 'Proteins'], breadth.map((n, i) => [i + 1, n]));
  const sizes = trimmed(RX.sizes.counts);
  const sizeBins = binsOf(RX.sizes.edges).slice(0, sizes.length);
  histo({
    mount: $('#pathway-size-chart'),
    bins: sizeBins,
    series: [{ key: 'pathways', label: 'Pathways', color: 'var(--data)', values: sizes }],
    tickEvery: 5,
    tickLabel: (b) => shortBases(b.lo),
    aria: 'Reactome pathways by number of proteins',
    describe: (i) => tipLines(`${exact(sizes[i])} pathways`, wholeRange(sizeBins[i], 'proteins')),
  });
  attachData('#pathway-size-chart', 'reactome-sizes', ['Proteins from', 'Proteins to', 'Pathways'], sizeBins.map((b, i) => [Math.ceil(b.lo), Math.ceil(b.hi) - 1, sizes[i]]));
  dotplot({
    mount: $('#largest-pathway-chart'),
    rows: RX.largest.map((p, i) => ({ key: p.id, label: p.name, value: p.proteins, hi: i === 0 })),
    scale: 'linear',
    labelShare: 0.5,
    labelMax: 340,
    aria: 'The largest Reactome pathways below the top level',
    href: (r) => reactomeLink(r.key),
    linkText: () => 'Open in Reactome',
    describe: (r) => tipLines(`${exact(r.value)} proteins`, r.label, r.key),
  });
  attachData('#largest-pathway-chart', 'reactome-largest', ['Pathway', 'Reactome ID', 'Proteins'], RX.largest.map((p) => [p.name, p.id, p.proteins]));
};

const ASPECT_NAMES = [['P', 'Biological process'], ['F', 'Molecular function'], ['C', 'Cellular component']];
const GO_GROUPS = [['experimental', 'Experimental', 'var(--s1)'], ['phylogenetic', 'Phylogenetic', 'var(--s2)'], ['computational', 'Computational', 'var(--s3)'], ['author', 'Author or curator', 'var(--s4)'], ['electronic', 'Electronic', 'var(--data-2)']];
const quickGo = (id) => `https://www.ebi.ac.uk/QuickGO/term/${id}`;

const renderOntology = () => {
  const ns = GO.namespaces;
  const others = Object.entries(GO.types).filter(([t]) => t !== 'protein');
  const rnas = sumOf(others.filter(([t]) => /RNA/.test(t)).map(([, n]) => n));
  $('#ontology-dek').textContent = `The ontology (release ${GO.release}) holds ${exact(ns.biological_process)} biological process, ${exact(ns.molecular_function)} molecular function and ${exact(ns.cellular_component)} cellular component terms, plus ${exact(GO.obsolete)} obsolete ones. Besides proteins, the human file annotates ${exact(rnas)} RNAs and ${exact(GO.types.protein_complex ?? 0)} protein complexes. ${exact(GO.not)} annotations say a product does NOT have a function. The most used evidence code is ${GO.codes[0].code}, ${exact(GO.codes[0].n)} times.`;
  hbars({
    mount: $('#go-evidence-chart'),
    rows: ASPECT_NAMES.map(([key, label]) => ({ key, label, values: GO.byAspect[key] })),
    segments: GO_GROUPS.map(([key, label, color]) => ({ key, label, color })),
    labelShare: 0.2,
    aria: 'GO annotations by aspect and evidence',
    describe: (r) => {
      const total = sumOf(Object.values(r.values));
      return tipLines(`${exact(total)} annotations`, r.label, ...GO_GROUPS.map(([k, label]) => `${label}: ${exact(r.values[k])}, ${pct(r.values[k] / total)}`));
    },
  });
  attachData('#go-evidence-chart', 'go-evidence', ['Aspect', ...GO_GROUPS.map(([, label]) => label)], ASPECT_NAMES.map(([key, label]) => [label, ...GO_GROUPS.map(([k]) => GO.byAspect[key][k])]));
  const counts = trimmed(GO.perProtein.counts);
  const bins = binsOf(GO.perProtein.edges).slice(0, counts.length);
  histo({
    mount: $('#go-per-protein-chart'),
    bins,
    series: [{ key: 'proteins', label: 'Proteins', color: 'var(--data)', values: counts }],
    tickEvery: 5,
    tickLabel: (b) => shortBases(b.lo),
    aria: 'Proteins by number of distinct GO terms',
    describe: (i) => tipLines(`${exact(counts[i])} proteins`, wholeRange(bins[i], 'GO terms')),
  });
  attachData('#go-per-protein-chart', 'go-terms-per-protein', ['Terms from', 'Terms to', 'Proteins'], bins.map((b, i) => [Math.ceil(b.lo), Math.ceil(b.hi) - 1, counts[i]]));
  dotplot({
    mount: $('#go-protein-chart'),
    rows: GO.topProteins.map((p, i) => ({ key: p.acc, label: p.symbol, value: p.terms, hi: i === 0 })),
    scale: 'linear',
    aria: 'Proteins with the most GO terms',
    href: (r) => `https://www.ebi.ac.uk/QuickGO/annotations?geneProductId=${r.key}`,
    linkText: () => 'Open in QuickGO',
    describe: (r) => tipLines(`${exact(r.value)} GO terms`, `${r.label}, UniProt ${r.key}`),
  });
  attachData('#go-protein-chart', 'go-top-proteins', ['Gene', 'UniProt', 'Distinct GO terms'], GO.topProteins.map((p) => [p.symbol, p.acc, p.terms]));
  for (const [sel, aspect, what] of [['#go-bp-chart', 'P', 'biological process'], ['#go-mf-chart', 'F', 'molecular function'], ['#go-cc-chart', 'C', 'cellular component']]) {
    const list = GO.topTerms[aspect];
    dotplot({
      mount: $(sel),
      rows: list.map((t, i) => ({ key: t.id, label: t.name ?? t.id, value: t.n, hi: i === 0 })),
      scale: 'linear',
      rowH: 22,
      aria: `The most annotated ${what} terms`,
      href: (r) => quickGo(r.key),
      linkText: () => 'Open in QuickGO',
      describe: (r) => tipLines(`${exact(r.value)} proteins`, r.label, r.key),
    });
    attachData(sel, `go-top-${what.replace(/ /g, '-')}`, ['GO term', 'Name', 'Proteins'], list.map((t) => [t.id, t.name, t.n]));
  }
  dotplot({
    mount: $('#go-groups-chart'),
    rows: GO.assignedBy.map((g, i) => ({ key: g.group, label: g.group, value: g.n, hi: i === 0 })),
    labelShare: 0.36,
    aria: 'GO annotations by contributing group',
    describe: (r) => tipLines(`${exact(r.value)} annotations`, r.label, pct(r.value / GO.annotations)),
  });
  attachData('#go-groups-chart', 'go-contributors', ['Group', 'Annotations'], GO.assignedBy.map((g) => [g.group, g.n]));
};

const REG = snapshot.regulation;
const TR = snapshot.trnas;
const ENC = snapshot.encode;
const span = (counts) => [counts.findIndex((n) => n > 0), counts.findLastIndex((n) => n > 0) + 1];
const sequencedMb = (chrom) => {
  const c = S.chromosomes.find((x) => x.chrom === chrom);
  return c ? acgtOf(c) / 1e6 : null;
};

const CCRE_CLASSES = [
  ['PLS', 'Promoter-like (PLS)', 'var(--s1)'], ['pELS', 'Proximal enhancer-like (pELS)', 'var(--s5)'], ['dELS', 'Distal enhancer-like (dELS)', 'var(--s3)'],
  ['CA-H3K4me3', 'Accessible with H3K4me3', 'var(--s2)'], ['CA-CTCF', 'Accessible with CTCF', 'var(--s4)'], ['CA-TF', 'Accessible with TF', 'var(--data)'],
  ['CA', 'Accessible only', 'var(--data-2)'], ['TF', 'TF only', 'var(--ink-2)'],
];

const renderCcres = () => {
  const { classes, byChrom, widths } = REG.ccre;
  const known = CCRE_CLASSES.filter(([key]) => classes[key]);
  const extra = Object.keys(classes).filter((k) => !CCRE_CLASSES.some(([key]) => key === k)).map((k) => [k, k, 'var(--muted)']);
  const list = [...known, ...extra];
  const bp = sumOf(Object.values(classes).map((c) => c.bp));
  $('#ccres-dek').textContent = `Together they span ${bases(bp)}. Every element is ${widths.summary.min} to ${widths.summary.max} bp wide, median ${exact(widths.summary.median)}. ${exact(classes.PLS?.count ?? 0)} are promoter-like, and ${exact(classes['CA-CTCF']?.count ?? 0)} are accessible sites whose only other mark is CTCF binding.`;
  dotplot({
    mount: $('#ccre-class-chart'),
    rows: list.map(([key, label], i) => ({ key, label, value: classes[key].count, hi: i === 2 })),
    labelShare: 0.48,
    aria: 'Candidate cis-regulatory elements by class',
    describe: (r) => {
      const c = classes[r.key];
      return tipLines(`${exact(c.count)} cCREs`, r.label, `${bases(c.bp)}, ${pct(c.bp / SN.acgt)} of GRCh38`, `mean width ${exact(c.bp / c.count)} bp`);
    },
  });
  attachData('#ccre-class-chart', 'ccre-classes', ['Class', 'cCREs', 'Bases'], list.map(([key]) => [key, classes[key].count, classes[key].bp]));
  const [a, b] = span(widths.counts);
  const wBins = binsOf(widths.edges).slice(a, b);
  const wCounts = widths.counts.slice(a, b);
  histo({
    mount: $('#ccre-width-chart'),
    bins: wBins,
    series: [{ key: 'ccres', label: 'cCREs', color: 'var(--data)', values: wCounts }],
    tickEvery: 2,
    tickLabel: (bin) => String(bin.lo),
    aria: 'cCREs by width',
    describe: (i) => {
      const top = Math.min(wBins[i].hi - 1, widths.summary.max);
      return tipLines(`${exact(wCounts[i])} cCREs`, top === wBins[i].lo ? `${top} bp wide` : `${wBins[i].lo} to ${top} bp wide`);
    },
  });
  attachData('#ccre-width-chart', 'ccre-widths', ['Width from', 'Width to', 'cCREs'], wBins.map((bin, i) => [bin.lo, Math.min(bin.hi - 1, widths.summary.max), wCounts[i]]));
  const rows = NUCLEAR.filter((c) => byChrom[c.chrom] && sequencedMb(c.chrom)).map((c) => ({
    key: c.chrom, label: c.chrom, values: Object.fromEntries(list.map(([k]) => [k, (byChrom[c.chrom][k] ?? 0) / sequencedMb(c.chrom)])),
  }));
  hbars({
    mount: $('#ccre-density-chart'),
    rows,
    segments: list.map(([key, label, color]) => ({ key, label, color })),
    labelShare: 0.05,
    rowH: 18,
    aria: 'cCREs per sequenced megabase by chromosome and class',
    describe: (r) => tipLines(`${exact(sumOf(Object.values(r.values)))} cCREs per Mb`, `chromosome ${r.key}`, ...list.map(([k, label]) => `${label}: ${r.values[k].toFixed(1)}`)),
  });
  attachData('#ccre-density-chart', 'ccre-density', ['Chromosome', ...list.map(([k]) => `${k} per Mb`)], rows.map((r) => [r.key, ...list.map(([k]) => Number(r.values[k].toFixed(2)))]));
};

const renderCpgIslands = () => {
  const cpg = REG.cpg;
  const p = cpg.promoters;
  const cats = Object.entries(p.bySpecificity).sort((x, y) => y[1].inIsland / y[1].genes - x[1].inIsland / x[1].genes);
  const rows = [['all', 'All MANE Select', { genes: p.tss, inIsland: p.inIsland }], ...cats.map(([k, v]) => [k, k, v])];
  const low = p.bySpecificity['Low tissue specificity'];
  const enriched = p.bySpecificity['Tissue enriched'];
  $('#cpgislands-dek').textContent = `Islands have a median length of ${exact(cpg.lengths.summary.median)} bp; the longest runs ${exact(cpg.lengths.summary.max)}.${low && enriched ? ` ${pct(low.inIsland / low.genes)} of genes with low tissue specificity start in an island, against ${pct(enriched.inIsland / enriched.genes)} of tissue-enriched genes.` : ''} Chromosome 19 carries ${(cpg.byChrom['19'] / sequencedMb('19')).toFixed(1)} islands per sequenced megabase.`;
  hbars({
    mount: $('#cpg-promoter-chart'),
    rows: rows.map(([key, label, v]) => ({
      key, label, v, values: { in: (v.inIsland / v.genes) * 100, out: 100 - (v.inIsland / v.genes) * 100 },
    })),
    segments: [{ key: 'in', label: 'Starts in a CpG island', color: 'var(--s2)' }, { key: 'out', label: 'Starts elsewhere', color: 'var(--data-2)' }],
    labelShare: 0.3,
    max: 100,
    R: 96,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => `${exact(r.v.genes)} genes`,
    aria: 'MANE Select transcripts starting in CpG islands, by tissue specificity',
    describe: (r) => tipLines(`${pct(r.v.inIsland / r.v.genes)} start in an island`, r.label, `${exact(r.v.inIsland)} of ${exact(r.v.genes)} transcripts`),
  });
  attachData('#cpg-promoter-chart', 'cpg-promoters', ['Group', 'Transcripts', 'Starting in a CpG island'], rows.map(([, label, v]) => [label, v.genes, v.inIsland]));
  const [a, b] = span(cpg.lengths.counts);
  const lBins = binsOf(cpg.lengths.edges).slice(a, b);
  const lCounts = cpg.lengths.counts.slice(a, b);
  histo({
    mount: $('#cpg-length-chart'),
    bins: lBins,
    series: [{ key: 'islands', label: 'CpG islands', color: 'var(--data)', values: lCounts }],
    tickEvery: 1,
    tickLabel: (bin, i) => (i === 0 || (i + a) % 10 === 0 ? shortBases(bin.lo) : null),
    aria: 'CpG islands by length',
    describe: (i) => tipLines(`${exact(lCounts[i])} islands`, `${exact(lBins[i].lo)} to ${exact(lBins[i].hi)} bp`),
  });
  attachData('#cpg-length-chart', 'cpg-island-lengths', ['Length from', 'Length to', 'Islands'], lBins.map((bin, i) => [Math.round(bin.lo), Math.round(bin.hi), lCounts[i]]));
  const dens = NUCLEAR.filter((c) => cpg.byChrom[c.chrom] && sequencedMb(c.chrom)).map((c) => ({ chrom: c.chrom, n: cpg.byChrom[c.chrom], perMb: cpg.byChrom[c.chrom] / sequencedMb(c.chrom) }));
  const top = dens.reduce((x, y) => (y.perMb > x.perMb ? y : x));
  dotplot({
    mount: $('#cpg-density-chart'),
    rows: dens.map((d) => ({
      key: d.chrom, label: d.chrom, value: d.perMb, hi: d === top,
    })),
    scale: 'linear',
    rowH: 18,
    labelShare: 0.05,
    format: (v) => v.toFixed(1),
    aria: 'CpG islands per sequenced megabase by chromosome',
    describe: (r) => {
      const d = dens.find((x) => x.chrom === r.key);
      return tipLines(`${d.perMb.toFixed(1)} islands per Mb`, `chromosome ${d.chrom}`, `${exact(d.n)} islands`);
    },
  });
  attachData('#cpg-density-chart', 'cpg-island-density', ['Chromosome', 'Islands', 'Per sequenced Mb'], dens.map((d) => [d.chrom, d.n, Number(d.perMb.toFixed(2))]));
  const ratioBins = binsOf(cpg.obsExp.edges);
  histo({
    mount: $('#cpg-oe-chart'),
    bins: ratioBins,
    series: [{ key: 'islands', label: 'CpG islands', color: 'var(--data)', values: cpg.obsExp.counts }],
    tickEvery: 4,
    tickLabel: (bin) => bin.lo.toFixed(1),
    aria: 'CpG islands by observed over expected CpG',
    describe: (i) => tipLines(`${exact(cpg.obsExp.counts[i])} islands`, i === ratioBins.length - 1 ? `observed over expected ${ratioBins[i].lo.toFixed(2)} or more` : `observed over expected ${ratioBins[i].lo.toFixed(2)} to ${ratioBins[i].hi.toFixed(2)}`),
  });
  attachData('#cpg-oe-chart', 'cpg-island-obs-exp', ['Ratio from', 'Ratio to', 'Islands'], ratioBins.map((bin, i) => [bin.lo, bin.hi, cpg.obsExp.counts[i]]));
  const gcBins = binsOf(cpg.gc.edges);
  histo({
    mount: $('#cpg-gc-chart'),
    bins: gcBins,
    series: [{ key: 'islands', label: 'CpG islands', color: 'var(--data)', values: cpg.gc.counts }],
    tickEvery: 5,
    tickLabel: (bin) => `${bin.lo}%`,
    aria: 'CpG islands by GC content',
    describe: (i) => tipLines(`${exact(cpg.gc.counts[i])} islands`, i === gcBins.length - 1 ? `${gcBins[i].lo}% GC or more` : `${gcBins[i].lo}% to ${gcBins[i].hi}% GC`),
  });
  attachData('#cpg-gc-chart', 'cpg-island-gc', ['GC from', 'GC to', 'Islands'], gcBins.map((bin, i) => [bin.lo, bin.hi, cpg.gc.counts[i]]));
};

const ISOTYPE_NAMES = { ...Object.fromEntries(Object.entries(AA3).filter(([k]) => k !== '*').map(([one, three]) => [three, AMINO_ACIDS[one]])), iMet: 'Initiator methionine', SeC: 'Selenocysteine' };

const renderTrna = () => {
  const pairs = trnaPairs(MN.codons, TR.anticodons);
  const fit = linearFit(pairs.map((x) => x.perK), pairs.map((x) => x.copies));
  const r = pearson(pairs.map((x) => x.perK), pairs.map((x) => x.copies));
  const zero = pairs.filter((x) => !x.copies);
  const [first, second] = Object.entries(TR.byChrom).sort((x, y) => y[1] - x[1]);
  $('#trna-dek').textContent = `Across the 61 sense codons, codon use and exact-match tRNA gene copies correlate with Pearson r = ${r.toFixed(2)}. The codons without an exact match are ${zero.map((x) => x.codon).join(', ')}. Chromosome ${first[0]} carries ${exact(first[1])} of the genes and chromosome ${second[0]} ${exact(second[1])}.`;
  scatter({
    mount: $('#trna-codon-chart'),
    points: pairs.map((x) => ({ ...x, key: x.codon, x: Number(x.perK.toFixed(2)), y: x.copies, hi: !x.copies })),
    xFormat: (v) => String(v),
    yFormat: (v) => String(v),
    xLabel: 'Codon use per thousand',
    yLabel: 'tRNA genes with the pairing anticodon',
    trend: fit,
    aria: 'tRNA gene copies against codon usage',
    describe: (pt) => tipLines(`${pt.copies} tRNA genes`, `${pt.codon}, ${pt.aa}; anticodon ${pt.anticodon}`, `${pt.perK.toFixed(1)} per thousand codons`),
  });
  attachData('#trna-codon-chart', 'trna-codons', ['Codon', 'Amino acid', 'Anticodon', 'Per thousand codons', 'tRNA genes'], pairs.map((x) => [x.codon, x.aa, x.anticodon, Number(x.perK.toFixed(3)), x.copies]));
  const iso = Object.entries(TR.isotypes).sort((x, y) => y[1] - x[1]);
  dotplot({
    mount: $('#trna-isotype-chart'),
    rows: iso.map(([k, n], i) => ({ key: k, label: `${ISOTYPE_NAMES[k] ?? k} (${k})`, value: n, hi: i === 0 })),
    scale: 'linear',
    rowH: 20,
    labelShare: 0.46,
    aria: 'tRNA genes by amino acid',
    describe: (row) => tipLines(`${exact(row.value)} tRNA genes`, row.label, Object.entries(TR.anticodons).filter(([k]) => k.startsWith(`${row.key}-`)).map(([k, n]) => `${k.split('-')[1]} ${n}`).join(', ')),
  });
  attachData('#trna-isotype-chart', 'trna-isotypes', ['Isotype', 'tRNA genes'], iso);
  const chroms = NUCLEAR.map((c) => c.chrom).filter((c) => TR.byChrom[c]);
  dotplot({
    mount: $('#trna-chrom-chart'),
    rows: chroms.map((c) => ({
      key: c, label: c, value: TR.byChrom[c], hi: c === first[0],
    })),
    scale: 'linear',
    rowH: 20,
    labelShare: 0.05,
    aria: 'tRNA genes by chromosome',
    describe: (row) => tipLines(`${exact(row.value)} tRNA genes`, `chromosome ${row.key}`),
  });
  attachData('#trna-chrom-chart', 'trna-chromosomes', ['Chromosome', 'tRNA genes'], chroms.map((c) => [c, TR.byChrom[c]]));
};

const renderEncode = () => {
  const share = (list, name) => (list.find((x) => x.name === name)?.n ?? 0) / ENC.total;
  $('#encode-dek').textContent = `${pct(share(ENC.biosamples, 'cell line'))} of the experiments use cell lines and ${pct(share(ENC.biosamples, 'tissue'))} tissues. The fourth phase, ENCODE4, released ${exact(ENC.phases.find((x) => x.name === 'ENCODE4')?.n ?? 0)} of them. Histone and transcription factor ChIP-seq together account for ${exact(sumOf(ENC.assays.filter((x) => /^(TF|Histone) ChIP-seq$/.test(x.name)).map((x) => x.n)))}.`;
  const ranked = (sel, name, list, aria, scale, labelShare) => {
    dotplot({
      mount: $(sel),
      rows: list.map((x, i) => ({ key: x.name, label: x.name, value: x.n, hi: i === 0 })),
      scale,
      labelShare,
      rowH: 20,
      aria,
      describe: (row) => tipLines(`${exact(row.value)} experiments`, row.label, pct(row.value / ENC.total)),
    });
    attachData(sel, name, ['Group', 'Experiments'], list.map((x) => [x.name, x.n]));
  };
  ranked('#encode-assay-chart', 'encode-assays', ENC.assays.slice(0, 25), 'ENCODE experiments by assay', 'log', 0.3);
  ranked('#encode-phase-chart', 'encode-phases', ENC.phases, 'ENCODE experiments by phase', 'linear', 0.44);
  ranked('#encode-sample-chart', 'encode-biosamples', ENC.biosamples, 'ENCODE experiments by biosample type', 'linear', 0.44);
  ranked('#encode-target-chart', 'encode-targets', ENC.targets, 'ENCODE experiments by target class', 'log', 0.44);
};

const ST = snapshot.structure;
const MT_KINDS = [
  { key: 'protein', label: 'Protein-coding', color: 'var(--s2)' },
  { key: 'rRNA', label: 'rRNA', color: 'var(--s4)' },
  { key: 'tRNA', label: 'tRNA', color: 'var(--s3)' },
  { key: 'dloop', label: 'D-loop', color: 'var(--data)' },
];

const renderMito = () => {
  const { genes, dloop, variants: v } = ST.mito;
  const mt = S.chromosomes.find((c) => c.chrom === 'MT');
  const length = mt.length;
  const minus = genes.filter((g) => g.strand === '-');
  const features = (strand) => genes.filter((g) => g.strand === strand).map((g) => ({ ...g, key: g.name, label: g.kind === 'tRNA' ? null : g.name }));
  const loop = dloop.map(([start, end]) => ({
    key: 'D-loop', name: 'D-loop', start, end, kind: 'dloop', strand: '+', product: 'control region', label: end - start > 400 ? 'D-loop' : null,
  }));
  $('#mito-dek').textContent = `The sequence is ${gcPct((mt.C + mt.G) / acgtOf(mt))} G and C. ${minus.length} genes, ${minus.map((g) => g.name).join(', ')}, lie on the minus strand of the reference; the rest on the plus strand. gnomAD v3.1 calls ${exact(v.variants)} variants across ${exact(v.samples)} samples; ${exact(v.homoplasmic)} are seen homoplasmic in someone and ${exact(v.heteroplasmicOnly)} only as heteroplasmy.`;
  genomeMap({
    mount: $('#mito-map-chart'),
    length,
    span: v.bins.length * 100,
    lanes: [{ key: '+', label: '+', features: [...features('+'), ...loop] }, { key: '-', label: '−', features: features('-') }],
    kinds: MT_KINDS,
    aria: 'Map of the mitochondrial genome',
    describe: (f) => tipLines(f.name, f.product ?? f.kind, `positions ${exact(f.start)} to ${exact(f.end)}, ${exact(f.end - f.start + 1)} bp`, `${f.strand === '-' ? 'minus' : 'plus'} strand of the reference`, f.kind === 'dloop' ? null : { href: `https://www.ncbi.nlm.nih.gov/gene/?term=MT-${f.name}`, text: 'Search NCBI Gene' }),
  });
  attachData('#mito-map-chart', 'mitochondrial-genes', ['Gene', 'Kind', 'Product', 'Start', 'End', 'Strand'], [...genes.map((g) => [g.name, g.kind, g.product, g.start, g.end, g.strand]), ...loop.map((l) => [l.name, 'control region', l.product, l.start, l.end, '+'])]);
  const totals = v.bins.map((b) => b.hom + b.het);
  histo({
    mount: $('#mito-variant-chart'),
    bins: v.bins.map((b, i) => ({ lo: i * 100 })),
    series: [{ key: 'variants', label: 'Variants', color: 'var(--data)', values: totals }],
    height: 140,
    tickEvery: 20,
    tickLabel: (b) => shortBases(b.lo),
    aria: 'gnomAD mitochondrial variants per 100 bp',
    describe: (i) => tipLines(`${exact(totals[i])} variants`, `positions ${exact(i * 100 + 1)} to ${exact(Math.min((i + 1) * 100, length))}`, `${exact(v.bins[i].hom)} seen homoplasmic`, `${exact(v.bins[i].het)} only heteroplasmic`),
  });
  attachData('#mito-variant-chart', 'mitochondrial-variants-per-100bp', ['From', 'To', 'Homoplasmic in someone', 'Only heteroplasmic'], v.bins.map((b, i) => [i * 100 + 1, Math.min((i + 1) * 100, length), b.hom, b.het]));
  const name = (t) => `m.${t.pos}${t.ref}>${t.alt}`;
  dotplot({
    mount: $('#mito-top-chart'),
    rows: v.top.map((t, i) => ({ key: name(t), label: name(t), value: t.afHom * 100, hi: i === 0 })),
    scale: 'linear',
    rowH: 20,
    format: (x) => `${x.toFixed(1)}%`,
    aria: 'Most frequent homoplasmic mitochondrial variants',
    describe: (r) => {
      const t = v.top.find((x) => name(x) === r.key);
      return tipLines(`${(t.afHom * 100).toFixed(2)}% homoplasmic`, r.label, `${(t.afHet * 100).toFixed(3)}% heteroplasmic`);
    },
  });
  attachData('#mito-top-chart', 'mitochondrial-top-variants', ['Variant', 'Homoplasmic share', 'Heteroplasmic share'], v.top.map((t) => [name(t), t.afHom, t.afHet]));
  const hBins = binsOf(v.heteroplasmy.edges);
  histo({
    mount: $('#mito-het-chart'),
    bins: hBins,
    series: [{ key: 'variants', label: 'Variants', color: 'var(--data)', values: v.heteroplasmy.counts }],
    yLog: true,
    tickEvery: 4,
    tickLabel: (b) => `${Math.round(b.lo * 100)}%`,
    aria: 'Variants by highest observed heteroplasmy',
    describe: (i) => tipLines(`${exact(v.heteroplasmy.counts[i])} variants`, `highest heteroplasmy ${Math.round(hBins[i].lo * 100)}% to ${Math.round(hBins[i].hi * 100)}%`),
  });
  attachData('#mito-het-chart', 'mitochondrial-heteroplasmy', ['Heteroplasmy from', 'Heteroplasmy to', 'Variants'], hBins.map((b, i) => [b.lo, b.hi, v.heteroplasmy.counts[i]]));
};

const SAT_GROUPS = [
  ['hor', 'Active α-satellite HOR', ['hor'], 'var(--s1)'],
  ['alpha', 'Other α-satellite', ['dhor', 'mon'], 'var(--s5)'],
  ['hsat1', 'HSat1A and HSat1B', ['hsat1A', 'hsat1B'], 'var(--s4)'],
  ['hsat2', 'HSat2', ['hsat2'], 'var(--s2)'],
  ['hsat3', 'HSat3', ['hsat3'], 'var(--s3)'],
  ['bsat', 'β-satellite', ['bsat'], 'var(--data)'],
  ['gsat', 'γ-satellite', ['gsat'], 'var(--ink-2)'],
  ['rDNA', 'rDNA arrays', ['rDNA'], 'var(--muted)'],
  ['other', 'Other satellites', [], 'var(--data-2)'],
];
const groupOf = (cls) => SAT_GROUPS.find(([, , members]) => members.includes(cls))?.[0] ?? 'other';
const satGroups = (byClass) => {
  const out = Object.fromEntries(SAT_GROUPS.map(([key]) => [key, 0]));
  for (const [cls, bp] of Object.entries(byClass ?? {})) out[groupOf(cls)] += bp;
  return out;
};

const renderCensat = () => {
  const { byChrom, classes, total, transition } = ST.censat;
  const t2tLength = A.t2t.chromosomes.filter((c) => c.chrom !== 'MT').reduce((n, c) => n + c.length, 0);
  const all = satGroups(classes);
  const hor = [...ST.censat.arrays.hor].sort((a, b) => b.bp - a.bp);
  $('#censat-dek').textContent = `Satellites make up ${pct(total / t2tLength)} of T2T-CHM13, and the transition regions between arrays another ${bases(transition)}. Active HOR arrays run from ${bases(hor.at(-1).bp)} on chromosome ${hor.at(-1).chrom} to ${bases(hor[0].bp)} on chromosome ${hor[0].chrom}. HSat3 is the largest family after α-satellite HORs, ${bases(all.hsat3)}.`;
  const chroms = A.t2t.chromosomes.filter((c) => c.chrom !== 'MT');
  hbars({
    mount: $('#censat-chart'),
    rows: chroms.map((c) => ({ key: c.chrom, label: c.chrom, values: Object.fromEntries(Object.entries(satGroups(byChrom[c.chrom])).map(([k, bp]) => [k, bp / 1e6])) })),
    segments: SAT_GROUPS.map(([key, label, , color]) => ({ key, label, color })),
    labelShare: 0.05,
    rowH: 18,
    format: (v) => `${v.toFixed(1)} Mb`,
    tickFormat: (t) => `${t} Mb`,
    aria: 'Satellite arrays per T2T-CHM13 chromosome by family',
    describe: (r) => tipLines(`${sumOf(Object.values(r.values)).toFixed(2)} Mb of satellites`, `chromosome ${r.key}`, ...SAT_GROUPS.filter(([key]) => r.values[key] > 0).map(([key, label]) => `${label}: ${r.values[key].toFixed(2)} Mb`)),
  });
  attachData('#censat-chart', 'satellites-by-chromosome', ['Chromosome', ...SAT_GROUPS.map(([key]) => `${key} bp`)], chroms.map((c) => {
    const g = satGroups(byChrom[c.chrom]);
    return [c.chrom, ...SAT_GROUPS.map(([key]) => g[key])];
  }));
  const groups = SAT_GROUPS.filter(([key]) => all[key] > 0).sort((a, b) => all[b[0]] - all[a[0]]);
  dotplot({
    mount: $('#censat-class-chart'),
    rows: groups.map(([key, label], i) => ({ key, label, value: all[key], hi: i === 0 })),
    format: bases,
    tickFormat: (t) => bases(t).replace(' ', ''),
    labelShare: 0.46,
    rowH: 20,
    aria: 'Satellite bases by family',
    describe: (r) => tipLines(bases(r.value), r.label, `${pct(r.value / t2tLength)} of T2T-CHM13`),
  });
  attachData('#censat-class-chart', 'satellite-families', ['Family', 'Bases'], groups.map(([key, label]) => [label, all[key]]));
  dumbbell({
    mount: $('#telomere-chart'),
    rows: ST.telomeres.map((t) => ({
      key: t.chrom, label: t.chrom, a: t.p, b: t.q,
    })),
    years: ['p arm', 'q arm'],
    format: (v) => `${(v / 1000).toFixed(1)} kb`,
    valueText: (r) => `${(r.a / 1000).toFixed(1)} and ${(r.b / 1000).toFixed(1)} kb`,
    tickFormat: (t) => `${t / 1000} kb`,
    aria: 'Telomere lengths at both ends of each chromosome',
    describe: (r) => tipLines(`p ${exact(r.a)} bp, q ${exact(r.b)} bp`, `chromosome ${r.key} telomeres`),
  });
  attachData('#telomere-chart', 't2t-telomeres', ['Chromosome', 'p-arm telomere bp', 'q-arm telomere bp'], ST.telomeres.map((t) => [t.chrom, t.p, t.q]));
};

const BIOTYPE_GROUPS = [
  ['protein', 'Protein-coding', 'var(--s1)'],
  ['lnc', 'lncRNA', 'var(--s2)'],
  ['pseudo', 'Pseudogenes', 'var(--s3)'],
  ['small', 'Small RNAs', 'var(--s4)'],
  ['other', 'Other', 'var(--data-2)'],
];
const SMALL_RNAS = new Set(['miRNA', 'snRNA', 'snoRNA', 'misc_RNA', 'scaRNA', 'rRNA', 'sRNA', 'vault_RNA', 'ribozyme', 'scRNA', 'Mt_tRNA', 'Mt_rRNA']);
const biotypeGroup = (type) => {
  if (type === 'protein_coding') return 'protein';
  if (type === 'lncRNA') return 'lnc';
  if (/pseudogene$/.test(type)) return 'pseudo';
  return SMALL_RNAS.has(type) ? 'small' : 'other';
};

const renderSexchrom = () => {
  const by = G.genes.byChrom;
  const grouped = (types) => {
    const out = Object.fromEntries(BIOTYPE_GROUPS.map(([key]) => [key, 0]));
    for (const [type, n] of Object.entries(types ?? {})) out[biotypeGroup(type)] += n;
    return out;
  };
  const autosomes = Array.from({ length: 22 }, (_, i) => grouped(by[String(i + 1)]));
  const mean = Object.fromEntries(BIOTYPE_GROUPS.map(([key]) => [key, sumOf(autosomes.map((a) => a[key])) / 22]));
  const rows = [
    { key: 'X', label: 'X', values: grouped(by.X) },
    { key: 'Y', label: 'Y', values: grouped(by.Y) },
    { key: 'autosome', label: 'Autosome mean', values: mean },
  ];
  const xy = NUCLEAR.filter((c) => c.chrom === 'X' || c.chrom === 'Y');
  const ySat = satGroups(ST.censat.byChrom.Y);
  const y2 = A.t2t.chromosomes.find((c) => c.chrom === 'Y').length;
  const par1 = PARS[0].x[1] - PARS[0].x[0] + 1;
  const par2 = PARS[1].x[1] - PARS[1].x[0] + 1;
  $('#sexchrom-dek').textContent = `PAR1 spans ${bases(par1)} at the tips of Xp and Yp and PAR2 ${bases(par2)} at the tips of Xq and Yq. Y is ${bases(xy[1].length)} long in GRCh38 and ${bases(y2)} in T2T-CHM13, which resolves its heterochromatin as ${bases(ySat.hsat3)} of HSat3 and ${bases(ySat.hsat1)} of HSat1. X carries ${exact(rows[0].values.protein)} protein-coding genes against an autosome mean of ${exact(mean.protein)}.`;
  ideogram({
    mount: $('#par-chart'),
    chromosomes: xy,
    bands: A.bands.filter((b) => b.chrom === 'X' || b.chrom === 'Y'),
    marks: PARS.flatMap((p) => [{ chrom: 'X', start: p.x[0], end: p.x[1], label: p.name }, { chrom: 'Y', start: p.y[0], end: p.y[1], label: p.name }]),
    aria: 'X and Y with the pseudoautosomal regions',
    describe: (band) => tipLines(`${band.chrom}${band.name}`, coords(band.chrom, band.start + 1, band.end), STAINS[band.stain] ?? band.stain),
  });
  attachData('#par-chart', 'pseudoautosomal-regions', ['Region', 'X start', 'X end', 'Y start', 'Y end'], PARS.map((p) => [p.name, ...p.x, ...p.y]));
  hbars({
    mount: $('#xy-genes-chart'),
    rows,
    segments: BIOTYPE_GROUPS.map(([key, label, color]) => ({ key, label, color })),
    labelShare: 0.26,
    format: (v) => exact(v),
    aria: 'Genes on X, Y and an average autosome by type',
    describe: (r) => tipLines(`${exact(sumOf(Object.values(r.values)))} genes`, r.label, ...BIOTYPE_GROUPS.map(([key, label]) => `${label}: ${exact(r.values[key])}`)),
  });
  attachData('#xy-genes-chart', 'x-y-genes', ['Chromosome', ...BIOTYPE_GROUPS.map(([, label]) => label)], rows.map((r) => [r.label, ...BIOTYPE_GROUPS.map(([key]) => Math.round(r.values[key]))]));
  const satTotal = sumOf(Object.values(ySat));
  const parts = [
    ...SAT_GROUPS.filter(([key]) => ySat[key] > 0).sort((a, b) => ySat[b[0]] - ySat[a[0]]).map(([key, label, , color]) => ({ label, value: ySat[key], color })),
    { label: 'Not satellite', value: y2 - satTotal, color: 'var(--g0)' },
  ];
  waffle({
    mount: $('#y-satellite-chart'), parts, aria: 'T2T-CHM13 Y chromosome by satellite family', format: bases,
  });
  attachData('#y-satellite-chart', 'y-satellites', ['Family', 'Bases'], parts.map((p) => [p.label, p.value]));
};

const VR = snapshot.variation;
const CN = snapshot.constraint;
const SV_CLASSES = [['deletion', 'Deletions', 'var(--s1)'], ['duplication', 'Duplications', 'var(--s2)'], ['cnv', 'Copy number variants', 'var(--s3)'], ['insertion', 'Insertions', 'var(--s4)']];

const renderDbsnp = () => {
  const d = VR.dbsnp;
  const rows = NUCLEAR.filter((c) => d.byChrom[c.chrom] && sequencedMb(c.chrom)).map((c) => ({ chrom: c.chrom, n: d.byChrom[c.chrom], perKb: d.byChrom[c.chrom] / (sequencedMb(c.chrom) * 1000) }));
  const mt = sequencedMb('MT');
  $('#dbsnp-dek').textContent = `Chromosome Y has ${(rows.find((r) => r.chrom === 'Y')?.perKb ?? 0).toFixed(1)} records per sequenced kilobase and the mitochondrial genome ${mt ? (d.byChrom.MT / (mt * 1000)).toFixed(0) : '—'}. ${exact(d.other.PAR ?? 0)} records sit in the pseudoautosomal regions and ${exact(d.other['Not On'] ?? 0)} on no assembled sequence.`;
  const top = rows.reduce((a, b) => (b.perKb > a.perKb ? b : a));
  dotplot({
    mount: $('#dbsnp-density-chart'),
    rows: rows.map((r) => ({
      key: r.chrom, label: r.chrom, value: r.perKb, hi: r === top,
    })),
    scale: 'linear',
    rowH: 18,
    labelShare: 0.05,
    format: (v) => v.toFixed(1),
    aria: 'dbSNP records per sequenced kilobase by chromosome',
    describe: (r) => {
      const x = rows.find((y) => y.chrom === r.key);
      return tipLines(`${x.perKb.toFixed(1)} per kb`, `chromosome ${x.chrom}`, `${exact(x.n)} reference SNP records`);
    },
  });
  attachData('#dbsnp-density-chart', 'dbsnp-by-chromosome', ['Chromosome', 'RS records', 'Per sequenced kb'], rows.map((r) => [r.chrom, r.n, Number(r.perKb.toFixed(3))]));
  const status = Object.entries(d.status);
  dotplot({
    mount: $('#dbsnp-status-chart'),
    rows: status.map(([k, n], i) => ({ key: k, label: k, value: n, hi: i === 0 })),
    labelShare: 0.4,
    aria: 'dbSNP records by status',
    describe: (r) => tipLines(exact(r.value), r.label, pct(r.value / d.total)),
  });
  attachData('#dbsnp-status-chart', 'dbsnp-status', ['Status', 'Records'], status);
};

const CLINICAL = [['pathogenic', 'Pathogenic or likely'], ['benign', 'Benign or likely'], ['uncertain', 'Uncertain'], ['conflicting', 'Conflicting'], ['other', 'Other assertions']];

const renderSv = () => {
  const s = VR.sv;
  const total = sumOf(Object.values(s.classes));
  const sized = SV_CLASSES.filter(([k]) => s.sizes[k]);
  const edges = s.sizes[sized[0][0]].edges;
  const bins = binsOf(edges);
  const gaps = SV_CLASSES.map(([k, label]) => [label, NUCLEAR.filter((c) => s.byChrom[c.chrom] && !s.byChrom[c.chrom][k]).map((c) => c.chrom)]).filter(([, list]) => list.length);
  const listed = (list) => (list.length > 1 ? `${list.slice(0, -1).join(', ')} and ${list.at(-1)}` : list[0]);
  const gapNote = gaps.map(([label, list]) => ` The published ${label.toLowerCase()} set has no records on chromosome${list.length > 1 ? 's' : ''} ${listed(list)}.`).join('');
  $('#sv-dek').textContent = `Median lengths: ${sized.map(([k, label]) => `${label.toLowerCase()} ${bases(s.sizes[k].summary.median)}`).join(', ')}. ${exact(s.mobile.Alu.insertion)} SVs are typed as Alu insertions and ${exact(s.mobile.L1.insertion)} as L1 insertions. ${exact(s.clinical.pathogenic)} carry a pathogenic or likely pathogenic assertion, ${pct(s.clinical.pathogenic / total)} of the ${exact(total)}.${gapNote}`;
  histo({
    mount: $('#sv-size-chart'),
    bins,
    series: sized.map(([k, label, color]) => ({ key: k, label, color, values: s.sizes[k].counts })),
    yLog: true,
    tickEvery: 2,
    tickLabel: (b) => shortBases(b.lo),
    aria: 'Structural variants by length and class, log scale',
    describe: (i) => tipLines(`${bases(bins[i].lo)} to ${bases(bins[i].hi)}`, 'length', ...sized.map(([k, label]) => `${label}: ${exact(s.sizes[k].counts[i])}`)),
  });
  attachData('#sv-size-chart', 'sv-lengths', ['Length from', 'Length to', ...sized.map(([, label]) => label)], bins.map((b, i) => [Math.round(b.lo), Math.round(b.hi), ...sized.map(([k]) => s.sizes[k].counts[i])]));
  const rows = NUCLEAR.filter((c) => s.byChrom[c.chrom] && sequencedMb(c.chrom)).map((c) => ({
    key: c.chrom, label: c.chrom, values: Object.fromEntries(SV_CLASSES.map(([k]) => [k, (s.byChrom[c.chrom][k] ?? 0) / sequencedMb(c.chrom)])),
  }));
  hbars({
    mount: $('#sv-chrom-chart'),
    rows,
    segments: SV_CLASSES.map(([key, label, color]) => ({ key, label, color })),
    labelShare: 0.05,
    rowH: 18,
    aria: 'Structural variants per sequenced megabase by chromosome and class',
    describe: (r) => tipLines(`${exact(sumOf(Object.values(r.values)))} SVs per Mb`, `chromosome ${r.key}`, ...SV_CLASSES.map(([k, label]) => `${label}: ${r.values[k].toFixed(0)}`)),
  });
  attachData('#sv-chrom-chart', 'sv-by-chromosome', ['Chromosome', ...SV_CLASSES.map(([k]) => k)], rows.map((r) => [r.key, ...SV_CLASSES.map(([k]) => s.byChrom[r.key][k] ?? 0)]));
  const mobile = Object.entries(s.mobile);
  hbars({
    mount: $('#sv-mobile-chart'),
    rows: mobile.map(([name, v]) => ({ key: name, label: name, values: v })),
    segments: [{ key: 'insertion', label: 'Insertions', color: 'var(--s4)' }, { key: 'deletion', label: 'Deletions', color: 'var(--s1)' }],
    labelShare: 0.14,
    aria: 'Mobile element SVs by family',
    describe: (r) => tipLines(`${exact(r.values.insertion + r.values.deletion)} SVs`, r.label, `${exact(r.values.insertion)} insertions`, `${exact(r.values.deletion)} deletions`),
  });
  attachData('#sv-mobile-chart', 'sv-mobile-elements', ['Family', 'Insertions', 'Deletions'], mobile.map(([name, v]) => [name, v.insertion, v.deletion]));
  dotplot({
    mount: $('#sv-study-chart'),
    rows: s.studies.map((x, i) => ({ key: x.study, label: x.study, value: x.n, hi: i === 0 })),
    rowH: 20,
    labelShare: 0.4,
    aria: 'Studies contributing the most structural variants',
    describe: (r) => tipLines(`${exact(r.value)} SVs`, r.label, pct(r.value / total)),
  });
  attachData('#sv-study-chart', 'sv-studies', ['Study', 'SVs'], s.studies.map((x) => [x.study, x.n]));
  dotplot({
    mount: $('#sv-clinical-chart'),
    rows: CLINICAL.map(([k, label], i) => ({ key: k, label, value: s.clinical[k], hi: i === 0 })),
    rowH: 20,
    labelShare: 0.46,
    aria: 'Structural variants by clinical assertion',
    describe: (r) => tipLines(`${exact(r.value)} SVs`, r.label),
  });
  attachData('#sv-clinical-chart', 'sv-clinical', ['Assertion', 'SVs'], [...CLINICAL.map(([k, label]) => [label, s.clinical[k]]), ['None', s.clinical.none]]);
};

const renderConstraint = () => {
  const g = CN.gnomad;
  const { sufficient, other } = CN.byDosage;
  const bins = binsOf(g.loeuf.edges);
  const share = (counts, n) => counts.map((c) => (c / n) * 100);
  const hiShare = share(sufficient.counts, sufficient.n);
  const otherShare = share(other.counts, other.n);
  $('#constraint-dek').textContent = `${exact(g.genes)} genes have a MANE Select LOEUF; the quartiles are ${g.loeuf.summary.q1.toFixed(2)}, ${g.loeuf.summary.median.toFixed(2)} and ${g.loeuf.summary.q3.toFixed(2)}. ${exact(g.pli09)} have a pLI of 0.9 or more. The ${exact(sufficient.n)} genes ClinGen finds haploinsufficient have a median LOEUF of ${sufficient.median.toFixed(2)}, against ${other.median.toFixed(2)} for the rest.`;
  histo({
    mount: $('#loeuf-chart'),
    bins,
    series: [
      { key: 'hi', label: `Haploinsufficient (${exact(sufficient.n)} genes)`, color: 'var(--s1)', values: hiShare },
      { key: 'other', label: `All other genes (${exact(other.n)})`, color: 'var(--data)', values: otherShare },
    ],
    yFormat: (v) => `${v}%`,
    tickEvery: 4,
    tickLabel: (b) => b.lo.toFixed(1),
    aria: 'LOEUF of haploinsufficient genes and all other genes',
    describe: (i) => tipLines(`LOEUF ${bins[i].lo.toFixed(2)} to ${bins[i].hi.toFixed(2)}`, 'share of each group', `Haploinsufficient: ${hiShare[i].toFixed(1)}%, ${exact(sufficient.counts[i])} genes`, `Others: ${otherShare[i].toFixed(1)}%, ${exact(other.counts[i])} genes`),
  });
  attachData('#loeuf-chart', 'loeuf', ['LOEUF from', 'LOEUF to', 'Haploinsufficient genes', 'Other genes', 'All genes'], bins.map((b, i) => [b.lo, b.hi, sufficient.counts[i], other.counts[i], g.loeuf.counts[i]]));
  const zBins = binsOf(g.misZ.edges);
  histo({
    mount: $('#misz-chart'),
    bins: zBins,
    series: [{ key: 'genes', label: 'Genes', color: 'var(--data)', values: g.misZ.counts }],
    tickEvery: 4,
    tickLabel: (b) => signed(b.lo, String),
    aria: 'Genes by missense Z score',
    describe: (i) => tipLines(`${exact(g.misZ.counts[i])} genes`, `missense Z ${zBins[i].lo} to ${zBins[i].hi}`),
  });
  attachData('#misz-chart', 'missense-z', ['Z from', 'Z to', 'Genes'], zBins.map((b, i) => [b.lo, b.hi, g.misZ.counts[i]]));
  dotplot({
    mount: $('#constrained-chart'),
    rows: g.top.map((t, i) => ({ key: t.gene, label: t.gene, value: t.loeuf, hi: i === 0 })),
    scale: 'linear',
    rowH: 20,
    format: (v) => v.toFixed(3),
    tickFormat: (t) => String(t),
    aria: 'The most loss-of-function constrained genes',
    href: (r) => `https://gnomad.broadinstitute.org/gene/${r.key}?dataset=gnomad_r4`,
    linkText: () => 'Open in gnomAD',
    describe: (r) => {
      const t = g.top.find((x) => x.gene === r.key);
      return tipLines(`LOEUF ${t.loeuf.toFixed(3)}`, t.gene, `pLI ${t.pli.toFixed(2)}`, `${t.expected.toFixed(1)} loss-of-function variants expected`);
    },
  });
  attachData('#constrained-chart', 'most-constrained-genes', ['Gene', 'LOEUF', 'pLI', 'Expected LoF'], g.top.map((t) => [t.gene, t.loeuf, t.pli, t.expected]));
};

const DOSAGE_SCORES = [['3', 'Sufficient (3)', 'var(--s1)'], ['2', 'Some (2)', 'var(--s5)'], ['1', 'Little (1)', 'var(--s3)'], ['0', 'None (0)', 'var(--data-2)'], ['30', 'Recessive (30)', 'var(--s2)'], ['40', 'Unlikely (40)', 'var(--s4)']];

const renderDosage = () => {
  const { hi, ts, date } = CN.dosage;
  const total = (x) => sumOf(Object.values(x));
  const m = /(\d{1,2}) (\w{3})\w*,\s*(\d{4})/.exec(date ?? '');
  const month = m ? ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'].indexOf(m[2]) + 1 : 0;
  const listed = month ? ` of ${day(`${m[3]}-${String(month).padStart(2, '0')}-${m[1].padStart(2, '0')}T00:00:00Z`)}` : '';
  $('#dosage-dek').textContent = `${exact(total(hi))} genes have a haploinsufficiency score and ${exact(total(ts))} a triplosensitivity score in the list${listed}. ${exact(hi['30'] ?? 0)} are scored 30, genes tied to a recessive disease, and ${exact(hi['40'] ?? 0)} are scored 40, dosage sensitivity unlikely.`;
  hbars({
    mount: $('#dosage-chart'),
    rows: [{ key: 'hi', label: 'Haploinsufficiency', values: hi }, { key: 'ts', label: 'Triplosensitivity', values: ts }],
    segments: DOSAGE_SCORES.map(([key, label, color]) => ({ key, label, color })),
    labelShare: 0.3,
    aria: 'ClinGen dosage scores',
    describe: (r) => tipLines(`${exact(total(r.values))} genes`, r.label, ...DOSAGE_SCORES.filter(([k]) => r.values[k]).map(([k, label]) => `${label}: ${exact(r.values[k])}`)),
  });
  attachData('#dosage-chart', 'clingen-dosage', ['Curation', ...DOSAGE_SCORES.map(([, label]) => label)], [['Haploinsufficiency', ...DOSAGE_SCORES.map(([k]) => hi[k] ?? 0)], ['Triplosensitivity', ...DOSAGE_SCORES.map(([k]) => ts[k] ?? 0)]]);
};

const yearOfDate = (iso) => {
  const [y, m, d] = iso.split('-').map(Number);
  return y + (Date.UTC(y, m - 1, d) - Date.UTC(y, 0, 1)) / (Date.UTC(y + 1, 0, 1) - Date.UTC(y, 0, 1));
};

const ANC = snapshot.ancestry;
const AAD = snapshot.ancient;
const KG = snapshot.igsr;
const ASM = snapshot.assemblies;

const ANCESTRY = [
  ['European', 'European', 'var(--s2)'],
  ['East Asian', 'East Asian', 'var(--s4)'],
  ['African American or Afro-Caribbean', 'African American or Afro-Caribbean', 'var(--s1)'],
  ['Hispanic or Latin American', 'Hispanic or Latin American', 'var(--s3)'],
  ['South Asian', 'South Asian', 'var(--s5)'],
  ['Multiple categories', 'Multiple categories', 'var(--ink-2)'],
  ['NR', 'NR (not reported)', 'var(--data-2)'],
  ['other', 'Other categories', 'var(--data)'],
];
const NAMED_ANCESTRY = ANCESTRY.map(([k]) => k).filter((k) => k !== 'other');
const ancestryOf = (counts) => ({
  ...Object.fromEntries(NAMED_ANCESTRY.map((k) => [k, counts[k] ?? 0])),
  other: sumOf(Object.entries(counts).filter(([k]) => !NAMED_ANCESTRY.includes(k)).map(([, n]) => n)),
});
const sharesOf = (values) => {
  const total = sumOf(Object.values(values));
  return Object.fromEntries(Object.entries(values).map(([k, v]) => [k, total ? (v / total) * 100 : 0]));
};
const pctOf = (v) => pct(v / 100);

const renderAncestry = () => {
  const years = Object.keys(ANC.byYear).map(Number).sort((a, b) => a - b);
  const counts = years.map((y) => ancestryOf(ANC.byYear[y]));
  const yearShares = counts.map(sharesOf);
  const totals = counts.map((c) => sumOf(Object.values(c)));
  const stages = [['initial', 'Discovery'], ['replication', 'Replication']].map(([key, label]) => {
    const values = ancestryOf(ANC.stages[key]);
    return {
      key, label, n: sumOf(Object.values(values)), counts: values, values: sharesOf(values),
    };
  });
  const replication = stages[1];
  const eu = (i) => yearShares[i].European;
  const large = years.map((_, i) => i).filter((i) => totals[i] >= 1e6);
  const low = large.reduce((a, b) => (eu(b) < eu(a) ? b : a), large[0]);
  const high = large.reduce((a, b) => (eu(b) > eu(a) ? b : a), large[0]);
  const last = years.length - 1;
  $('#ancestry-dek').textContent = `${exact(ANC.studies)} studies list ${wordy(ANC.participants)} participant entries; a cohort counts again in every study that uses it. At replication, ${pctOf(replication.values.European)} are listed as European and ${pctOf(replication.values.NR)} have no ancestry reported. Among years with a million or more entries, the European share of discovery participants ranged from ${pctOf(eu(low))} in ${years[low]} to ${pctOf(eu(high))} in ${years[high]}; it is ${pctOf(eu(last))} for studies published in ${years[last]}.`;
  stackedArea({
    mount: $('#ancestry-years-chart'),
    years,
    series: ANCESTRY.map(([key, label, color]) => ({
      key, label, color, values: yearShares.map((s) => s[key]),
    })),
    max: 100,
    format: pctOf,
    yFormat: (v) => `${v}%`,
    head: (i) => `${years[i]}: ${exact(ANC.studiesByYear[years[i]] ?? 0)} studies, ${compact(totals[i])} entries`,
    aria: 'Share of GWAS discovery participants by broad ancestral category and publication year',
  });
  const ranked = (tables) => {
    const sums = {};
    for (const t of tables) for (const [c, n] of Object.entries(t)) sums[c] = (sums[c] ?? 0) + n;
    return Object.keys(sums).sort((a, b) => sums[b] - sums[a]);
  };
  const yearCats = ranked(years.map((y) => ANC.byYear[y]));
  attachData('#ancestry-years-chart', 'gwas-ancestry-by-year', ['Year', 'Studies', 'Discovery participant entries', ...yearCats], years.map((y, i) => [y, ANC.studiesByYear[y] ?? 0, totals[i], ...yearCats.map((c) => ANC.byYear[y][c] ?? 0)]));
  hbars({
    mount: $('#ancestry-stages-chart'),
    rows: stages,
    segments: ANCESTRY.map(([key, label, color]) => ({ key, label, color })),
    max: 100,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => compact(r.n),
    labelShare: 0.2,
    aria: 'Share of GWAS participant entries by broad ancestral category at discovery and replication',
    describe: (r) => tipLines(`${exact(r.n)} participant entries`, `${r.label} stage`, ...ANCESTRY.map(([k, label]) => `${label}: ${pctOf(r.values[k])}`)),
  });
  const stageCats = ranked(stages.map((r) => ANC.stages[r.key]));
  attachData('#ancestry-stages-chart', 'gwas-ancestry-by-stage', ['Stage', ...stageCats], stages.map((r) => [r.label, ...stageCats.map((c) => ANC.stages[r.key][c] ?? 0)]));
  dotplot({
    mount: $('#ancestry-countries-chart'),
    rows: ANC.countries.map((c, i) => ({
      key: c.country, label: c.country, value: c.participants, hi: i === 0,
    })),
    format: compact,
    rowH: 20,
    aria: 'Countries of recruitment with the most GWAS participant entries, log scale',
    describe: (r) => tipLines(`${exact(r.value)} participant entries`, r.key, `${pct(r.value / ANC.participants)} of all entries`),
  });
  attachData('#ancestry-countries-chart', 'gwas-countries', ['Country of recruitment', 'Participant entries'], ANC.countries.map((c) => [c.country, c.participants]));
};

const binRange = (edges, i, unit, format = exact) => {
  if (i === 0) return `under ${format(edges[1])} ${unit}`;
  if (i === edges.length - 2) return `${format(edges[i])} ${unit} or more`;
  return `${format(edges[i])} to ${format(edges[i + 1])} ${unit}`;
};
const cellSpan = (v, size, pos, neg) => (v >= 0 ? `${v}–${v + size}°${pos}` : `${-v - size}–${-v}°${neg}`);
const cellName = (c) => `${cellSpan(c.lat, 2, 'N', 'S')}, ${cellSpan(c.lon, 2, 'E', 'W')}`;

const renderAncient = () => {
  const {
    ages, snps, dataTypes, sex, grid,
  } = AAD;
  const old = sumOf(ages.counts.filter((_, i) => ages.edges[i] >= 1e4));
  const densest = grid.reduce((a, b) => (b.n > a.n ? b : a));
  const release = AAD.released ? ` of ${longDate(AAD.released)}` : '';
  $('#ancient-dek').textContent = `Ages are mean dates in years before 1950 CE; the median is ${exact(ages.median)}, and ${exact(old)} individuals date to 10,000 or more years before present. ${exact(dataTypes.Capture ?? 0)} were genotyped by targeted capture, ${exact(dataTypes.Shotgun ?? 0)} by shotgun sequencing and ${exact(dataTypes['Capture and shotgun'] ?? 0)} by both. Molecular sex: ${exact(sex.M ?? 0)} male, ${exact(sex.F ?? 0)} female, ${exact(sex.U ?? 0)} not determined. The densest 2° cell, ${cellName(densest)}, holds ${exact(densest.n)}. Release ${AAD.version}${release} also lists ${exact(AAD.present)} present-day people, left out here.`;
  histo({
    mount: $('#ancient-ages-chart'),
    bins: ages.counts.map((_, i) => i),
    series: [{
      key: 'ages', label: 'Individuals', color: 'var(--data)', values: ages.counts,
    }],
    tickEvery: 4,
    tickLabel: (b) => exact(ages.edges[b]),
    height: 200,
    aria: 'Ancient individuals by mean date, log scale of time',
    describe: (i) => tipLines(`${exact(ages.counts[i])} individuals`, binRange(ages.edges, i, 'years before present')),
  });
  attachData('#ancient-ages-chart', 'aadr-ages', ['Years BP from', 'Years BP to', 'Individuals'], ages.counts.map((n, i) => [Math.round(ages.edges[i]), Math.round(ages.edges[i + 1]), n]));
  const listed = Object.keys(AAD.published).map(Number);
  const first = Math.min(...listed);
  const releaseYear = AAD.released ? Number(AAD.released.slice(0, 4)) : null;
  const published = Array.from({ length: Math.max(...listed) - first + 1 }, (_, i) => ({
    year: first + i, n: AAD.published[first + i] ?? 0, partial: first + i === releaseYear,
  }));
  columns({
    mount: $('#ancient-published-chart'),
    rows: published,
    key: (r) => r.year,
    value: (r) => r.n,
    height: 200,
    aria: 'Ancient individuals by the year they were first published',
    describe: (r) => tipLines(`${exact(r.n)} individuals`, `first published in ${r.year}${r.partial ? `, up to the release${release}` : ''}`),
  });
  attachData('#ancient-published-chart', 'aadr-first-published', ['Year', 'Individuals'], published.map((r) => [r.year, r.n]));
  dotmap({
    mount: $('#ancient-map-chart'),
    cells: grid,
    size: 2,
    top: Math.max(84, ...grid.map((c) => c.lat + 2)),
    bottom: Math.min(-56, ...grid.map((c) => c.lat)),
    aria: 'Ancient individuals per 2 by 2 degree cell of latitude and longitude',
    describe: (c) => tipLines(`${exact(c.n)} individuals`, cellName(c), `median ${exact(c.age)} years before present`),
  });
  attachData('#ancient-map-chart', 'aadr-grid', ['Latitude from', 'Longitude from', 'Individuals', 'Median years BP'], grid.map((c) => [c.lat, c.lon, c.n, c.age]));
  dotplot({
    mount: $('#ancient-countries-chart'),
    rows: AAD.countries.map((c, i) => ({
      key: c.country, label: c.country, value: c.n, hi: i === 0,
    })),
    scale: 'linear',
    rowH: 20,
    aria: 'Countries with the most ancient individuals',
    describe: (r) => tipLines(`${exact(r.value)} individuals`, r.key, `${pct(r.value / AAD.individuals)} of all`),
  });
  attachData('#ancient-countries-chart', 'aadr-countries', ['Country', 'Individuals'], AAD.countries.map((c) => [c.country, c.n]));
  histo({
    mount: $('#ancient-snps-chart'),
    bins: snps.counts.map((_, i) => i),
    series: [{
      key: 'snps', label: 'Individuals', color: 'var(--data)', values: snps.counts,
    }],
    tickEvery: 4,
    tickLabel: (b) => compact(snps.edges[b]),
    height: 200,
    aria: 'Ancient individuals by the number of 1240k sites covered, log scale',
    describe: (i) => tipLines(`${exact(snps.counts[i])} individuals`, binRange(snps.edges, i, 'sites covered')),
  });
  attachData('#ancient-snps-chart', 'aadr-sites-covered', ['Sites from', 'Sites to', 'Individuals'], snps.counts.map((n, i) => [Math.round(snps.edges[i]), Math.round(snps.edges[i + 1]), n]));
  for (const [kind, list, label] of [['y', AAD.y, 'Y'], ['mt', AAD.mt, 'mtDNA']]) {
    dotplot({
      mount: $(`#ancient-${kind}-chart`),
      rows: list.map((g, i) => ({
        key: g.group, label: g.group, value: g.n, hi: i === 0,
      })),
      scale: 'linear',
      rowH: 18,
      labelShare: 0.1,
      aria: `Ancient individuals by major ${label} haplogroup`,
      describe: (r) => tipLines(`${exact(r.value)} individuals`, `${label} haplogroup ${r.key}`),
    });
    attachData(`#ancient-${kind}-chart`, `aadr-${kind}-haplogroups`, [`${label} haplogroup`, 'Individuals'], list.map((g) => [g.group, g.n]));
  }
};

const SUPER_COLORS = {
  AFR: 'var(--s1)', AMR: 'var(--s3)', EAS: 'var(--s4)', EUR: 'var(--s2)', SAS: 'var(--s5)',
};
const QUALIFIER = { about: '~', 'more than': '>' };
const shortCite = (key) => SOURCE_TITLES[key]?.[1].split(', ')[0] ?? key;

const renderPanels = () => {
  const { panel } = KG;
  const supers = [...new Set(panel.populations.map((p) => p.superpopulation))].sort();
  const superName = Object.fromEntries(panel.populations.map((p) => [p.superpopulation, p.superpopulationName]));
  const rows = supers.flatMap((sp) => panel.populations.filter((p) => p.superpopulation === sp).map((p) => ({
    key: p.code, label: `${p.name} (${p.code})`, values: { [sp]: p.n }, pop: p,
  })));
  const bySuper = (sp) => sumOf(panel.populations.filter((p) => p.superpopulation === sp).map((p) => p.n));
  $('#panels-dek').textContent = `The 30x set has ${exact(panel.sex.M ?? 0)} male and ${exact(panel.sex.F ?? 0)} female samples: ${supers.map((sp) => `${exact(bySuper(sp))} ${sp}`).join(', ')}. IGSR holds ${exact(KG.total)} samples across ${exact(KG.collections.length)} data collections. gnomAD v4 pools exomes from ${exact(GNOMAD_V4.exomes)} people, ${exact(GNOMAD_V4.ukBiobank)} of them from UK Biobank, and genomes from ${exact(GNOMAD_V4.genomes)}; about ${exact(GNOMAD_V4.nonEuropean)} are of non-European genetic ancestry.`;
  hbars({
    mount: $('#panels-1kg-chart'),
    rows,
    segments: supers.map((sp) => ({ key: sp, label: `${superName[sp]} (${sp})`, color: SUPER_COLORS[sp] ?? 'var(--data)' })),
    labelShare: 0.34,
    labelMax: 180,
    rowH: 18,
    aria: 'Samples per population in the 1000 Genomes 30x set, by superpopulation',
    describe: (r) => tipLines(`${exact(r.pop.n)} samples`, `${r.pop.name} (${r.pop.code})`, `${r.pop.superpopulationName} (${r.pop.superpopulation})`),
  });
  attachData('#panels-1kg-chart', '1000-genomes-30x-populations', ['Code', 'Population', 'Superpopulation', 'Samples'], rows.map((r) => [r.pop.code, r.pop.name, r.pop.superpopulationName, r.pop.n]));
  dotplot({
    mount: $('#panels-biobanks-chart'),
    rows: BIOBANKS.map((b, i) => ({
      key: b.name, label: b.name, value: b.people, hi: i === 0, bank: b,
    })),
    scale: 'linear',
    format: (v, r) => `${QUALIFIER[r.bank.qualifier] ?? ''}${exact(v)}`,
    aria: 'Participants in eight large biobanks, as each cohort paper reports them',
    describe: (r) => tipLines(`${r.bank.qualifier ? `${r.bank.qualifier} ` : ''}${exact(r.value)} ${r.bank.measure}`, r.bank.name, { href: SOURCES[r.bank.source], text: shortCite(r.bank.source) }),
  });
  attachData('#panels-biobanks-chart', 'biobank-sizes', ['Biobank', 'Qualifier', 'Count', 'Counted as', 'Source'], BIOBANKS.map((b) => [b.name, b.qualifier ?? '', b.people, b.measure, SOURCES[b.source]]));
};

const ASSEMBLY_LEVELS = [['Contig', 'var(--s4)'], ['Scaffold', 'var(--s3)'], ['Chromosome', 'var(--s2)'], ['Complete Genome', 'var(--s5)']];

const renderAssemblies = () => {
  const known = ASSEMBLY_LEVELS.map(([l]) => l);
  const levels = [...ASSEMBLY_LEVELS.filter(([l]) => l in ASM.levels), ...Object.keys(ASM.levels).filter((l) => !known.includes(l)).map((l) => [l, 'var(--data)'])];
  const colorOf = Object.fromEntries(levels);
  const listed = Object.keys(ASM.byYear).map(Number);
  const first = Math.min(...listed);
  const years = Array.from({ length: Math.max(...listed) - first + 1 }, (_, i) => first + i);
  const perYear = years.map((y) => ASM.byYear[y] ?? {});
  const n50Of = (y) => ASM.n50.filter((p) => p[0].startsWith(String(y))).map((p) => p[1]);
  const busy = years.filter((y) => n50Of(y).length >= 10);
  const [y0, y1] = [busy[0], busy.at(-1)];
  const over100 = ASM.n50.filter((p) => p[1] >= 1e8).length;
  const best = ASM.best;
  $('#assemblies-dek').textContent = `By assembly level: ${Object.entries(ASM.levels).sort((a, b) => b[1] - a[1]).map(([l, n]) => `${l.toLowerCase()} ${exact(n)}`).join(', ')}. The median contig N50 was ${bases(summary(n50Of(y0)).median)} for assemblies released in ${y0} and ${bases(summary(n50Of(y1)).median)} for ${y1}. ${exact(over100)} have a contig N50 of 100 Mb or more; the longest is ${best.name} (${best.accession}) from ${best.submitter}.`;
  const toDate = perYear.map((_, i) => Object.fromEntries(levels.map(([l]) => [l, sumOf(perYear.slice(0, i + 1).map((x) => x[l] ?? 0))])));
  stackedArea({
    mount: $('#assemblies-years-chart'),
    years,
    series: levels.map(([l, color]) => ({
      key: l, label: l, color, values: toDate.map((x) => x[l]),
    })),
    head: (i) => `${years[i]}: ${exact(sumOf(Object.values(toDate[i])))} by year end, ${exact(sumOf(Object.values(perYear[i])))} released that year`,
    aria: 'Human genome assemblies at NCBI released by the end of each year, by assembly level',
  });
  attachData('#assemblies-years-chart', 'assemblies-by-year', ['Year', ...levels.map(([l]) => `${l} released`), 'All released by year end'], years.map((y, i) => [y, ...levels.map(([l]) => perYear[i][l] ?? 0), sumOf(Object.values(toDate[i]))]));
  const points = ASM.n50.filter(([date]) => /^\d{4}-\d{2}-\d{2}/.test(date)).map(([date, n50, level, accession]) => ({
    key: accession, x: yearOfDate(date), y: n50, color: colorOf[level] ?? 'var(--data)', hi: accession === best.accession, date, level, accession,
  }));
  scatter({
    mount: $('#assemblies-n50-chart'),
    points,
    yLog: true,
    xFormat: (v) => String(v),
    yFormat: (v) => (v >= 1e9 ? `${v / 1e9} Gb` : v >= 1e6 ? `${v / 1e6} Mb` : v >= 1e3 ? `${v / 1e3} kb` : `${v} bp`),
    xLabel: 'release date',
    yLabel: 'contig N50, log scale',
    legend: levels.map(([label, color]) => ({ label, color })),
    aria: 'Contig N50 of each human genome assembly at NCBI by release date, log scale',
    describe: (p) => tipLines(bases(p.y), `contig N50 of ${p.accession}`, `${p.level}, released ${longDate(p.date)}`),
  });
  attachData('#assemblies-n50-chart', 'assemblies-contig-n50', ['Accession', 'Release date', 'Level', 'Contig N50'], ASM.n50.map(([date, n50, level, accession]) => [accession, date, level, n50]));
  dotplot({
    mount: $('#assemblies-submitters-chart'),
    rows: ASM.submitters.map((s, i) => ({
      key: s.submitter, label: s.submitter, value: s.n, hi: i === 0,
    })),
    scale: 'linear',
    rowH: 20,
    labelShare: 0.48,
    labelMax: 530,
    aria: 'Submitters of the most human genome assemblies at NCBI',
    describe: (r) => tipLines(`${exact(r.value)} assemblies`, r.key, pct(r.value / ASM.total)),
  });
  attachData('#assemblies-submitters-chart', 'assembly-submitters', ['Submitter', 'Assemblies'], ASM.submitters.map((s) => [s.submitter, s.n]));
};

const CV = snapshot.clinvar;
const CT = snapshot.clinvarTable;
const CLASS_KEYS = ['Pathogenic', 'Pathogenic/Likely pathogenic', 'Likely pathogenic', 'Uncertain significance', 'Likely benign', 'Benign/Likely benign', 'Benign', 'Conflicting', 'Other'];
const CLASS_SHORT = ['P', 'P/LP', 'LP', 'VUS', 'LB', 'B/LB', 'B', 'Confl.', 'Other'];
const GROUPS = [
  ['P/LP', 'Pathogenic or likely pathogenic', 'var(--accent)'],
  ['VUS', 'Uncertain significance', 'var(--s5)'],
  ['B/LB', 'Benign or likely benign', 'var(--s2)'],
  ['Conflicting', 'Conflicting', 'var(--s3)'],
  ['Other', 'Other classifications', 'var(--data-2)'],
];
const GROUP_OF = {
  Pathogenic: 'P/LP', 'Pathogenic/Likely pathogenic': 'P/LP', 'Likely pathogenic': 'P/LP', 'Uncertain significance': 'VUS', 'Likely benign': 'B/LB', 'Benign/Likely benign': 'B/LB', Benign: 'B/LB', Conflicting: 'Conflicting', Other: 'Other',
};
const byGroup = (classes) => Object.fromEntries(GROUPS.map(([g]) => [g, sumOf(CLASS_KEYS.filter((c) => GROUP_OF[c] === g).map((c) => classes[c] ?? 0))]));

const renderClinvar = () => {
  const cons = Object.keys(CV.consequences);
  const counts = cons.map((k) => CLASS_KEYS.map((c) => CV.consequences[k][c] ?? 0));
  const totals = counts.map(sumOf);
  const rowShares = counts.map((r, i) => r.map((v) => (totals[i] ? v / totals[i] : 0)));
  const rowShare = (k, group) => {
    const i = cons.indexOf(k);
    return sumOf(CLASS_KEYS.map((c, j) => (GROUP_OF[c] === group ? counts[i][j] : 0))) / totals[i];
  };
  $('#clinvar-dek').textContent = `${pct(rowShare('nonsense', 'P/LP'))} of nonsense and ${pct(rowShare('frameshift', 'P/LP'))} of frameshift variants are pathogenic or likely pathogenic, against ${pct(rowShare('missense', 'P/LP'))} of the ${exact(totals[cons.indexOf('missense')])} missense variants, ${pct(rowShare('missense', 'VUS'))} of which are of uncertain significance. ${exact(CV.classes.none)} records carry no aggregate germline classification and are left out.`;
  shareGrid({
    mount: $('#clinvar-grid-chart'),
    rows: cons,
    cols: CLASS_SHORT,
    values: rowShares,
    totals,
    aria: 'ClinVar germline classes as a share of each molecular consequence',
    describe: (i, j) => tipLines(pct(rowShares[i][j]), `of ${cons[i]} variants are ${CLASS_KEYS[j].toLowerCase()}`, `${exact(counts[i][j])} of ${exact(totals[i])}`),
  });
  attachData('#clinvar-grid-chart', 'clinvar-classes-by-consequence', ['Consequence', ...CLASS_KEYS, 'Total'], cons.map((k, i) => [k, ...counts[i], totals[i]]));
  const types = Object.entries(CT.types).sort((a, b) => b[1] - a[1]);
  dotplot({
    mount: $('#clinvar-types-chart'),
    rows: types.map(([t, n], i) => ({
      key: t, label: t, value: n, hi: i === 0,
    })),
    aria: 'ClinVar variants on GRCh38 by type, log scale',
    describe: (r) => tipLines(`${exact(r.value)} variants`, r.label, pct(r.value / CT.records)),
  });
  attachData('#clinvar-types-chart', 'clinvar-types', ['Type', 'Variants'], types);
  const GN = CV.genes;
  dumbbell({
    mount: $('#clinvar-genes-chart'),
    rows: GN.top.map((g) => ({
      key: g.gene, label: g.gene, a: g.plp, b: g.vus,
    })),
    years: ['P/LP alleles', 'VUS alleles'],
    valueText: (r) => `${exact(r.a)} · ${exact(r.b)}`,
    aria: 'Protein-coding genes with the most pathogenic or likely pathogenic ClinVar alleles, with their VUS, log scale',
    describe: (r) => {
      const g = GN.top.find((x) => x.gene === r.key);
      return tipLines(`${exact(g.plp)} P/LP alleles`, g.gene, `${exact(g.vus)} of uncertain significance`, `${exact(g.conflicts)} with conflicting classifications`, `${exact(g.alleles)} alleles in all`, { href: `https://www.ncbi.nlm.nih.gov/clinvar/?term=${encodeURIComponent(`${g.gene}[gene]`)}`, text: 'Open in ClinVar' });
    },
  });
  attachData('#clinvar-genes-chart', 'clinvar-genes', ['Gene', 'P/LP alleles', 'VUS alleles', 'Conflicting alleles', 'All alleles'], GN.top.map((g) => [g.gene, g.plp, g.vus, g.conflicts, g.alleles]));
  lines({
    mount: $('#clinvar-lorenz-chart'),
    series: [{ label: 'P/LP alleles', rows: GN.lorenz.map((v, x) => ({ x, y: v * 100 })), cls: 'alt' }],
    scale: 'linear',
    height: 240,
    yFormat: (v) => `${v}%`,
    labelW: 84,
    aria: 'Cumulative share of P/LP alleles held by the genes with the most',
    describe: (x) => tipLines(`${(GN.lorenz[x] * 100).toFixed(1)}%`, `of P/LP alleles sit in the top ${x}% of the ${exact(GN.withPlp)} genes with any`),
  });
  attachData('#clinvar-lorenz-chart', 'clinvar-lorenz', ['Top share of genes (%)', 'Share of P/LP alleles (%)'], GN.lorenz.map((v, x) => [x, Number((v * 100).toFixed(2))]));
};

const STARS = [['0 stars', 'var(--data-2)'], ['1 star', 'var(--s3)'], ['2 stars', 'var(--s2)'], ['3 stars', 'var(--s4)'], ['4 stars', 'var(--ink)']];

const renderReview = () => {
  const rows = CLASS_KEYS.map((c, i) => {
    const s = CV.stars[c];
    const n = sumOf(s);
    return {
      key: c, label: CLASS_SHORT[i], n, counts: s, values: Object.fromEntries(s.map((v, k) => [k, n ? (v / n) * 100 : 0])),
    };
  });
  const sub = CT.submitters;
  const groupNames = ['P/LP', 'VUS', 'B/LB', 'Conflicting'];
  const several = (g) => 1 - sub[g][0] / sumOf(sub[g]);
  const O = CV.organizations;
  const lead = O.top[0];
  $('#review-dek').textContent = `Practice guidelines cover ${exact(sumOf(Object.values(CV.stars).map((s) => s[4])))} variants. ${pct(several('P/LP'))} of P/LP variants have two or more submitters, against ${pct(several('VUS'))} of VUS. ${exact(O.organizations)} organizations have made ${exact(O.submissions)} submissions; ${lead.name.split(';')[0]} alone made ${exact(lead.submissions)}, ${pct(lead.submissions / O.submissions)}.`;
  hbars({
    mount: $('#review-stars-chart'),
    rows,
    segments: STARS.map(([label, color], k) => ({ key: String(k), label, color })),
    max: 100,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => compact(r.n),
    labelShare: 0.1,
    aria: 'ClinVar review stars as a share of each germline class',
    describe: (r) => tipLines(`${exact(r.n)} variants`, r.key, ...r.counts.map((v, k) => `${STARS[k][0]}: ${exact(v)}, ${r.n ? pct(v / r.n) : '0%'}`)),
  });
  attachData('#review-stars-chart', 'clinvar-stars', ['Class', ...STARS.map(([label]) => label)], rows.map((r) => [r.key, ...r.counts]));
  const share = (g) => sub[g].map((v) => (v / sumOf(sub[g])) * 100);
  histo({
    mount: $('#review-submitters-chart'),
    bins: sub.bins,
    series: groupNames.map((g) => ({
      key: g, label: g, color: GROUPS.find((x) => x[0] === g)[2], values: share(g),
    })),
    centered: true,
    tickEvery: 1,
    tickLabel: (b) => b,
    yFormat: (v) => `${v}%`,
    aria: 'Submitters per ClinVar variant, as a share of each group',
    describe: (i) => tipLines(`${sub.bins[i]} submitter${sub.bins[i] === '1' ? '' : 's'}`, 'share of each group', ...groupNames.map((g) => `${g}: ${share(g)[i].toFixed(1)}%, ${exact(sub[g][i])}`)),
  });
  attachData('#review-submitters-chart', 'clinvar-submitters', ['Submitters', ...groupNames], sub.bins.map((b, i) => [b, ...groupNames.map((g) => sub[g][i])]));
  const years = new Map();
  for (const g of groupNames) {
    for (const [y, n] of Object.entries(CT.evaluated[g])) if (y !== 'none') years.set(Number(y), (years.get(Number(y)) ?? 0) + n);
  }
  const last = Math.max(...years.keys());
  const early = sumOf([...years].filter(([y]) => y < 2005).map(([, n]) => n));
  const evalRows = [{ key: 'early', label: '<2005', value: early }, ...Array.from({ length: last - 2004 }, (_, i) => ({ key: 2005 + i, label: String(2005 + i), value: years.get(2005 + i) ?? 0 }))];
  const undated = sumOf(groupNames.map((g) => CT.evaluated[g].none ?? 0));
  columns({
    mount: $('#review-evaluated-chart'),
    rows: evalRows,
    key: (r) => r.key,
    value: (r) => r.value,
    label: (r) => r.label,
    aria: 'ClinVar classified variants by the year their classification was last evaluated',
    describe: (r) => tipLines(`${exact(r.value)} variants`, r.key === 'early' ? 'last evaluated before 2005' : `last evaluated in ${r.label}${r.key === last ? ', to date' : ''}`, `${exact(undated)} have no date`),
  });
  attachData('#review-evaluated-chart', 'clinvar-last-evaluated', ['Year', ...groupNames], [['before 2005', ...groupNames.map((g) => sumOf(Object.entries(CT.evaluated[g]).filter(([y]) => y !== 'none' && Number(y) < 2005).map(([, n]) => n)))], ...evalRows.slice(1).map((r) => [r.label, ...groupNames.map((g) => CT.evaluated[g][r.label] ?? 0)]), ['no date', ...groupNames.map((g) => CT.evaluated[g].none ?? 0)]]);
  dotplot({
    mount: $('#review-orgs-chart'),
    rows: O.top.map((o, i) => ({
      key: o.name, label: o.name, value: o.submissions, hi: i === 0,
    })),
    labelShare: 0.42,
    labelMax: 300,
    aria: 'Organizations with the most ClinVar submissions, log scale',
    describe: (r) => {
      const o = O.top.find((x) => x.name === r.key);
      return tipLines(`${exact(o.submissions)} submissions`, o.name, `${o.type}, ${o.country}`, pct(o.submissions / O.submissions));
    },
  });
  attachData('#review-orgs-chart', 'clinvar-submitters-top', ['Organization', 'Type', 'Country', 'Submissions'], O.top.map((o) => [o.name, o.type, o.country, o.submissions]));
  dotplot({
    mount: $('#review-types-chart'),
    rows: O.types.map((t, i) => ({
      key: t.type, label: `${t.type} (${exact(t.orgs)})`, value: t.submissions, hi: i === 0,
    })),
    aria: 'ClinVar submissions by type of organization, log scale',
    describe: (r) => {
      const t = O.types.find((x) => x.type === r.key);
      return tipLines(`${exact(t.submissions)} submissions`, `${t.type}: ${exact(t.orgs)} organizations`, pct(t.submissions / O.submissions));
    },
  });
  attachData('#review-types-chart', 'clinvar-submitter-types', ['Type', 'Organizations', 'Submissions'], O.types.map((t) => [t.type, t.orgs, t.submissions]));
};

const renderGrowth = () => {
  const { points } = snapshot.clinvarGrowth;
  const xs = points.map((p) => yearOfDate(p.release));
  const groups = points.map((p) => byGroup(p.classes));
  const first = points[0];
  const last = points.at(-1);
  const share = (p, g) => byGroup(p.classes)[g] / sumOf(Object.values(byGroup(p.classes)));
  $('#growth-dek').textContent = `Counts of the first GRCh38 VCF in each year’s archive folder, written in the last days of the year before, and the current file of ${longDate(last.release)}. The share of uncertain significance went from ${pct(share(first, 'VUS'))} to ${pct(share(last, 'VUS'))}; pathogenic or likely pathogenic from ${pct(share(first, 'P/LP'))} to ${pct(share(last, 'P/LP'))}.`;
  stackedArea({
    mount: $('#growth-chart'),
    years: xs,
    series: GROUPS.map(([g, label, color]) => ({
      key: g, label, color, values: groups.map((x) => x[g]),
    })),
    xTicks: xs,
    xTick: (v) => String(Math.floor(v)),
    xFormat: (v) => longDate(points[xs.indexOf(v)].release),
    aria: 'ClinVar GRCh38 VCF records by germline class at the end of each year',
  });
  attachData('#growth-chart', 'clinvar-growth', ['File date', 'Records', ...GROUPS.map(([g]) => g), 'No germline class'], points.map((p, i) => [p.release, p.records, ...GROUPS.map(([g]) => groups[i][g]), p.classes.none]));
};

const VARIANT_KINDS = [
  { key: 'missense', label: 'Missense', color: 'var(--s2)' },
  { key: 'truncating', label: 'Nonsense or frameshift', color: 'var(--accent)' },
  { key: 'splice', label: 'Splice site, at the exon edge', color: 'var(--s4)' },
  { key: 'inframe', label: 'In-frame indel', color: 'var(--s3)' },
  { key: 'synonymous', label: 'Synonymous', color: 'var(--data-2)' },
  { key: 'other', label: 'Start or stop codon', color: 'var(--s5)' },
];
const KIND_LABEL = Object.fromEntries(VARIANT_KINDS.map((k) => [k.key, k.label.toLowerCase()]));
const AA3_TO_1 = Object.fromEntries(Object.entries(AA3).map(([one, three]) => [three, one]));

const sitesOf = (lollipop) => {
  const map = new Map();
  for (const side of ['up', 'down']) {
    for (const [pos, kind, n, aa] of lollipop[side]) {
      const key = `${side}:${pos}`;
      const s = map.get(key) ?? {
        key, side, pos, n: 0, kinds: {}, aa: '',
      };
      s.n += n;
      s.kinds[kind] = (s.kinds[kind] ?? 0) + n;
      s.aa ||= aa;
      map.set(key, s);
    }
  }
  return [...map.values()].map((s) => ({ ...s, kind: Object.entries(s.kinds).sort((a, b) => b[1] - a[1])[0][0] }));
};

const shortDomain = (name) => name.replace(/^Zinc finger, C3HC4 type \(RING finger\)$/, 'RING').replace(/^BRCA1 C Terminus \(BRCT\) domain$/, 'BRCT').replace(/^ABC transporter transmembrane region$/, 'TMD').replace(/^ABC transporter$/, 'NBD').replace(/^Cystic fibrosis TM conductance regulator \(CFTR\), regulator domain$/, 'R domain')
  .replace(/^Myosin head \(motor domain\)$/, 'Motor domain')
  .replace(/^Myosin N-terminal SH3-like domain$/, 'SH3')
  .replace(/^Myosin tail$/, 'Coiled-coil tail')
  .replace(/^P53 DNA-binding domain$/, 'DNA-binding')
  .replace(/^P53 tetramerisation motif$/, 'Tet')
  .replace(/^P53 transactivation motif$/, 'TAD1')
  .replace(/^Transactivation domain 2$/, 'TAD2')
  .replace(/^Serine-rich domain associated with BRCT$/, 'Ser-rich');

const renderLollipops = () => {
  const genes = ['TP53', 'BRCA1', 'CFTR', 'MYH7'];
  const placed = (list) => sumOf(list.map((r) => r[2]));
  $('#lollipops-dek').textContent = `${genes.map((g) => `${g} ${exact(placed(CT.lollipops[g].up))} P/LP and ${exact(placed(CT.lollipops[g].down))} B/LB`).join('; ')}, placed on ${genes.map((g) => CT.lollipops[g].transcript).join(', ')}. Large deletions, deep intronic and untranslated variants have no residue and are left out.`;
  for (const gene of genes) {
    const l = CT.lollipops[gene];
    const p = snapshot.domains[gene];
    const sites = sitesOf(l);
    const ids = [...new Set(p.domains.map((d) => d.id))];
    const domains = p.domains.map((d) => ({ ...d, short: shortDomain(d.name), color: CAT[ids.indexOf(d.id) % CAT.length] }));
    const labels = sites.filter((s) => s.side === 'up' && s.kind !== 'splice').sort((a, b) => b.n - a.n || a.pos - b.pos).slice(0, 3)
      .map((s) => ({ ...s, label: `${AA3_TO_1[s.aa] ?? ''}${s.pos}` }));
    lollipopChart({
      mount: $(`#lollipop-${gene}-chart`),
      length: p.length ?? Math.max(...sites.map((s) => s.pos)),
      domains,
      sites,
      kinds: VARIANT_KINDS,
      labels,
      sideLabel: `P/LP above, ${exact(placed(l.up))} · B/LB below, ${exact(placed(l.down))}`,
      aria: `${gene} protein with ClinVar P/LP variants above and B/LB variants below, by residue`,
      describe: (s) => tipLines(`${s.aa ? `${s.aa}${s.pos}` : `codon ${s.pos}`}: ${exact(s.n)} variant${s.n === 1 ? '' : 's'}`, `${gene}, ${s.side === 'up' ? 'pathogenic or likely pathogenic' : 'benign or likely benign'}`, ...Object.entries(s.kinds).map(([k, n]) => `${KIND_LABEL[k]}: ${exact(n)}`)),
      describeDomain: (d) => tipLines(d.name, `${d.id}, residues ${exact(d.start)} to ${exact(d.end)}`, { href: `https://www.ebi.ac.uk/interpro/entry/pfam/${d.id}/`, text: 'Open in InterPro' }),
    });
    attachData(`#lollipop-${gene}-chart`, `lollipop-${gene.toLowerCase()}`, ['Side', 'Residue', 'Amino acid', 'Kind', 'Variants'], [...l.up.map((r) => ['P/LP', r[0], r[3] || null, r[1], r[2]]), ...l.down.map((r) => ['B/LB', r[0], r[3] || null, r[1], r[2]])]);
  }
};

const SUB_COLORS = ['var(--s2)', 'var(--ink)', 'var(--accent)', 'var(--data-2)', 'var(--s3)', 'var(--s4)'];

const renderSpectrum = () => {
  const s = CV.spectrum;
  const nUp = sumOf(s.plp);
  const nDown = sumOf(s.blb);
  const up = s.plp.map((v) => (v / nUp) * 100);
  const down = s.blb.map((v) => (v / nDown) * 100);
  const cpg = (list) => sumOf([2, 6, 10, 14].map((k) => list[32 + k])) / sumOf(list);
  $('#spectrum-dek').textContent = `${exact(nUp)} P/LP and ${exact(nDown)} B/LB single-nucleotide variants, each read with its two neighbours from the hg38 sequence and turned to the pyrimidine strand. C>T at CpG makes up ${pct(cpg(s.plp))} of P/LP and ${pct(cpg(s.blb))} of B/LB SNVs; transitions per transversion are ${titv(s.plp).toFixed(2)} and ${titv(s.blb).toFixed(2)}. ${s.mismatches ? `${exact(s.mismatches)} reference alleles disagree with hg38.` : 'Every reference allele matches hg38.'}`;
  spectrum({
    mount: $('#spectrum-chart'),
    channels: CHANNELS,
    up,
    down,
    colors: SUB_COLORS,
    groups: SUBSTITUTIONS,
    upLabel: `P/LP, ${exact(nUp)} SNVs`,
    downLabel: `B/LB, ${exact(nDown)} SNVs`,
    aria: 'Trinucleotide substitution profile of P/LP SNVs above and B/LB SNVs below',
    describe: (i) => tipLines(CHANNELS[i], `${SUBSTITUTIONS[Math.floor(i / 16)]} with ${CHANNELS[i][0]} before and ${CHANNELS[i][6]} after`, `P/LP: ${up[i].toFixed(2)}%, ${exact(s.plp[i])}`, `B/LB: ${down[i].toFixed(2)}%, ${exact(s.blb[i])}`),
  });
  attachData('#spectrum-chart', 'clinvar-snv-spectrum', ['Channel', 'P/LP SNVs', 'B/LB SNVs'], CHANNELS.map((c, i) => [c, s.plp[i], s.blb[i]]));
  const sixUp = SUBSTITUTIONS.map((_, g) => sumOf(up.slice(g * 16, g * 16 + 16)));
  const sixDown = SUBSTITUTIONS.map((_, g) => sumOf(down.slice(g * 16, g * 16 + 16)));
  histo({
    mount: $('#spectrum-classes-chart'),
    bins: SUBSTITUTIONS,
    series: [
      { key: 'plp', label: 'P/LP', color: 'var(--accent)', values: sixUp },
      { key: 'blb', label: 'B/LB', color: 'var(--s2)', values: sixDown },
    ],
    centered: true,
    tickEvery: 1,
    tickLabel: (b) => b,
    yFormat: (v) => `${v}%`,
    height: 200,
    aria: 'Six substitution classes as a share of P/LP and B/LB SNVs',
    describe: (i) => tipLines(SUBSTITUTIONS[i], 'share of SNVs', `P/LP: ${sixUp[i].toFixed(1)}%`, `B/LB: ${sixDown[i].toFixed(1)}%`),
  });
  attachData('#spectrum-classes-chart', 'clinvar-substitution-classes', ['Substitution', 'P/LP share (%)', 'B/LB share (%)'], SUBSTITUTIONS.map((x, i) => [x, Number(sixUp[i].toFixed(2)), Number(sixDown[i].toFixed(2))]));
};

const VALIDITY = [['Definitive', 'var(--accent)'], ['Strong', 'var(--s4)'], ['Moderate', 'var(--s5)'], ['Limited', 'var(--s3)'], ['Disputed', 'var(--s2)'], ['Refuted', 'var(--ink)'], ['No Known Disease Relationship', 'var(--data-2)']];
const MOI = {
  AD: 'autosomal dominant', AR: 'autosomal recessive', XL: 'X-linked', SD: 'semidominant', MT: 'mitochondrial', UD: 'undetermined',
};

const renderValidity = () => {
  const V = snapshot.validity;
  const mois = Object.keys(V.byMoi).sort((a, b) => sumOf(Object.values(V.byMoi[b])) - sumOf(Object.values(V.byMoi[a])));
  const { AD, AR } = V.loeufByMoi;
  $('#validity-dek').textContent = `${exact(V.diseases)} diseases, as MONDO terms, from the curation file of ${longDate(V.date)}. Genes with a definitive, strong or moderate dominant link have a median LOEUF of ${AD.median.toFixed(2)}, recessive ones ${AR.median.toFixed(2)}. Of the ${V.gceps.length} expert panels, the ${V.gceps[0].name.replace(/ Gene Curation Expert Panel$/, '')} panel has made the most curations, ${exact(V.gceps[0].n)}.`;
  hbars({
    mount: $('#validity-chart'),
    rows: mois.map((m) => ({ key: m, label: m, values: V.byMoi[m] })),
    segments: VALIDITY.map(([key, color]) => ({ key, label: key, color })),
    labelShare: 0.08,
    aria: 'ClinGen gene-disease curations by mode of inheritance and classification',
    describe: (r) => tipLines(`${exact(sumOf(Object.values(r.values)))} curations`, `${r.key}, ${MOI[r.key] ?? r.key}`, ...VALIDITY.filter(([k]) => r.values[k]).map(([k]) => `${k}: ${exact(r.values[k])}`)),
  });
  attachData('#validity-chart', 'clingen-validity', ['Inheritance', ...VALIDITY.map(([k]) => k)], mois.map((m) => [m, ...VALIDITY.map(([k]) => V.byMoi[m][k] ?? 0)]));
  const years = Object.keys(V.byYear).map(Number).sort((a, b) => a - b);
  stackedArea({
    mount: $('#validity-years-chart'),
    years,
    series: VALIDITY.map(([key, color]) => ({
      key, label: key, color, values: years.map((y) => V.byYear[y][key] ?? 0),
    })),
    aria: 'ClinGen gene-disease curations by year of classification',
  });
  attachData('#validity-years-chart', 'clingen-validity-years', ['Year', ...VALIDITY.map(([k]) => k)], years.map((y) => [y, ...VALIDITY.map(([k]) => V.byYear[y][k] ?? 0)]));
  const bins = binsOf(linearEdges(0, 2, 0.05));
  const share = (g) => g.counts.map((c) => (c / g.n) * 100);
  histo({
    mount: $('#validity-loeuf-chart'),
    bins,
    series: [
      { key: 'AD', label: `Dominant (${exact(AD.n)} genes)`, color: 'var(--accent)', values: share(AD) },
      { key: 'AR', label: `Recessive (${exact(AR.n)} genes)`, color: 'var(--s2)', values: share(AR) },
    ],
    tickEvery: 4,
    tickLabel: (b) => b.lo.toFixed(1),
    yFormat: (v) => `${v}%`,
    height: 200,
    aria: 'LOEUF of genes with a supported dominant or recessive disease link',
    describe: (i) => tipLines(`LOEUF ${bins[i].lo.toFixed(2)} to ${bins[i].hi.toFixed(2)}`, 'share of each group', `Dominant: ${share(AD)[i].toFixed(1)}%, ${exact(AD.counts[i])}`, `Recessive: ${share(AR)[i].toFixed(1)}%, ${exact(AR.counts[i])}`),
  });
  attachData('#validity-loeuf-chart', 'clingen-validity-loeuf', ['LOEUF from', 'LOEUF to', 'Dominant genes', 'Recessive genes'], bins.map((b, i) => [b.lo, b.hi, AD.counts[i], AR.counts[i]]));
  const top = V.gceps.slice(0, 12);
  dotplot({
    mount: $('#validity-gcep-chart'),
    rows: top.map((g, i) => ({
      key: g.name, label: g.name.replace(/ Gene Curation Expert Panel$/, ''), value: g.n, hi: i === 0,
    })),
    scale: 'linear',
    labelShare: 0.4,
    aria: 'ClinGen gene curation expert panels with the most curations',
    describe: (r) => tipLines(`${exact(r.value)} curations`, r.key),
  });
  attachData('#validity-gcep-chart', 'clingen-expert-panels', ['Expert panel', 'Curations'], V.gceps.map((g) => [g.name, g.n]));
};

const renderCnvRegions = () => {
  const RG = snapshot.regions;
  const sizes = RG.sensitive.map((r) => r.end - r.start + 1);
  const both = RG.sensitive.filter((r) => r.hi === 3 && r.ts === 3).length;
  $('#cnvregions-dek').textContent = `${exact(RG.sensitive.length)} regions have sufficient evidence for at least one direction, ${exact(both)} for both a loss and a gain; ${exact(RG.sensitive.filter((r) => r.recurrent).length)} carry ClinGen’s recurrent region name. They span ${bases(Math.min(...sizes))} to ${bases(Math.max(...sizes))}. ${exact(RG.hi[40] ?? 0)} regions are scored as unlikely to be haploinsufficient.`;
  ideogram({
    mount: $('#regions-ideogram-chart'),
    chromosomes: NUCLEAR,
    bands: A.bands,
    marks: RG.sensitive.map((r) => ({
      chrom: r.chrom, start: r.start - 1, end: r.end, label: r.id,
    })),
    aria: 'ClinGen dosage-sensitive regions on the GRCh38 chromosomes',
    describe: (band) => tipLines(`${band.chrom}${band.name}`, coords(band.chrom, band.start + 1, band.end), ...RG.sensitive.filter((r) => r.chrom === band.chrom && r.start <= band.end && r.end > band.start).map((r) => r.name)),
  });
  attachData('#regions-ideogram-chart', 'clingen-dosage-regions', ['ID', 'Region', 'Band', 'Chromosome', 'Start', 'End', 'Haploinsufficiency', 'Triplosensitivity'], RG.sensitive.map((r) => [r.id, r.name, r.band, r.chrom, r.start, r.end, r.hi, r.ts]));
  const rows = RG.sensitive.map((r) => ({
    key: r.id, label: r.name, value: r.end - r.start + 1, hi: r.hi === 3 && r.ts === 3,
  })).sort((a, b) => b.value - a.value);
  dotplot({
    mount: $('#regions-size-chart'),
    rows,
    format: bases,
    tickFormat: bases,
    rowH: 20,
    labelShare: 0.5,
    labelMax: 380,
    aria: 'Size of each dosage-sensitive region, log scale',
    describe: (row) => {
      const r = RG.sensitive.find((x) => x.id === row.key);
      return tipLines(bases(row.value), r.name, `${r.band}, ${coords(r.chrom, r.start, r.end)}`, `Haploinsufficiency ${r.hi ?? 'not scored'}, triplosensitivity ${r.ts ?? 'not scored'}`, { href: `https://search.clinicalgenome.org/kb/regions/${r.id}`, text: 'Open in ClinGen' });
    },
  });
  attachData('#regions-size-chart', 'clingen-region-sizes', ['ID', 'Region', 'Bases', 'Haploinsufficiency', 'Triplosensitivity', 'Recurrent'], rows.map((row) => {
    const r = RG.sensitive.find((x) => x.id === row.key);
    return [r.id, r.name, row.value, r.hi, r.ts, r.recurrent ? 'yes' : 'no'];
  }));
};

const RANGE_KINDS = [
  { key: 'benign', label: 'Benign', color: 'var(--data-2)' },
  { key: 'intermediate', label: 'Intermediate', color: 'var(--s5)' },
  { key: 'pathogenic', label: 'Pathogenic', color: 'var(--accent)' },
];
const LOCATIONS = [["5' UTR", "5′ UTR or promoter", 'var(--s5)'], ['coding exon', 'Coding exon', 'var(--accent)'], ['exon', 'Unspecified exon', 'var(--s4)'], ['intron', 'Intron', 'var(--s2)'], ["3' UTR", "3′ UTR", 'var(--s3)'], ['non-coding', 'Non-coding RNA', 'var(--ink)'], ['unknown', 'Not stated', 'var(--data-2)']];

const renderExpansions = () => {
  const E = snapshot.expansions;
  const ranged = E.loci.filter((l) => ['Definitive', 'Strong'].includes(l.evidence) && l.pathogenic).sort((a, b) => a.pathogenic[0] - b.pathogenic[0] || a.gene.localeCompare(b.gene));
  const rows = ranged.map((l) => ({
    key: l.id,
    label: `${l.gene} ${l.motif ?? ''}`.trim(),
    segments: [['benign', l.benign], ['intermediate', l.intermediate], ['pathogenic', l.pathogenic]].filter(([, r]) => r).map(([kind, [lo, hi]]) => ({ kind, lo, hi })),
    locus: l,
  }));
  const contracts = (l) => l.benign && l.pathogenic[1] < l.benign[0];
  const grows = ranged.filter((l) => !contracts(l));
  const shrink = ranged.filter(contracts);
  const lows = grows[0];
  const highs = grows.at(-1);
  $('#expansions-dek').textContent = `${exact(ranged.length)} loci with definitive or strong evidence, from STRchive’s data${E.updated ? ` of ${longDate(E.updated)}` : ''}. Pathogenic expansions start at ${exact(lows.pathogenic[0])} repeats in ${lows.gene} and at ${exact(highs.pathogenic[0])} in ${highs.gene}${shrink.length ? `; at ${shrink.map((l) => l.gene).join(' and ')} the pathogenic allele is shorter than the benign one` : ''}. ${exact(E.loci.filter((l) => l.location === 'coding exon').length)} of the ${exact(E.loci.length)} loci sit in coding exons, ${exact(E.loci.filter((l) => l.location === 'intron').length)} in introns.`;
  ranges({
    mount: $('#expansions-chart'),
    rows,
    kinds: RANGE_KINDS,
    valueText: (r) => (contracts(r.locus) ? `≤${exact(r.locus.pathogenic[1])}` : `≥${exact(r.locus.pathogenic[0])}`),
    aria: 'Benign, intermediate and pathogenic repeat counts at each locus, log scale',
    describe: (r) => {
      const l = r.locus;
      return tipLines(`${l.gene} ${l.motif ?? ''}`, `${l.name} (${l.disease}), ${l.inheritance}`, ...r.segments.map((s) => `${s.kind}: ${exact(s.lo)} to ${exact(s.hi)} repeats`), l.onset ? `Typical onset at ${l.onset[0]} to ${l.onset[1]} years` : null, l.year ? `Found in ${l.year}` : null, l.genereviews ? { href: `https://www.ncbi.nlm.nih.gov/books/${l.genereviews}/`, text: 'GeneReviews chapter' } : null);
    },
  });
  attachData('#expansions-chart', 'repeat-expansion-ranges', ['Locus', 'Gene', 'Motif', 'Disease', 'Inheritance', 'Benign from', 'Benign to', 'Intermediate from', 'Intermediate to', 'Pathogenic from', 'Pathogenic to', 'Year'], ranged.map((l) => [l.id, l.gene, l.motif, l.name, l.inheritance, l.benign?.[0] ?? null, l.benign?.[1] ?? null, l.intermediate?.[0] ?? null, l.intermediate?.[1] ?? null, l.pathogenic[0], l.pathogenic[1], l.year]));
  const dated = E.loci.filter((l) => l.year);
  const years = Array.from({ length: Math.max(...dated.map((l) => l.year)) - Math.min(...dated.map((l) => l.year)) + 1 }, (_, i) => Math.min(...dated.map((l) => l.year)) + i);
  const present = LOCATIONS.filter(([key]) => dated.some((l) => l.location === key));
  stackedArea({
    mount: $('#expansions-years-chart'),
    years,
    series: present.map(([key, label, color]) => ({
      key, label, color, values: years.map((y) => dated.filter((l) => l.location === key && l.year <= y).length),
    })),
    aria: 'Repeat expansion loci found by each year, by where they sit in the gene',
  });
  attachData('#expansions-years-chart', 'repeat-expansion-years', ['Year', ...present.map(([, label]) => label)], years.map((y) => [y, ...present.map(([key]) => dated.filter((l) => l.location === key && l.year === y).length)]));
};

const systemName = (name) => {
  const s = name.replace(/^Abnormality of (the )?/, '');
  return `${s[0].toUpperCase()}${s.slice(1)}`;
};

const renderHpo = () => {
  const H = snapshot.hpo;
  const { genes } = H;
  $('#hpo-dek').textContent = `${exact(H.sources.OMIM ?? 0)} OMIM, ${exact(H.sources.ORPHA ?? 0)} Orphanet and ${exact(H.sources.DECIPHER ?? 0)} DECIPHER diseases carry ${exact(H.annotations)} phenotype annotations in the release of ${longDate(H.version)}; ${pct((H.evidence.IEA ?? 0) / H.annotations)} are inferred electronically. ${exact(genes.genes)} genes have a Mendelian link to ${exact(genes.diseases)} diseases; ${genes.top[0].gene} has the most, ${exact(genes.top[0].diseases)}.`;
  dotplot({
    mount: $('#hpo-systems-chart'),
    rows: H.systems.map((s, i) => ({
      key: s.id, label: systemName(s.name), value: s.diseases, hi: i === 0,
    })),
    scale: 'linear',
    rowH: 20,
    labelShare: 0.4,
    aria: 'Rare diseases with at least one phenotype in each organ system',
    describe: (r) => tipLines(`${exact(r.value)} diseases`, r.label, pct(r.value / H.diseases), { href: `https://hpo.jax.org/browse/term/${r.key}`, text: 'Open in HPO' }),
  });
  attachData('#hpo-systems-chart', 'hpo-organ-systems', ['HPO term', 'Organ system', 'Diseases'], H.systems.map((s) => [s.id, s.name, s.diseases]));
  dotplot({
    mount: $('#hpo-terms-chart'),
    rows: H.topTerms.map((t, i) => ({
      key: t.id, label: t.name, value: t.diseases, hi: i === 0,
    })),
    scale: 'linear',
    rowH: 20,
    labelShare: 0.4,
    aria: 'Phenotype terms annotated to the most rare diseases',
    describe: (r) => tipLines(`${exact(r.value)} diseases`, `${r.label}, ${r.key}`, pct(r.value / H.diseases), { href: `https://hpo.jax.org/browse/term/${r.key}`, text: 'Open in HPO' }),
  });
  attachData('#hpo-terms-chart', 'hpo-top-terms', ['HPO term', 'Name', 'Diseases'], H.topTerms.map((t) => [t.id, t.name, t.diseases]));
  const edges = H.perDisease.edges;
  const label = (i) => (edges[i + 1] - edges[i] === 1 ? String(edges[i]) : `${exact(edges[i])}–${exact(edges[i + 1] - 1)}`);
  histo({
    mount: $('#hpo-perdisease-chart'),
    bins: H.perDisease.counts.map((_, i) => i),
    series: [{ key: 'diseases', label: 'Diseases', color: 'var(--data)', values: H.perDisease.counts }],
    tickEvery: 1,
    tickLabel: (b) => exact(edges[b]),
    height: 200,
    aria: 'Rare diseases by the number of phenotype terms annotated',
    describe: (i) => tipLines(`${exact(H.perDisease.counts[i])} diseases`, `${label(i)} phenotype terms`),
  });
  attachData('#hpo-perdisease-chart', 'hpo-terms-per-disease', ['Terms from', 'Terms to', 'Diseases'], H.perDisease.counts.map((n, i) => [edges[i], edges[i + 1] - 1, n]));
  dotplot({
    mount: $('#hpo-genes-chart'),
    rows: genes.top.map((g, i) => ({
      key: g.gene, label: g.gene, value: g.diseases, hi: i === 0,
    })),
    scale: 'linear',
    rowH: 20,
    aria: 'Genes linked to the most Mendelian diseases',
    describe: (r) => tipLines(`${exact(r.value)} diseases`, r.key),
  });
  attachData('#hpo-genes-chart', 'hpo-genes', ['Gene', 'Mendelian diseases'], genes.top.map((g) => [g.gene, g.diseases]));
};

const PREVALENCE_ORDER = ['<1 / 1 000 000', '1-9 / 1 000 000', '1-9 / 100 000', '1-5 / 10 000', '6-9 / 10 000', '>1 / 1000', 'Unknown', 'Not yet documented', 'Case reports only', 'No point prevalence'];
const ONSET_ORDER = ['Antenatal', 'Neonatal', 'Infancy', 'Childhood', 'Adolescent', 'Adult', 'Elderly', 'All ages', 'No data available'];
const ordered = (counts, order) => [...order.filter((k) => k in counts), ...Object.keys(counts).filter((k) => !order.includes(k)).sort((a, b) => counts[b] - counts[a])];

const renderOrphanet = () => {
  const { genes, prevalence, ages } = snapshot.orphanet;
  $('#orphanet-dek').textContent = `Orphadata products of ${longDate(prevalence.date)}, counted at the level of single disorders, not groups or subtypes. ${exact(genes.genes)} genes are tied to ${exact(genes.disorders)} disorders, ${pct(genes.perDisorder.counts[0] / genes.disorders)} of them to one gene. ${exact(ages.noOnset)} disorders give no age of onset and ${exact(ages.noInheritance)} no inheritance; one disorder can list several of each.`;
  const plot = (sel, counts, order, aria, noun, name, head) => {
    const keys = ordered(counts, order);
    dotplot({
      mount: $(sel),
      rows: keys.map((k) => ({ key: k, label: k, value: counts[k] })),
      rowH: 20,
      labelShare: 0.42,
      aria,
      describe: (r) => tipLines(`${exact(r.value)} ${noun}`, r.label),
    });
    attachData(sel, name, [head, noun[0].toUpperCase() + noun.slice(1)], keys.map((k) => [k, counts[k]]));
  };
  plot('#orpha-prevalence-chart', prevalence.classes, PREVALENCE_ORDER, 'Orphanet disorders by point prevalence class, log scale', 'disorders', 'orphanet-prevalence', 'Point prevalence');
  plot('#orpha-onset-chart', ages.onset, ONSET_ORDER, 'Orphanet disorders by age of onset, log scale', 'disorders', 'orphanet-onset', 'Age of onset');
  plot('#orpha-inheritance-chart', ages.inheritance, [], 'Orphanet disorders by type of inheritance, log scale', 'disorders', 'orphanet-inheritance', 'Inheritance');
  plot('#orpha-genes-chart', genes.associations, [], 'Orphanet gene-disorder associations by type, log scale', 'associations', 'orphanet-gene-associations', 'Association');
};

const GENE_COUNT_LABELS = ['No gene listed', 'One gene', '2 to 10 genes', '11 to 100 genes', '101 to 1,000 genes', 'More than 1,000 genes'];

const renderGtr = () => {
  const GTRS = snapshot.gtr;
  const dev = GTRS.development;
  const lab = dev['Test developed by laboratory (no manufacturer test name)'] ?? 0;
  const [first, second] = GTRS.countries;
  $('#gtr-dek').textContent = `${exact(dev['FDA-reviewed (has FDA test name)'] ?? 0)} tests are FDA-reviewed and ${exact(lab)} were developed by the laboratory; ${exact(dev['not stated'] ?? 0)} do not say. ${first.country} hosts ${exact(first.tests)} tests from ${exact(first.labs)} laboratories, ${second.country} ${exact(second.tests)} from ${exact(second.labs)}.`;
  const methods = GTRS.methods.slice(0, 15);
  dotplot({
    mount: $('#gtr-methods-chart'),
    rows: methods.map((m, i) => ({
      key: m.method, label: m.method, value: m.tests, hi: i === 0,
    })),
    rowH: 20,
    labelShare: 0.45,
    aria: 'Current GTR tests by method, log scale',
    describe: (r) => tipLines(`${exact(r.value)} tests`, r.label, pct(r.value / GTRS.tests)),
  });
  attachData('#gtr-methods-chart', 'gtr-methods', ['Method', 'Tests'], GTRS.methods.map((m) => [m.method, m.tests]));
  dotplot({
    mount: $('#gtr-genecount-chart'),
    rows: GTRS.genesPerTest.counts.map((n, i) => ({ key: GTRS.genesPerTest.bins[i], label: GENE_COUNT_LABELS[i], value: n })),
    rowH: 20,
    labelShare: 0.36,
    aria: 'Current GTR tests by the number of genes they cover, log scale',
    describe: (r) => tipLines(`${exact(r.value)} tests`, r.label, pct(r.value / GTRS.tests)),
  });
  attachData('#gtr-genecount-chart', 'gtr-genes-per-test', ['Genes', 'Tests'], GTRS.genesPerTest.counts.map((n, i) => [GENE_COUNT_LABELS[i], n]));
  dotplot({
    mount: $('#gtr-countries-chart'),
    rows: GTRS.countries.map((c, i) => ({
      key: c.country, label: `${c.country} (${exact(c.labs)})`, value: c.tests, hi: i === 0,
    })),
    rowH: 20,
    labelShare: 0.36,
    aria: 'Current GTR tests by the country of the laboratory, log scale',
    describe: (r) => {
      const c = GTRS.countries.find((x) => x.country === r.key);
      return tipLines(`${exact(c.tests)} tests`, `${c.country}, ${exact(c.labs)} laboratories`);
    },
  });
  attachData('#gtr-countries-chart', 'gtr-countries', ['Country', 'Laboratories', 'Tests'], GTRS.countries.map((c) => [c.country, c.labs, c.tests]));
  dotplot({
    mount: $('#gtr-genes-chart'),
    rows: GTRS.topGenes.map((g, i) => ({
      key: g.gene, label: g.gene, value: g.tests, hi: i === 0,
    })),
    scale: 'linear',
    rowH: 20,
    aria: 'Genes covered by the most current GTR tests',
    describe: (r) => tipLines(`${exact(r.value)} tests`, r.key, { href: `https://www.ncbi.nlm.nih.gov/gtr/all/tests/?term=${encodeURIComponent(r.key)}`, text: 'Open in GTR' }),
  });
  attachData('#gtr-genes-chart', 'gtr-genes', ['Gene', 'Tests'], GTRS.topGenes.map((g) => [g.gene, g.tests]));
};

const THERAPY_GROUPS = [['AAV', 'var(--accent)'], ['Other in vivo vectors', 'var(--s5)'], ['Blood stem cells', 'var(--s2)'], ['CAR-T and TCR-T cells', 'var(--s3)'], ['Other gene-modified cells', 'var(--s4)']];
const TRIAL_COLORS = ['var(--accent)', 'var(--s2)', 'var(--s5)', 'var(--s3)', 'var(--s4)', 'var(--ink)', 'var(--data-2)'];
const PHASE_ORDER = [['Early phase 1', 'var(--data-2)'], ['Phase 1', 'var(--s3)'], ['Phase 2', 'var(--s2)'], ['Phase 3', 'var(--accent)'], ['Phase 4', 'var(--ink)'], ['Not applicable', 'var(--s5)'], ['Not stated', 'var(--s4)']];

const renderTherapy = () => {
  const T = { fda: snapshot.fda, trials: snapshot.trials };
  const listed = new Set(T.fda.map((p) => p.trade));
  const present = GENE_THERAPIES.filter((g) => listed.has(g.trade));
  const others = T.fda.length - present.length;
  const now = Number(snapshot.generated.slice(0, 4));
  const first = Math.min(...present.map((g) => g.year));
  const years = Array.from({ length: now - first + 1 }, (_, i) => first + i);
  const cart = T.trials.cart;
  const later = sumOf(Object.values(T.trials).map((t) => sumOf(Object.entries(t.byYear).filter(([y]) => Number(y) > now).map(([, n]) => n))));
  $('#therapy-dek').textContent = `The ${exact(present.length)} products here are FDA’s gene therapies and gene-modified cells; the first approval year and vector of each were read from its package insert on ${longDate(GENE_THERAPIES_CHECKED)}. FDA’s list holds ${exact(others)} more: cord blood, tissue and cell products without genetic modification, and any approved since then. ${exact(cart.total)} interventional studies mention CAR T cells, ${exact(cart.recruiting)} of them recruiting.`;
  stackedArea({
    mount: $('#therapy-approvals-chart'),
    years,
    series: THERAPY_GROUPS.map(([key, color]) => ({
      key, label: key, color, values: years.map((y) => present.filter((g) => g.group === key && g.year <= y).length),
    })),
    aria: 'FDA-approved gene therapies and gene-modified cell products, cumulative, by kind',
  });
  attachData('#therapy-approvals-chart', 'fda-gene-therapies', ['Product', 'Proper name', 'First approval', 'Vector', 'Kind', 'Condition'], present.map((g) => [g.trade, g.proper, g.year, g.vector, g.group, g.condition]));
  const vectors = new Map();
  for (const g of present) vectors.set(g.vector, (vectors.get(g.vector) ?? 0) + 1);
  const vectorRows = [...vectors].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]));
  dotplot({
    mount: $('#therapy-vectors-chart'),
    rows: vectorRows.map(([v, n]) => ({ key: v, label: v, value: n })),
    scale: 'linear',
    rowH: 20,
    labelShare: 0.45,
    format: (v) => String(v),
    tickFormat: (t) => String(t),
    aria: 'FDA-approved gene therapies by vector or editing method',
    describe: (r) => tipLines(`${exact(r.value)} products`, r.label, present.filter((g) => g.vector === r.key).map((g) => `${g.trade} ${g.year}, ${g.condition}`).join('; ')),
  });
  attachData('#therapy-vectors-chart', 'fda-gene-therapy-vectors', ['Vector', 'Products'], vectorRows);
  const startYears = Array.from({ length: now - 2000 + 1 }, (_, i) => 2000 + i);
  lines({
    mount: $('#therapy-trials-chart'),
    series: TRIAL_MODALITIES.map((m, i) => ({
      label: m.label, color: TRIAL_COLORS[i], rows: startYears.map((y) => ({ x: y, y: T.trials[m.key].byYear[y] ?? 0 })),
    })),
    scale: 'log',
    logStep: 1,
    height: 320,
    labelW: 170,
    aria: 'Interventional studies starting each year, by the gene therapy modality their intervention names, log scale',
    describe: (y) => tipLines(String(y), 'studies starting that year', ...TRIAL_MODALITIES.map((m) => `${m.label}: ${exact(T.trials[m.key].byYear[y] ?? 0)}`)),
  });
  attachData('#therapy-trials-chart', 'trials-by-modality', ['Start year', ...TRIAL_MODALITIES.map((m) => m.label)], startYears.map((y) => [y, ...TRIAL_MODALITIES.map((m) => T.trials[m.key].byYear[y] ?? 0)]));
  $('#therapy-trials-note').textContent = `${exact(later)} studies give a planned start after ${now} and are left out.`;
  hbars({
    mount: $('#therapy-phases-chart'),
    rows: TRIAL_MODALITIES.map((m) => {
      const t = T.trials[m.key];
      return {
        key: m.key, label: m.label, n: t.total, values: Object.fromEntries(PHASE_ORDER.map(([p]) => [p, ((t.phases[p] ?? 0) / t.total) * 100])),
      };
    }),
    segments: PHASE_ORDER.map(([key, color]) => ({ key, label: key, color })),
    max: 100,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => exact(r.n),
    labelShare: 0.3,
    aria: 'Interventional studies by their highest phase, as a share of each modality',
    describe: (r) => {
      const t = T.trials[r.key];
      return tipLines(`${exact(t.total)} studies`, r.label, ...PHASE_ORDER.filter(([p]) => t.phases[p]).map(([p]) => `${p}: ${exact(t.phases[p])}`), `${exact(t.recruiting)} recruiting`);
    },
  });
  attachData('#therapy-phases-chart', 'trials-by-phase', ['Modality', ...PHASE_ORDER.map(([p]) => p), 'Recruiting'], TRIAL_MODALITIES.map((m) => [m.label, ...PHASE_ORDER.map(([p]) => T.trials[m.key].phases[p] ?? 0), T.trials[m.key].recruiting]));
};

const STUDY = { miller2010: 'Miller 2010', clark2018: 'Clark 2018', srivastava2019: 'Srivastava 2019' };
const SHORT_TEST = {
  'G-banded karyotype': 'Karyotype', 'Chromosomal microarray': 'Microarray', 'Exome sequencing': 'Exome', 'Genome sequencing': 'Genome',
};
const shortTest = (test) => test.replace(/^(G-banded karyotype|Chromosomal microarray|Exome sequencing|Genome sequencing)/, (m) => SHORT_TEST[m]);
const MEASURES = [{ key: 'diagnostic', label: 'Diagnostic yield', color: 'var(--accent)' }, { key: 'clinical', label: 'Changed clinical management', color: 'var(--s2)' }];
const pctRange = (lo, hi) => `${lo}–${hi}%`;

const renderYield = () => {
  $('#yield-dek').textContent = 'Clark and colleagues pooled 37 studies of 20,068 children; Miller and colleagues 33 studies of 21,698 people with developmental delay, intellectual disability, autism or congenital anomalies; Srivastava and colleagues 30 articles on neurodevelopmental disorders (NDD). The karyotype figure leaves out Down syndrome and other recognizable chromosomal syndromes.';
  intervals({
    mount: $('#yield-chart'),
    rows: DIAGNOSTIC_YIELD.map((y) => ({
      ...y, kind: y.measure, label: `${shortTest(y.test)} (${STUDY[y.source].slice(-4)})`,
    })),
    domain: [0, 60],
    ticks: [0, 10, 20, 30, 40, 50, 60],
    kinds: MEASURES,
    valueText: (r) => (r.value === null ? pctRange(r.lo, r.hi) : `${r.value}%${r.lo === null ? '' : ` (${pctRange(r.lo, r.hi)})`}`),
    aria: 'Diagnostic yield and clinical utility of genetic tests in meta-analyses, with 95% confidence intervals',
    describe: (r) => tipLines(r.value === null ? pctRange(r.lo, r.hi) : `${r.value}%`, `${r.test}: ${r.measure === 'diagnostic' ? 'diagnostic yield' : 'changed clinical management'}`, r.lo !== null && r.value !== null ? `95% CI ${pctRange(r.lo, r.hi)}` : null, { href: SOURCES[r.source], text: STUDY[r.source] }),
  });
  attachData('#yield-chart', 'diagnostic-yield', ['Test', 'Measure', 'Percent', 'Low', 'High', 'Source'], DIAGNOSTIC_YIELD.map((y) => [y.test, y.measure, y.value, y.lo, y.hi, STUDY[y.source]]));
};

const renderPrenatal = () => {
  $('#prenatal-dek').textContent = `Gil and colleagues pooled 35 studies of singleton pregnancies. Each false-positive rate is shown beside its detection rate. ACMG’s 2021 guidance offers carrier screening to every pregnancy for conditions with a carrier frequency of at least ${CARRIER_TIER3.frequency}: ${exact(CARRIER_TIER3.genes)} autosomal recessive genes, ${exact(CARRIER_TIER3.definitive)} of them definitive in ClinGen when published, plus X-linked genes.`;
  intervals({
    mount: $('#prenatal-chart'),
    rows: CFDNA_SCREENING.map((s) => ({
      key: s.condition, label: s.condition, value: s.detection, lo: s.lo, hi: s.hi, kind: 'dr', screen: s,
    })),
    domain: [60, 100],
    ticks: [60, 70, 80, 90, 100],
    kinds: [{ key: 'dr', label: 'Detection rate', color: 'var(--accent)' }],
    valueText: (r) => `${r.value}%, FPR ${r.screen.fpr}%`,
    aria: 'Detection rates of cell-free DNA screening with 95% confidence intervals',
    describe: (r) => {
      const s = r.screen;
      return tipLines(`${s.detection}% detected`, s.condition, `95% CI ${pctRange(s.lo, s.hi)}`, `False-positive rate ${s.fpr}% (${pctRange(s.fprLo, s.fprHi)})`, `${exact(s.cases)} affected and ${exact(s.unaffected)} unaffected pregnancies`, { href: SOURCES.gil2017, text: 'Gil 2017' });
    },
  });
  attachData('#prenatal-chart', 'cfdna-screening', ['Condition', 'Detection rate', 'Low', 'High', 'False-positive rate', 'FPR low', 'FPR high', 'Affected', 'Unaffected'], CFDNA_SCREENING.map((s) => [s.condition, s.detection, s.lo, s.hi, s.fpr, s.fprLo, s.fprHi, s.cases, s.unaffected]));
};

const DRIVER_ROLES = [['Act', 'Activating', 'var(--s1)'], ['LoF', 'Loss of function', 'var(--s2)'], ['ambiguous', 'No majority role', 'var(--data-2)']];
const AGE_GROUPS = [{ key: 'Adult', label: 'Adult cohorts', color: 'var(--s2)' }, { key: 'Pediatric', label: 'Pediatric cohorts', color: 'var(--s1)' }];
const perSampleText = (v) => (v < 10 ? v.toFixed(1) : exact(v));

const renderDrivers = () => {
  const { drivers: drv, cohorts: coh, release } = snapshot.intogen;
  const nameOf = (code) => coh.names[code] ?? code;
  const byCancer = Object.fromEntries(coh.burden.map((b) => [b.cancer, { cohorts: b.cohorts.length, samples: sumOf(b.cohorts.map((c) => c.samples)) }]));
  const lead = drv.byType[0];
  $('#drivers-dek').textContent = `IntOGen’s release of ${longDate(release)} covers ${exact(coh.cohorts)} cohorts and ${exact(coh.samples)} tumour samples. Of its ${exact(drv.genes)} driver genes, ${exact(drv.roles.Act ?? 0)} act mostly through activating mutations, ${exact(drv.roles.LoF ?? 0)} through loss of function, and ${exact(drv.roles.ambiguous ?? 0)} have no majority role. ${nameOf(lead.type)} has the most drivers, ${exact(lead.genes)}.`;
  hbars({
    mount: $('#drivers-types-chart'),
    rows: drv.byType.slice(0, 30).map((t) => ({ key: t.type, label: nameOf(t.type), values: { d: t.genes } })),
    segments: [{ key: 'd', label: 'Driver genes', color: 'var(--data)' }],
    labelShare: 0.4,
    labelMax: 380,
    rowH: 18,
    aria: 'Cancer types with the most IntOGen driver genes',
    describe: (r) => tipLines(`${exact(r.values.d)} driver genes`, r.label, `${exact(byCancer[r.key]?.cohorts ?? 0)} cohorts, ${exact(byCancer[r.key]?.samples ?? 0)} samples`),
  });
  attachData('#drivers-types-chart', 'intogen-drivers-by-cancer', ['Code', 'Cancer type', 'Driver genes', 'Cohorts', 'Samples'], drv.byType.map((t) => [t.type, nameOf(t.type), t.genes, byCancer[t.type]?.cohorts ?? 0, byCancer[t.type]?.samples ?? 0]));
  hbars({
    mount: $('#drivers-genes-chart'),
    rows: drv.top.map((g) => ({
      key: g.gene, label: g.gene, values: { Act: g.Act, LoF: g.LoF, ambiguous: g.ambiguous }, gene: g,
    })),
    segments: DRIVER_ROLES.map(([key, label, color]) => ({ key, label, color })),
    labelShare: 0.12,
    rowH: 18,
    aria: 'Genes that IntOGen calls drivers in the most cancer types, by role',
    describe: (r) => tipLines(`${exact(r.gene.types)} cancer types`, r.gene.gene, ...DRIVER_ROLES.map(([key, label]) => `${label}: ${exact(r.gene[key])}`)),
  });
  attachData('#drivers-genes-chart', 'intogen-top-drivers', ['Gene', 'Cancer types', ...DRIVER_ROLES.map(([, label]) => label)], drv.top.map((g) => [g.gene, g.types, g.Act, g.LoF, g.ambiguous]));
};

const renderTcga = () => {
  const gd = snapshot.gdc;
  const projects = gd.tcga.projects.filter((p) => p.tested > 0);
  const { genes } = gd.tcga;
  const casesOf = (gene, p) => gd.tcga.mutated[gene]?.[p.project] ?? 0;
  const values = genes.map((gene) => projects.map((p) => casesOf(gene, p) / p.tested));
  const tested = sumOf(projects.map((p) => p.tested));
  const across = (gene) => sumOf(projects.map((p) => casesOf(gene, p))) / tested;
  $('#tcga-dek').textContent = `GDC Data Release ${gd.release} of ${longDate(gd.date)}: ${exact(tested)} cases in ${exact(projects.length)} TCGA projects have exome data. Across them, ${genes[0]} carries a somatic mutation in ${pct(across(genes[0]))} of cases and ${genes[1]} in ${pct(across(genes[1]))}. The rows are the ${exact(genes.length)} genes IntOGen calls drivers in the most cancer types.`;
  shareGrid({
    mount: $('#tcga-chart'),
    rows: genes,
    cols: projects.map((p) => p.project.replace(/^TCGA-/, '')),
    values,
    aria: 'Share of TCGA cases with a somatic mutation in each driver gene, by project',
    describe: (i, j) => tipLines(`${genes[i]}: ${pct(values[i][j])}`, `of ${projects[j].name} cases, ${projects[j].project}`, `${exact(casesOf(genes[i], projects[j]))} of ${exact(projects[j].tested)} cases with exome data`),
  });
  attachData('#tcga-chart', 'tcga-driver-mutations', ['Gene', 'Project', 'Cases with a mutation', 'Cases with exome data'], genes.flatMap((gene) => projects.map((p) => [gene, p.project, casesOf(gene, p), p.tested])));
};

const renderBurden = () => {
  const coh = snapshot.intogen.cohorts;
  const rowsFor = (platform) => coh.burden.map((b) => {
    const list = b.cohorts.filter((c) => c.platform === platform);
    return list.length ? {
      key: b.cancer, label: b.name, median: summary(list.map((c) => c.perSample)).median, points: list.map((c) => ({ value: Math.max(c.perSample, 0.1), kind: c.age, cohort: c })),
    } : null;
  }).filter(Boolean).sort((x, y) => y.median - x.median || x.label.localeCompare(y.label));
  const exomes = rowsFor('WXS');
  const genomes = rowsFor('WGS');
  $('#burden-dek').textContent = `Each dot is one cohort’s mutations divided by its samples, and the bar marks the median cohort; rows run from the highest median down. Exome cohorts count coding mutations and genome cohorts count mutations anywhere in the genome, so the two are ranked apart: ${exact(exomes.length)} cancer types have exome cohorts and ${exact(genomes.length)} have genome cohorts. ${exact(coh.ages.Pediatric ?? 0)} of the ${exact(coh.cohorts)} cohorts are pediatric and ${exact(coh.kinds.Metastasis ?? 0)} are metastases.`;
  for (const [id, rows, platform, word] of [['exome', exomes, 'WXS', 'exome'], ['genome', genomes, 'WGS', 'genome']]) {
    dotRows({
      mount: $(`#burden-${id}-chart`),
      rows,
      kinds: AGE_GROUPS,
      labelShare: 0.38,
      labelMax: 420,
      rowH: 14,
      aria: `Mutations per sample in each IntOGen ${word} cohort, by cancer type, log scale`,
      describe: (row, p) => tipLines(`${perSampleText(p.cohort.perSample)} mutations per sample`, p.cohort.cohort, `${row.label}; ${word}, ${exact(p.cohort.samples)} samples, ${p.cohort.age.toLowerCase()}`),
    });
    attachData(`#burden-${id}-chart`, `intogen-${word}-burden`, ['Cohort', 'Cancer type', 'Age group', 'Samples', 'Mutations per sample'], coh.burden.flatMap((b) => b.cohorts.filter((c) => c.platform === platform).map((c) => [c.cohort, b.name, c.age, c.samples, c.perSample])));
  }
};

const ONCOGENICITY = ['Oncogenic', 'Oncogenic/Likely oncogenic', 'Likely oncogenic', 'Uncertain significance', 'Likely benign', 'Benign', 'Conflicting classifications of oncogenicity'];

const renderSomatic = () => {
  const gd = snapshot.gdc;
  const onco = CV.oncogenicity;
  const impact = CV.clinicalImpact;
  const largest = gd.programs[0];
  $('#somatic-dek').textContent = `GDC Data Release ${gd.release}: ${exact(gd.cases)} cases in ${exact(gd.projects)} projects from ${exact(gd.programs.length)} programs; ${largest.program} is the largest, with ${exact(largest.cases)}. In ClinVar, ${exact(CV.origin.somatic ?? 0)} variants have only a somatic origin. Its oncogenicity classes follow the ClinGen, CGC and VICC standard, and its somatic clinical impact uses the AMP/ASCO/CAP tiers; ${exact(impact['Tier I - Strong'] ?? 0)} variants reach Tier I.`;
  dotplot({
    mount: $('#somatic-programs-chart'),
    rows: gd.programs.map((p, i) => ({
      key: p.program, label: p.program, value: p.cases, hi: i === 0,
    })),
    labelShare: 0.3,
    rowH: 20,
    aria: 'Cases at the NCI Genomic Data Commons by program, log scale',
    describe: (r) => {
      const p = gd.programs.find((x) => x.program === r.key);
      return tipLines(`${exact(p.cases)} cases`, `${p.program}, ${exact(p.projects)} project${p.projects === 1 ? '' : 's'}`, pct(p.cases / gd.cases));
    },
  });
  attachData('#somatic-programs-chart', 'gdc-programs', ['Program', 'Projects', 'Cases'], gd.programs.map((p) => [p.program, p.projects, p.cases]));
  dotplot({
    mount: $('#somatic-sites-chart'),
    rows: gd.sites.slice(0, 20).map((s, i) => ({
      key: s.key, label: s.key, value: s.n, hi: i === 0,
    })),
    labelShare: 0.42,
    labelMax: 300,
    rowH: 20,
    aria: 'Cases at the NCI Genomic Data Commons by primary site, log scale',
    describe: (r) => tipLines(`${exact(r.value)} cases`, r.key, pct(r.value / gd.cases)),
  });
  attachData('#somatic-sites-chart', 'gdc-primary-sites', ['Primary site', 'Cases'], gd.sites.map((s) => [s.key, s.n]));
  const oncoKeys = [...ONCOGENICITY.filter((k) => onco[k]), ...Object.keys(onco).filter((k) => !ONCOGENICITY.includes(k))];
  dotplot({
    mount: $('#somatic-oncogenicity-chart'),
    rows: oncoKeys.map((k) => ({ key: k, label: k, value: onco[k] })),
    labelShare: 0.45,
    labelMax: 300,
    rowH: 20,
    aria: 'ClinVar somatic variants by oncogenicity classification, log scale',
    describe: (r) => tipLines(`${exact(r.value)} variants`, r.key, pct(r.value / sumOf(Object.values(onco)))),
  });
  attachData('#somatic-oncogenicity-chart', 'clinvar-oncogenicity', ['Oncogenicity', 'Variants'], oncoKeys.map((k) => [k, onco[k]]));
  const tiers = Object.keys(impact).sort();
  dotplot({
    mount: $('#somatic-impact-chart'),
    rows: tiers.map((k) => ({ key: k, label: k, value: impact[k] })),
    labelShare: 0.45,
    labelMax: 300,
    rowH: 20,
    aria: 'ClinVar somatic variants by clinical impact tier, log scale',
    describe: (r) => tipLines(`${exact(r.value)} variants`, r.key, pct(r.value / sumOf(Object.values(impact)))),
  });
  attachData('#somatic-impact-chart', 'clinvar-somatic-impact', ['Clinical impact', 'Variants'], tiers.map((k) => [k, impact[k]]));
};

const CPIC_LEVEL_COLORS = {
  A: 'var(--s1)', 'A/B': 'var(--s5)', B: 'var(--s3)', 'B/C': 'var(--s4)', C: 'var(--s2)', 'C/D': 'var(--data)', D: 'var(--data-2)', Retired: 'var(--ink-2)',
};
const CPIC_ORDER = Object.keys(CPIC_LEVEL_COLORS);
const CLINPGX_ORDER = ['1A', '1B', '2A', '2B', '3', '4'];
const CLINPGX_COLORS = {
  '1A': 'var(--s1)', '1B': 'var(--s5)', '2A': 'var(--s3)', '2B': 'var(--s4)', 3: 'var(--s2)', 4: 'var(--data-2)',
};
const FUNCTION_COLORS = [['Normal function', 'var(--s4)'], ['Increased function', 'var(--s2)'], ['Decreased function', 'var(--s3)'], ['No function', 'var(--s1)'], ['Uncertain function', 'var(--data)'], ['Unknown function', 'var(--data-2)'], ['Gene-specific class', 'var(--s5)'], ['Not assigned', 'var(--ink-2)']];
const METABOLIZERS = [['Poor Metabolizer', 'var(--s1)'], ['Intermediate Metabolizer', 'var(--s3)'], ['Normal Metabolizer', 'var(--s4)'], ['Ultrarapid Metabolizer', 'var(--s2)'], ['Indeterminate', 'var(--data-2)']];
const TESTING_LEVELS = [['Testing Required', 'var(--s1)'], ['Testing Recommended', 'var(--s5)'], ['Actionable PGx', 'var(--s3)'], ['Informative PGx', 'var(--s2)'], ['No Clinical PGx', 'var(--data)'], ['No level given', 'var(--data-2)']];

const renderCpic = () => {
  const { pairs, alleles } = snapshot.cpic;
  const levels = ordered(pairs.levels, CPIC_ORDER);
  $('#cpic-dek').textContent = `CPIC lists ${exact(pairs.pairs)} gene–drug pairs for ${exact(pairs.genes)} genes and ${exact(pairs.drugs)} drugs, not counting ${exact(pairs.removed)} removed pairs. At level A genotype should change prescribing and at level B it could; at levels C and D no prescribing action is recommended. FDA label annotations call testing required for ${exact(pairs.testing['Testing Required'] ?? 0)} of the pairs and recommended for ${exact(pairs.testing['Testing Recommended'] ?? 0)}. CPIC curates ${exact(alleles.alleles)} alleles in ${exact(alleles.genes.length)} genes.`;
  hbars({
    mount: $('#cpic-genes-chart'),
    rows: pairs.byGene.slice(0, 30).map((g) => ({
      key: g.gene, label: g.gene, values: g.levels, gene: g,
    })),
    segments: levels.map((l) => ({ key: l, label: `Level ${l}`, color: CPIC_LEVEL_COLORS[l] ?? 'var(--data)' })),
    labelShare: 0.12,
    rowH: 18,
    aria: 'CPIC gene–drug pairs per gene by CPIC level',
    describe: (r) => tipLines(`${exact(r.gene.total)} pairs`, r.gene.gene, ...levels.filter((l) => r.gene.levels[l]).map((l) => `Level ${l}: ${exact(r.gene.levels[l])}`)),
  });
  attachData('#cpic-genes-chart', 'cpic-pairs-by-gene', ['Gene', ...levels.map((l) => `Level ${l}`), 'Total'], pairs.byGene.map((g) => [g.gene, ...levels.map((l) => g.levels[l] ?? 0), g.total]));
  const cols = [...CLINPGX_ORDER, 'none'];
  const rowTotals = levels.map((l) => sumOf(Object.values(pairs.clinpgx[l] ?? {})));
  const counts = levels.map((l) => cols.map((c) => pairs.clinpgx[l]?.[c] ?? 0));
  shareGrid({
    mount: $('#cpic-levels-chart'),
    rows: levels.map((l) => `Level ${l}`),
    cols: [...CLINPGX_ORDER, 'None'],
    values: counts.map((row, i) => row.map((n) => (rowTotals[i] ? n / rowTotals[i] : 0))),
    totals: rowTotals,
    aria: 'CPIC levels against ClinPGx levels of evidence, as a share of each CPIC level',
    describe: (i, j) => tipLines(pct(rowTotals[i] ? counts[i][j] / rowTotals[i] : 0), `of level ${levels[i]} pairs have ${cols[j] === 'none' ? 'no ClinPGx level' : `ClinPGx level ${cols[j]}`}`, `${exact(counts[i][j])} of ${exact(rowTotals[i])}`),
  });
  attachData('#cpic-levels-chart', 'cpic-and-clinpgx-levels', ['CPIC level', ...cols.map((c) => (c === 'none' ? 'No ClinPGx level' : `ClinPGx ${c}`)), 'Pairs'], levels.map((l, i) => [l, ...counts[i], rowTotals[i]]));
  const fnKeys = FUNCTION_COLORS.filter(([f]) => alleles.genes.some((g) => g.functions[f]));
  hbars({
    mount: $('#cpic-alleles-chart'),
    rows: alleles.genes.map((g) => ({
      key: g.gene, label: g.gene, values: g.functions, gene: g,
    })),
    segments: fnKeys.map(([key, color]) => ({ key, label: key, color })),
    labelShare: 0.12,
    rowH: 18,
    aria: 'Alleles CPIC curates per gene, by clinical function',
    describe: (r) => tipLines(`${exact(r.gene.total)} alleles`, r.gene.gene, ...fnKeys.filter(([f]) => r.gene.functions[f]).map(([f]) => `${f}: ${exact(r.gene.functions[f])}`), r.gene.classes.length ? `Gene-specific classes: ${r.gene.classes.join('; ')}` : null),
  });
  attachData('#cpic-alleles-chart', 'cpic-alleles-by-function', ['Gene', ...fnKeys.map(([f]) => f), 'Gene-specific classes', 'Total'], alleles.genes.map((g) => [g.gene, ...fnKeys.map(([f]) => g.functions[f] ?? 0), g.classes.join('; '), g.total]));
};

const renderPhenotypes = () => {
  const { groups, genes } = snapshot.cpic.phenotypes;
  const values = groups.map((g) => genes.map((x) => x.shares[g] ?? null));
  const allHigh = genes.filter((x) => Object.values(x.shares).every((v) => v > 0.9));
  const flagged = allHigh.map((x) => `${x.gene} is above 90% in every group listed; CPIC flags ${exact(x.highRisk.length)} of its ${exact(x.results)} phenotypes.`).join(' ');
  $('#phenotypes-dek').textContent = `CPIC flags each predicted phenotype for electronic health records; each cell adds the frequencies of the phenotypes flagged as abnormal, priority or high risk. The ${exact(groups.length)} biogeographic groups are CPIC’s. ${flagged}${flagged ? ' ' : ''}Dashed cells have no frequency in CPIC’s tables.`;
  shareGrid({
    mount: $('#phenotypes-chart'),
    rows: groups,
    cols: genes.map((x) => x.gene),
    values,
    color: 'var(--accent)',
    aria: 'Share of each CPIC biogeographic group predicted to have a high-risk phenotype, by gene',
    describe: (i, j) => {
      const x = genes[j];
      const v = values[i][j];
      return v === null ? tipLines(x.gene, `No ${groups[i]} frequency in CPIC’s tables`) : tipLines(`${x.gene}: ${pct(v)}`, `of ${groups[i]} people predicted high risk`, `Counted: ${x.highRisk.join(', ')}`);
    },
  });
  attachData('#phenotypes-chart', 'cpic-high-risk-phenotypes', ['Group', ...genes.map((x) => x.gene)], groups.map((g, i) => [g, ...values[i].map((v) => v ?? '')]));
};

const renderCyp2d6 = () => {
  const { groups } = snapshot.cpic.phenotypes;
  const cy = snapshot.cpic.cyp2d6;
  const functionColor = Object.fromEntries(FUNCTION_COLORS);
  const known = (freqs) => Object.entries(freqs).filter(([, f]) => f !== null);
  const peaks = cy.top.filter((a) => a.function !== 'Normal function' && known(a.freqs).length).map((a) => {
    const [group, v] = known(a.freqs).reduce((x, y) => (y[1] > x[1] ? y : x));
    return { name: a.name, group, v };
  }).sort((x, y) => y.v - x.v).slice(0, 3);
  $('#cyp2d6-dek').textContent = `Each of the ${exact(cy.alleles)} CYP2D6 alleles in CPIC’s tables has a function and an activity value; a person’s activity score adds the values of both alleles and sets the predicted phenotype: 0 poor, 0.25 to 1 intermediate, 1.25 to 2.25 normal, above 2.25 ultrarapid. ${exact(cy.functions['No function'] ?? 0)} alleles have no function. ${peaks.map((p) => `${p.name} reaches ${pct(p.v)} in the ${p.group} group`).join('; ')}.`;
  const alleleValues = groups.map((g) => cy.top.map((a) => a.freqs[g] ?? null));
  shareGrid({
    mount: $('#cyp2d6-alleles-chart'),
    rows: groups,
    cols: cy.top.map((a) => a.name),
    values: alleleValues,
    max: Math.max(...alleleValues.flat().filter((v) => v !== null)),
    colKeys: cy.top.map((a) => functionColor[a.function] ?? 'var(--data)'),
    legend: FUNCTION_COLORS.filter(([f]) => cy.top.some((a) => a.function === f)).map(([label, color]) => ({ label, color })),
    aria: 'Frequency of the commonest CYP2D6 alleles in each CPIC biogeographic group',
    describe: (i, j) => {
      const a = cy.top[j];
      const v = alleleValues[i][j];
      return tipLines(`${a.name} in ${groups[i]}: ${v === null ? 'no frequency' : pct(v)}`, `${a.function}, activity value ${a.activity ?? 'n/a'}`);
    },
  });
  attachData('#cyp2d6-alleles-chart', 'cyp2d6-allele-frequencies', ['Allele', 'Function', 'Activity value', ...groups], cy.top.map((a) => [a.name, a.function, a.activity ?? '', ...groups.map((g) => a.freqs[g] ?? '')]));
  hbars({
    mount: $('#cyp2d6-phenotypes-chart'),
    rows: groups.map((g) => ({ key: g, label: g, values: Object.fromEntries(METABOLIZERS.map(([m]) => [m, (cy.phenotypes[g]?.[m] ?? 0) * 100])) })),
    segments: METABOLIZERS.map(([key, color]) => ({ key, label: key, color })),
    max: 100,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => `PM ${pctOf(r.values['Poor Metabolizer'])}`,
    labelShare: 0.3,
    labelMax: 220,
    aria: 'Predicted CYP2D6 phenotypes in each CPIC biogeographic group',
    describe: (r) => tipLines(r.key, 'predicted CYP2D6 phenotypes', ...METABOLIZERS.map(([m]) => `${m}: ${pctOf(r.values[m])}`)),
  });
  attachData('#cyp2d6-phenotypes-chart', 'cyp2d6-phenotypes', ['Group', ...METABOLIZERS.map(([m]) => m)], groups.map((g) => [g, ...METABOLIZERS.map(([m]) => cy.phenotypes[g]?.[m] ?? 0)]));
  const scoreValues = groups.map((g) => cy.scores.groups[g] ?? cy.scores.bins.map(() => null));
  shareGrid({
    mount: $('#cyp2d6-scores-chart'),
    rows: groups,
    cols: cy.scores.bins,
    values: scoreValues,
    max: Math.max(...scoreValues.flat().filter((v) => v !== null)),
    color: 'var(--s5)',
    aria: 'Share of each CPIC biogeographic group at each CYP2D6 activity score',
    describe: (i, j) => tipLines(scoreValues[i][j] === null ? 'No frequency' : pct(scoreValues[i][j]), `of ${groups[i]} people have activity score ${cy.scores.bins[j]}`),
  });
  attachData('#cyp2d6-scores-chart', 'cyp2d6-activity-scores', ['Group', ...cy.scores.bins], groups.map((g, i) => [g, ...scoreValues[i].map((v) => v ?? '')]));
  const mixKeys = FUNCTION_COLORS.filter(([f]) => groups.some((g) => cy.mix[g]?.[f]));
  hbars({
    mount: $('#cyp2d6-function-chart'),
    rows: groups.map((g) => ({ key: g, label: g, values: Object.fromEntries(mixKeys.map(([f]) => [f, (cy.mix[g]?.[f] ?? 0) * 100])) })),
    segments: mixKeys.map(([key, color]) => ({ key, label: key, color })),
    max: 100,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => `none ${pctOf(r.values['No function'] ?? 0)}`,
    labelShare: 0.3,
    labelMax: 220,
    aria: 'CYP2D6 allele frequency by clinical function in each CPIC biogeographic group',
    describe: (r) => tipLines(r.key, 'CYP2D6 alleles by function', ...mixKeys.map(([f]) => `${f}: ${pctOf(r.values[f])}`)),
  });
  attachData('#cyp2d6-function-chart', 'cyp2d6-function-mix', ['Group', ...mixKeys.map(([f]) => f)], groups.map((g) => [g, ...mixKeys.map(([f]) => cy.mix[g]?.[f] ?? 0)]));
};

const renderFdaLabels = () => {
  const fda = snapshot.fdaPgx;
  const boxed = fda.sections.find((s) => s.section === 'Boxed Warning')?.rows ?? 0;
  const typo = fda.folded ? ` ${fda.folded === 1 ? 'One row' : `${exact(fda.folded)} rows`} spelled Onoclogy ${fda.folded === 1 ? 'is' : 'are'} counted as Oncology.` : '';
  $('#fdalabels-dek').textContent = `FDA’s table, current as of ${longDate(fda.updated)}, gives each drug a row per biomarker, ${exact(fda.rows)} rows in all. ${fda.top[0].biomarker} appears in the most rows, ${exact(fda.top[0].rows)}. ${exact(boxed)} rows put the biomarker in a boxed warning.${typo}`;
  dotplot({
    mount: $('#fdalabels-areas-chart'),
    rows: fda.areas.map((a, i) => ({
      key: a.area, label: a.area, value: a.rows, hi: i === 0,
    })),
    labelShare: 0.4,
    rowH: 20,
    aria: 'Rows of FDA’s pharmacogenomic biomarker table by therapeutic area, log scale',
    describe: (r) => tipLines(`${exact(r.value)} rows`, r.key, pct(r.value / fda.rows)),
  });
  attachData('#fdalabels-areas-chart', 'fda-pgx-areas', ['Therapeutic area', 'Rows'], fda.areas.map((a) => [a.area, a.rows]));
  dotplot({
    mount: $('#fdalabels-biomarkers-chart'),
    rows: fda.top.map((b, i) => ({
      key: b.biomarker, label: b.biomarker, value: b.rows, hi: i === 0,
    })),
    scale: 'linear',
    labelShare: 0.4,
    labelMax: 320,
    rowH: 20,
    aria: 'Biomarkers in the most rows of FDA’s pharmacogenomic biomarker table',
    describe: (r) => tipLines(`${exact(r.value)} rows`, r.key, pct(r.value / fda.rows)),
  });
  attachData('#fdalabels-biomarkers-chart', 'fda-pgx-biomarkers', ['Biomarker', 'Rows'], fda.top.map((b) => [b.biomarker, b.rows]));
  dotplot({
    mount: $('#fdalabels-sections-chart'),
    rows: fda.sections.map((s, i) => ({
      key: s.section, label: s.section, value: s.rows, hi: i === 0,
    })),
    scale: 'linear',
    labelShare: 0.4,
    rowH: 20,
    aria: 'Label sections that carry pharmacogenomic information in FDA’s table',
    describe: (r) => tipLines(`${exact(r.value)} rows`, r.key, pct(r.value / fda.rows)),
  });
  attachData('#fdalabels-sections-chart', 'fda-pgx-label-sections', ['Label section', 'Rows'], fda.sections.map((s) => [s.section, s.rows]));
};

const renderClinpgx = () => {
  const ann = snapshot.clinpgx.annotations;
  const lab = snapshot.clinpgx.labels;
  const levels = ordered(ann.levels, CLINPGX_ORDER);
  const cats = Object.keys(ann.categories).sort((a, b) => sumOf(Object.values(ann.categories[b])) - sumOf(Object.values(ann.categories[a])));
  const agencies = Object.keys(lab.sources).sort((a, b) => sumOf(Object.values(lab.sources[b])) - sumOf(Object.values(lab.sources[a])));
  $('#clinpgx-dek').textContent = `Clinical annotations file of ${longDate(ann.created)}; drug label annotations of ${longDate(lab.created)}. A level 1A annotation is a variant–drug pair in a guideline or in clinical use; level 4 rests on case reports, non-significant studies or lab assays. ${exact(ann.pediatric)} annotations concern children. ${exact(lab.labels)} label annotations come from ${exact(agencies.length)} agencies; ${exact(lab.onFdaList)} are on FDA’s biomarker list and ${exact(lab.cancer)} concern tumour genomes.`;
  const levelTotals = levels.map((l) => ann.levels[l]);
  const counts = levels.map((l) => cats.map((c) => ann.categories[c][l] ?? 0));
  shareGrid({
    mount: $('#clinpgx-levels-chart'),
    rows: levels.map((l) => `Level ${l}`),
    cols: cats,
    values: counts.map((row, i) => row.map((n) => n / levelTotals[i])),
    totals: levelTotals,
    aria: 'ClinPGx clinical annotations by level of evidence and phenotype category, as a share of each level',
    describe: (i, j) => tipLines(pct(counts[i][j] / levelTotals[i]), `of level ${levels[i]} annotations are in ${cats[j]}`, `${exact(counts[i][j])} of ${exact(levelTotals[i])}`),
  });
  attachData('#clinpgx-levels-chart', 'clinpgx-annotations-by-level', ['Level', ...cats, 'Annotations'], levels.map((l, i) => [l, ...counts[i], levelTotals[i]]));
  hbars({
    mount: $('#clinpgx-genes-chart'),
    rows: ann.genes.map((g) => ({
      key: g.gene, label: g.gene, values: g.levels, gene: g,
    })),
    segments: levels.map((l) => ({ key: l, label: `Level ${l}`, color: CLINPGX_COLORS[l] ?? 'var(--data)' })),
    labelShare: 0.12,
    rowH: 18,
    aria: 'Genes with the most ClinPGx clinical annotations, by level of evidence',
    describe: (r) => tipLines(`${exact(r.gene.total)} annotations`, r.gene.gene, ...levels.filter((l) => r.gene.levels[l]).map((l) => `Level ${l}: ${exact(r.gene.levels[l])}`)),
  });
  attachData('#clinpgx-genes-chart', 'clinpgx-annotations-by-gene', ['Gene', ...levels.map((l) => `Level ${l}`), 'Total'], ann.genes.map((g) => [g.gene, ...levels.map((l) => g.levels[l] ?? 0), g.total]));
  const testing = TESTING_LEVELS.filter(([t]) => agencies.some((a) => lab.sources[a][t]));
  hbars({
    mount: $('#clinpgx-labels-chart'),
    rows: agencies.map((a) => ({ key: a, label: a, values: lab.sources[a] })),
    segments: testing.map(([key, color]) => ({ key, label: key, color })),
    labelShare: 0.14,
    aria: 'ClinPGx drug label annotations by agency and pharmacogenomic testing level',
    describe: (r) => tipLines(`${exact(sumOf(Object.values(r.values)))} label annotations`, r.key, ...testing.filter(([t]) => r.values[t]).map(([t]) => `${t}: ${exact(r.values[t])}`)),
  });
  attachData('#clinpgx-labels-chart', 'clinpgx-label-annotations', ['Agency', ...testing.map(([t]) => t)], agencies.map((a) => [a, ...testing.map(([t]) => lab.sources[a][t] ?? 0)]));
};

const MLOG_5E8 = -Math.log10(5e-8);
const SUPERSCRIPT = '⁰¹²³⁴⁵⁶⁷⁸⁹';
const pValue = (mlog) => {
  let e = Math.ceil(mlog);
  let m = 10 ** (e - mlog);
  if (m >= 9.95) {
    m = 1;
    e -= 1;
  }
  return `p = ${m.toFixed(1)}×10⁻${String(e).replace(/\d/g, (d) => SUPERSCRIPT[d])}`;
};
const capitalized = (s) => `${s[0].toUpperCase()}${s.slice(1)}`;

const renderManhattan = () => {
  const gw = snapshot.gwas;
  const catalogBins = Object.entries(gw.manhattan).flatMap(([chrom, list]) => list.map((b) => [chrom, b]));
  const overCap = catalogBins.filter(([, b]) => b[1] > 300).length;
  const [chrom, top] = catalogBins.reduce((x, y) => (y[1][1] > x[1][1] ? y : x));
  $('#manhattan-dek').textContent = `${exact(gw.placed)} of the ${exact(gw.associations)} associations carry a GRCh38 position, and ${exact(catalogBins.length)} megabases hold at least one. ${exact(overCap)} megabases pass −log10 p of 300 and are drawn at the cap; the strongest, ${top[4]} for ${top[3]} on chromosome ${chrom}, reaches ${pValue(top[1])}. The catalog’s latest addition is dated ${longDate(gw.latestAdded)}.`;
  manhattan({
    mount: $('#manhattan-chart'),
    chromosomes: NUCLEAR.map((c) => ({ chrom: c.chrom, length: c.length, bins: gw.manhattan[c.chrom] ?? [] })),
    cap: 300,
    threshold: MLOG_5E8,
    aria: 'Strongest GWAS Catalog association in each megabase of GRCh38, capped at -log10 p of 300',
    describe: (c, b) => tipLines(pValue(b[1]), `${b[3]}${b[4] ? `, ${b[4]}` : ''}`, `${b[5]}, ${megabase(c, b[0])}`, `${exact(b[2])} associations in this megabase`),
  });
  attachData('#manhattan-chart', 'gwas-catalog-manhattan', ['Chromosome', 'Megabase', 'Strongest -log10 p', 'Associations', 'Trait', 'Mapped gene', 'Variant'], catalogBins.map(([c, b]) => [c, b[0], b[1], b[2], b[3], b[4], b[5]]));
  const contexts = Object.entries(gw.context).sort((x, y) => y[1] - x[1]);
  dotplot({
    mount: $('#manhattan-context-chart'),
    rows: contexts.map(([k, n], i) => ({
      key: k, label: k, value: n, hi: i === 0,
    })),
    labelShare: 0.4,
    rowH: 20,
    aria: 'GWAS Catalog associations by the functional context of the lead variant, log scale',
    describe: (r) => tipLines(`${exact(r.value)} associations`, r.key, pct(r.value / gw.associations)),
  });
  attachData('#manhattan-context-chart', 'gwas-catalog-contexts', ['Context', 'Associations'], contexts);
};

const renderGwasGrowth = () => {
  const gw = snapshot.gwas;
  const latestYear = gw.latestAdded.slice(0, 4);
  const rows = Object.keys(gw.byYear).sort().map((y) => ({ year: y, partial: y === latestYear, ...gw.byYear[y] }));
  const { sizes } = gw;
  const none = sumOf(sizes.counts[0]);
  $('#gwasgrowth-dek').textContent = `The catalog records associations with p below 1×10⁻⁵ from each study it curates. ${exact(sizes.studies)} studies list a discovery sample; ${exact(none)} of them, ${pct(none / sizes.studies)}, have no association in the catalog.`;
  columns({
    mount: $('#gwasgrowth-years-chart'),
    rows,
    key: (r) => r.year,
    value: (r) => r.associations,
    overlay: (r) => r.studies,
    label: (r) => r.year,
    aria: 'GWAS Catalog associations and studies by publication year',
    describe: (r) => tipLines(`${exact(r.associations)} associations`, `from studies published in ${r.year}${r.partial ? ', to date' : ''}`, `${exact(r.studies)} studies with at least one`),
  });
  attachData('#gwasgrowth-years-chart', 'gwas-catalog-by-year', ['Year', 'Associations', 'Studies'], rows.map((r) => [r.year, r.associations, r.studies]));
  const most = Math.max(...sizes.counts.flat());
  shareGrid({
    mount: $('#gwasgrowth-sizes-chart'),
    rows: sizes.assocLabels,
    cols: sizes.sizeLabels,
    values: sizes.counts.map((row) => row.map((n) => n / most)),
    format: (v, i, j) => compact(sizes.counts[i][j]),
    color: 'var(--s4)',
    aria: 'GWAS Catalog studies by discovery sample size and number of associations',
    describe: (i, j) => tipLines(`${exact(sizes.counts[i][j])} studies`, `${sizes.assocLabels[i]} associations`, `discovery sample ${sizes.sizeLabels[j]}`),
  });
  attachData('#gwasgrowth-sizes-chart', 'gwas-catalog-study-sizes', ['Associations', ...sizes.sizeLabels], sizes.assocLabels.map((label, i) => [label, ...sizes.counts[i]]));
};

const TRAIT_COLORS = ['var(--s2)', 'var(--s4)', 'var(--s1)', 'var(--s3)', 'var(--s5)', 'var(--ink-2)', 'var(--w3)', 'var(--data-2)', 'var(--data)'];

const renderGwasTraits = () => {
  const gw = snapshot.gwas;
  const { totals, byYear } = gw.categories;
  const rankedCats = Object.keys(totals).sort((x, y) => totals[y] - totals[x]);
  const shown = rankedCats.slice(0, 8);
  const rest = rankedCats.slice(8);
  const series = [...shown, ...(rest.length ? ['other'] : [])];
  const years = Object.keys(byYear).sort();
  const valueOf = (y, key) => (key === 'other' ? sumOf(rest.map((c) => byYear[y][c] ?? 0)) : byYear[y][key] ?? 0);
  const yearTotal = (y) => sumOf(Object.values(byYear[y]));
  $('#gwastraits-dek').textContent = `Each association takes the parent category of its first mapped trait, matched by ontology term or by trait name. ${shown[0]} holds ${pct(totals[shown[0]] / gw.associations)} of all associations and ${shown[1]} ${pct(totals[shown[1]] / gw.associations)}; ${pct((totals['Not mapped'] ?? 0) / gw.associations)} could not be placed in a category.`;
  stackedArea({
    mount: $('#gwastraits-years-chart'),
    years: years.map(Number),
    series: series.map((key, i) => ({
      key, label: key === 'other' ? 'Other categories' : key, color: TRAIT_COLORS[i], values: years.map((y) => (valueOf(y, key) / yearTotal(y)) * 100),
    })),
    max: 100,
    format: pctOf,
    yFormat: (v) => `${v}%`,
    head: (i) => `${years[i]}: ${exact(yearTotal(years[i]))} associations`,
    aria: 'Share of GWAS Catalog associations by parent trait category and publication year',
  });
  attachData('#gwastraits-years-chart', 'gwas-catalog-categories-by-year', ['Year', ...rankedCats], years.map((y) => [y, ...rankedCats.map((c) => byYear[y][c] ?? 0)]));
  dotplot({
    mount: $('#gwastraits-traits-chart'),
    rows: gw.topTraits.map((t, i) => ({
      key: t.trait, label: t.trait, value: t.n, hi: i === 0,
    })),
    scale: 'linear',
    labelShare: 0.42,
    rowH: 20,
    aria: 'Mapped traits with the most GWAS Catalog associations',
    describe: (r) => {
      const t = gw.topTraits.find((x) => x.trait === r.key);
      return tipLines(`${exact(t.n)} associations`, capitalized(t.trait), t.parent);
    },
  });
  attachData('#gwastraits-traits-chart', 'gwas-catalog-top-traits', ['Trait', 'Parent category', 'Associations'], gw.topTraits.map((t) => [t.trait, t.parent, t.n]));
  dotplot({
    mount: $('#gwastraits-genes-chart'),
    rows: gw.topGenes.map((g, i) => ({
      key: g.gene, label: g.gene, value: g.n, hi: i === 0,
    })),
    scale: 'linear',
    labelShare: 0.2,
    rowH: 20,
    aria: 'Mapped genes with the most GWAS Catalog associations',
    describe: (r) => tipLines(`${exact(r.value)} associations`, r.key),
  });
  attachData('#gwastraits-genes-chart', 'gwas-catalog-top-genes', ['Mapped gene', 'Associations'], gw.topGenes.map((g) => [g.gene, g.n]));
};

const renderEffects = () => {
  const e = snapshot.gwas.effects;
  const orLabel = (i) => (i === e.orEdges.length - 1 ? `≥${e.orEdges[i]}` : `${e.orEdges[i]}–${e.orEdges[i + 1]}`);
  const rafLabels = e.rafEdges.slice(0, -1).map((lo, j) => `${lo.toFixed(1)}–${e.rafEdges[j + 1].toFixed(1)}`);
  const order = e.orEdges.map((_, i) => i).reverse();
  const most = Math.max(...e.counts.flat());
  const large = e.counts.slice(4);
  const small = e.counts[0];
  $('#effects-dek').textContent = `Disease traits only, where the catalog’s effect column holds odds ratios above 1. Of the associations with an odds ratio of 2 or more, ${pct(sumOf(large.map((r) => r[0])) / sumOf(large.flat()))} have a risk allele frequency below 0.1; of those below 1.1, ${pct(small[0] / sumOf(small))} do.`;
  shareGrid({
    mount: $('#effects-chart'),
    rows: order.map((i) => `OR ${orLabel(i)}`),
    cols: rafLabels,
    values: order.map((i) => e.counts[i].map((n) => n / most)),
    format: (v, r, c) => compact(e.counts[order[r]][c]),
    color: 'var(--s5)',
    aria: 'Disease associations in the GWAS Catalog by odds ratio and risk allele frequency',
    describe: (r, c) => tipLines(`${exact(e.counts[order[r]][c])} associations`, `odds ratio ${orLabel(order[r])}`, `risk allele frequency ${rafLabels[c]}`),
  });
  attachData('#effects-chart', 'gwas-catalog-effect-sizes', ['Odds ratio', ...rafLabels], order.map((i) => [orLabel(i), ...e.counts[i]]));
};

const FAMILY_ORDER = ['Penalized regression', 'LDpred', 'PRS-CS', 'Pruning and thresholding', 'Genome-wide significant variants', 'SBayes', 'MegaPRS and LDAK', 'Other'];

const renderPgs = () => {
  const pg = snapshot.pgs;
  const fams = ordered(pg.families, FAMILY_ORDER);
  const years = Object.keys(pg.byYear).sort();
  const toDate = years.map((_, i) => Object.fromEntries(fams.map((f) => [f, sumOf(years.slice(0, i + 1).map((y) => pg.byYear[y][f] ?? 0))])));
  $('#pgs-dek').textContent = `Release of ${longDate(pg.catalog.release)}: ${exact(pg.catalog.scores)} scores from ${exact(pg.catalog.publications)} publications. Their ${exact(pg.methods)} method names are grouped here into ${exact(fams.length)} families; penalized regression built ${exact(pg.families['Penalized regression'] ?? 0)} scores and genome-wide significant variants alone ${exact(pg.families['Genome-wide significant variants'] ?? 0)}.`;
  stackedArea({
    mount: $('#pgs-years-chart'),
    years: years.map(Number),
    series: fams.map((f, i) => ({
      key: f, label: f, color: TRAIT_COLORS[i], values: toDate.map((x) => x[f]),
    })),
    head: (i) => `${years[i]}: ${exact(sumOf(Object.values(toDate[i])))} scores by year end, ${exact(sumOf(Object.values(pg.byYear[years[i]])))} released that year`,
    aria: 'Polygenic scores in the PGS Catalog by release year and method family',
  });
  attachData('#pgs-years-chart', 'pgs-catalog-by-year', ['Year', ...fams], years.map((y) => [y, ...fams.map((f) => pg.byYear[y][f] ?? 0)]));
  const v = pg.variants;
  histo({
    mount: $('#pgs-variants-chart'),
    bins: v.counts.map((_, i) => i),
    series: [{
      key: 'scores', label: 'Scores', color: 'var(--data)', values: v.counts,
    }],
    tickEvery: 2,
    tickLabel: (b) => compact(v.edges[b]),
    height: 200,
    aria: 'Polygenic scores by number of variants, log scale',
    describe: (i) => tipLines(`${exact(v.counts[i])} scores`, `${exact(v.edges[i])} to ${exact(v.edges[i + 1])} variants`),
  });
  attachData('#pgs-variants-chart', 'pgs-catalog-variants', ['Variants from', 'Variants to', 'Scores'], v.counts.map((n, i) => [Math.round(v.edges[i]), Math.round(v.edges[i + 1]), n]));
  dotplot({
    mount: $('#pgs-traits-chart'),
    rows: pg.traits.map((t, i) => ({
      key: t.trait, label: t.trait, value: t.n, hi: i === 0,
    })),
    scale: 'linear',
    labelShare: 0.42,
    rowH: 20,
    aria: 'Traits with the most polygenic scores',
    describe: (r) => tipLines(`${exact(r.value)} scores`, capitalized(r.key)),
  });
  attachData('#pgs-traits-chart', 'pgs-catalog-traits', ['Trait', 'Scores'], pg.traits.map((t) => [t.trait, t.n]));
};

const PGS_ANCESTRY_COLORS = {
  European: 'var(--s2)', 'East Asian': 'var(--s4)', 'African American or Afro-Caribbean': 'var(--s1)', 'Hispanic or Latin American': 'var(--s3)', 'South Asian': 'var(--s5)', 'Multiple categories': 'var(--ink-2)', 'Not reported': 'var(--data-2)',
};

const renderPgsAncestry = () => {
  const anc = snapshot.pgs.ancestry;
  const stages = [['gwas', 'GWAS behind the scores'], ['development', 'Score development'], ['evaluation', 'Score evaluation']];
  const sums = {};
  for (const [key] of stages) for (const [c, n] of Object.entries(anc[key])) sums[c] = (sums[c] ?? 0) + n;
  const cats = Object.keys(sums).sort((x, y) => sums[y] - sums[x]);
  const shown = cats.slice(0, 7);
  const segments = [...shown.map((c) => ({ key: c, label: c, color: PGS_ANCESTRY_COLORS[c] ?? 'var(--data)' })), ...(cats.length > 7 ? [{ key: 'other', label: 'Other categories', color: 'var(--w3)' }] : [])];
  const rows = stages.map(([key, label]) => {
    const total = sumOf(Object.values(anc[key]));
    const share = (c) => ((c === 'other' ? sumOf(cats.slice(7).map((x) => anc[key][x] ?? 0)) : anc[key][c] ?? 0) / total) * 100;
    return {
      key, label, total, values: Object.fromEntries(segments.map((s) => [s.key, share(s.key)])),
    };
  });
  $('#pgsancestry-dek').textContent = `Individuals are summed over the samples listed for each score, so a cohort used for many scores counts each time: ${compact(rows[0].total)} entries behind the GWAS, ${compact(rows[1].total)} in development and ${compact(rows[2].total)} in evaluation. Comma-joined categories count as Multiple categories.`;
  hbars({
    mount: $('#pgsancestry-chart'),
    rows,
    segments,
    max: 100,
    tickFormat: (t) => `${t}%`,
    valueText: (r) => compact(r.total),
    labelShare: 0.26,
    aria: 'Ancestry of the samples behind polygenic scores, by stage',
    describe: (r) => tipLines(`${exact(r.total)} individual entries`, r.label, ...segments.filter((s) => r.values[s.key]).map((s) => `${s.label}: ${pctOf(r.values[s.key])}`)),
  });
  attachData('#pgsancestry-chart', 'pgs-catalog-ancestry', ['Stage', ...cats], stages.map(([key, label]) => [label, ...cats.map((c) => anc[key][c] ?? 0)]));
};

const renderPgsPerformance = () => {
  const perf = snapshot.pgs.performance;
  const au = perf.auroc;
  $('#pgsperformance-dek').textContent = `${exact(perf.metrics)} performance metrics, each as its paper reported it: ${exact(au.n)} AUROCs, ${exact(perf.cindex.n)} C-indices with a median of ${perf.cindex.median.toFixed(2)}, and ${exact(perf.r2.n)} R² values with a median of ${perf.r2.median.toFixed(2)}. Covariates differ between papers, so values compare loosely.`;
  histo({
    mount: $('#pgsperformance-auroc-chart'),
    bins: au.counts.map((_, i) => i),
    series: [{
      key: 'auroc', label: 'Evaluations', color: 'var(--data)', values: au.counts,
    }],
    tickEvery: 4,
    tickLabel: (b) => au.edges[b].toFixed(1),
    height: 200,
    aria: 'Reported AUROC of polygenic scores',
    describe: (i) => tipLines(`${exact(au.counts[i])} evaluations`, `AUROC ${au.edges[i].toFixed(3)} to ${au.edges[i + 1].toFixed(3)}`),
  });
  attachData('#pgsperformance-auroc-chart', 'pgs-catalog-auroc', ['AUROC from', 'AUROC to', 'Evaluations'], au.counts.map((n, i) => [au.edges[i], au.edges[i + 1], n]));
  intervals({
    mount: $('#pgsperformance-categories-chart'),
    rows: au.byCategory.map((c) => ({
      key: c.category, label: c.category, value: c.median, lo: c.q1, hi: c.q3, kind: 'auc', cat: c,
    })),
    domain: [0.5, 0.9],
    ticks: [0.5, 0.6, 0.7, 0.8, 0.9],
    kinds: [{ key: 'auc', label: 'Median AUROC, interquartile range', color: 'var(--accent)' }],
    format: (v) => v.toFixed(2),
    valueText: (r) => `${r.value.toFixed(3)}, n ${exact(r.cat.n)}`,
    aria: 'Median reported AUROC of polygenic scores by trait category',
    describe: (r) => tipLines(`AUROC ${r.value.toFixed(3)}`, r.label, `middle half ${r.lo.toFixed(3)} to ${r.hi.toFixed(3)}`, `${exact(r.cat.n)} evaluations`),
  });
  attachData('#pgsperformance-categories-chart', 'pgs-catalog-auroc-by-category', ['Category', 'Evaluations', 'Lower quartile', 'Median', 'Upper quartile'], au.byCategory.map((c) => [c.category, c.n, c.q1, c.median, c.q3]));
};

const TYPE_LABEL = {
  continuous: 'Continuous measures', biomarkers: 'Biomarkers', categorical: 'Questionnaire answers', icd10: 'ICD-10 diagnoses', phecode: 'Phecodes', prescriptions: 'Prescriptions',
};
const POP_NAMES = {
  AFR: 'African', AMR: 'Admixed American', CSA: 'Central/South Asian', EAS: 'East Asian', EUR: 'European', MID: 'Middle Eastern',
};

const renderHeritability = () => {
  const hp = snapshot.panukb;
  const { eur } = hp;
  const pops = Object.keys(hp.byPop).sort((x, y) => hp.byPop[y].passing - hp.byPop[x].passing || x.localeCompare(y));
  $('#heritability-dek').textContent = `Pan-UKB estimated SNP heritability for ${exact(hp.phenotypes)} phenotypes in six genetic ancestry groups; ${exact(eur.passing)} European-ancestry estimates pass all of its quality flags. Estimates are on the liability scale, which equals the observed scale for continuous traits.`;
  histo({
    mount: $('#heritability-hist-chart'),
    bins: eur.counts.map((_, i) => i),
    series: [{
      key: 'h2', label: 'Traits', color: 'var(--data)', values: eur.counts,
    }],
    tickEvery: 4,
    tickLabel: (b) => eur.edges[b].toFixed(1),
    height: 200,
    aria: 'SNP heritability of UK Biobank traits in European-ancestry participants',
    describe: (i) => tipLines(`${exact(eur.counts[i])} traits`, `heritability ${eur.edges[i].toFixed(3)} to ${eur.edges[i + 1].toFixed(3)}`),
  });
  attachData('#heritability-hist-chart', 'panukb-heritability', ['Heritability from', 'Heritability to', 'Traits'], eur.counts.map((n, i) => [eur.edges[i], eur.edges[i + 1], n]));
  intervals({
    mount: $('#heritability-types-chart'),
    rows: eur.byType.map((t) => ({
      key: t.type, label: TYPE_LABEL[t.type] ?? t.type, value: t.median, lo: t.q1, hi: t.q3, kind: 'h2', t,
    })),
    domain: [0, 0.4],
    ticks: [0, 0.1, 0.2, 0.3, 0.4],
    kinds: [{ key: 'h2', label: 'Median heritability, interquartile range', color: 'var(--accent)' }],
    format: (v) => v.toFixed(1),
    valueText: (r) => `${r.value.toFixed(2)}, n ${exact(r.t.n)}`,
    aria: 'SNP heritability by Pan-UKB trait type',
    describe: (r) => tipLines(`median ${r.value.toFixed(3)}`, r.label, `middle half ${r.lo.toFixed(3)} to ${r.hi.toFixed(3)}`, `${exact(r.t.n)} traits`),
  });
  attachData('#heritability-types-chart', 'panukb-heritability-by-type', ['Trait type', 'Traits', 'Lower quartile', 'Median', 'Upper quartile'], eur.byType.map((t) => [t.type, t.n, t.q1, t.median, t.q3]));
  dotplot({
    mount: $('#heritability-top-chart'),
    rows: eur.top.map((t, i) => ({
      key: `${i}`, label: t.description, value: t.h2, hi: i === 0, t,
    })),
    scale: 'linear',
    format: (v) => v.toFixed(2),
    tickFormat: (v) => v.toFixed(1),
    labelShare: 0.5,
    labelMax: 420,
    rowH: 20,
    aria: 'UK Biobank traits with the highest SNP heritability in European-ancestry participants',
    describe: (r) => tipLines(`h² ${r.t.h2.toFixed(3)} ± ${r.t.se.toFixed(3)}`, r.t.description, TYPE_LABEL[r.t.type] ?? r.t.type),
  });
  attachData('#heritability-top-chart', 'panukb-top-heritability', ['Trait', 'Type', 'Heritability', 'Standard error'], eur.top.map((t) => [t.description, t.type, t.h2, t.se]));
  hbars({
    mount: $('#heritability-pops-chart'),
    rows: pops.map((p) => ({ key: p, label: `${p}, ${POP_NAMES[p] ?? p}`, values: { pass: hp.byPop[p].passing, rest: hp.byPop[p].estimates - hp.byPop[p].passing } })),
    segments: [{ key: 'pass', label: 'Pass all quality flags', color: 'var(--accent)' }, { key: 'rest', label: 'Other estimates', color: 'var(--data-2)' }],
    labelShare: 0.3,
    aria: 'Pan-UKB heritability estimates by ancestry group and quality',
    describe: (r) => tipLines(`${exact(hp.byPop[r.key].estimates)} estimates`, r.label, `${exact(hp.byPop[r.key].passing)} pass all quality flags`),
  });
  attachData('#heritability-pops-chart', 'panukb-estimates-by-group', ['Group', 'Estimates', 'Passing all flags'], pops.map((p) => [p, hp.byPop[p].estimates, hp.byPop[p].passing]));
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
  ['STRING version API, for the current human network files', FILES.stringVersion],
  ['Reactome database version', FILES.reactomeVersion],
  ['Reactome human pathways', FILES.reactomePathways],
  ['Reactome pathway hierarchy', FILES.reactomeRelations],
  ['Reactome UniProt mapping, all levels', FILES.reactomeUniprot],
  ['GO annotations, goa_human', FILES.goaHuman],
  ['Gene Ontology, go-basic', FILES.goBasic],
  ['ENCODE SCREEN cCRE registry V4, GRCh38', FILES.ccre],
  ['UCSC hg38 CpG islands', FILES.cpgIslands],
  ['GtRNAdb hg38 high-confidence mature tRNAs', FILES.trnas],
  ['ENCODE portal, released human experiments', FILES.encode],
  ['NCBI feature table for NC_012920.1, the rCRS', FILES.mitoFeatures],
  ['gnomAD v3.1 mitochondrial variants', FILES.gnomadMito],
  ['T2T-CHM13 CenSat satellite annotation', FILES.censat],
  ['T2T-CHM13 telomere annotation', FILES.telomeres],
  ['dbSNP latest release notes', FILES.dbsnpNotes],
  ['dbVar non-redundant deletions, GRCh38', FILES.nr_deletions],
  ['dbVar non-redundant duplications, GRCh38', FILES.nr_duplications],
  ['dbVar non-redundant insertions, GRCh38', FILES.nr_insertions],
  ['gnomAD v4.1 constraint metrics', FILES.gnomadConstraint],
  ['ClinGen gene curation list, GRCh38', FILES.clingenGenes],
  ['GWAS Catalog ancestry download', FILES.gwasAncestry],
  ['Allen Ancient DNA Resource on Harvard Dataverse', FILES.aadrDataset],
  ['IGSR sample search API', FILES.igsrSamples],
  ['NCBI Datasets human genome assembly reports', assemblyPage()],
  ['ClinVar GRCh38 VCF', FILES.clinvarVcf],
  ['ClinVar GRCh38 VCF archive', FILES.clinvarArchive],
  ['ClinVar variant_summary', FILES.clinvarTable],
  ['ClinVar gene_specific_summary', FILES.clinvarGenes],
  ['ClinVar organization_summary', FILES.clinvarOrgs],
  ['ClinGen gene-disease validity download', FILES.clingenValidity],
  ['ClinGen region curation list, GRCh38', FILES.clingenRegions],
  ['InterPro Pfam matches by UniProt protein', FILES.interproProtein],
  ['STRchive loci, GitHub', FILES.strchive],
  ['HPO disease annotations, phenotype.hpoa', FILES.hpoAnnotations],
  ['HPO ontology, hp.obo', FILES.hpoOntology],
  ['HPO genes to diseases', FILES.hpoGenes],
  ['Orphadata genes and disorders', FILES.orphaGenes],
  ['Orphadata prevalence', FILES.orphaPrevalence],
  ['Orphadata age of onset and inheritance', FILES.orphaAges],
  ['GTR test versions', FILES.gtrTests],
  ['FDA approved cellular and gene therapy products', FILES.fdaCgt],
  ['ClinicalTrials.gov API, studies', FILES.trials],
  ['IntOGen downloads, driver compendium and cohorts', FILES.intogenDownload],
  ['GDC API, data release status', FILES.gdcStatus],
  ['GDC API, projects', FILES.gdcProjects],
  ['GDC API, cases by primary site', FILES.gdcSites],
  ['GDC API, TCGA cases with exome data', FILES.gdcTested],
  ['GDC API, simple somatic mutation occurrences, one gene at a time', gdcSsmPage('TP53')],
  ['CPIC API, gene-drug pairs', FILES.cpicPairs],
  ['CPIC API, alleles with functions and frequencies', FILES.cpicAlleles],
  ['CPIC API, phenotype results with frequencies', FILES.cpicResults],
  ['FDA Table of Pharmacogenomic Biomarkers in Drug Labeling', FILES.fdaPgx],
  ['ClinPGx clinical annotations download', FILES.clinpgxAnnotations],
  ['ClinPGx drug label annotations download', FILES.clinpgxLabels],
  ['GWAS Catalog associations with ontology annotations', FILES.gwasAssociations],
  ['GWAS Catalog EFO trait mappings', FILES.gwasTraitMappings],
  ['PGS Catalog score metadata', FILES.pgsScores],
  ['PGS Catalog performance metrics', FILES.pgsPerformance],
  ['PGS Catalog score development samples', FILES.pgsDevelopment],
  ['PGS Catalog evaluation samples', FILES.pgsEvaluation],
  ['PGS Catalog REST API, release', FILES.pgsInfo],
  ['PGS Catalog REST API, trait categories', FILES.pgsCategories],
  ['Pan-UKB heritability manifest', FILES.panukbH2],
  ['Pan-UKB phenotype manifest', FILES.panukbPhenotypes],
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
    h('p', {}, `Built from public data. Snapshot ${snapshot.generated.slice(0, 16).replace('T', ' ')} UTC. Not an official site of NIH, NCBI, EMBL-EBI, UCSC, ClinGen, FDA or any other source.`),
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

const RENDERERS = [
  renderChromosomes, renderGaps, renderT2t, renderComposition, renderTransposons, renderLandscape, renderSegdups, renderTandem,
  renderGenecount, renderDensity, renderGenesize, renderStructure, renderIsoforms, renderNames, renderPseudogenes, renderImmune,
  renderSplicing, renderTranslation, renderCodons, renderProteome, renderStructures, renderMembrane, renderDomains, renderModifications,
  renderTargets, renderSpecificity, renderGtex, renderPlasma, renderInteractions, renderPathways, renderOntology, renderCcres,
  renderCpgIslands, renderTrna, renderEncode, renderMito, renderCensat, renderSexchrom, renderDbsnp, renderSv, renderConstraint,
  renderDosage, renderAncestry, renderAncient, renderPanels, renderAssemblies, renderClinvar, renderReview, renderGrowth, renderLollipops,
  renderSpectrum, renderValidity, renderCnvRegions, renderExpansions, renderHpo, renderOrphanet, renderGtr, renderTherapy, renderYield,
  renderPrenatal, renderDrivers, renderTcga, renderBurden, renderSomatic, renderCpic, renderPhenotypes, renderCyp2d6, renderFdaLabels,
  renderClinpgx, renderManhattan, renderGwasGrowth, renderGwasTraits, renderEffects, renderPgs, renderPgsAncestry, renderPgsPerformance,
  renderHeritability,
];

let drawRun = 0;

const drawCharts = () => {
  drawRun += 1;
  const run = drawRun;
  delete document.documentElement.dataset.drawn;
  measureWidths($$('.chart'));
  const queue = [...RENDERERS];
  const step = () => {
    if (run !== drawRun) return;
    const until = performance.now() + 32;
    while (queue.length && performance.now() < until) {
      const render = queue.shift();
      try {
        render();
      } catch (error) {
        reportError(error);
      }
    }
    if (queue.length) {
      setTimeout(step, 0);
      return;
    }
    forgetWidths();
    applySearch();
    document.documentElement.dataset.drawn = 'true';
  };
  step();
};

const paintTheme = () => {
  $('#theme').textContent = theme() === 'dark' ? 'Light' : 'Dark';
  $('#theme').setAttribute('aria-pressed', String(theme() === 'dark'));
};

const initTheme = () => {
  const saved = store.get(THEME_KEY);
  if (saved) document.documentElement.dataset.theme = saved;
  paintTheme();
  darkQuery.addEventListener('change', paintTheme);
  $('#theme').addEventListener('click', () => {
    const next = theme() === 'dark' ? 'light' : 'dark';
    document.documentElement.dataset.theme = next;
    store.set(THEME_KEY, next);
    paintTheme();
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
