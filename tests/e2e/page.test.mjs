import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { chromium, devices } from 'playwright';
import { start } from '../../scripts/serve.mjs';
import snapshot from '../../src/snapshot.js';
import { EXAMPLES, DIAGNOSTIC_YIELD, CFDNA_SCREENING } from '../../src/facts.js';
import { exact } from '../../src/lib.js';
import { findings } from '../../src/findings.js';

let server;
let browser;
let base;

const open = async ({ width = 1280, height = 900, search = '', motion = 'reduce' } = {}) => {
  const page = await browser.newPage({ viewport: { width, height }, reducedMotion: motion });
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  page.on('console', (m) => m.type() === 'error' && !/net::ERR_FAILED|fonts\.g/.test(m.text()) && errors.push(m.text()));
  await page.route(/^https?:\/\/(?!127\.0\.0\.1)/, (route) => route.abort());
  await page.goto(`${base}${search}`);
  await page.waitForSelector('html[data-drawn]', { state: 'attached' });
  return { page, errors };
};

const openPhone = async () => {
  const context = await browser.newContext({ ...devices['iPhone 13'], reducedMotion: 'reduce' });
  const page = await context.newPage();
  await page.route(/^https?:\/\/(?!127\.0\.0\.1)/, (route) => route.abort());
  await page.goto(base);
  await page.waitForSelector('html[data-drawn]', { state: 'attached' });
  return { page, context };
};

const tapAt = async (page, selector, fx = 0.5, fy = 0.5) => {
  await page.locator(selector).scrollIntoViewIfNeeded();
  const box = await page.locator(selector).boundingBox();
  await page.touchscreen.tap(box.x + box.width * fx, box.y + box.height * fy);
  await page.waitForTimeout(250);
};

const tipLinesOf = (page) => page.$$eval('#tip > *', (els) => els.map((e) => e.textContent));

const CHAPTERS = {
  assembly: ['chromosomes', 'gaps', 't2t', 'composition'],
  repeats: ['transposons', 'landscape', 'segdups', 'tandem'],
  genes: ['genecount', 'density', 'genesize', 'structure', 'isoforms', 'names', 'pseudogenes', 'immune'],
  transcripts: ['splicing', 'translation', 'codons'],
  proteins: ['proteome', 'structures', 'membrane', 'domains', 'modifications', 'targets'],
  expression: ['specificity', 'gtex', 'plasma'],
  networks: ['interactions', 'pathways', 'ontology'],
  regulation: ['ccres', 'cpgislands', 'trna', 'encode'],
  structure: ['mito', 'censat', 'sexchrom'],
  variation: ['dbsnp', 'sv'],
  constraint: ['constraint', 'dosage'],
  clinical: ['clinvar', 'review', 'growth', 'lollipops', 'spectrum', 'validity', 'cnvregions', 'expansions'],
  rare: ['hpo', 'orphanet', 'gtr', 'therapy', 'yield', 'prenatal'],
  reference: ['decoder', 'sources'],
};

const lastBin = (counts) => counts.findLastIndex((n) => n > 0) + 1;

before(async () => {
  server = await start(0);
  base = `http://127.0.0.1:${server.address().port}/`;
  browser = await chromium.launch();
});

after(async () => {
  await browser.close();
  server.close();
});

test('loads without errors, leads with the bases of a complete genome and states every finding', async () => {
  const { page, errors } = await open();
  assert.equal(await page.title(), 'Human Genome in Numbers');
  const t2t = snapshot.assembly.t2t.chromosomes.reduce((n, c) => n + c.length, 0);
  assert.equal((await page.textContent('[data-hero="bases"]')).trim(), exact(t2t));
  const f = findings(snapshot);
  for (const node of await page.$$('[data-finding]')) {
    const key = await node.getAttribute('data-finding');
    assert.equal((await node.textContent()).trim(), f[key], key);
  }
  assert.equal(await page.locator('[data-finding]').count(), Object.keys(f).length);
  assert.deepEqual(errors, []);
  await page.close();
});

test('every chart draws an svg', async () => {
  const { page } = await open();
  const empty = await page.$$eval('.chart', (els) => els.filter((e) => !e.querySelector('svg')).map((e) => e.id));
  assert.deepEqual(empty, []);
  assert.ok(await page.locator('.chart svg').count() >= 4);
  await page.close();
});

test('the ideogram draws every chromosome and band to scale and names a band on hover', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#ideogram-chart [data-chrom]').count(), 24);
  assert.equal(await page.locator('#ideogram-chart [data-band]').count(), snapshot.assembly.bands.length);
  const width = (c) => page.$eval(`#ideogram-chart [data-chrom="${c}"] .body`, (n) => n.getBBox().width);
  const ratio = (await width('21')) / (await width('1'));
  const real = snapshot.assembly.grch38.chromosomes.find((c) => c.chrom === '21').length / snapshot.assembly.grch38.chromosomes.find((c) => c.chrom === '1').length;
  assert.ok(Math.abs(ratio - real) < 0.01, `${ratio} vs ${real}`);
  await page.hover('#ideogram-chart [data-band="17p13.1"]');
  assert.match(await page.textContent('#tip'), /17p13\.1[\s\S]*chr17:/);
  await page.close();
});

test('the gap chart splits each chromosome\'s Ns by gap type', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#gaps-chart [data-row]').count(), 24);
  const types = new Set(Object.values(snapshot.assembly.gaps).flatMap((t) => Object.keys(t)));
  assert.equal(await page.locator('#gaps-chart .legend li').count(), types.size);
  await page.hover('#gaps-chart [data-row="Y"]');
  assert.match(await page.textContent('#tip'), /chromosome Y[\s\S]*heterochromatin/);
  await page.close();
});

test('the T2T chart shows the bases each chromosome gains', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#t2t-chart [data-row]').count(), 24);
  const y = snapshot.assembly.t2t.chromosomes.find((c) => c.chrom === 'Y').length - snapshot.assembly.grch38.chromosomes.find((c) => c.chrom === 'Y').ungapped;
  await page.hover('#t2t-chart [data-row="Y"]');
  assert.match(await page.textContent('#tip'), new RegExp(`\\+${exact(y)}`));
  await page.close();
});

test('the GC landscape draws one strip per chromosome with a cell per megabase', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#gc-chart [data-chrom]').count(), 24);
  const bins = Object.entries(snapshot.sequence.bins).filter(([c]) => c !== 'MT').reduce((n, [, b]) => n + b.filter((v) => v !== null).length, 0);
  assert.equal(await page.locator('#gc-chart [data-bin]').count(), bins);
  await page.hover('#gc-chart [data-bin="19:10"]');
  assert.match(await page.textContent('#tip'), /% GC[\s\S]*chr19:10,000,001-11,000,000/);
  await page.close();
});

test('the dinucleotide matrix shows CpG far below expectation', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#dinuc-chart [data-cell]').count(), 16);
  await page.hover('#dinuc-chart [data-cell="CG"]');
  const lines = await tipLinesOf(page);
  assert.match(lines[0], /^0\.2\d observed\/expected$/);
  assert.ok(lines.some((l) => /CpG/.test(l)));
  await page.close();
});

test('CpG and G-quadruplex charts list every chromosome', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#cpg-chart [data-row]').count(), 25);
  assert.equal(await page.locator('#g4-chart [data-row]').count(), 24);
  assert.equal(await page.locator('#g4-chart .legend li').count(), 2);
  await page.close();
});

test('the repeat treemap tiles the largest families and names them on hover', async () => {
  const { page } = await open();
  assert.ok(await page.locator('#repeat-chart [data-key]').count() >= 12);
  await page.hover('#repeat-chart [data-key="LINE/L1"]');
  assert.match(await page.textContent('#tip'), /L1[\s\S]*LINE[\s\S]*copies/);
  await page.close();
});

test('the repeat landscape stacks four transposon classes by divergence', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#landscape-chart path[data-key]').count(), 4);
  assert.equal(await page.locator('#landscape-chart .legend li').count(), 4);
  assert.equal(await page.locator('#alu-chart .legend li').count(), 3);
  await page.close();
});

test('segmental duplication and tandem repeat charts cover every bin, chromosome and unit', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#segdup-chart [data-bin]').count(), snapshot.repeats.superDups.identity.counts.length);
  assert.equal(await page.locator('#segdup-chrom-chart [data-row]').count(), 24);
  assert.equal(await page.locator('#tandem-chart [data-row]').count(), 9);
  assert.equal(await page.locator('#motif-chart [data-row]').count(), 15);
  await page.hover('#motif-chart [data-row="AC"]');
  assert.match(await page.textContent('#tip'), /\(AC\)n/);
  await page.close();
});

test('gene count charts compare sources and biotypes', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#source-count-chart [data-row]').count(), 4);
  assert.ok(await page.locator('#biotype-chart [data-row]').count() >= 12);
  await page.hover('#biotype-chart [data-row="protein_coding"]');
  assert.match(await page.textContent('#tip'), /protein-coding[\s\S]*genes/i);
  await page.close();
});

test('gene density pairs megabase strips with a GC scatter', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#density-chart [data-chrom]').count(), 24);
  const points = await page.locator('#density-gc-chart [data-point]').count();
  assert.ok(points > 2500, `${points}`);
  await page.hover('#density-gc-chart [data-point="19:10"]');
  assert.match(await page.textContent('#tip'), /protein-coding genes[\s\S]*GC/);
  await page.close();
});

test('gene size, structure and isoform histograms draw their bins', async () => {
  const { page } = await open();
  for (const id of ['span-chart', 'intron-chart', 'exon-count-chart', 'isoform-chart']) {
    assert.ok(await page.locator(`#${id} [data-bin]`).count() >= 10, id);
  }
  assert.equal(await page.locator('#longest-chart [data-row]').count(), 15);
  assert.equal(await page.locator('#territory-chart .legend li').count(), 5);
  await page.close();
});

test('naming, pseudogene and immune locus charts list their rows', async () => {
  const { page } = await open();
  assert.ok(await page.locator('#approvals-chart [data-bin]').count() >= 30);
  assert.equal(await page.locator('#groups-chart [data-row]').count(), 15);
  assert.equal(await page.locator('#parents-chart [data-row]').count(), 15);
  assert.equal(await page.locator('#immune-chart [data-row]').count(), 7);
  assert.equal(await page.locator('#immune-chart .legend li').count(), 5);
  await page.close();
});

test('sequence logos draw one stack of four letters per position', async () => {
  const { page } = await open();
  for (const [id, n] of [['donor-chart', 9], ['acceptor-chart', 23], ['kozak-chart', 10]]) {
    assert.equal(await page.locator(`#${id} [data-pos]`).count(), n, id);
    assert.equal(await page.locator(`#${id} text.glyph`).count(), n * 4, id);
  }
  await page.hover('#donor-chart [data-pos="3"]');
  assert.match(await page.textContent('#tip'), /bits[\s\S]*\+1[\s\S]*G/);
  assert.equal(await page.locator('#intron-class-chart [data-row]').count(), 4);
  await page.close();
});

test('the codon table shows all 64 codons with their amino acids', async () => {
  const { page } = await open();
  assert.equal(await page.locator('#codon-chart [data-codon]').count(), 64);
  await page.hover('#codon-chart [data-codon="ATG"]');
  assert.match(await page.textContent('#tip'), /per thousand[\s\S]*ATG[\s\S]*Methionine/);
  assert.equal(await page.locator('#amino-chart [data-row]').count(), 20);
  await page.close();
});

test('protein charts draw lengths, evidence, structures, membranes, domains and modifications', async () => {
  const { page } = await open();
  const P = snapshot.proteins;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#protein-length-chart [data-bin]'), P.lengths.counts.length);
  assert.equal(await count('#existence-chart [data-row]'), Object.keys(P.existence).length);
  assert.equal(await count('#longest-protein-chart [data-row]'), 10);
  assert.equal(await count('#pdb-count-chart [data-bin]'), 50);
  assert.equal(await count('#pdb-top-chart [data-row]'), 15);
  assert.equal(await count('#tm-chart [data-bin]'), 25);
  assert.equal(await count('#location-chart [data-row]'), 20);
  assert.equal(await count('#hpa-location-chart [data-row]'), Math.min(20, snapshot.expression.hpa.locations.length));
  assert.equal(await count('#pfam-chart [data-row]'), 25);
  await page.hover('#pfam-chart [data-row="PF00096"]');
  assert.match(await page.textContent('#tip'), /PF00096[\s\S]*Zinc finger, C2H2 type/);
  assert.equal(await count('#pfam-per-protein-chart [data-bin]'), 16);
  assert.equal(await count('#disorder-chart [data-bin]'), 10);
  assert.equal(await count('#ptm-chart [data-row]'), 9);
  await page.close();
});

test('the druggable genome splits Pharos development levels by family', async () => {
  const { page } = await open();
  const T = snapshot.targets;
  assert.equal(await page.locator('#tdl-chart .legend li').count(), 4);
  assert.equal(await page.locator('#family-chart [data-row]').count(), T.families.length);
  assert.equal(await page.locator('#family-chart .legend li').count(), 4);
  await page.hover('#family-chart [data-row="Kinase"]');
  assert.match(await page.textContent('#tip'), /Kinase[\s\S]*Tclin/);
  assert.equal(await page.locator('#novelty-chart [data-bin]').count(), T.novelty.length);
  assert.equal(await page.locator('#pubmed-score-chart [data-bin]').count(), T.pubmed.length);
  await page.close();
});

test('expression charts compare six HPA atlases, GTEx tissues and plasma levels', async () => {
  const { page } = await open();
  const { gtex, hpa } = snapshot.expression;
  const n = gtex.tissues.length;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#contexts-chart [data-row]'), 6);
  assert.equal(await count('#contexts-chart .legend li'), 5);
  assert.equal(await count('#enriched-chart [data-row]'), 20);
  assert.equal(await count('#celltype-chart [data-row]'), 20);
  assert.equal(await count('#tissue-genes-chart [data-row]'), n);
  assert.equal(await count('#tissue-corr-chart [data-cell]'), n * n);
  assert.equal(await count('#tissue-corr-chart text.lab'), n);
  await page.hover('#tissue-corr-chart [data-cell="0:1"]');
  const lines = await tipLinesOf(page);
  assert.match(lines[0], /^r = 0\.\d{3}$/);
  assert.ok(lines[1].includes(gtex.sites[0].name) && lines[1].includes(gtex.sites[1].name), lines[1]);
  assert.equal(await count('#egene-chart [data-point]'), gtex.sites.filter((s) => s.eGenes && s.samples).length);
  assert.equal(await count('#plasma-chart [data-bin]'), hpa.blood.edges.length - 1);
  assert.equal(await count('#plasma-chart .legend li'), 2);
  assert.equal(await count('#plasma-top-chart [data-row]'), 15);
  assert.equal(await count('#secretome-chart [data-row]'), Object.keys(hpa.secretome).length);
  await page.close();
});

test('axes speak the units of their data and count bins label their own bars', async () => {
  const { page } = await open();
  const plasmaTicks = await page.$$eval('#plasma-top-chart text.tick', (ts) => ts.map((t) => t.textContent));
  assert.ok(plasmaTicks.length >= 2);
  for (const t of plasmaTicks) assert.match(t, /g\/L$/);
  const offsets = await page.$$eval('#pfam-per-protein-chart svg', ([svg]) => [...svg.querySelectorAll('text.tick')]
    .filter((t) => t.getAttribute('text-anchor') === 'middle' && /^\d+\+?$/.test(t.textContent))
    .map((t) => {
      const i = Number(t.textContent.replace('+', ''));
      const bar = svg.querySelector(`[data-bin="${i}"] rect.bar`);
      return Math.abs(Number(t.getAttribute('x')) - (Number(bar.getAttribute('x')) + Number(bar.getAttribute('width')) / 2));
    }));
  assert.ok(offsets.length >= 3);
  for (const d of offsets) assert.ok(d < 1, `${d}`);
  await page.close();
});

test('the family chart compares shares, with each family\'s size at the end', async () => {
  const { page } = await open();
  const rows = await page.$$eval('#family-chart [data-row]', (gs) => gs.map((g) => ({
    width: [...g.querySelectorAll('rect.seg')].reduce((n, r) => n + Number(r.getAttribute('width')), 0),
    value: g.querySelector('text.val').textContent,
  })));
  const widths = rows.map((r) => r.width);
  assert.ok(Math.max(...widths) - Math.min(...widths) < 2, widths.join(', '));
  for (const r of rows) assert.match(r.value, /^[\d,]+ targets$/);
  await page.close();
});

test('network charts draw degrees on log axes, hubs, evidence channels and scores', async () => {
  const { page } = await open();
  const N = snapshot.networks;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#degree-chart [data-bin]'), lastBin(N.degree.counts));
  const yTicks = await page.$$eval('#degree-chart text.tick[text-anchor="end"]', (ts) => ts.map((t) => t.textContent));
  assert.ok(yTicks.length >= 3, yTicks.join(','));
  for (const t of yTicks) assert.match(t, /^(1|10|100|1K|10K)$/);
  assert.equal(await count('#hub-chart [data-row]'), 20);
  assert.equal(await count('#channel-chart [data-row]'), 7);
  assert.equal(await count('#channel-chart .legend li'), 2);
  assert.equal(await count('#score-chart [data-bin]'), N.scores.counts.length);
  await page.close();
});

test('pathway charts rank, overlap and size Reactome pathways', async () => {
  const { page } = await open();
  const R = snapshot.pathways;
  const n = R.topLevel.length;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#toplevel-chart [data-row]'), n);
  assert.equal(await count('#pathway-overlap-chart [data-cell]'), n * n);
  await page.hover('#pathway-overlap-chart [data-cell="0:1"]');
  const lines = await tipLinesOf(page);
  assert.match(lines[0], /^Jaccard index \d\.\d{3}$/);
  assert.ok(lines[1].includes(R.topLevel[0].name) && lines[1].includes(R.topLevel[1].name), lines[1]);
  assert.equal(await count('#breadth-chart [data-bin]'), R.breadth.findLastIndex((v) => v > 0));
  assert.equal(await count('#pathway-size-chart [data-bin]'), lastBin(R.sizes.counts));
  assert.equal(await count('#largest-pathway-chart [data-row]'), 15);
  await page.close();
});

test('ontology charts split GO annotations by evidence and rank proteins, terms and curators', async () => {
  const { page } = await open();
  const O = snapshot.ontology;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#go-evidence-chart [data-row]'), 3);
  assert.equal(await count('#go-evidence-chart .legend li'), 5);
  assert.equal(await count('#go-per-protein-chart [data-bin]'), lastBin(O.perProtein.counts));
  assert.equal(await count('#go-protein-chart [data-row]'), 15);
  for (const [id, aspect] of [['go-bp-chart', 'P'], ['go-mf-chart', 'F'], ['go-cc-chart', 'C']]) assert.equal(await count(`#${id} [data-row]`), O.topTerms[aspect].length, id);
  assert.equal(await count('#go-groups-chart [data-row]'), O.assignedBy.length);
  await page.close();
});

test('regulation charts cover cCRE classes and density, CpG islands and promoters', async () => {
  const { page } = await open();
  const R = snapshot.regulation;
  const classes = Object.keys(R.ccre.classes).length;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#ccre-class-chart [data-row]'), classes);
  assert.equal(await count('#ccre-density-chart [data-row]'), 24);
  assert.equal(await count('#ccre-density-chart .legend li'), classes);
  assert.equal(await count('#cpg-promoter-chart [data-row]'), Object.keys(R.cpg.promoters.bySpecificity).length + 1);
  assert.ok(await count('#cpg-length-chart [data-bin]') >= 10);
  assert.equal(await count('#cpg-oe-chart [data-bin]'), R.cpg.obsExp.counts.length);
  assert.equal(await count('#cpg-gc-chart [data-bin]'), R.cpg.gc.counts.length);
  assert.equal(await count('#cpg-density-chart [data-row]'), 24);
  await page.close();
});

test('tRNA genes are set against codon usage, isotype and chromosome', async () => {
  const { page } = await open();
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#trna-codon-chart [data-point]'), 61);
  await page.hover('#trna-codon-chart [data-point="GCT"]');
  assert.match(await page.textContent('#tip'), /GCT[\s\S]*AGC/);
  assert.equal(await count('#trna-isotype-chart [data-row]'), Object.keys(snapshot.trnas.isotypes).length);
  assert.equal(await count('#trna-chrom-chart [data-row]'), Object.keys(snapshot.trnas.byChrom).length);
  await page.close();
});

test('ENCODE charts count experiments by assay, phase, sample type and target', async () => {
  const { page } = await open();
  const E = snapshot.encode;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#encode-assay-chart [data-row]'), Math.min(25, E.assays.length));
  assert.equal(await count('#encode-phase-chart [data-row]'), E.phases.length);
  assert.equal(await count('#encode-sample-chart [data-row]'), E.biosamples.length);
  assert.equal(await count('#encode-target-chart [data-row]'), E.targets.length);
  await page.close();
});

test('the mitochondrial genome map draws every gene on its strand and the variation along it', async () => {
  const { page } = await open();
  const M = snapshot.structure.mito;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#mito-map-chart [data-feature]'), M.genes.length + M.dloop.length);
  assert.equal(await count('#mito-map-chart .legend li'), 4);
  await page.hover('#mito-map-chart [data-feature="ND6"]');
  assert.match(await page.textContent('#tip'), /ND6[\s\S]*14,149[\s\S]*minus/);
  assert.equal(await count('#mito-variant-chart [data-bin]'), M.variants.bins.length);
  assert.equal(await count('#mito-top-chart [data-row]'), M.variants.top.length);
  assert.equal(await count('#mito-het-chart [data-bin]'), M.variants.heteroplasmy.counts.length);
  await page.close();
});

test('satellite, telomere and sex chromosome charts cover every chromosome and both PARs', async () => {
  const { page } = await open();
  const T = snapshot.structure;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#censat-chart [data-row]'), 24);
  assert.ok(await count('#censat-chart .legend li') >= 6);
  assert.ok(await count('#censat-class-chart [data-row]') >= 6);
  assert.equal(await count('#telomere-chart [data-row]'), T.telomeres.length);
  assert.equal(await count('#xy-genes-chart [data-row]'), 3);
  assert.equal(await count('#par-chart [data-mark]'), 4);
  assert.equal(await count('#par-chart [data-chrom]'), 2);
  assert.ok(await count('#y-satellite-chart .legend li') >= 3);
  await page.close();
});

test('variation charts cover dbSNP density and status, and SV sizes, chromosomes, mobile elements, studies and clinical calls', async () => {
  const { page } = await open();
  const V = snapshot.variation;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#dbsnp-density-chart [data-row]'), 24);
  assert.equal(await count('#dbsnp-status-chart [data-row]'), Object.keys(V.dbsnp.status).length);
  assert.equal(await count('#sv-size-chart [data-bin]'), V.sv.sizes.deletion.counts.length);
  assert.equal(await count('#sv-size-chart .legend li'), Object.keys(V.sv.sizes).length);
  assert.equal(await count('#sv-chrom-chart [data-row]'), 24);
  assert.equal(await count('#sv-mobile-chart [data-row]'), 4);
  assert.equal(await count('#sv-study-chart [data-row]'), V.sv.studies.length);
  assert.equal(await count('#sv-clinical-chart [data-row]'), 5);
  await page.close();
});

test('constraint charts compare LOEUF for haploinsufficient genes and all others', async () => {
  const { page } = await open();
  const C = snapshot.constraint;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#loeuf-chart [data-bin]'), C.gnomad.loeuf.counts.length);
  assert.equal(await count('#loeuf-chart .legend li'), 2);
  assert.equal(await count('#misz-chart [data-bin]'), C.gnomad.misZ.counts.length);
  assert.equal(await count('#constrained-chart [data-row]'), C.gnomad.top.length);
  assert.equal(await count('#dosage-chart [data-row]'), 2);
  await page.close();
});

test('ClinVar charts cross classes with consequences and cover types, genes, review, submitters, growth and the spectrum', async () => {
  const { page } = await open();
  const C = snapshot.clinvar;
  const T = snapshot.clinvarTable;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#clinvar-grid-chart [data-cell]'), Object.keys(C.consequences).length * 9);
  assert.equal(await count('#clinvar-grid-chart text.lab'), Object.keys(C.consequences).length);
  assert.equal(await count('#clinvar-types-chart [data-row]'), Object.keys(T.types).length);
  assert.equal(await count('#clinvar-genes-chart [data-row]'), C.genes.top.length);
  assert.equal(await count('#clinvar-lorenz-chart path.line'), 1);
  assert.equal(await count('#review-stars-chart [data-row]'), 9);
  assert.equal(await count('#review-stars-chart .legend li'), 5);
  assert.equal(await count('#review-submitters-chart [data-bin]'), T.submitters.bins.length);
  assert.equal(await count('#review-submitters-chart .legend li'), 4);
  const years = Object.keys(T.evaluated['P/LP']).filter((y) => y !== 'none').map(Number);
  assert.equal(await count('#review-evaluated-chart [data-year]'), Math.max(...years) - 2005 + 2);
  assert.equal(await count('#review-orgs-chart [data-row]'), C.organizations.top.length);
  assert.equal(await count('#review-types-chart [data-row]'), C.organizations.types.length);
  assert.equal(await count('#growth-chart path.area'), 5);
  assert.equal(await count('#growth-chart text.tick[text-anchor="middle"]'), snapshot.clinvarGrowth.points.length);
  assert.equal(await count('#spectrum-chart [data-channel]'), 96);
  assert.equal(await count('#spectrum-chart text.sub-lab'), 6);
  assert.equal(await count('#spectrum-classes-chart [data-bin]'), 6);
  await page.close();
});

test('lollipops mirror P/LP and B/LB variants over each protein\'s Pfam domains and name a residue on hover', async () => {
  const { page } = await open();
  const T = snapshot.clinvarTable;
  for (const gene of ['TP53', 'BRCA1', 'CFTR', 'MYH7']) {
    const l = T.lollipops[gene];
    const sites = new Set([...l.up.map((r) => `up:${r[0]}`), ...l.down.map((r) => `down:${r[0]}`)]);
    assert.equal(await page.locator(`#lollipop-${gene}-chart [data-site]`).count(), sites.size, gene);
    assert.equal(await page.locator(`#lollipop-${gene}-chart [data-domain]`).count(), snapshot.domains[gene].domains.length, gene);
  }
  const top = [...T.lollipops.TP53.up].sort((a, b) => b[2] - a[2] || a[0] - b[0])[0];
  await page.locator(`#lollipop-TP53-chart [data-site="up:${top[0]}"]`).first().hover();
  assert.match((await tipLinesOf(page)).join(' '), new RegExp(`(?<!\\d)${top[0]}(?!\\d)`));
  await page.close();
});

test('ClinGen and STRchive charts split validity by inheritance, map dosage regions and range repeat thresholds', async () => {
  const { page } = await open();
  const V = snapshot.validity;
  const R = snapshot.regions;
  const E = snapshot.expansions;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#validity-chart [data-row]'), Object.keys(V.byMoi).length);
  assert.equal(await count('#validity-chart .legend li'), 7);
  assert.equal(await count('#validity-years-chart path.area'), 7);
  assert.equal(await count('#validity-loeuf-chart .legend li'), 2);
  assert.equal(await count('#validity-gcep-chart [data-row]'), 12);
  assert.equal(await count('#regions-ideogram-chart [data-mark]'), R.sensitive.length);
  assert.equal(await count('#regions-size-chart [data-row]'), R.sensitive.length);
  const ranged = E.loci.filter((l) => ['Definitive', 'Strong'].includes(l.evidence) && l.pathogenic);
  assert.equal(await count('#expansions-chart [data-row]'), ranged.length);
  assert.equal(await count('#expansions-chart rect.range'), ranged.reduce((n, l) => n + [l.benign, l.intermediate, l.pathogenic].filter(Boolean).length, 0));
  assert.ok(await count('#expansions-years-chart path.area') >= 4);
  await page.close();
});

test('rare disease charts cover HPO systems and terms, Orphanet prevalence, onset and inheritance, and GTR tests', async () => {
  const { page } = await open();
  const H = snapshot.hpo;
  const O = snapshot.orphanet;
  const G = snapshot.gtr;
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#hpo-systems-chart [data-row]'), H.systems.length);
  assert.equal(await count('#hpo-terms-chart [data-row]'), H.topTerms.length);
  assert.equal(await count('#hpo-perdisease-chart [data-bin]'), H.perDisease.counts.length);
  assert.equal(await count('#hpo-genes-chart [data-row]'), H.genes.top.length);
  assert.equal(await count('#orpha-prevalence-chart [data-row]'), Object.keys(O.prevalence.classes).length);
  assert.equal(await count('#orpha-onset-chart [data-row]'), Object.keys(O.ages.onset).length);
  assert.equal(await count('#orpha-inheritance-chart [data-row]'), Object.keys(O.ages.inheritance).length);
  assert.equal(await count('#orpha-genes-chart [data-row]'), Object.keys(O.genes.associations).length);
  assert.equal(await count('#gtr-methods-chart [data-row]'), Math.min(15, G.methods.length));
  assert.equal(await count('#gtr-genecount-chart [data-row]'), G.genesPerTest.bins.length);
  assert.equal(await count('#gtr-countries-chart [data-row]'), G.countries.length);
  assert.equal(await count('#gtr-genes-chart [data-row]'), G.topGenes.length);
  await page.close();
});

test('therapy and testing charts cover approvals, vectors, trials, diagnostic yield and cfDNA screening with intervals', async () => {
  const { page } = await open();
  const T = { trials: snapshot.trials };
  const count = (sel) => page.locator(sel).count();
  assert.equal(await count('#therapy-approvals-chart path.area'), 5);
  assert.ok(await count('#therapy-vectors-chart [data-row]') >= 8);
  assert.equal(await count('#therapy-trials-chart path.line'), Object.keys(T.trials).length);
  assert.equal(await count('#therapy-phases-chart [data-row]'), Object.keys(T.trials).length);
  assert.equal(await count('#yield-chart [data-row]'), DIAGNOSTIC_YIELD.length);
  assert.equal(await count('#yield-chart line.ci'), DIAGNOSTIC_YIELD.filter((y) => y.lo !== null).length);
  assert.equal(await count('#prenatal-chart [data-row]'), CFDNA_SCREENING.length);
  assert.equal(await count('#prenatal-chart line.ci'), CFDNA_SCREENING.length);
  await page.close();
});

test('on a phone the tissue heatmap drops its labels and still fits', async () => {
  const { page, context } = await openPhone();
  const n = snapshot.expression.gtex.tissues.length;
  assert.equal(await page.locator('#tissue-corr-chart [data-cell]').count(), n * n);
  assert.equal(await page.locator('#tissue-corr-chart text.lab').count(), 0);
  const box = await page.locator('#tissue-corr-chart svg').boundingBox();
  assert.ok(box.x >= 0 && box.x + box.width <= page.viewportSize().width + 0.5, JSON.stringify(box));
  await context.close();
});

test('figure tables and CSV files are built only when asked for', async () => {
  const { page } = await open();
  assert.equal(await page.locator('.fig figcaption details.data table').count(), 0);
  assert.ok(await page.locator('.fig figcaption details.data').count() >= 4);
  const csv = page.locator('#ideogram-chart').locator('xpath=ancestor::figure').locator('a.csv');
  assert.doesNotMatch(await csv.getAttribute('href') ?? '', /^data:/);
  await page.locator('#ideogram-chart').locator('xpath=ancestor::figure').locator('details.data summary').click();
  assert.equal(await page.locator('#ideogram-chart').locator('xpath=ancestor::figure').locator('details.data table').count(), 1);
  await csv.dispatchEvent('pointerdown');
  assert.match(await csv.getAttribute('href'), /^data:text\/csv;charset=utf-8,Chromosome/);
  await page.close();
});

test('one chart that fails to draw does not stop the others', async () => {
  const page = await browser.newPage({ viewport: { width: 1280, height: 900 }, reducedMotion: 'reduce' });
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  await page.addInitScript(() => document.addEventListener('DOMContentLoaded', () => document.querySelector('#gaps-chart').remove()));
  await page.route(/^https?:\/\/(?!127\.0\.0\.1)/, (route) => route.abort());
  await page.goto(base);
  await page.waitForSelector('html[data-drawn]', { state: 'attached' });
  assert.equal(await page.locator('#prenatal-chart svg').count(), 1);
  assert.equal(errors.length, 1);
  await page.close();
});

test('figure tables keep the decimals of their numbers', async () => {
  const { page } = await open();
  const fig = page.locator('#loeuf-chart').locator('xpath=ancestor::figure');
  await fig.locator('details.data summary').click();
  assert.deepEqual((await fig.locator('details.data tbody tr').first().locator('td').allTextContents()).slice(0, 2), ['0', '0.05']);
  await page.close();
});

test('switching the theme keeps the drawn charts', async () => {
  const { page } = await open();
  await page.evaluate(() => { document.querySelector('#ideogram-chart svg').dataset.mark = 'kept'; });
  await page.click('#theme');
  assert.equal(await page.locator('#ideogram-chart svg[data-mark="kept"]').count(), 1);
  await page.close();
});

test('the page draws without long main-thread stalls and keeps the DOM small', async () => {
  const page = await browser.newPage({ viewport: { width: 1280, height: 900 }, reducedMotion: 'reduce' });
  await page.addInitScript(() => {
    window.stalls = [];
    new PerformanceObserver((list) => list.getEntries().forEach((e) => window.stalls.push(e.duration))).observe({ type: 'longtask', buffered: true });
  });
  await page.route(/^https?:\/\/(?!127\.0\.0\.1)/, (route) => route.abort());
  await page.goto(base);
  await page.waitForSelector('html[data-drawn]', { state: 'attached' });
  const { longest, nodes } = await page.evaluate(() => ({ longest: Math.max(0, ...window.stalls), nodes: document.getElementsByTagName('*').length }));
  assert.ok(longest < 1000, `longest task ${Math.round(longest)} ms`);
  assert.ok(nodes < 100000, `${nodes} nodes`);
  await page.close();
});

test('every figure offers its numbers as a table and as CSV', async () => {
  const { page } = await open();
  await page.evaluate(() => document.querySelectorAll('.fig figcaption details.data').forEach((d) => { d.open = true; }));
  await page.waitForFunction(() => [...document.querySelectorAll('.fig figcaption details.data')].every((d) => d.querySelector('table')));
  await page.evaluate(() => document.querySelectorAll('.fig figcaption a.csv').forEach((a) => a.dispatchEvent(new PointerEvent('pointerdown'))));
  const figs = await page.$$eval('.fig', (els) => els.map((f) => ({
    id: f.querySelector('.chart')?.id,
    table: !!f.querySelector('figcaption details.data table'),
    head: [...f.querySelectorAll('figcaption details.data thead th')].map((th) => th.textContent).join(','),
    csv: f.querySelector('figcaption a.csv')?.getAttribute('href') ?? '',
    name: f.querySelector('figcaption a.csv')?.getAttribute('download') ?? '',
  })));
  assert.ok(figs.length >= 4);
  for (const f of figs) {
    assert.ok(f.table, `${f.id} table`);
    assert.match(f.csv, /^data:text\/csv;charset=utf-8,/, f.id);
    assert.equal(decodeURIComponent(f.csv.split(',').slice(1).join(',')).split('\r\n')[0], f.head, f.id);
    assert.match(f.name, /^[a-z0-9-]+\.csv$/, f.id);
  }
  await page.close();
});

test('sections sit in their parts with a reading guide up front', async () => {
  const { page } = await open();
  const layout = await page.evaluate(() => {
    const out = {};
    let current = null;
    for (const n of document.querySelectorAll('main > .chapter, main > section.story')) {
      if (n.classList.contains('chapter')) {
        current = n.dataset.chapter;
        out[current] = [];
      } else out[current].push(n.id);
    }
    return out;
  });
  assert.deepEqual(layout, CHAPTERS);
  assert.ok(await page.locator('#readme li').count() >= 5);
  assert.equal(await page.locator('#jump a.jump-ch').count(), Object.keys(CHAPTERS).length);
  await page.close();
});

test('lens buttons hide sections and update the URL', async () => {
  const { page } = await open();
  await page.click('#lenses button[data-lens="clinic"]');
  assert.equal(await page.isVisible('#chromosomes'), false);
  assert.equal(await page.isVisible('#decoder'), true);
  assert.match(page.url(), /\?lens=clinic$/);
  await page.click('#lenses button[data-lens="all"]');
  assert.equal(await page.isVisible('#chromosomes'), true);
  assert.doesNotMatch(page.url(), /lens=/);
  await page.close();
});

test('a lens in the URL is applied on load and hides empty parts', async () => {
  const { page } = await open({ search: '?lens=clinic' });
  assert.equal(await page.isVisible('[data-chapter="assembly"]'), false);
  assert.equal(await page.getAttribute('#lenses button[data-lens="clinic"]', 'aria-pressed'), 'true');
  await page.close();
});

test('theme toggle flips the color scheme', async () => {
  const { page } = await open();
  const initial = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
  await page.click('#theme');
  assert.notEqual(await page.evaluate(() => getComputedStyle(document.body).backgroundColor), initial);
  await page.close();
});

test('decoder identifies typed and example identifiers', async () => {
  const { page } = await open();
  await page.fill('#decoder-input', 'NM_000546.6:c.215C>G');
  await page.waitForSelector('#decoder-out a');
  assert.match(await page.textContent('#decoder-out'), /HGVS variant[\s\S]*RefSeq mRNA[\s\S]*substitution/);
  assert.equal(await page.getAttribute('#decoder-out a', 'href'), 'https://mutalyzer.nl/normalizer/NM_000546.6%3Ac.215C%3EG');
  await page.click('#examples button[data-id="17p13.1"]');
  assert.equal(await page.inputValue('#decoder-input'), '17p13.1');
  assert.match(await page.textContent('#decoder-out'), /Cytogenetic band[\s\S]*chr17:/);
  await page.fill('#decoder-input', 'not an id');
  assert.match(await page.textContent('#decoder-out'), /No match/);
  assert.equal(await page.locator('#examples button').count(), EXAMPLES.length);
  await page.close();
});

test('on desktop hover previews, click pins and Escape closes', async () => {
  const { page } = await open();
  await page.hover('#t2t-chart [data-row="9"]');
  assert.equal(await page.isVisible('#tip'), true);
  assert.notEqual(await page.getAttribute('#tip', 'data-pinned'), 'true');
  await page.click('#t2t-chart [data-row="9"]');
  assert.equal(await page.getAttribute('#tip', 'data-pinned'), 'true');
  await page.keyboard.press('Escape');
  assert.equal(await page.isVisible('#tip'), false);
  await page.close();
});

test('each chart is one tab stop and arrow keys step through its marks', async () => {
  const { page } = await open();
  const stops = await page.evaluate(() => ({
    marks: document.querySelectorAll('.chart [data-k][tabindex="0"]').length,
    charts: document.querySelectorAll('.chart svg[tabindex="0"]').length,
  }));
  assert.equal(stops.marks, 0);
  assert.ok(stops.charts >= 4, `${stops.charts}`);
  await page.focus('#t2t-chart svg');
  await page.keyboard.press('Home');
  assert.ok((await tipLinesOf(page)).includes('chromosome 1'));
  await page.keyboard.press('ArrowDown');
  assert.ok((await tipLinesOf(page)).includes('chromosome 2'));
  await page.close();
});

test('the page never scrolls sideways on a phone', async () => {
  const { page, context } = await openPhone();
  const over = await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
  assert.ok(over <= 0, `${over}px`);
  await context.close();
});

test('on a phone a tap pins the popup and a tap outside closes it', async () => {
  const { page, context } = await openPhone();
  await tapAt(page, '#t2t-chart [data-row="13"]', 0.6, 0.5);
  assert.equal(await page.getAttribute('#tip', 'data-pinned'), 'true');
  const tip = await page.locator('#tip').boundingBox();
  assert.ok(tip.x >= 0 && tip.x + tip.width <= page.viewportSize().width);
  await tapAt(page, '#t2t .dek', 0.5, 0.5);
  assert.equal(await page.isVisible('#tip'), false);
  await context.close();
});

test('chart labels are never cut off and axis labels never overlap', async () => {
  for (const width of [1280, 375]) {
    const { page } = await open({ width });
    const cut = await page.$$eval('.chart text.lab', (ts) => ts.map((t) => t.textContent).filter((s) => s.endsWith('…')));
    assert.deepEqual(cut, [], `${width}px`);
    const clashes = await page.$$eval('.chart svg', (svgs) => svgs.flatMap((svg) => {
      const ticks = [...svg.querySelectorAll('text.tick')].map((t) => ({ t: t.textContent, b: t.getBBox() }));
      return ticks.flatMap((a, i) => ticks.slice(i + 1).filter((c) => Math.abs(a.b.y - c.b.y) < 2 && a.b.x < c.b.x + c.b.width && c.b.x < a.b.x + a.b.width)
        .map((c) => `${svg.closest('[id]')?.id}: ${a.t}/${c.t}`));
    }));
    assert.deepEqual(clashes, [], `${width}px`);
    await page.close();
  }
});

test('scatter points sit inside their plot', async () => {
  const { page } = await open();
  const out = await page.$$eval('.chart svg', (svgs) => svgs.flatMap((svg) => {
    const box = svg.viewBox.baseVal;
    return [...svg.querySelectorAll('circle.pt')]
      .filter((c) => Number(c.getAttribute('cx')) > box.width || Number(c.getAttribute('cx')) < 0 || Number(c.getAttribute('cy')) < 0 || Number(c.getAttribute('cy')) > box.height)
      .map((c) => `${svg.closest('[id]')?.id}: ${c.dataset.point}`);
  }));
  assert.deepEqual(out, []);
  const right = await page.$$eval('#egene-chart svg', ([svg]) => {
    const ticks = [...svg.querySelectorAll('text.tick[text-anchor="middle"]')].map((t) => Number(t.getAttribute('x')));
    const dots = [...svg.querySelectorAll('circle.pt')].map((c) => Number(c.getAttribute('cx')));
    return { tick: Math.max(...ticks), dot: Math.max(...dots) };
  });
  assert.ok(right.dot <= right.tick + 0.5, JSON.stringify(right));
  await page.close();
});

test('chart labels stay inside their chart', async () => {
  for (const width of [1280, 375]) {
    const { page } = await open({ width });
    const out = await page.$$eval('.chart svg', (svgs) => svgs.flatMap((svg) => {
      const box = svg.viewBox.baseVal;
      return [...svg.querySelectorAll('text.lab')].map((t) => ({ t: t.textContent, b: t.getBBox() }))
        .filter(({ b }) => b.x < -1 || b.x + b.width > box.width + 1)
        .map(({ t }) => `${svg.closest('[id]')?.id}: ${t}`);
    }));
    assert.deepEqual(out, [], `${width}px`);
    await page.close();
  }
});

test('the page states numbers, not opinions', async () => {
  const { page } = await open();
  const text = await page.evaluate(() => document.querySelector('main').innerText);
  assert.doesNotMatch(text, /\b(amazing|astonishing|stunning|shocking|staggering|incredible|remarkabl\w*|surprising\w*|merely|sadly|interestingly)\b/i);
  await page.close();
});

test('sources are titled and grouped, not bare URLs', async () => {
  const { page } = await open();
  assert.ok(await page.locator('#sources .srcgroup').count() >= 2);
  const texts = await page.$$eval('#sources .srcgroup a', (as) => as.map((a) => a.textContent.trim()));
  assert.ok(texts.length >= 8);
  assert.deepEqual(texts.filter((t) => /^(https?:|www\.|[a-z0-9.-]+\.(gov|org|com|edu)\/)/.test(t)), []);
  await page.close();
});

test('the data behind the page can be downloaded and cited', async () => {
  const { page } = await open();
  const href = await page.getAttribute('#data-use a[href$="snapshot.json"]', 'href');
  const json = await (await page.request.get(new URL(href, page.url()).href)).json();
  assert.equal(json.generated, snapshot.generated);
  assert.ok(await page.locator('#data-use a[href="https://github.com/collide2954/human-genome"]').count());
  assert.match(await page.textContent('#data-use'), /Cite as:.*https:\/\/collide2954\.github\.io\/human-genome\//);
  await page.close();
});

test('shared links get a title, description and preview image', async () => {
  const { page } = await open();
  const meta = (sel) => page.getAttribute(sel, 'content');
  assert.equal(await meta('meta[property="og:title"]'), 'Human Genome in Numbers');
  assert.ok((await meta('meta[property="og:description"]')).length > 40);
  assert.equal(await meta('meta[property="og:url"]'), 'https://collide2954.github.io/human-genome/');
  assert.equal(await meta('meta[property="og:image"]'), 'https://collide2954.github.io/human-genome/og.png');
  assert.equal(await meta('meta[name="twitter:card"]'), 'summary_large_image');
  const res = await page.request.get(new URL('og.png', page.url()).href);
  assert.equal(res.status(), 200);
  const png = await res.body();
  assert.deepEqual([png.readUInt32BE(16), png.readUInt32BE(20)], [1200, 630]);
  await page.close();
});

test('the footer says AI helped and that a person reviewed the work', async () => {
  const { page } = await open();
  const note = await page.textContent('#foot .made');
  assert.match(note, /^Made by a person, with help from AI tools\./);
  assert.match(note, /reviewed before publishing/);
  assert.match(await page.textContent('#foot'), /Not an official site of NIH, NCBI, EMBL-EBI, UCSC, ClinGen, FDA or any other source/);
  await page.close();
});

test('the reading guide names the weekly sources and the fixed values from papers', async () => {
  const { page } = await open();
  const guide = await page.textContent('#readme');
  assert.match(guide, /ClinGen/);
  assert.match(guide, /Fixed values/);
  assert.doesNotMatch(await page.textContent('#stamp'), /from NCBI, UCSC, GENCODE and EBI files/);
  await page.close();
});

test('figures animate in when scrolled into view', async () => {
  const { page } = await open({ motion: 'no-preference' });
  assert.equal(await page.locator('#t2t .fig.in').count(), 0);
  await page.locator('#t2t-chart').scrollIntoViewIfNeeded();
  await page.waitForSelector('#t2t .fig.in');
  await page.close();
});
