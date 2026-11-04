import { test, before, after } from 'node:test';
import assert from 'node:assert/strict';
import { chromium, devices } from 'playwright';
import { start } from '../../scripts/serve.mjs';
import snapshot from '../../src/snapshot.js';
import { EXAMPLES } from '../../src/facts.js';
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
  await page.waitForSelector('#ideogram-chart svg', { state: 'attached' });
  return { page, errors };
};

const openPhone = async () => {
  const context = await browser.newContext({ ...devices['iPhone 13'], reducedMotion: 'reduce' });
  const page = await context.newPage();
  await page.route(/^https?:\/\/(?!127\.0\.0\.1)/, (route) => route.abort());
  await page.goto(base);
  await page.waitForSelector('#ideogram-chart svg', { state: 'attached' });
  return { page, context };
};

const CHAPTERS = {
  assembly: ['chromosomes', 'gaps', 't2t'],
  reference: ['decoder', 'sources'],
};

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

test('every figure offers its numbers as a table and as CSV', async () => {
  const { page } = await open();
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
  assert.match(await page.textContent('#tip'), /chromosome 1\b/);
  await page.keyboard.press('ArrowDown');
  assert.match(await page.textContent('#tip'), /chromosome 2\b/);
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
  const row = page.locator('#t2t-chart [data-row="13"]');
  await row.scrollIntoViewIfNeeded();
  const box = await row.boundingBox();
  await page.touchscreen.tap(box.x + box.width * 0.6, box.y + box.height / 2);
  await page.waitForTimeout(250);
  assert.equal(await page.getAttribute('#tip', 'data-pinned'), 'true');
  const tip = await page.locator('#tip').boundingBox();
  assert.ok(tip.x >= 0 && tip.x + tip.width <= page.viewportSize().width);
  const dek = page.locator('#t2t .dek');
  const d = await dek.boundingBox();
  await page.touchscreen.tap(d.x + 5, d.y + 5);
  await page.waitForTimeout(250);
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
  assert.match(await page.textContent('#foot'), /Not an official NIH, NCBI, EMBL-EBI or UCSC site/);
  await page.close();
});

test('figures animate in when scrolled into view', async () => {
  const { page } = await open({ motion: 'no-preference' });
  assert.equal(await page.locator('#t2t .fig.in').count(), 0);
  await page.locator('#t2t-chart').scrollIntoViewIfNeeded();
  await page.waitForSelector('#t2t .fig.in');
  await page.close();
});
