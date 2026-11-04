import { existsSync } from 'node:fs';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { gunzipSync } from 'node:zlib';
import { retry } from '../src/lib.js';
import { FILES, USER_AGENT } from '../src/queries.js';
import { buildAssembly, toModule } from './parsers.mjs';
import { cacheName, downloadTo } from './streams.mjs';
import { runSections } from './sections.mjs';

const OUT = new URL('../src/snapshot.js', import.meta.url);
const HEADERS = { 'user-agent': USER_AGENT };
const CACHE = process.env.SNAPSHOT_CACHE ?? await mkdtemp(join(tmpdir(), 'human-genome-'));
await mkdir(CACHE, { recursive: true });

const file = (url) => retry(async () => {
  const path = join(CACHE, cacheName(url));
  return existsSync(path) ? path : downloadTo(url, path, HEADERS);
}, { tries: 5, wait: 3000 });

const text = async (url) => {
  const buffer = await readFile(await file(url));
  return (url.endsWith('.gz') ? gunzipSync(buffer) : buffer).toString('utf8');
};

const assembly = async () => buildAssembly({
  report38: await text(FILES.grch38Report),
  stats38: await text(FILES.grch38Stats),
  reportT2t: await text(FILES.t2tReport),
  statsT2t: await text(FILES.t2tStats),
  cyto: await text(FILES.cytobands),
  gaps: await text(FILES.gaps),
  cen: await text(FILES.centromeres),
});

const SECTIONS = { assembly };

const now = new Date();
const only = process.argv.slice(2).filter((a) => SECTIONS[a]);
const previous = existsSync(OUT) ? (await import(`${OUT.href}?t=${now.getTime()}`)).default : {};
const snapshot = await runSections(SECTIONS, { now, previous, only });
await writeFile(OUT, toModule(snapshot));
process.stderr.write(`wrote ${OUT.pathname}\n`);
