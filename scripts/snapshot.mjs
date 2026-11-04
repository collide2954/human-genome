import { existsSync } from 'node:fs';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { gunzipSync } from 'node:zlib';
import { retry } from '../src/lib.js';
import { FILES, USER_AGENT } from '../src/queries.js';
import {
  buildAssembly, latestFile, toModule, rmskSummary, superDupsSummary, tandemSummary, parseAssemblyReport,
} from './parsers.mjs';
import { gencodeSummary } from './genes.mjs';
import { cacheName, downloadTo, gzLines } from './streams.mjs';
import { runSections } from './sections.mjs';
import { readTwoBit } from './twobit.mjs';
import { parseManeGtf, sequenceSection } from './sequence.mjs';

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

const sequence = async () => {
  const gtf = latestFile(await text(FILES.maneDir), /MANE\.GRCh38\.v([\d.]+)\.ensembl_genomic\.gtf\.gz/);
  if (!gtf) throw new Error('No MANE GTF in the current release folder');
  const transcripts = parseManeGtf((await text(`${FILES.maneDir}${gtf.file}`)).split('\n'));
  const twoBit = readTwoBit(await readFile(await file(FILES.hg38TwoBit)));
  return { maneRelease: gtf.version, ...sequenceSection(twoBit, transcripts) };
};

const summarize = async (url, accumulator) => {
  for await (const line of gzLines(await file(url))) if (line) accumulator.add(line);
  return accumulator.result();
};

const repeats = async () => ({
  rmsk: await summarize(FILES.rmsk, rmskSummary()),
  superDups: await summarize(FILES.superDups, superDupsSummary()),
  tandem: await summarize(FILES.simpleRepeat, tandemSummary()),
});

const genes = async () => {
  const gtf = latestFile(await text(FILES.gencodeDir), /gencode\.v(\d+)\.annotation\.gtf\.gz/);
  if (!gtf) throw new Error('No GENCODE annotation in the latest release folder');
  const lengths = Object.fromEntries(parseAssemblyReport(await text(FILES.grch38Report)).sequences
    .filter((s) => s.role === 'assembled-molecule').map((s) => [s.molecule, s.length]));
  return summarize(`${FILES.gencodeDir}${gtf.file}`, gencodeSummary({ lengths }));
};

const SECTIONS = {
  assembly, sequence, repeats, genes,
};

const now = new Date();
const only = process.argv.slice(2).filter((a) => SECTIONS[a]);
const previous = existsSync(OUT) ? (await import(`${OUT.href}?t=${now.getTime()}`)).default : {};
const snapshot = await runSections(SECTIONS, { now, previous, only });
await writeFile(OUT, toModule(snapshot));
process.stderr.write(`wrote ${OUT.pathname}\n`);
