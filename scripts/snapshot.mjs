import { existsSync } from 'node:fs';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { gunzipSync } from 'node:zlib';
import { retry } from '../src/lib.js';
import {
  FILES, USER_AGENT, PHAROS_QUERY, pharosLevelQuery,
} from '../src/queries.js';
import {
  buildAssembly, latestFile, toModule, rmskSummary, superDupsSummary, tandemSummary, parseAssemblyReport,
} from './parsers.mjs';
import { gencodeSummary } from './genes.mjs';
import { hgncSummary, geneInfoSummary, maneSummary } from './catalogs.mjs';
import { cacheName, downloadTo, gzLines, unzip } from './streams.mjs';
import {
  uniprotSummary, hpaSummary, gtexSummary, pharosSummary, tissueSites, pfamNames, alignSites,
} from './proteins.mjs';
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

const json = async (url) => JSON.parse(await text(url));

const headersOf = (url) => retry(async () => {
  const res = await fetch(url, { method: 'HEAD', headers: HEADERS });
  if (!res.ok) throw new Error(`HTTP ${res.status} ${url}`);
  return res.headers;
}, { tries: 5, wait: 3000 });

const graphql = (url, query) => retry(async () => {
  const res = await fetch(url, { method: 'POST', headers: { ...HEADERS, 'content-type': 'application/json' }, body: JSON.stringify({ query }) });
  if (!res.ok) throw new Error(`HTTP ${res.status} ${url}`);
  const out = await res.json();
  if (out.errors?.length) throw new Error(out.errors[0].message);
  return out;
}, { tries: 5, wait: 3000 });

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

const catalogs = async () => {
  const summary = latestFile(await text(FILES.maneDir), /MANE\.GRCh38\.v([\d.]+)\.summary\.txt\.gz/);
  if (!summary) throw new Error('No MANE summary in the current release folder');
  return {
    hgnc: hgncSummary(await text(FILES.hgnc)),
    geneInfo: geneInfoSummary(await text(FILES.geneInfo)),
    mane: { release: summary.version, ...maneSummary(await text(`${FILES.maneDir}${summary.file}`)) },
  };
};

const proteins = async () => {
  const head = await headersOf(new URL('/uniprotkb/search?query=accession:P04637&fields=accession&format=tsv', FILES.uniprot).href);
  const uniprot = await summarize(FILES.uniprot, uniprotSummary());
  const entries = {};
  for (const { id } of uniprot.pfam.top) entries[id] = await json(`${FILES.interproPfam}${id}/`);
  return {
    release: head.get('x-uniprot-release'),
    ...uniprot,
    pfam: { ...uniprot.pfam, top: pfamNames(uniprot.pfam.top, entries) },
  };
};

const targets = async () => {
  const byLevel = {};
  for (const level of ['Tclin', 'Tchem', 'Tbio', 'Tdark']) byLevel[level] = await graphql(FILES.pharos, pharosLevelQuery(level));
  return pharosSummary(await graphql(FILES.pharos, PHAROS_QUERY), byLevel);
};

const expression = async () => {
  const zip = unzip(await readFile(await file(FILES.hpa)));
  const tsv = Object.keys(zip).find((name) => name.endsWith('.tsv'));
  if (!tsv) throw new Error('No TSV inside the HPA download');
  const modified = (await headersOf(FILES.hpa)).get('last-modified');
  const gtex = gtexSummary((await text(FILES.gtexMedian)).split('\n'));
  return {
    hpa: { file: modified ? new Date(modified).toISOString().slice(0, 10) : null, ...hpaSummary(zip[tsv].toString('utf8')) },
    gtex: { release: 'v10', ...gtex, sites: alignSites(gtex.tissues, tissueSites((await json(FILES.gtexTissues)).data)) },
  };
};

const SECTIONS = {
  assembly, sequence, repeats, genes, catalogs, proteins, targets, expression,
};

const now = new Date();
const only = process.argv.slice(2).filter((a) => SECTIONS[a]);
const previous = existsSync(OUT) ? (await import(`${OUT.href}?t=${now.getTime()}`)).default : {};
const snapshot = await runSections(SECTIONS, { now, previous, only });
await writeFile(OUT, toModule(snapshot));
process.stderr.write(`wrote ${OUT.pathname}\n`);
