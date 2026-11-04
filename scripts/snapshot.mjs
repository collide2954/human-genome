import { existsSync } from 'node:fs';
import { mkdir, mkdtemp, readFile, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { gunzipSync } from 'node:zlib';
import { retry } from '../src/lib.js';
import {
  FILES, USER_AGENT, PHAROS_QUERY, pharosLevelQuery, stringFiles, LOLLIPOPS, TRIAL_MODALITIES, trialPage, assemblyPage, intogenFile, gdcSsmPage,
} from '../src/queries.js';
import {
  latestIntogen, intogenDrivers, intogenCohorts, gdcRelease, gdcPrograms, facetCounts, casesByProject,
} from './cancer.mjs';
import {
  cpicPairs, cpicAlleles, cpicPhenotypes, cyp2d6, fdaBiomarkers, clinpgxAnnotations, clinpgxLabels, createdDate,
} from './pgx.mjs';
import {
  traitParents, gwasAssociations, studySizes, pgsScores, pgsAncestry, pgsCategories, pgsPerformance, heritability,
} from './traits.mjs';
import {
  gwasAncestry, aadrSummary, aadrAnnoFile, igsrSamples, assemblySummary,
} from './populations.mjs';
import {
  hpoTerms, hpoAnnotations, geneDiseases, orphaGenes, orphaPrevalence, orphaAges, gtrTests, fdaProducts, trialYears,
} from './rare.mjs';
import { linearEdges } from '../src/stats.js';
import {
  buildAssembly, latestFile, toModule, rmskSummary, superDupsSummary, tandemSummary, parseAssemblyReport,
} from './parsers.mjs';
import { gencodeSummary } from './genes.mjs';
import {
  hgncSummary, geneInfoSummary, maneSummary, codingSymbols,
} from './catalogs.mjs';
import {
  clinvarVcf, vcfClassCounts, earliestFile, geneSummary, organizations, variantTable, pfamDomains, validitySummary, regionSummary, repeatLoci,
} from './clinical.mjs';
import {
  cacheName, downloadTo, anyLines, unzip, zipLines,
} from './streams.mjs';
import {
  stringNames, stringSummary, physicalSummary, reactomeSummary, oboTerms, goSummary,
} from './networks.mjs';
import {
  ccreSummary, cpgIslands, tssInIslands, specificityByGene, trnaSummary, encodeSummary,
} from './regulation.mjs';
import {
  mitoFeatures, mitoVariants, censatSummary, telomeres,
} from './structure.mjs';
import {
  dbsnpNotes, svSummary, constraintSummary, clingenDosage, groupedValues,
} from './variation.mjs';
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
  return (/\.b?gz$/.test(url) ? gunzipSync(buffer) : buffer).toString('utf8');
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

const maneGtf = async () => {
  const gtf = latestFile(await text(FILES.maneDir), /MANE\.GRCh38\.v([\d.]+)\.ensembl_genomic\.gtf\.gz/);
  if (!gtf) throw new Error('No MANE GTF in the current release folder');
  return { version: gtf.version, transcripts: parseManeGtf((await text(`${FILES.maneDir}${gtf.file}`)).split('\n')) };
};

const hpaTable = async () => {
  const zip = unzip(await readFile(await file(FILES.hpa)));
  const tsv = Object.keys(zip).find((name) => name.endsWith('.tsv'));
  if (!tsv) throw new Error('No TSV inside the HPA download');
  return zip[tsv].toString('utf8');
};

const sequence = async () => {
  const { version, transcripts } = await maneGtf();
  const twoBit = readTwoBit(await readFile(await file(FILES.hg38TwoBit)));
  return { maneRelease: version, ...sequenceSection(twoBit, transcripts) };
};

const summarize = async (url, accumulator) => {
  for await (const line of await anyLines(await file(url))) if (line) accumulator.add(line);
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
  const table = await hpaTable();
  const modified = (await headersOf(FILES.hpa)).get('last-modified');
  const gtex = gtexSummary((await text(FILES.gtexMedian)).split('\n'));
  return {
    hpa: { file: modified ? new Date(modified).toISOString().slice(0, 10) : null, ...hpaSummary(table) },
    gtex: { release: 'v10', ...gtex, sites: alignSites(gtex.tissues, tissueSites((await json(FILES.gtexTissues)).data)) },
  };
};

const networks = async () => {
  const [{ string_version: version }] = await json(FILES.stringVersion);
  const files = stringFiles(version);
  const names = stringNames(await text(files.info));
  return {
    version,
    ...await summarize(files.detailed, stringSummary({ names })),
    physical: await summarize(files.physical, physicalSummary()),
  };
};

const pathways = async () => {
  const accumulator = reactomeSummary({ pathways: await text(FILES.reactomePathways), relations: await text(FILES.reactomeRelations) });
  return { release: (await text(FILES.reactomeVersion)).trim(), ...await summarize(FILES.reactomeUniprot, accumulator) };
};

const ontology = async () => {
  const obo = oboTerms(await text(FILES.goBasic));
  return {
    release: obo.version?.replace(/^releases\//, '') ?? null,
    namespaces: obo.counts,
    obsolete: obo.obsolete,
    ...await summarize(FILES.goaHuman, goSummary(obo.terms)),
  };
};

const regulation = async () => {
  const { intervals, ...islands } = cpgIslands(await text(FILES.cpgIslands));
  const promoters = tssInIslands((await maneGtf()).transcripts, intervals, specificityByGene(await hpaTable()));
  return {
    ccre: await summarize(FILES.ccre, ccreSummary()),
    cpg: { ...islands, promoters },
  };
};

const trnas = async () => trnaSummary(await text(FILES.trnas));

const encode = async () => encodeSummary(await json(FILES.encode));

const structure = async () => ({
  mito: {
    ...mitoFeatures(await text(FILES.mitoFeatures)),
    variants: await summarize(FILES.gnomadMito, mitoVariants({ length: 16569 })),
  },
  censat: censatSummary(await text(FILES.censat)),
  telomeres: telomeres(await text(FILES.telomeres)),
});

const feed = async (url, accumulator) => {
  for await (const line of await anyLines(await file(url))) if (line) accumulator.add(line);
  return accumulator;
};

const variation = async () => {
  const sv = svSummary();
  for (const key of ['nr_deletions', 'nr_duplications', 'nr_insertions']) await feed(FILES[key], sv);
  return { dbsnp: dbsnpNotes(await text(FILES.dbsnpNotes)), sv: sv.result() };
};

let constraintRun;
const constraintFile = () => {
  constraintRun ??= feed(FILES.gnomadConstraint, constraintSummary());
  return constraintRun;
};

const constraint = async () => {
  const acc = await constraintFile();
  const gnomad = acc.result();
  const { hiGenes, ...dosage } = clingenDosage(await text(FILES.clingenGenes));
  return {
    gnomad,
    dosage,
    byDosage: groupedValues(acc.values(), new Map([...hiGenes].map((g) => [g, 'sufficient'])), gnomad.loeuf.edges),
  };
};

const genomeBase = (twoBit) => {
  const seqs = new Map();
  return (chrom, pos) => {
    const name = chrom === 'MT' ? 'chrM' : `chr${chrom}`;
    if (!seqs.has(name)) seqs.set(name, twoBit.sequence(name));
    const s = seqs.get(name);
    return s ? 'ACGTN'[s.at(pos - 1)] : null;
  };
};

const clinvar = async () => {
  const twoBit = readTwoBit(await readFile(await file(FILES.hg38TwoBit)));
  return {
    ...await summarize(FILES.clinvarVcf, clinvarVcf({ base: genomeBase(twoBit) })),
    genes: geneSummary(await text(FILES.clinvarGenes), codingSymbols(await text(FILES.hgnc))),
    organizations: organizations(await text(FILES.clinvarOrgs)),
  };
};

const clinvarTable = async () => summarize(FILES.clinvarTable, variantTable({ genes: LOLLIPOPS.map((l) => l.gene) }));

const domains = async () => {
  const byGene = {};
  for (const { gene, accession } of LOLLIPOPS) {
    const pfam = await json(`${FILES.interproProtein}${accession}/`);
    byGene[gene] = { accession, length: pfam.results?.[0]?.proteins?.[0]?.protein_length ?? null, domains: pfamDomains(pfam) };
  }
  return byGene;
};

const clinvarGrowth = async () => {
  const years = [...(await text(FILES.clinvarArchive)).matchAll(/href="(\d{4})\/"/g)].map((m) => Number(m[1])).filter((y) => y >= 2018).sort((a, b) => a - b);
  const points = [];
  for (const year of years) {
    const first = earliestFile(await text(`${FILES.clinvarArchive}${year}/`), /clinvar_(\d{8})\.vcf\.gz/);
    if (first) points.push(await summarize(`${FILES.clinvarArchive}${year}/${first.file}`, vcfClassCounts()));
  }
  points.push(await summarize(FILES.clinvarVcf, vcfClassCounts()));
  return { points };
};

const validity = async () => {
  const { moiGenes, ...summary } = validitySummary(await text(FILES.clingenValidity));
  return { ...summary, loeufByMoi: groupedValues((await constraintFile()).values(), moiGenes, linearEdges(0, 2, 0.05)) };
};

const regions = async () => regionSummary(await text(FILES.clingenRegions));

const expansions = async () => {
  const [last] = await json(FILES.strchiveCommits).catch(() => []);
  return { updated: last?.commit?.committer?.date?.slice(0, 10) ?? null, ...repeatLoci(await json(FILES.strchive)) };
};

const hpo = async () => {
  const onto = hpoTerms(await text(FILES.hpoOntology));
  return {
    ontology: onto.version, terms: onto.terms.size, ...hpoAnnotations(await text(FILES.hpoAnnotations), onto), genes: geneDiseases(await text(FILES.hpoGenes)),
  };
};

const orphanet = async () => ({
  genes: orphaGenes(await text(FILES.orphaGenes)),
  prevalence: orphaPrevalence(await text(FILES.orphaPrevalence)),
  ages: orphaAges(await text(FILES.orphaAges)),
});

const gtr = async () => gtrTests(await text(FILES.gtrTests));

const allStudies = async (query) => {
  const studies = [];
  let token = null;
  do {
    const page = await json(trialPage(query, token));
    studies.push(...(page.studies ?? []));
    token = page.nextPageToken ?? null;
  } while (token);
  return studies;
};

const fda = async () => fdaProducts(await text(FILES.fdaCgt));

const trials = async () => {
  const out = {};
  for (const m of TRIAL_MODALITIES) out[m.key] = trialYears(await allStudies(m.query));
  return out;
};

const postJson = (url, body) => retry(async () => {
  const res = await fetch(url, { method: 'POST', headers: { ...HEADERS, 'content-type': 'application/json' }, body: JSON.stringify(body) });
  if (!res.ok) throw new Error(`HTTP ${res.status} ${url}`);
  return res.json();
}, { tries: 5, wait: 3000 });

const ancestry = async () => gwasAncestry(await text(FILES.gwasAncestry));

const ancient = async () => {
  const anno = aadrAnnoFile(await json(FILES.aadrDataset));
  if (!anno) throw new Error('No 1240K annotation file in the latest AADR release');
  return { version: anno.version, released: anno.released, ...aadrSummary(await text(`${FILES.dataverseFile}${anno.id}`)) };
};

const igsr = async () => {
  const byName = new Map();
  let total = Infinity;
  for (let from = 0; from < total; from += 500) {
    const page = await postJson(FILES.igsrSamples, {
      size: 500, from, _source: ['name', 'sex', 'populations.code', 'populations.name', 'populations.superpopulationCode', 'populations.superpopulationName', 'dataCollections.title'],
    });
    const { hits } = page.hits;
    total = typeof page.hits.total === 'object' ? page.hits.total.value : page.hits.total;
    if (!hits.length) break;
    for (const hit of hits) byName.set(hit._source.name, hit);
  }
  return igsrSamples({ hits: { total, hits: [...byName.values()] } }, '1000 Genomes 30x on GRCh38');
};

const assemblies = async () => {
  const reports = [];
  let token = null;
  do {
    const page = await json(assemblyPage(token));
    reports.push(...(page.reports ?? []));
    token = page.next_page_token ?? null;
  } while (token);
  return assemblySummary(reports);
};

const zipFiles = async (url) => unzip(await readFile(await file(url)));

const zipEntry = (files, pattern) => {
  const name = Object.keys(files).find((n) => pattern.test(n));
  if (!name) throw new Error(`No file matching ${pattern} in the download`);
  return files[name].toString('utf8');
};

const intogenCompendium = async () => {
  const latest = latestIntogen(await text(FILES.intogenDownload));
  if (!latest) throw new Error('No IntOGen drivers and cohorts files on the download page');
  return { latest, drivers: intogenDrivers(zipEntry(await zipFiles(intogenFile(latest.drivers)), /Compendium_Cancer_Genes\.tsv$/)) };
};

const intogen = async () => {
  const { latest, drivers } = await intogenCompendium();
  return { release: latest.release, drivers, cohorts: intogenCohorts(zipEntry(await zipFiles(intogenFile(latest.cohorts)), /cohorts\.tsv$/)) };
};

const TCGA_GENES = 20;

const mutatedCases = async (gene) => {
  const hits = [];
  let total = Infinity;
  for (let from = 0; from < total; from += 5000) {
    const page = await json(gdcSsmPage(gene, from));
    total = page.data.pagination.total;
    if (!page.data.hits.length) break;
    hits.push(...page.data.hits);
  }
  return casesByProject(hits);
};

const gdc = async () => {
  const programs = gdcPrograms(await json(FILES.gdcProjects));
  const tested = Object.fromEntries(facetCounts(await json(FILES.gdcTested), 'project.project_id').map((b) => [b.key, b.n]));
  const drivers = (await intogenCompendium()).drivers.top.slice(0, TCGA_GENES).map((g) => g.gene);
  const mutated = {};
  for (const gene of drivers) mutated[gene] = await mutatedCases(gene);
  return {
    ...gdcRelease(await json(FILES.gdcStatus)),
    projects: programs.projects,
    cases: programs.cases,
    programs: programs.programs,
    sites: facetCounts(await json(FILES.gdcSites), 'primary_site'),
    tcga: { genes: drivers, projects: programs.tcga.map((p) => ({ ...p, tested: tested[p.project] ?? 0 })), mutated },
  };
};

const cpic = async () => {
  const alleles = await json(FILES.cpicAlleles);
  const results = await json(FILES.cpicResults);
  return {
    pairs: cpicPairs(await json(FILES.cpicPairs)),
    alleles: cpicAlleles(alleles),
    phenotypes: cpicPhenotypes(results),
    cyp2d6: cyp2d6(alleles.filter((a) => a.genesymbol === 'CYP2D6'), results.filter((r) => r.genesymbol === 'CYP2D6')),
  };
};

const fdaPgx = async () => {
  const table = fdaBiomarkers(await text(FILES.fdaPgx));
  if (!table?.rows) throw new Error('No rows in the FDA table of pharmacogenomic biomarkers');
  return table;
};

const clinpgx = async () => {
  const annotations = await zipFiles(FILES.clinpgxAnnotations);
  const labels = await zipFiles(FILES.clinpgxLabels);
  return {
    annotations: { created: createdDate(Object.keys(annotations)), ...clinpgxAnnotations(zipEntry(annotations, /^clinical_annotations\.tsv$/)) },
    labels: { created: createdDate(Object.keys(labels)), ...clinpgxLabels(zipEntry(labels, /^drugLabels\.tsv$/)) },
  };
};

const gwas = async () => {
  const acc = gwasAssociations({ parents: traitParents(await text(FILES.gwasTraitMappings)) });
  for await (const line of await zipLines(await file(FILES.gwasAssociations), /\.tsv$/)) if (line) acc.add(line);
  return { ...acc.result(), sizes: studySizes(acc.perStudy, await text(FILES.gwasAncestry)) };
};

const pgs = async () => {
  const info = (await json(FILES.pgsInfo)).latest_release ?? {};
  const scores = await text(FILES.pgsScores);
  return {
    catalog: {
      release: info.date ?? null, scores: info.scores ?? null, publications: info.publications ?? null, traits: info.traits ?? null,
    },
    ...pgsScores(scores),
    ancestry: pgsAncestry(await text(FILES.pgsDevelopment), await text(FILES.pgsEvaluation)),
    performance: pgsPerformance(await text(FILES.pgsPerformance), scores, pgsCategories(await json(FILES.pgsCategories))),
  };
};

const panukb = async () => heritability(await text(FILES.panukbH2), await text(FILES.panukbPhenotypes));

const SECTIONS = {
  assembly, sequence, repeats, genes, catalogs, proteins, targets, expression, networks, pathways, ontology, regulation, trnas, encode, structure, variation, constraint, ancestry, ancient, igsr, assemblies, clinvar, clinvarTable, domains, clinvarGrowth, validity, regions, expansions, hpo, orphanet, gtr, fda, trials, intogen, gdc, cpic, fdaPgx, clinpgx, gwas, pgs, panukb,
};

const now = new Date();
const only = process.argv.slice(2).filter((a) => SECTIONS[a]);
const previous = existsSync(OUT) ? (await import(`${OUT.href}?t=${now.getTime()}`)).default : {};
const snapshot = await runSections(SECTIONS, { now, previous, only });
await writeFile(OUT, toModule(snapshot));
process.stderr.write(`wrote ${OUT.pathname}\n`);
