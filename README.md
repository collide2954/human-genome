# Human Genome in Numbers

One static page about the human genome: the GRCh38 and T2T-CHM13 assemblies, their bands, gaps and alternate sequence, and, part by part, genes, proteins, variation, clinical genetics and the research behind them. Hand-built SVG charts with no runtime dependencies, each with its numbers as a table and a CSV.

Lenses at the top filter the page: `?lens=sequence`, `?lens=genes`, `?lens=variation` or `?lens=clinic`.

## Run it

| Command | Does |
| --- | --- |
| `npm run serve` | Serves the page at http://127.0.0.1:8080/ |
| `npm run lint` | ESLint, Stylelint and html-validate, each with a no-comments rule |
| `npm test` | Lint, unit tests, and Playwright end-to-end tests at desktop and phone widths |
| `npm run snapshot` | Rebuilds `src/snapshot.js` from the files and APIs listed in `src/queries.js` |
| `npm run snapshot -- assembly` | Rebuilds only the named sections |

Set `SNAPSHOT_CACHE` to a folder to keep the downloaded files between runs. The end-to-end tests need a Chromium for Playwright: `npx playwright install chromium`.

## Publish it

The `Publish site` workflow deploys to `https://collide2954.github.io/human-genome/` on every push to `master` and after each weekly snapshot (Settings → Pages → Source: GitHub Actions). `.github/workflows/snapshot.yml` rebuilds the data every Sunday and commits it when numbers change. A section whose source fails keeps last week's numbers and its own date, and the run log names it.

## Where things live

- `src/queries.js`: every file and API behind the numbers.
- `scripts/parsers.mjs`: readers for NCBI assembly reports and statistics and the UCSC tables. `scripts/streams.mjs` streams and unzips the large files; `scripts/sections.mjs` runs each section and keeps the last good values when one fails.
- `scripts/genes.mjs`, `catalogs.mjs`, `sequence.mjs`, `proteins.mjs`, `networks.mjs`, `regulation.mjs`, `structure.mjs` and `variation.mjs`: readers for GENCODE, HGNC, MANE, the hg38 2bit file, UniProt, HPA, GTEx, STRING, Reactome, GO, ENCODE, dbSNP, dbVar and gnomAD.
- `scripts/clinical.mjs`: ClinVar's VCF, archive and tab-delimited files, ClinGen's validity and dosage region curations, InterPro domains and STRchive's repeat loci. `scripts/rare.mjs`: HPO, Orphanet, the Genetic Testing Registry, FDA's list of cell and gene therapies and ClinicalTrials.gov.
- `src/findings.js`: the section headlines, computed from the data so they stay true after each update.
- `src/facts.js`: sources, decoder examples and the few fixed values taken from cited papers: gene therapy approvals checked against FDA package inserts, diagnostic yields and prenatal screening rates.
- `src/charts.js`: the SVG chart engine and animations. `src/stats.js`: quantiles, binning, Lorenz curves and interval merging.
- `scripts/snapshot.mjs`: the weekly data pull. `scripts/export.mjs` writes the snapshot out as `snapshot.json`, which the published site links for download.
- `og.png`: the preview image shown when the link is shared.

## Reading the numbers

Positions are GRCh38.p14 coordinates counted from 1 with both ends included; UCSC tables count from 0 and are converted. Chromosome totals cover 1 to 22, X and Y unless the mitochondrial genome is named. "Sequenced bases" means bases that are not N. T2T-CHM13v2.0 is the CHM13 cell line with chromosome Y from HG002, so chromosome-by-chromosome differences mix assembly gains with differences between people.

ClinVar counts use the aggregate germline classification of each record in the GRCh38 VCF, grouped as P, P/LP, LP, VUS, LB, B/LB, B, conflicting and other; records without a germline classification are counted apart. Captions that say "fixed values" come from the cited paper and change only when the page is edited.
