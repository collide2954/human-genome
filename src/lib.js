export const toNumber = (v) => {
  if (v === null || v === undefined || v === '') return NaN;
  return typeof v === 'number' ? v : Number(String(v).replace(/,/g, ''));
};

export const exact = (v) => {
  const n = toNumber(v);
  return Number.isFinite(n) ? Math.round(n).toLocaleString('en-US') : '—';
};

const UNITS = [
  [1e12, 'T'],
  [1e9, 'B'],
  [1e6, 'M'],
  [1e3, 'K'],
];

export const compact = (v) => {
  const n = toNumber(v);
  if (!Number.isFinite(n)) return '—';
  const scaledBy = (size) => Math.round((n / size) * 10) / 10;
  for (const [i, [size, unit]] of UNITS.entries()) {
    if (Math.abs(n) >= size) {
      const scaled = scaledBy(size);
      return Math.abs(scaled) >= 1000 && i > 0 ? `${scaledBy(UNITS[i - 1][0])}${UNITS[i - 1][1]}` : `${scaled}${unit}`;
    }
  }
  return Math.abs(Math.round(n)) >= 1000 ? `${scaledBy(1e3)}K` : String(Math.round(n));
};

export const pct = (ratio) => {
  const v = ratio * 100;
  const places = v !== 0 && Math.abs(v) < 0.1 ? 2 : 1;
  return `${Number(v.toFixed(places))}%`;
};

export const usd = (v) => `$${compact(v)}`;

export const decimal = (v, places = 1) => {
  const n = toNumber(v);
  return Number.isFinite(n) ? n.toFixed(places) : '—';
};

const WORD_UNITS = [[1e12, 'trillion'], [1e9, 'billion'], [1e6, 'million']];

export const wordy = (n) => {
  const unit = WORD_UNITS.find(([size]) => n >= size);
  if (!unit) return exact(n);
  const v = n / unit[0];
  const digits = v >= 100 ? 0 : v >= 10 ? 1 : 2;
  return `${Number(v.toFixed(digits))} ${unit[1]}`;
};

export const ordinal = (n) => {
  const suffix = ['th', 'st', 'nd', 'rd'];
  const v = n % 100;
  return `${n}${suffix[(v - 20) % 10] || suffix[v] || suffix[0]}`;
};

const BASE_UNITS = [[1e9, 'Gb', 2], [1e6, 'Mb', 1], [1e3, 'kb', 1]];

export const bases = (v) => {
  const n = toNumber(v);
  if (!Number.isFinite(n)) return '—';
  if (Math.abs(n) < 1e4) return `${exact(n)} bp`;
  const at = BASE_UNITS.findIndex(([size]) => Math.abs(n) >= size);
  const [size, unit, places] = BASE_UNITS[at];
  const scaled = Number((n / size).toFixed(places));
  if (Math.abs(scaled) >= 1000 && at > 0) {
    const [upSize, upUnit, upPlaces] = BASE_UNITS[at - 1];
    return `${Number((n / upSize).toFixed(upPlaces))} ${upUnit}`;
  }
  return `${scaled} ${unit}`;
};

const SUPERSCRIPT = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };

export const sci = (p) => {
  if (p === 0) return '0';
  if (p >= 0.001) return String(Number(p.toPrecision(2)));
  const power = Math.floor(Math.log10(p));
  let mantissa = Number((p / 10 ** power).toFixed(1));
  let exponent = power;
  if (mantissa >= 10) {
    mantissa = 1;
    exponent += 1;
  }
  return `${mantissa}×10${[...String(exponent)].map((c) => SUPERSCRIPT[c]).join('')}`;
};

export const mlog = (p) => -Math.log10(p);

const chromName = (chrom) => {
  const bare = String(chrom).replace(/^chr/i, '');
  return `chr${bare === 'MT' ? 'M' : bare}`;
};

export const coords = (chrom, start, end) => `${chromName(chrom)}:${exact(start)}-${exact(end)}`;

export const density = (count, length) => `${Number(((count / length) * 1e6).toFixed(1))} per Mb`;

export const signed = (v, format = exact) => {
  if (v === 0) return '0';
  return v > 0 ? `+${format(v)}` : `−${format(-v)}`;
};

const csvCell = (v) => {
  if (v === null || v === undefined) return '';
  const s = String(v);
  return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
};

export const toCsv = (columns, rows) => `${[columns, ...rows].map((r) => r.map(csvCell).join(',')).join('\r\n')}\r\n`;

export const retry = async (fn, { tries = 3, wait = 1000 } = {}) => {
  for (let attempt = 1; ; attempt += 1) {
    try {
      return await fn();
    } catch (error) {
      if (attempt >= tries) throw error;
      await new Promise((r) => setTimeout(r, wait * 2 ** (attempt - 1)));
    }
  }
};

const ENSEMBL_FEATURES = { G: 'gene', T: 'transcript', P: 'protein', E: 'exon', R: 'regulatory feature' };

const REFSEQ_PREFIXES = {
  NC: ['RefSeq chromosome', 'complete genomic molecule', 'nuccore'],
  NG: ['RefSeqGene region', 'genomic region kept as a gene reference', 'nuccore'],
  NT: ['RefSeq contig', 'genomic contig', 'nuccore'],
  NW: ['RefSeq contig', 'genomic contig from whole-genome shotgun data', 'nuccore'],
  NM: ['RefSeq mRNA', 'curated mRNA', 'nuccore'],
  NR: ['RefSeq noncoding RNA', 'curated noncoding RNA', 'nuccore'],
  XM: ['RefSeq model mRNA', 'predicted mRNA', 'nuccore'],
  XR: ['RefSeq model noncoding RNA', 'predicted noncoding RNA', 'nuccore'],
  NP: ['RefSeq protein', 'curated protein', 'protein'],
  XP: ['RefSeq model protein', 'predicted protein', 'protein'],
  YP: ['RefSeq protein', 'protein with no transcript record', 'protein'],
};

const GRCH38_VERSIONS = [11, 12, 12, 12, 10, 12, 14, 11, 12, 11, 10, 12, 11, 9, 10, 10, 11, 10, 10, 11, 9, 11, 11, 10];
const CHROMOSOMES = [...Array.from({ length: 22 }, (_, i) => String(i + 1)), 'X', 'Y'];

export const refseqMolecule = (accession) => {
  const m = String(accession).match(/^NC_(\d{6})\.(\d+)$/);
  if (!m) return null;
  const n = Number(m[1]);
  const version = Number(m[2]);
  if (n === 12920 && version === 1) return { molecule: 'Mitochondrion (rCRS)', chrom: 'MT', assembly: 'GRCh38' };
  if (n >= 1 && n <= 24) {
    const chrom = CHROMOSOMES[n - 1];
    const v38 = GRCH38_VERSIONS[n - 1];
    const assembly = version === v38 ? 'GRCh38' : version === v38 - 1 ? 'GRCh37' : null;
    return { molecule: `Chromosome ${chrom}`, chrom, assembly };
  }
  if (n >= 60925 && n <= 60948 && version === 1) {
    const chrom = CHROMOSOMES[n - 60925];
    return { molecule: `Chromosome ${chrom}`, chrom, assembly: 'T2T-CHM13v2.0' };
  }
  return null;
};

const ASSEMBLIES = {
  'GCF_000001405.40': 'GRCh38.p14',
  'GCA_000001405.29': 'GRCh38.p14',
  'GCF_000001405.39': 'GRCh38.p13',
  'GCF_000001405.26': 'GRCh38',
  'GCF_000001405.25': 'GRCh37.p13',
  'GCF_009914755.1': 'T2T-CHM13v2.0',
  'GCA_009914755.4': 'T2T-CHM13v2.0',
};

const MIM_PREFIXES = {
  '*': 'gene',
  '#': 'phenotype, molecular basis known',
  '+': 'gene and phenotype, combined',
  '%': 'phenotype or locus, molecular basis unknown',
  '^': 'entry removed or moved',
};

const MIM_RANGES = {
  1: 'Autosomal, entered before May 1994',
  2: 'Autosomal, entered before May 1994',
  3: 'X-linked',
  4: 'Y-linked',
  5: 'Mitochondrial',
  6: 'Autosomal, entered since May 1994',
};

const HLA_EXPRESSION = {
  N: 'null, not expressed',
  L: 'low cell-surface expression',
  S: 'secreted, not on the cell surface',
  C: 'in the cytoplasm, not on the cell surface',
  A: 'aberrant expression',
  Q: 'expression questionable',
};

const HGVS_SYSTEMS = {
  c: 'c. coding DNA, numbered from the A of the ATG start codon',
  g: 'g. linear genomic DNA',
  m: 'm. mitochondrial DNA',
  n: 'n. noncoding DNA',
  r: 'r. RNA',
  p: 'p. protein',
};

const ONTOLOGIES = {
  GO: ['Gene Ontology term', 'Gene Ontology', (id) => `https://amigo.geneontology.org/amigo/term/${id}`],
  HP: ['Human Phenotype Ontology term', 'Human Phenotype Ontology', (id) => `https://hpo.jax.org/browse/term/${id}`],
  MONDO: ['MONDO disease', 'Mondo Disease Ontology', (id) => `https://monarchinitiative.org/${id}`],
  EFO: ['Experimental Factor Ontology term', 'Experimental Factor Ontology', (id) => `https://www.ebi.ac.uk/ols4/ontologies/efo/classes?obo_id=${id}`],
};

const ENCODE_OBJECTS = {
  SR: 'ENCODE experiment', FF: 'ENCODE file', BS: 'ENCODE biosample', DO: 'ENCODE donor', AB: 'ENCODE antibody', AN: 'ENCODE annotation',
};

const regionOf = (change) => {
  const m = change.match(/^([-*]?)(\d+)([+-]\d+)?/);
  if (!m) return 'coding sequence';
  if (m[1] === '-') return '5′ UTR, before the start codon';
  if (m[1] === '*') return '3′ UTR, after the stop codon';
  if (m[3]) return Math.abs(Number(m[3])) <= 10 ? 'intron, next to an exon' : 'intron';
  return Number(m[2]) <= 3 ? 'start codon' : 'coding sequence';
};

export const hgvsChange = (system, change) => {
  const body = change.replace(/^\((.*)\)$/, '$1');
  if (system === 'p') {
    if (/fs/.test(body)) return { kind: 'frameshift' };
    if (/=$/.test(body)) return { kind: 'synonymous' };
    if (/(Ter|\*)$/.test(body)) return { kind: 'nonsense' };
    if (/delins/.test(body)) return { kind: 'deletion-insertion' };
    if (/del/.test(body)) return { kind: 'deletion' };
    if (/dup/.test(body)) return { kind: 'duplication' };
    if (/ins/.test(body)) return { kind: 'insertion' };
    if (/ext/.test(body)) return { kind: 'extension' };
    return { kind: 'missense' };
  }
  const kinds = [[/delins/, 'deletion-insertion'], [/del/, 'deletion'], [/dup/, 'duplication'], [/ins/, 'insertion'], [/inv/, 'inversion'], [/\[/, 'repeated sequence'], [/>/, 'substitution'], [/=/, 'no change']];
  const kind = kinds.find(([re]) => re.test(body))?.[1] ?? 'other';
  return system === 'c' ? { kind, region: regionOf(body) } : { kind };
};

const describeReference = (ref) => {
  const mol = refseqMolecule(ref);
  if (mol?.assembly) return `${ref}, ${mol.molecule.toLowerCase()}, ${mol.assembly}`;
  const prefix = REFSEQ_PREFIXES[ref.slice(0, 2)];
  if (prefix && ref[2] === '_') return `${ref}, ${prefix[0]}`;
  if (/^ENS/.test(ref)) return `${ref}, Ensembl ${ENSEMBL_FEATURES[ref[3]]}`;
  return ref;
};

const ucsc = (chrom, start, end) => `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${encodeURIComponent(`chr${chrom === 'MT' ? 'M' : chrom}:${start}-${end}`)}`;

const validChrom = (c) => CHROMOSOMES.includes(c) || c === 'M' || c === 'MT';

export const bandInterval = (bands, chrom, name) => {
  const hits = (bands ?? []).filter((b) => b.chrom === chrom && (b.name === name || b.name.startsWith(name)));
  if (!hits.length) return null;
  return { chrom, start: Math.min(...hits.map((b) => b.start)), end: Math.max(...hits.map((b) => b.end)) };
};

const rule = (re, build) => ({ re, build });

const RULES = [
  rule(/^([A-Z]{2}_\d+\.\d+|LRG_\d+(?:t\d+|p\d+)?|ENS[GTP]\d{11}\.\d+)(?:\([A-Za-z0-9-]+\))?:([cgmnrp])\.(.+)$/, (m) => {
    const [, ref, system, change] = m;
    const found = hgvsChange(system, change);
    const predicted = system === 'p' && /^\(.*\)$/.test(change);
    return {
      type: 'HGVS variant', id: m[0], db: 'HGVS sequence variant nomenclature', url: `https://mutalyzer.nl/normalizer/${encodeURIComponent(m[0])}`,
      note: 'Reference sequence, coordinate system, position and change',
      parts: [
        ['Reference', describeReference(ref)],
        ['Coordinates', predicted ? 'p. protein, predicted' : HGVS_SYSTEMS[system]],
        ['Change', found.kind],
        found.region ? ['Region', found.region] : null,
      ].filter(Boolean),
    };
  }),
  rule(/^([A-Z]{2}_\d+\.\d+):(\d+):([ACGTN]*):([ACGTN]*)$/i, (m) => {
    const ref = m[1].toUpperCase();
    const pos = Number(m[2]);
    const id = `${ref}:${pos}:${m[3].toUpperCase()}:${m[4].toUpperCase()}`;
    return {
      type: 'SPDI variant', id, db: 'NCBI Variation Services', url: `https://api.ncbi.nlm.nih.gov/variation/v0/spdi/${id}/hgvs`,
      note: 'Sequence, position, deleted and inserted bases',
      parts: [['Sequence', describeReference(ref)], ['Position', `${exact(pos)} counted from 0, so ${exact(pos + 1)} counted from 1`], ['Deleted', m[3].toUpperCase() || 'nothing'], ['Inserted', m[4].toUpperCase() || 'nothing']],
    };
  }),
  rule(/^(?:chr)?(\d{1,2}|X|Y|M|MT):([\d,]+)(?:-([\d,]+))?$/i, (m) => {
    const chrom = m[1].toUpperCase();
    if (!validChrom(chrom)) return null;
    const start = Number(m[2].replace(/,/g, ''));
    const end = m[3] ? Number(m[3].replace(/,/g, '')) : start;
    if (end < start || start < 1) return null;
    return {
      type: 'GRCh38 interval', id: coords(chrom, start, end), db: 'GRCh38 reference coordinates', url: ucsc(chrom, start, end),
      note: 'Chromosome, then first and last base, counted from 1',
      parts: [['Chromosome', chrom === 'MT' ? 'M' : chrom], ['Length', bases(end - start + 1)]],
    };
  }),
  rule(/^(\d{1,2}|X|Y|M|MT)-(\d+)-([ACGT]+)-([ACGT]+)$/i, (m) => {
    const chrom = m[1].toUpperCase();
    if (!validChrom(chrom)) return null;
    const id = `${chrom}-${m[2]}-${m[3].toUpperCase()}-${m[4].toUpperCase()}`;
    const snv = m[3].length === 1 && m[4].length === 1;
    return {
      type: 'gnomAD variant', id, db: 'gnomAD', url: `https://gnomad.broadinstitute.org/variant/${id}?dataset=gnomad_r4`,
      note: 'Chromosome, GRCh38 position, reference and alternate alleles',
      parts: [['Chromosome', chrom], ['Position', exact(Number(m[2]))], ['Change', `${m[3].toUpperCase()}>${m[4].toUpperCase()}`], ['Kind', snv ? 'single-nucleotide variant' : 'insertion or deletion']],
    };
  }),
  rule(/^HLA-([A-Z0-9]+)\*(\d{2,3})(?::(\d{2,3}))?(?::(\d{2,3}))?(?::(\d{2,3}))?([NLSCAQ])?$/i, (m) => {
    const gene = `HLA-${m[1].toUpperCase()}`;
    const suffix = m[6]?.toUpperCase();
    const id = `${gene}*${[m[2], m[3], m[4], m[5]].filter(Boolean).join(':')}${suffix ?? ''}`;
    return {
      type: 'HLA allele', id, db: 'IPD-IMGT/HLA', url: 'https://hla.alleles.org/nomenclature/naming.html',
      note: 'Gene, allele group, protein, synonymous change, noncoding change, expression',
      parts: [
        ['Gene', gene], ['Allele group', m[2]], m[3] ? ['Protein', m[3]] : null, m[4] ? ['Synonymous', m[4]] : null,
        m[5] ? ['Noncoding', m[5]] : null, suffix ? ['Expression', `${suffix}: ${HLA_EXPRESSION[suffix]}`] : null,
      ].filter(Boolean),
    };
  }),
  rule(/^([A-Z][A-Z0-9]{1,9})\*(\d+[A-Z]?)(?:\.(\d{3}))?$/i, (m) => {
    const gene = m[1].toUpperCase();
    const allele = `*${m[2].toUpperCase()}`;
    return {
      type: 'PharmVar star allele', id: `${gene}${allele}${m[3] ? `.${m[3]}` : ''}`, db: 'Pharmacogene Variation Consortium', url: `https://www.pharmvar.org/gene/${gene}`,
      note: 'Gene, star allele, suballele', parts: [['Gene', gene], ['Allele', allele], m[3] ? ['Suballele', m[3]] : null].filter(Boolean),
    };
  }),
  rule(/^ENS([GTPER])(\d{11})(?:\.(\d+))?$/i, (m) => {
    const letter = m[1].toUpperCase();
    const stable = `ENS${letter}${m[2]}`;
    const feature = ENSEMBL_FEATURES[letter];
    return {
      type: `Ensembl ${feature}`, id: m[3] ? `${stable}.${m[3]}` : stable, db: 'Ensembl and GENCODE',
      url: 'GTP'.includes(letter) ? `https://www.ensembl.org/id/${stable}` : `https://www.ensembl.org/Homo_sapiens/Search/Results?q=${stable}`,
      note: 'ENS, feature letter, eleven digits, version', parts: [['Feature', feature], m[3] ? ['Version', m[3]] : null].filter(Boolean),
    };
  }),
  rule(/^(N[CGTWMRP]|X[MRP]|YP)_(\d+)(?:\.(\d+))?$/i, (m) => {
    const prefix = m[1].toUpperCase();
    const id = `${prefix}_${m[2]}${m[3] ? `.${m[3]}` : ''}`;
    const [type, meaning, db] = REFSEQ_PREFIXES[prefix];
    const mol = refseqMolecule(id);
    return {
      type, id, db: 'NCBI RefSeq', url: `https://www.ncbi.nlm.nih.gov/${db}/${id}`, note: 'Prefix, number, version',
      parts: [
        ['Prefix', `${prefix}_: ${meaning}`], m[3] ? ['Version', m[3]] : null,
        mol ? ['Molecule', mol.molecule] : null, mol?.assembly ? ['Assembly', mol.assembly] : null,
      ].filter(Boolean),
    };
  }),
  rule(/^(GC[AF])_(\d{9})\.(\d+)$/i, (m) => {
    const id = `${m[1].toUpperCase()}_${m[2]}.${m[3]}`;
    return {
      type: 'NCBI assembly', id, db: 'NCBI Datasets', url: `https://www.ncbi.nlm.nih.gov/datasets/genome/${id}/`, note: 'GCA for GenBank, GCF for RefSeq, then number and version',
      parts: [['Database', m[1].toUpperCase() === 'GCF' ? 'RefSeq' : 'GenBank'], ASSEMBLIES[id] ? ['Name', ASSEMBLIES[id]] : null].filter(Boolean),
    };
  }),
  rule(/^LRG_(\d+)(t\d+|p\d+)?$/i, (m) => ({
    type: 'Locus Reference Genomic', id: `LRG_${m[1]}${m[2] ?? ''}`, db: 'LRG', url: `https://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_${m[1]}.xml`,
    note: 'A fixed genomic reference for reporting variants', parts: m[2] ? [['Part', m[2][0] === 't' ? `transcript ${m[2].slice(1)}` : `protein ${m[2].slice(1)}`]] : [],
  })),
  rule(/^(rs|ss)(\d+)$/i, (m) => {
    const kind = m[1].toLowerCase();
    const id = `${kind}${m[2]}`;
    return kind === 'rs'
      ? { type: 'dbSNP RefSNP', id, db: 'dbSNP', url: `https://www.ncbi.nlm.nih.gov/snp/${id}`, note: 'One ID per variant position, merging every submission', parts: [] }
      : { type: 'dbSNP submitted SNP', id, db: 'dbSNP', url: `https://www.ncbi.nlm.nih.gov/snp/?term=${id}`, note: 'One ID per submission', parts: [] };
  }),
  rule(/^(VCV|RCV|SCV)(\d{9})(?:\.(\d+))?$/i, (m) => {
    const kind = m[1].toUpperCase();
    const id = `${kind}${m[2]}${m[3] ? `.${m[3]}` : ''}`;
    const types = {
      VCV: ['ClinVar variation', 'Every submission about one variant', `https://www.ncbi.nlm.nih.gov/clinvar/variation/${Number(m[2])}/`],
      RCV: ['ClinVar record', 'One variant and one condition', `https://www.ncbi.nlm.nih.gov/clinvar/${kind}${m[2]}/`],
      SCV: ['ClinVar submission', 'One submitter’s classification', `https://www.ncbi.nlm.nih.gov/clinvar/?term=${kind}${m[2]}`],
    };
    const [type, note, url] = types[kind];
    return { type, id, db: 'ClinVar', url, note, parts: m[3] ? [['Version', m[3]]] : [] };
  }),
  rule(/^(nsv|esv|nssv|essv)(\d+)$/i, (m) => {
    const kind = m[1].toLowerCase();
    const id = `${kind}${m[2]}`;
    return {
      type: kind.endsWith('ssv') ? 'dbVar variant call' : 'dbVar variant region', id, db: 'dbVar', url: `https://www.ncbi.nlm.nih.gov/dbvar/variants/${id}/`,
      note: kind.startsWith('n') ? 'Accessioned by NCBI' : 'Accessioned by EBI', parts: [],
    };
  }),
  rule(/^(nstd|estd)(\d+)$/i, (m) => {
    const id = `${m[1].toLowerCase()}${m[2]}`;
    return { type: 'dbVar study', id, db: 'dbVar', url: `https://www.ncbi.nlm.nih.gov/dbvar/studies/${id}/`, note: 'A structural variant study', parts: [] };
  }),
  rule(/^CA(\d+)$/i, (m) => ({
    type: 'ClinGen canonical allele', id: `CA${m[1]}`, db: 'ClinGen Allele Registry',
    url: `https://reg.clinicalgenome.org/redmine/projects/registry/genboree_registry/by_caid?caid=CA${m[1]}`, note: 'One ID for the same allele across assemblies and transcripts', parts: [],
  })),
  rule(/^COS([VM])(\d+)$/i, (m) => ({
    type: m[1].toUpperCase() === 'V' ? 'COSMIC genomic variant' : 'COSMIC legacy mutation', id: `COS${m[1].toUpperCase()}${m[2]}`, db: 'COSMIC',
    url: `https://cancer.sanger.ac.uk/cosmic/search?q=COS${m[1].toUpperCase()}${m[2]}`, note: 'Somatic variants seen in cancers', parts: [],
  })),
  rule(/^HGNC:(\d+)$/i, (m) => ({
    type: 'HGNC gene', id: `HGNC:${m[1]}`, db: 'HUGO Gene Nomenclature Committee', url: `https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:${m[1]}`,
    note: 'Stable ID behind an approved gene symbol', parts: [],
  })),
  rule(/^(?:GeneID|NCBIGene|EntrezGene|Entrez|Gene ?ID)[: ]?(\d+)$/i, (m) => ({
    type: 'NCBI Gene', id: `GeneID:${m[1]}`, db: 'NCBI Gene', url: `https://www.ncbi.nlm.nih.gov/gene/${m[1]}`, note: 'Integer ID for a gene record', parts: [],
  })),
  rule(/^(?:(?:MIM|OMIM)[:\s#*]*|([*#+%^]))(\d{6})$/i, (m) => {
    const range = MIM_RANGES[m[2][0]];
    if (!range) return null;
    return {
      type: 'OMIM entry', id: `MIM:${m[2]}`, db: 'Online Mendelian Inheritance in Man', url: `https://www.omim.org/entry/${m[2]}`,
      note: 'Six digits; the first gives the inheritance class', parts: [m[1] ? ['Prefix', `${m[1]}: ${MIM_PREFIXES[m[1]]}`] : null, ['Number', range]].filter(Boolean),
    };
  }),
  rule(/^(GO|HP|MONDO|EFO)[:_](\d{7})$/i, (m) => {
    const prefix = m[1].toUpperCase();
    const id = `${prefix}:${m[2]}`;
    const [type, db, url] = ONTOLOGIES[prefix];
    return { type, id, db, url: url(id), note: 'Ontology prefix and seven digits', parts: [] };
  }),
  rule(/^ORPHA(?:NET)?[:_ ]?(\d+)$/i, (m) => ({
    type: 'Orphanet disease', id: `ORPHA:${m[1]}`, db: 'Orphanet', url: `https://www.orpha.net/en/disease/detail/${m[1]}`, note: 'Rare disease number', parts: [],
  })),
  rule(/^GCST(\d{6,9})$/i, (m) => ({
    type: 'GWAS Catalog study', id: `GCST${m[1]}`, db: 'NHGRI-EBI GWAS Catalog', url: `https://www.ebi.ac.uk/gwas/studies/GCST${m[1]}`, note: 'One study accession per trait analysed', parts: [],
  })),
  rule(/^PMID:? ?(\d{1,9})$/i, (m) => ({
    type: 'PubMed record', id: `PMID:${m[1]}`, db: 'PubMed', url: `https://pubmed.ncbi.nlm.nih.gov/${m[1]}/`, note: 'PubMed identifier', parts: [],
  })),
  rule(/^ENC(SR|FF|BS|DO|AB|AN)(\d{3})([A-Z]{3})$/i, (m) => {
    const id = `ENC${m[1].toUpperCase()}${m[2]}${m[3].toUpperCase()}`;
    return { type: ENCODE_OBJECTS[m[1].toUpperCase()], id, db: 'ENCODE', url: `https://www.encodeproject.org/${id}/`, note: 'ENC, object type, three digits, three letters', parts: [] };
  }),
  rule(/^EH38E(\d{7})$/i, (m) => ({
    type: 'ENCODE candidate cis-regulatory element', id: `EH38E${m[1]}`, db: 'ENCODE Registry of cCREs', url: `https://screen.wenglab.org/search?assembly=GRCh38&q=EH38E${m[1]}`,
    note: 'EH38 for GRCh38, E for element', parts: [],
  })),
  rule(/^R-HSA-(\d+)$/i, (m) => ({
    type: 'Reactome pathway', id: `R-HSA-${m[1]}`, db: 'Reactome', url: `https://reactome.org/content/detail/R-HSA-${m[1]}`, note: 'R, species code HSA, number', parts: [],
  })),
  rule(/^CCDS(\d+)(?:\.(\d+))?$/i, (m) => ({
    type: 'Consensus CDS', id: `CCDS${m[1]}${m[2] ? `.${m[2]}` : ''}`, db: 'CCDS project', url: `https://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi?REQUEST=CCDS&DATA=CCDS${m[1]}`,
    note: 'A coding sequence agreed by NCBI, Ensembl, UCSC and HGNC', parts: m[2] ? [['Version', m[2]]] : [],
  })),
  rule(/^DB(\d{5})$/i, (m) => ({
    type: 'DrugBank drug', id: `DB${m[1]}`, db: 'DrugBank', url: `https://go.drugbank.com/drugs/DB${m[1]}`, note: 'DB and five digits', parts: [],
  })),
  rule(/^NBK(\d+)$/i, (m) => ({
    type: 'NCBI Bookshelf chapter', id: `NBK${m[1]}`, db: 'NCBI Bookshelf, including GeneReviews', url: `https://www.ncbi.nlm.nih.gov/books/NBK${m[1]}/`, note: 'Bookshelf accession', parts: [],
  })),
  rule(/^C(\d{7})$/, (m) => ({
    type: 'MedGen concept', id: `C${m[1]}`, db: 'MedGen, using UMLS concept IDs', url: `https://www.ncbi.nlm.nih.gov/medgen/C${m[1]}`, note: 'C and seven digits', parts: [],
  })),
  rule(/^(?:chr)?(\d{1,2}|X|Y)([pq])(\d)(\d)?(?:\.(\d+))?$/i, (m, ctx) => {
    const chrom = m[1].toUpperCase();
    if (!CHROMOSOMES.includes(chrom)) return null;
    const arm = m[2].toLowerCase();
    const name = `${arm}${m[3]}${m[4] ?? ''}${m[5] ? `.${m[5]}` : ''}`;
    const span = bandInterval(ctx.bands, chrom, name);
    return {
      type: 'Cytogenetic band', id: `${chrom}${name}`, db: 'ISCN banding on GRCh38',
      url: span ? ucsc(chrom, span.start + 1, span.end) : `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${encodeURIComponent(`${chrom}${name}`)}`,
      note: 'Chromosome, arm, region, band, sub-band',
      parts: [
        ['Chromosome', chrom], ['Arm', arm === 'p' ? 'p, short arm' : 'q, long arm'], ['Region', m[3]], m[4] ? ['Band', m[4]] : null, m[5] ? ['Sub-band', m[5]] : null,
        span ? ['GRCh38', coords(chrom, span.start + 1, span.end)] : null,
      ].filter(Boolean),
    };
  }),
  rule(/^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-(\d+))?$/, (m) => ({
    type: 'UniProtKB accession', id: m[0], db: 'UniProtKB', url: `https://www.uniprot.org/uniprotkb/${m[1]}/entry`, note: 'Six or ten characters, then an optional isoform number',
    parts: m[2] ? [['Isoform', m[2]]] : [],
  })),
  rule(/^([A-Z0-9]{1,5})_HUMAN$/i, (m) => ({
    type: 'UniProtKB entry name', id: `${m[1].toUpperCase()}_HUMAN`, db: 'UniProtKB', url: `https://www.uniprot.org/uniprotkb?query=${m[1].toUpperCase()}_HUMAN`,
    note: 'Mnemonic, then species code', parts: [],
  })),
  rule(/^(?:pdb_0000)?([1-9](?=[A-Z0-9]*[A-Z])[A-Z0-9]{3})$/i, (m) => ({
    type: 'PDB entry', id: m[0].toLowerCase().startsWith('pdb_') ? m[0].toLowerCase() : m[1].toUpperCase(), db: 'Protein Data Bank', url: `https://www.rcsb.org/structure/${m[1].toUpperCase()}`,
    note: 'A digit then three letters or digits', parts: [],
  })),
];

export const identify = (input, ctx = {}) => {
  const value = String(input ?? '').trim();
  if (!value) return null;
  for (const { re, build } of RULES) {
    const m = value.match(re);
    if (m) {
      const hit = build(m, ctx);
      if (hit) return hit;
    }
  }
  return null;
};
