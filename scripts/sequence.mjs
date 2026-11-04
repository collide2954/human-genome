const CODE = { A: 0, C: 1, G: 2, T: 3 };
const LETTER = 'ACGTN';

export const encode = (text) => Uint8Array.from(text.toUpperCase(), (ch) => CODE[ch] ?? 4);

export const decodeText = (codes) => Array.from(codes, (c) => LETTER[c]).join('');

export const composition = (codes, binSize = 1e6) => {
  const counts = [0, 0, 0, 0, 0];
  const dinucs = new Array(16).fill(0);
  const bins = [];
  let prev = 4;
  let gc = 0;
  let acgt = 0;
  let binStart = 0;
  const close = (end) => {
    bins.push(acgt >= (end - binStart) / 2 ? gc / acgt : null);
    gc = 0;
    acgt = 0;
    binStart = end;
  };
  for (let i = 0; i < codes.length; i += 1) {
    if (i - binStart === binSize) close(i);
    const c = codes[i];
    counts[c] += 1;
    if (c < 4) {
      acgt += 1;
      if (c === 1 || c === 2) gc += 1;
      if (prev < 4) dinucs[prev * 4 + c] += 1;
    }
    prev = c;
  }
  if (codes.length > binStart) close(codes.length);
  return {
    A: counts[0], C: counts[1], G: counts[2], T: counts[3], N: counts[4], cpg: dinucs[6], dinucs, bins,
  };
};

const COMP = [3, 2, 1, 0, 4];

export const gtfRecord = (line) => {
  const c = line.split('\t');
  const attrs = {};
  for (const [, key, quoted, plain] of c[8].matchAll(/(\w+) (?:"([^"]*)"|(\S+?));/g)) {
    const value = quoted ?? plain;
    if (key === 'tag') (attrs.tag ??= []).push(value);
    else attrs[key] = value;
  }
  const chrom = c[0].replace(/^chr/, '');
  return {
    chrom: chrom === 'M' ? 'MT' : chrom, source: c[1], feature: c[2], start: Number(c[3]), end: Number(c[4]), strand: c[6], frame: c[7], attrs,
  };
};

const FEATURES = { exon: 'exons', CDS: 'cds', start_codon: 'start', stop_codon: 'stop' };

export const parseManeGtf = (lines) => {
  const txs = new Map();
  for (const line of lines) {
    if (!line || line.startsWith('#')) continue;
    const r = gtfRecord(line);
    const id = r.attrs.transcript_id;
    if (!id) continue;
    if (!txs.has(id)) {
      txs.set(id, {
        id, gene: r.attrs.gene_name, chrom: r.chrom, strand: r.strand, tag: (r.attrs.tag ?? []).find((t) => t.startsWith('MANE')) ?? null, exons: [], cds: [], start: [], stop: [],
      });
    }
    const list = FEATURES[r.feature];
    if (list) txs.get(id)[list].push([r.start, r.end]);
  }
  return [...txs.values()].map((t) => {
    const order = (list) => list.sort((a, b) => (t.strand === '-' ? b[0] - a[0] : a[0] - b[0]));
    return { ...t, exons: order(t.exons), cds: order(t.cds), start: order(t.start), stop: order(t.stop) };
  });
};

export const introns = (tx) => tx.exons.slice(1).map((next, k) => {
  const prev = tx.exons[k];
  return tx.strand === '-' ? [next[1] + 1, prev[0] - 1] : [prev[1] + 1, next[0] - 1];
});

const windowOf = (codes, from, to, strand) => {
  if (from < 1 || to > codes.length || to < from) return null;
  const out = codes.slice(from - 1, to);
  if (strand !== '-') return out;
  out.reverse();
  for (let i = 0; i < out.length; i += 1) out[i] = COMP[out[i]];
  return out;
};

export const spliceSites = (codes, tx) => {
  const donors = [];
  const acceptors = [];
  for (const [gs, ge] of introns(tx)) {
    const donor = tx.strand === '-' ? windowOf(codes, ge - 5, ge + 3, '-') : windowOf(codes, gs - 3, gs + 5, '+');
    const acceptor = tx.strand === '-' ? windowOf(codes, gs - 3, gs + 19, '-') : windowOf(codes, ge - 19, ge + 3, '+');
    if (donor && acceptor) {
      donors.push(donor);
      acceptors.push(acceptor);
    }
  }
  return { donors, acceptors };
};

export const intronClass = (donor, acceptor) => {
  const pair = `${decodeText(donor.slice(3, 5))}-${decodeText(acceptor.slice(18, 20))}`;
  return ['GT-AG', 'GC-AG', 'AT-AC'].includes(pair) ? pair : 'other';
};

const spliced = (codes, segments, strand) => {
  const parts = segments.map(([s, e]) => windowOf(codes, s, e, strand));
  if (parts.some((p) => !p)) return null;
  const out = new Uint8Array(parts.reduce((n, p) => n + p.length, 0));
  let at = 0;
  for (const p of parts) {
    out.set(p, at);
    at += p.length;
  }
  return out;
};

export const kozakWindow = (codes, tx) => {
  if (!tx.start.length) return null;
  const pos = tx.strand === '-' ? tx.start[0][1] : tx.start[0][0];
  let acc = 0;
  let offset = -1;
  for (const [s, e] of tx.exons) {
    if (pos >= s && pos <= e) {
      offset = acc + (tx.strand === '-' ? e - pos : pos - s);
      break;
    }
    acc += e - s + 1;
  }
  const seq = spliced(codes, tx.exons, tx.strand);
  if (!seq || offset < 6 || offset + 4 > seq.length) return null;
  return seq.slice(offset - 6, offset + 4);
};

export const cdsSequence = (codes, tx) => {
  const segments = [...tx.cds, ...tx.stop].sort((a, b) => (tx.strand === '-' ? b[0] - a[0] : a[0] - b[0]));
  return spliced(codes, segments, tx.strand);
};

export const pwm = (windows, length) => {
  const counts = Array.from({ length }, () => [0, 0, 0, 0]);
  for (const w of windows) {
    for (let i = 0; i < length; i += 1) if (w[i] < 4) counts[i][w[i]] += 1;
  }
  return counts;
};

const STOPS = { 48: 'TAA', 50: 'TAG', 56: 'TGA' };

export const codonCounts = (seqs) => {
  const codons = new Array(64).fill(0);
  const stops = { TAA: 0, TAG: 0, TGA: 0 };
  let coding = 0;
  let broken = 0;
  let nonAtg = 0;
  for (const seq of seqs) {
    const n = seq.length / 3;
    const last = n >= 1 ? seq[seq.length - 3] * 16 + seq[seq.length - 2] * 4 + seq[seq.length - 1] : -1;
    if (!Number.isInteger(n) || !STOPS[last] || seq.some((c) => c > 3)) {
      broken += 1;
      continue;
    }
    coding += 1;
    stops[STOPS[last]] += 1;
    if (seq[0] * 16 + seq[1] * 4 + seq[2] !== 14) nonAtg += 1;
    for (let i = 0; i < seq.length - 3; i += 3) codons[seq[i] * 16 + seq[i + 1] * 4 + seq[i + 2]] += 1;
  }
  return {
    codons, stops, coding, broken, nonAtg,
  };
};

export const homopolymers = (codes) => {
  const best = [0, 1, 2, 3].map(() => ({ length: 0, start: -1 }));
  let base = 4;
  let start = 0;
  for (let i = 0; i <= codes.length; i += 1) {
    const c = i < codes.length ? codes[i] : 5;
    if (c === base) continue;
    if (base < 4 && i - start > best[base].length) best[base] = { length: i - start, start };
    base = c;
    start = i;
  }
  return { A: best[0], C: best[1], G: best[2], T: best[3] };
};

export const g4Motifs = (codes) => {
  let gRun = 0;
  let gEnd = -1e9;
  let gChain = 0;
  let gCount = 0;
  let cRun = 0;
  let cEnd = -1e9;
  let cChain = 0;
  let cCount = 0;
  for (let i = 0; i <= codes.length; i += 1) {
    const c = i < codes.length ? codes[i] : 4;
    if (c === 2) gRun += 1;
    else {
      if (gRun >= 3) {
        const loop = i - gRun - gEnd;
        gChain = gChain > 0 && loop >= 1 && loop <= 7 ? gChain + 1 : 1;
        gEnd = i;
        if (gChain === 4) {
          gCount += 1;
          gChain = 0;
        }
      }
      gRun = 0;
    }
    if (c === 1) cRun += 1;
    else {
      if (cRun >= 3) {
        const loop = i - cRun - cEnd;
        cChain = cChain > 0 && loop >= 1 && loop <= 7 ? cChain + 1 : 1;
        cEnd = i;
        if (cChain === 4) {
          cCount += 1;
          cChain = 0;
        }
      }
      cRun = 0;
    }
    if (c === 4) {
      gChain = 0;
      cChain = 0;
    }
  }
  return { plus: gCount, minus: cCount };
};
