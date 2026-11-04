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
