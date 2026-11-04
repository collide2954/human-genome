export const quantile = (sorted, q) => {
  if (!sorted.length) return NaN;
  const pos = (sorted.length - 1) * q;
  const lo = Math.floor(pos);
  const hi = Math.ceil(pos);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo);
};

export const summary = (values) => {
  if (!values.length) return null;
  const s = [...values].sort((a, b) => a - b);
  return {
    n: s.length,
    min: s[0],
    q1: quantile(s, 0.25),
    median: quantile(s, 0.5),
    q3: quantile(s, 0.75),
    max: s.at(-1),
    mean: s.reduce((a, b) => a + b, 0) / s.length,
  };
};

const tidy = (v) => Number(v.toPrecision(12));

export const logEdges = (lo, hi, perDecade) => Array.from({ length: (hi - lo) * perDecade + 1 }, (_, i) => tidy(10 ** (lo + i / perDecade)));

export const linearEdges = (lo, hi, step) => Array.from({ length: Math.round((hi - lo) / step) + 1 }, (_, i) => tidy(lo + i * step));

export const binIndex = (edges, v) => {
  const last = edges.length - 2;
  if (v < edges[1]) return 0;
  if (v >= edges[last]) return last;
  let lo = 1;
  let hi = last;
  while (lo < hi) {
    const mid = (lo + hi + 1) >> 1;
    if (edges[mid] <= v) lo = mid;
    else hi = mid - 1;
  }
  return lo;
};

export const histogram = (values, edges) => {
  const counts = new Array(edges.length - 1).fill(0);
  for (const v of values) counts[binIndex(edges, v)] += 1;
  return counts;
};

export const lorenz = (values, points = 100) => {
  const s = [...values].sort((a, b) => a - b);
  const total = s.reduce((a, b) => a + b, 0);
  const prefix = [0];
  for (const v of s) prefix.push(prefix.at(-1) + v);
  return Array.from({ length: points + 1 }, (_, k) => {
    const x = k / points;
    return [x, total ? tidy(prefix[Math.round(x * s.length)] / total) : 0];
  });
};

export const gini = (values) => {
  const s = [...values].sort((a, b) => a - b);
  const total = s.reduce((a, b) => a + b, 0);
  if (!total) return 0;
  const weighted = s.reduce((acc, v, i) => acc + (i + 1) * v, 0);
  return tidy((2 * weighted) / (s.length * total) - (s.length + 1) / s.length);
};

export const mergeIntervals = (intervals) => {
  const sorted = [...intervals].sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const out = [];
  for (const [start, end] of sorted) {
    const last = out.at(-1);
    if (last && start <= last[1]) last[1] = Math.max(last[1], end);
    else out.push([start, end]);
  }
  return out;
};

export const coverage = (intervals) => mergeIntervals(intervals).reduce((n, [s, e]) => n + e - s, 0);

const COMPLEMENT = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' };

export const reverseComplement = (seq) => [...seq].reverse().map((c) => COMPLEMENT[c] ?? c).join('');

const rotations = (s) => Array.from({ length: s.length }, (_, i) => s.slice(i) + s.slice(0, i));

export const canonicalMotif = (motif) => {
  const m = motif.toUpperCase();
  return [...rotations(m), ...rotations(reverseComplement(m))].sort()[0];
};

export const topBy = (rows, key, n) => [...rows].sort((a, b) => b[key] - a[key]).slice(0, n);
