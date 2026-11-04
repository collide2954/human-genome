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
