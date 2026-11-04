export const tsvRecords = (text) => {
  const lines = text.split('\n').map((line) => line.replace(/\r$/, '')).filter(Boolean);
  if (!lines.length) return [];
  const head = lines[0].split('\t').map((h) => h.trim());
  return lines.slice(1).map((line) => {
    const cells = line.split('\t');
    return Object.fromEntries(head.map((h, i) => [h, (cells[i] ?? '').trim()]));
  });
};

export const bump = (obj, key, by = 1) => {
  obj[key] = (obj[key] ?? 0) + by;
};

export const round = (v, places) => Number(v.toFixed(places));
