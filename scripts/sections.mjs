const stderr = (message) => process.stderr.write(`${message}\n`);

export const runSections = async (sections, { now = new Date(), previous = {}, only = [], log = stderr } = {}) => {
  const stamp = now.toISOString();
  const out = { generated: stamp, asOf: {}, failed: [] };
  const keep = (key) => {
    out[key] = previous[key];
    out.asOf[key] = previous.asOf?.[key] ?? previous.generated;
  };
  for (const [key, build] of Object.entries(sections)) {
    if (only.length && !only.includes(key)) {
      if (key in previous) keep(key);
      continue;
    }
    log(key);
    try {
      out[key] = await build();
      out.asOf[key] = stamp;
    } catch (error) {
      if (!(key in previous)) throw new Error(`${key}: ${error.message}`, { cause: error });
      keep(key);
      out.failed.push({ section: key, error: error.message });
      log(`${key} failed: ${error.message}; keeping values from ${out.asOf[key]}`);
    }
  }
  return out;
};

export const diffCounts = (before, after) => {
  if (!before) return [];
  return Object.keys(after)
    .filter((key) => typeof before[key] === 'number' && before[key] !== after[key])
    .map((key) => ({ key, before: before[key], after: after[key], delta: after[key] - before[key] }))
    .sort((a, b) => Math.abs(b.delta / b.before) - Math.abs(a.delta / a.before));
};
