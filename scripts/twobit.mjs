const SIGNATURE = 0x1a412743;
const MAP = [3, 1, 0, 2];
const TABLE = new Uint8Array(1024);
for (let b = 0; b < 256; b += 1) {
  for (let j = 0; j < 4; j += 1) TABLE[b * 4 + j] = MAP[(b >> (6 - 2 * j)) & 3];
}

const decode = (buffer, at, size, nBlocks) => {
  const out = new Uint8Array(size);
  const full = size >> 2;
  for (let i = 0; i < full; i += 1) {
    const t = buffer[at + i] * 4;
    const k = i * 4;
    out[k] = TABLE[t];
    out[k + 1] = TABLE[t + 1];
    out[k + 2] = TABLE[t + 2];
    out[k + 3] = TABLE[t + 3];
  }
  for (let k = full * 4; k < size; k += 1) out[k] = TABLE[buffer[at + (k >> 2)] * 4 + (k & 3)];
  for (const [start, n] of nBlocks) out.fill(4, start, start + n);
  return out;
};

export const readTwoBit = (buffer) => {
  const le = buffer.readUInt32LE(0) === SIGNATURE;
  if (!le && buffer.readUInt32BE(0) !== SIGNATURE) throw new Error('Not a 2bit file');
  const u32 = (o) => (le ? buffer.readUInt32LE(o) : buffer.readUInt32BE(o));
  if (u32(4) !== 0) throw new Error('Only version 0 2bit files are read');
  const offsets = new Map();
  let o = 16;
  for (let k = 0; k < u32(8); k += 1) {
    const len = buffer.readUInt8(o);
    offsets.set(buffer.toString('latin1', o + 1, o + 1 + len), u32(o + 1 + len));
    o += 1 + len + 4;
  }
  const sequence = (name) => {
    let p = offsets.get(name);
    if (p === undefined) return null;
    const size = u32(p);
    const nCount = u32(p + 4);
    p += 8;
    const nBlocks = Array.from({ length: nCount }, (_, i) => [u32(p + 4 * i), u32(p + 4 * (nCount + i))]);
    p += 8 * nCount;
    const mCount = u32(p);
    p += 4;
    let masked = 0;
    for (let i = 0; i < mCount; i += 1) masked += u32(p + 4 * (mCount + i));
    p += 8 * mCount + 4;
    const dna = p;
    return { size, nBlocks, masked, codes: () => decode(buffer, dna, size, nBlocks) };
  };
  return { names: [...offsets.keys()], sequence };
};
