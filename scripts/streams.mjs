import { createHash } from 'node:crypto';
import { createReadStream, createWriteStream } from 'node:fs';
import { open, rename } from 'node:fs/promises';
import { Readable } from 'node:stream';
import { pipeline } from 'node:stream/promises';
import { StringDecoder } from 'node:string_decoder';
import { createGunzip, inflateRawSync } from 'node:zlib';

const trimCr = (s) => (s.endsWith('\r') ? s.slice(0, -1) : s);

export async function* lines(chunks) {
  const decoder = new StringDecoder('utf8');
  let rest = '';
  for await (const chunk of chunks) {
    const parts = (rest + (typeof chunk === 'string' ? chunk : decoder.write(chunk))).split('\n');
    rest = parts.pop();
    for (const part of parts) yield trimCr(part);
  }
  rest += decoder.end();
  if (rest) yield trimCr(rest);
}

export function* bufferLines(buffer) {
  let start = 0;
  while (start < buffer.length) {
    let end = buffer.indexOf(10, start);
    if (end < 0) end = buffer.length;
    yield buffer.toString('utf8', start, end > start && buffer[end - 1] === 13 ? end - 1 : end);
    start = end + 1;
  }
}

export const unzip = (buffer) => {
  let end = buffer.length - 22;
  while (end >= 0 && buffer.readUInt32LE(end) !== 0x06054b50) end -= 1;
  const count = buffer.readUInt16LE(end + 10);
  let at = buffer.readUInt32LE(end + 16);
  const files = {};
  for (let k = 0; k < count; k += 1) {
    const method = buffer.readUInt16LE(at + 10);
    const size = buffer.readUInt32LE(at + 20);
    const nameLen = buffer.readUInt16LE(at + 28);
    const extraLen = buffer.readUInt16LE(at + 30);
    const commentLen = buffer.readUInt16LE(at + 32);
    const offset = buffer.readUInt32LE(at + 42);
    const name = buffer.toString('utf8', at + 46, at + 46 + nameLen);
    const start = offset + 30 + buffer.readUInt16LE(offset + 26) + buffer.readUInt16LE(offset + 28);
    const data = buffer.subarray(start, start + size);
    files[name] = method === 8 ? inflateRawSync(data) : data;
    at += 46 + nameLen + extraLen + commentLen;
  }
  return files;
};

export const gzLines = (path) => lines(createReadStream(path).pipe(createGunzip()));

export const fileLines = (path) => lines(createReadStream(path));

export const anyLines = async (path) => {
  const handle = await open(path, 'r');
  try {
    const { buffer, bytesRead } = await handle.read(Buffer.alloc(2), 0, 2, 0);
    return bytesRead === 2 && buffer[0] === 0x1f && buffer[1] === 0x8b ? gzLines(path) : fileLines(path);
  } finally {
    await handle.close();
  }
};

export const cacheName = (url) => {
  const u = new URL(url);
  const name = `${u.hostname}${u.pathname}${u.search}`.replace(/[^A-Za-z0-9.-]+/g, '_').replace(/_+$/, '');
  return name.length <= 160 ? name : `${name.slice(0, 140)}_${createHash('sha256').update(url).digest('hex').slice(0, 16)}`;
};

export const downloadTo = async (url, file, headers = {}) => {
  const res = await fetch(url, { headers });
  if (!res.ok) throw new Error(`HTTP ${res.status} ${url}`);
  const part = `${file}.part`;
  await pipeline(Readable.fromWeb(res.body), createWriteStream(part));
  await rename(part, file);
  return file;
};
