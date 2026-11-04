import { test } from 'node:test';
import assert from 'node:assert/strict';
import { gzipSync, deflateRawSync, crc32 } from 'node:zlib';
import { mkdtemp, writeFile, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import {
  lines, bufferLines, unzip, gzLines, fileLines, cacheName,
} from '../../scripts/streams.mjs';
import { retry } from '../../src/lib.js';

const collect = async (iterable) => {
  const out = [];
  for await (const x of iterable) out.push(x);
  return out;
};

test('lines joins chunks split mid-line and drops carriage returns', async () => {
  const chunks = ['ab\ncd', '\r\nef', '\ngh'].map((s) => Buffer.from(s));
  assert.deepEqual(await collect(lines(chunks)), ['ab', 'cd', 'ef', 'gh']);
  assert.deepEqual(await collect(lines([Buffer.from('x\n')])), ['x']);
});

test('lines keeps multibyte characters split across chunks', async () => {
  const bytes = Buffer.from('5′ UTR\nok');
  const chunks = [bytes.subarray(0, 2), bytes.subarray(2)];
  assert.deepEqual(await collect(lines(chunks)), ['5′ UTR', 'ok']);
});

test('bufferLines walks a buffer line by line', () => {
  assert.deepEqual([...bufferLines(Buffer.from('a\r\nb\n\nc'))], ['a', 'b', '', 'c']);
});

const zipOf = (files) => {
  const locals = [];
  const centrals = [];
  let offset = 0;
  for (const [name, text] of Object.entries(files)) {
    const raw = Buffer.from(text);
    const packed = deflateRawSync(raw);
    const nameBytes = Buffer.from(name);
    const local = Buffer.alloc(30);
    local.writeUInt32LE(0x04034b50, 0);
    local.writeUInt16LE(8, 8);
    local.writeUInt32LE(crc32(raw), 14);
    local.writeUInt32LE(packed.length, 18);
    local.writeUInt32LE(raw.length, 22);
    local.writeUInt16LE(nameBytes.length, 26);
    const central = Buffer.alloc(46);
    central.writeUInt32LE(0x02014b50, 0);
    central.writeUInt16LE(8, 10);
    central.writeUInt32LE(crc32(raw), 16);
    central.writeUInt32LE(packed.length, 20);
    central.writeUInt32LE(raw.length, 24);
    central.writeUInt16LE(nameBytes.length, 28);
    central.writeUInt32LE(offset, 42);
    locals.push(local, nameBytes, packed);
    centrals.push(central, nameBytes);
    offset += 30 + nameBytes.length + packed.length;
  }
  const dir = Buffer.concat(centrals);
  const end = Buffer.alloc(22);
  end.writeUInt32LE(0x06054b50, 0);
  end.writeUInt16LE(Object.keys(files).length, 8);
  end.writeUInt16LE(Object.keys(files).length, 10);
  end.writeUInt32LE(dir.length, 12);
  end.writeUInt32LE(offset, 16);
  return Buffer.concat([...locals, dir, end]);
};

test('unzip inflates every entry by name', () => {
  const files = unzip(zipOf({ 'a.tsv': 'x\ty\n1\t2\n', 'b.txt': 'hello' }));
  assert.deepEqual(Object.keys(files), ['a.tsv', 'b.txt']);
  assert.equal(files['a.tsv'].toString(), 'x\ty\n1\t2\n');
  assert.equal(files['b.txt'].toString(), 'hello');
});

test('gzLines and fileLines stream a file from disk', async () => {
  const dir = await mkdtemp(join(tmpdir(), 'hg-'));
  try {
    await writeFile(join(dir, 'a.gz'), gzipSync('one\ntwo\nthree'));
    await writeFile(join(dir, 'b.txt'), 'four\r\nfive\n');
    assert.deepEqual(await collect(gzLines(join(dir, 'a.gz'))), ['one', 'two', 'three']);
    assert.deepEqual(await collect(fileLines(join(dir, 'b.txt'))), ['four', 'five']);
  } finally {
    await rm(dir, { recursive: true });
  }
});

test('cacheName turns a URL into a safe file name', () => {
  assert.equal(cacheName('https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz'), 'ftp.ncbi.nlm.nih.gov_gene_DATA_gene2pubmed.gz');
  assert.equal(cacheName('https://rest.uniprot.org/uniprotkb/stream?query=a:b&format=tsv'), 'rest.uniprot.org_uniprotkb_stream_query_a_b_format_tsv');
});

test('cacheName shortens long URLs and keeps them apart', () => {
  const long = (field) => `https://rest.uniprot.org/uniprotkb/stream?fields=${Array.from({ length: 40 }, (_, i) => `${field}_${i}`).join(',')}`;
  const a = cacheName(long('ft_region'));
  const b = cacheName(long('ft_domain'));
  assert.ok(a.length <= 160, `${a.length}`);
  assert.ok(b.length <= 160, `${b.length}`);
  assert.notEqual(a, b);
  assert.match(a, /^rest\.uniprot\.org_uniprotkb_stream_fields_ft_region/);
  assert.equal(cacheName(long('ft_region')), a);
});

test('retry returns the first success and rethrows after the last try', async () => {
  let calls = 0;
  const ok = await retry(async () => {
    calls += 1;
    if (calls < 3) throw new Error('flaky');
    return 'done';
  }, { tries: 3, wait: 1 });
  assert.equal(ok, 'done');
  assert.equal(calls, 3);
  await assert.rejects(retry(async () => { throw new Error('down'); }, { tries: 2, wait: 1 }), /down/);
});
