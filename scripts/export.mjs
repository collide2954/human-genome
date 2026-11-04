import { writeFile } from 'node:fs/promises';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';

const SNAPSHOT = new URL('../src/snapshot.js', import.meta.url);

export const snapshotJson = async () => JSON.stringify((await import(`${SNAPSHOT.href}?t=${Date.now()}`)).default);

if (process.argv[1] === fileURLToPath(import.meta.url)) {
  await writeFile(join(process.argv[2] ?? '.', 'snapshot.json'), await snapshotJson());
}
