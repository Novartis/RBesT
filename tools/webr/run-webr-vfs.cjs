// Mount a prebuilt Emscripten VFS library image into webR and run an R script
// against it, under Node.
//
// This is the check that `.github/workflows/webr-image.yaml` runs against the
// image it just built: it proves the delivery model end to end -- no package
// installation, no repository, everything from the mounted image.
//
// Usage:
//   WEBR_VFS=<base> node tools/webr/run-webr-vfs.cjs <script.R>
//
// <base> is the image basename; the driver accepts either
//   <base>.data + <base>.js.metadata     (uncompressed) or
//   <base>.data.gz + <base>.js.metadata  (as produced with `compress: true`)
// and decompresses in memory when needed. Emscripten's `gzip` metadata flag
// is a hint for the web server's Content-Encoding; the mount itself always
// wants the plain bytes.
//
// Exits non-zero if the R code signals an error, so CI fails loudly.

const { WebR } = require('webr');
const { readFileSync, existsSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');

function readImageBlob(base) {
  const plain = `${base}.data`;
  if (existsSync(plain)) return readFileSync(plain);
  const gz = `${base}.data.gz`;
  if (existsSync(gz)) return gunzipSync(readFileSync(gz));
  throw new Error(`neither ${plain} nor ${gz} exists`);
}

async function main() {
  const scriptPath = process.argv[2];
  if (!scriptPath) throw new Error('usage: run-webr-vfs.cjs <script.R>');
  const code = readFileSync(scriptPath, 'utf8');
  const base = process.env.WEBR_VFS;

  const t0 = Date.now();
  const webR = new WebR({ interactive: false });
  await webR.init();
  console.error(`[driver] webR ready (${Date.now() - t0} ms)`);

  // Deliberately point the package repo at a dead address: nothing may be
  // downloaded. Everything must come from the mounted image.
  await webR.evalRVoid('options(webr.pkg.urls = "http://127.0.0.1:1/dead")');

  if (base) {
    const t1 = Date.now();
    const metadata = JSON.parse(readFileSync(`${base}.js.metadata`, 'utf8'));
    const blob = readImageBlob(base);
    await webR.FS.mkdir('/vfslib');
    await webR.FS.mount('WORKERFS', { packages: [{ metadata, blob }] }, '/vfslib');
    console.error(`[driver] mounted ${base} ` +
      `(${(blob.length / 1e6).toFixed(1)} MB, ${metadata.files.length} files) ` +
      `in ${Date.now() - t1} ms`);
    await webR.evalRVoid('.libPaths(c("/vfslib", .libPaths()))');
  } else {
    console.error('[driver] WEBR_VFS unset -- running against the stock webR library only');
  }

  const shelter = await new webR.Shelter();
  try {
    // `captureConditions: false` is what keeps the buffered output available
    // when the script errors -- with conditions captured, `captureR()` throws
    // and everything printed up to that point is lost, which is exactly the
    // output needed to diagnose a failure. The cost is that an R error is no
    // longer signalled to the caller, so completion is detected explicitly:
    // the sentinel below is the last expression evaluated, and an error
    // anywhere above it aborts the remaining expressions and leaves it unset.
    const result = await shelter.captureR(`${code}\n.webr_run_completed <- TRUE\n`, {
      withAutoprint: true, captureStreams: true, captureConditions: false,
    });
    for (const o of result.output) {
      process.stdout.write(`${o.type === 'stderr' ? '! ' : ''}${o.data}\n`);
    }
    const completed = await webR.evalRBoolean(
      'isTRUE(mget(".webr_run_completed", envir = globalenv(), ifnotfound = FALSE)[[1]])'
    );
    if (!completed) {
      console.error(`[driver] ${scriptPath} did not run to completion`);
      process.exitCode = 1;
    }
  } catch (e) {
    console.error('[driver] R error:', (e && e.message) || e);
    process.exitCode = 1;
  } finally {
    await shelter.purge();
    await webR.close();
  }
}

main().catch((e) => {
  console.error('[driver] fatal:', (e && e.stack) || e);
  process.exit(1);
});
