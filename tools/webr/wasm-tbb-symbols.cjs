// Report Intel TBB symbols that a WebAssembly module imports but does not
// itself export -- i.e. the ones that must be resolved by something else at
// `dlopen()` time, and that make `library(rstan)` fail in stock webR when
// nothing provides them.
//
// This is the mechanical check behind the TBB fix described in
// design/howto-build-rbest-webr.md section 9. Do not eyeball `wasm-objdump`
// output instead: nearly every `tbb::detail::d1::*` import is a header-defined
// weak symbol that the same module also exports and therefore self-resolves,
// so the raw import list is misleading. Only the difference matters.
//
// Usage:
//   node tools/webr/wasm-tbb-symbols.cjs <image-base>      # VFS library image
//   node tools/webr/wasm-tbb-symbols.cjs <directory>       # tree of .so files
//   node tools/webr/wasm-tbb-symbols.cjs <file.so>         # a single module
//
// <image-base> is the basename `rwasm::make_vfs_library()` produced, i.e. the
// path without the extension; both `<base>.data` and `<base>.data.gz` are
// accepted, alongside `<base>.js.metadata`.
//
// Options:
//   --pattern <re>  symbols to consider           (default: tbb)
//   --quiet         print only the summary line
//
// Exit status is 0 when nothing is unresolved and 1 otherwise, so it can be
// used directly as a CI gate.

const fs = require('node:fs');
const path = require('node:path');
const { gunzipSync } = require('node:zlib');

// --- minimal wasm binary reader ------------------------------------------
// Only the import (id 2) and export (id 7) sections are decoded; every other
// section is skipped by its declared length, so this stays independent of the
// rest of the binary format.

function readVarUint(buf, cursor) {
  let result = 0;
  let shift = 0;
  let byte;
  do {
    if (cursor.p >= buf.length) throw new Error('truncated LEB128');
    byte = buf[cursor.p++];
    result += (byte & 0x7f) * 2 ** shift;
    shift += 7;
  } while (byte & 0x80);
  return result;
}

function readName(buf, cursor) {
  const len = readVarUint(buf, cursor);
  const str = buf.toString('utf8', cursor.p, cursor.p + len);
  cursor.p += len;
  return str;
}

// A limits record: flags, minimum, and a maximum only when the flag is set.
function skipLimits(buf, cursor) {
  const flags = readVarUint(buf, cursor);
  readVarUint(buf, cursor);
  if (flags & 0x01) readVarUint(buf, cursor);
}

function readImportsExports(buf) {
  if (buf.length < 8 || buf.readUInt32LE(0) !== 0x6d736100) {
    throw new Error('not a WebAssembly module');
  }
  const cursor = { p: 8 };
  const imports = new Set();
  const exports = new Set();

  while (cursor.p < buf.length) {
    const id = buf[cursor.p++];
    const size = readVarUint(buf, cursor);
    const end = cursor.p + size;

    if (id === 2) {
      const count = readVarUint(buf, cursor);
      for (let i = 0; i < count; i++) {
        readName(buf, cursor); // module
        const field = readName(buf, cursor);
        const kind = buf[cursor.p++];
        if (kind === 0) {
          readVarUint(buf, cursor); // type index
        } else if (kind === 1) {
          cursor.p++; // reftype
          skipLimits(buf, cursor);
        } else if (kind === 2) {
          skipLimits(buf, cursor);
        } else if (kind === 3) {
          cursor.p += 2; // valtype, mutability
        } else if (kind === 4) {
          cursor.p++; // attribute
          readVarUint(buf, cursor); // type index
        } else {
          throw new Error(`unknown import kind ${kind}`);
        }
        imports.add(field);
      }
    } else if (id === 7) {
      const count = readVarUint(buf, cursor);
      for (let i = 0; i < count; i++) {
        const field = readName(buf, cursor);
        cursor.p++; // kind
        readVarUint(buf, cursor); // index
        exports.add(field);
      }
    }

    if (end > buf.length) throw new Error('section runs past end of module');
    cursor.p = end;
  }
  return { imports, exports };
}

// --- module sources -------------------------------------------------------

function modulesFromImage(base) {
  const metadata = JSON.parse(fs.readFileSync(`${base}.js.metadata`, 'utf8'));
  const plain = `${base}.data`;
  const blob = fs.existsSync(plain)
    ? fs.readFileSync(plain)
    : gunzipSync(fs.readFileSync(`${base}.data.gz`));
  return metadata.files
    .filter((f) => f.filename.endsWith('.so'))
    .map((f) => ({ name: f.filename, buf: blob.subarray(f.start, f.end) }));
}

function modulesFromDir(dir) {
  const out = [];
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) out.push(...modulesFromDir(full));
    else if (entry.name.endsWith('.so')) out.push({ name: full, buf: fs.readFileSync(full) });
  }
  return out;
}

function collectModules(target) {
  if (fs.existsSync(`${target}.js.metadata`)) return modulesFromImage(target);
  if (fs.existsSync(target) && fs.statSync(target).isDirectory()) return modulesFromDir(target);
  if (fs.existsSync(target)) return [{ name: target, buf: fs.readFileSync(target) }];
  throw new Error(
    `no such target: ${target} (expected a VFS image base, a directory or a .so file)`
  );
}

// --- main -----------------------------------------------------------------

function main(argv) {
  let pattern = 'tbb';
  let quiet = false;
  const positional = [];
  for (let i = 0; i < argv.length; i++) {
    if (argv[i] === '--pattern') pattern = argv[++i];
    else if (argv[i] === '--quiet') quiet = true;
    else positional.push(argv[i]);
  }
  if (positional.length !== 1) {
    throw new Error('usage: wasm-tbb-symbols.cjs [--pattern <re>] [--quiet] <image-base|dir|so>');
  }
  const re = new RegExp(pattern);

  const modules = collectModules(positional[0]);
  if (!modules.length) throw new Error(`no .so modules found under ${positional[0]}`);

  let unresolved = 0;
  let skipped = 0;
  for (const mod of modules) {
    let sections;
    try {
      sections = readImportsExports(mod.buf);
    } catch (e) {
      // Not fatal on its own, but never silent: a module that cannot be parsed
      // is a module that was not checked.
      console.error(`[tbb] SKIPPED ${mod.name}: ${e.message}`);
      skipped++;
      continue;
    }
    const missing = [...sections.imports]
      .filter((name) => re.test(name) && !sections.exports.has(name))
      .sort();
    if (missing.length && !quiet) {
      console.log(`[tbb] ${mod.name}`);
      for (const name of missing) console.log(`[tbb]     ${name}`);
    }
    unresolved += missing.length;
  }

  console.log(
    `[tbb] ${modules.length - skipped} module(s) checked, ` +
      `${unresolved} unresolved /${pattern}/ import(s)` +
      (skipped ? `, ${skipped} skipped` : '')
  );
  return unresolved === 0 && skipped === 0 ? 0 : 1;
}

try {
  process.exitCode = main(process.argv.slice(2));
} catch (e) {
  console.error(`[tbb] fatal: ${(e && e.message) || e}`);
  process.exitCode = 2;
}
