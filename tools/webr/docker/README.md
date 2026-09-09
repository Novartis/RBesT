# Native-architecture webR build for RBesT (option C)

Cross-compile environment for RBesT's webR/WebAssembly binary, option C from
`../../../design/howto-build-rbest-webr.md` §4.5 (paths there are repo-root
relative): it runs on the host's **native** architecture, so neither Rosetta
(EOL after macOS 27) nor qemu (unstable for this workload, §4.4) is involved.

**Result: it works.** RBesT compiles to WebAssembly on arm64 and a real
`gMAP()` fit runs in webR under node with posteriors matching a native run.

This directory is tracked in git and is **the** webR build: the top-level
`Makefile` target `r-binary-webr` (manual only: `make r-binary-webr`) runs it
locally, extracting `build/RBesT_<version>.tar.gz` (the `r-source-release`
artifact) for the build's `/src`, and `.github/workflows/webr-image.yaml` runs
the same setup on a GitHub runner with the checkout as `/src`. It started as a
prototype at `design/wasm-build-native/`, which still holds `e2e/`, a manual
browser-verification harness left behind (still git-ignored, disposable);
`design/wasm-build/` (the amd64/qemu harness for option B) is untouched.

## Status

| step | status |
| --- | --- |
| arch-agnostic image builds natively | ✅ verified on arm64 (aarch64, colima) |
| toolchain self-checks (native LLVM, emcc 5.0.7, hello-world link) | ✅ |
| prebuilt wasm dependency closure | ✅ 64 prebuilt, only RBesT and rstan compiled |
| RBesT `.so` cross-compiled | ✅ `RBesT_1.12-0.tgz`, 1.36 MB (`notbb` variant) |
| VFS library image | ✅ 213.9 MB raw / 46.6 MB gz, 17 251 files, mounts in 1.0 s — **projected 62.3 MB raw / ~28.8 MB gz** once `include,lib` stripping is rebuilt (see below) |
| unresolved TBB symbols | ✅ 33 modules scanned, **0** unresolved `tbb::` imports |
| `gMAP()` in webR under node | ✅ 5.4 s, `Rhat` = 1, agrees with native |
| amd64 | designed for, **not** verified locally (no amd64 host); the GitHub Actions workflow is what exercises it |

**The CRAN-only Stan source is verified end to end** (2026-09-08, arm64): rstan
2.32.7 (Stan 2.32.2) and StanHeaders 2.32.10, loaded from the mounted image in
**stock** webR 0.6.0 (npm install diffed against the published tarball, byte
for byte) with the package repository pointed at a dead address, so nothing
could be silently downloaded. `gMAP()` on `AS` gave
`theta_resp_pred` mean 0.250, sd 0.074 against 0.257, sd 0.077 for the same
model run natively — agreement within MCMC noise — and `automixfit()` on the
predictive draws converged.

Two things are worth recording because they were the open risks of the switch:

- **rstan 2.32.7 cross-compiles to wasm cleanly.** It had never been done
  publicly. No change to `patch-rstan-tarball.R` was needed.
- **The existing two TBB no-op stubs are sufficient for it.** The symbol set
  rstan 2.32 imports is the same one rstan 2.39 imported; the mechanical scan
  reports zero unresolved `tbb::` imports across all 33 modules.

The image is also smaller than the r-universe/rstan 2.39 one it replaces:
46.6 MB gzipped against 54.7 MB, i.e. about 15% less to download.

**Header stripping.** `include` and `lib` are now in the default `RBEST_STRIP`
set. They hold the C++ headers and static archives that the `LinkingTo`
packages exist for, and nothing in a browser compiles: measured on the 204 MB
image above they are 141.7 MB of it — 122 MB Boost headers from `BH`, 8.4 MB
`RcppEigen`, 7.2 MB `StanHeaders` — taking the image to a projected 62.3 MB
unpacked and 28.8 MB gzipped. That is 69% off what the browser holds in memory
and 37% off the download. `StanHeaders/stanc.js` must *not* be removed:
`rstan:::.onLoad()` sources it with `mustWork = TRUE`. It is a file, and
`strip` removes only directories, so it survives. See the how-to, §7.

The rows above are all from branch `issue-webr-cran-stan`, whose versions were
still pins at the time; the derivation now in `install-host-deps.R` resolves the
same pair. Supporting checks made on the host side before the build (R 4.6.1,
aarch64, with `rstan` 2.32.7 / `StanHeaders` 2.32.10 / `RcppParallel` 6.2.1
installed):

- **Dependency resolution.** The closure walk against `repo.r-wasm.org`
  resolves 65 packages: 64 prebuilt, plus `rstan`, built here. The wasm
  `StanHeaders` is 2.32.10, equal to the host version, so the equality check
  passes; a deliberate mismatch was confirmed to fail the build.
- **Both tarballs download and validate**, `StanHeaders` via the
  `src/contrib/Archive/` fallback, `rstan` from `src/contrib/`.
- **The oneTBB path is the one taken.** `StanHeaders:::CxxFlags()` emits
  `-DTBB_INTERFACE_NEW`, so `init_threadpool_tbb.hpp` uses
  `tbb/global_control.h` + `tbb/task_arena.h`, not the removed
  `tbb/task_scheduler_init.h`.
- **`NEW_RSTAN` is *not* defined** with rstan 2.32.7 (< 2.36), so
  `stan::rng_t` is `boost::ecuyer1988` — matching rstan 2.32.7's own
  `src/Module.cpp`. This is the compatibility argument, confirmed at runtime
  rather than assumed.
- **Failure class (d) does not arise.** `rstantools::rstan_config()` run with
  the resolved rstan reproduces the tracked `src/stanExports_gMAP.h` byte for
  byte; the `.cc` differs only in spelling the same type as
  `boost::ecuyer1988` rather than `boost::random::ecuyer1988` (Boost declares
  `using random::ecuyer1988;`). No `rng_t`/`mixmax` mismatch anywhere.
- **The native TBB linkage hazard is still live and still handled.**
  `StanHeaders:::LdFlags()` returns `-L…/RcppParallel/lib -ltbb -ltbbmalloc`,
  which is exactly what the `src/Makevars.webr` overrides strip.

How the end-to-end verification was run, for reproduction:

```sh
make r-binary-webr

# mechanical: no tbb:: import may be left unresolved in any module
node tools/webr/wasm-tbb-symbols.cjs build/RBesT-webr-library_1.12-0

# stock webR, installed outside the checkout and diffed against the
# published tarball so nothing in this repo can have patched it
webr_dir="$(mktemp -d)"; cd "$webr_dir"
npm install --no-package-lock --no-audit --no-fund webr@0.6.0
npm pack webr@0.6.0 >/dev/null && tar xzf webr-0.6.0.tgz
diff -r package/dist node_modules/webr/dist
cd -

NODE_PATH="$webr_dir/node_modules" \
WEBR_VFS="$PWD/build/RBesT-webr-library_1.12-0" \
  node tools/webr/run-webr-vfs.cjs tools/webr/verify-rbest.R
```

The npm `webr` version must match `WEBR_TAG` (the container the image was
built against), or the image is loaded by a different R.

For a readable fit rather than the CI assertions, use `make webr-demo` (see
[End-to-end check](#end-to-end-check)), which does the npm setup itself.

## How it is put together

Three stages, only one of which ever executes anything:

1. `FROM --platform=linux/amd64 ghcr.io/r-wasm/webr:v0.6.0 AS webr` — used
   **only** as a `COPY --from` source (`wasm/`, `R/R-VERSION`, `R/build`,
   `host/bin`, `tools/`, `libs/`, `packages/`, and the emscripten port cache).
   Nothing is `RUN` there, so the amd64 pin costs no emulation.
2. `FROM emscripten/emsdk:5.0.7 AS emsdk` — **multi-arch**, so it resolves to
   the build platform automatically; this is where the native toolchain comes
   from. Relocated to `/opt/emsdk` because webR's wasm `Makeconf` bakes that
   absolute path in.
3. `FROM rocker/r-ver:4.6.0` — multi-arch, provides native host R 4.6.0, the
   native `LinkingTo` closure, `rwasm`, and the corporate CA/proxy plumbing.

Dependencies are **never cross-compiled**: their wasm binaries are downloaded
prebuilt from `repo.r-wasm.org`. That is what removes the amd64-only `flang`
requirement — RBesT and `rstan` are the only things compiled.

## Where the Stan packages come from

**CRAN only.** The `stan-dev` r-universe is not consulted by any route: the
GitHub Actions workflow (`.github/workflows/webr-image.yaml`) runs *this same
Docker setup* on a GitHub runner, so there is one build to reason about, not
two.

| | host (native) | wasm |
| --- | --- | --- |
| `StanHeaders` | CRAN `src/contrib/`, `src/contrib/Archive/<pkg>/` fallback | prebuilt from `repo.r-wasm.org` |
| `rstan` | CRAN `src/contrib/` | built here from the *same* tarball, patched |

Neither version is pinned. `install-host-deps.R` **derives** both:

- `StanHeaders` is whatever `repo.r-wasm.org` advertises for the target R minor
  in its emscripten `PACKAGES` index, because the wasm binary is taken prebuilt
  and the host headers must be the same release. CRAN keeps only the current
  release under `src/contrib/`, so a superseded one is fetched from
  `src/contrib/Archive/<pkg>/`.
- `rstan` is CRAN's current release. Nothing constrains it — no wasm `rstan` is
  published anywhere, so it is cross-compiled from the CRAN tarball regardless.

As of 2026-09-09 that resolves to **StanHeaders 2.32.10** (Archive) and
**rstan 2.32.7** (current), i.e. exactly the pair this route was verified
against. Deriving rather than pinning means no edit is needed when
`repo.r-wasm.org` catches up with CRAN.

Why this works:

- **It needs no extra cross-compile.** `repo.r-wasm.org` publishes a wasm
  `StanHeaders` binary, so the host source and the wasm binary are the same
  release. The full closure was walked against that index: 65 packages, and
  `rstan` is the only one with no wasm binary — and no loadable wasm `rstan` is
  published anywhere, so it is compiled here either way.
- **CRAN tarballs are immutable.** This removes the rolling-development-version
  hazard (`howto-build-rbest-webr.md` §7): the same
  `StanHeaders_2.39.0.9000.tar.gz` URL served different bytes hours apart.
  A CRAN version identifies fixed content; an r-universe URL 404s within days.
- **The `boost::ecuyer1988` incompatibility does not arise.** rstan 2.32.7 and
  StanHeaders 2.32.10 are a matched release pair. (For the record, it would not
  arise with CRAN's *current* StanHeaders 2.39.1 either: it ships
  `inst/include/overrides/create_rng.hpp` and guards
  `src/stan/services/util/create_rng.hpp` with `#ifndef NEW_RSTAN` →
  `using rng_t = boost::ecuyer1988`, and `StanHeaders::CxxFlags()` emits
  `-DNEW_RSTAN` only when rstan >= 2.36 is installed. The incompatibility is a
  property of the r-universe *development* headers, not of CRAN.)
- **`src/stanExports_gMAP.*` needs no regeneration.** The tracked sources were
  generated by stanc 2.32.2, which is rstan 2.32.7's stanc, and
  `inst/stan/gMAP.stan` documents that RBesT targets Stan 2.32.
  `rstantools::rstan_config()` is still run during staging — it also generates
  `src/Makevars` — but it is currently a no-op for the exported C++.
- **The `std::from_chars` patch is inert.** StanHeaders 2.32.10 uses
  `boost::lexical_cast`, a form `patch-stanheaders.R` already recognises and
  passes through unchanged.

Two versions derived independently can in principle drift apart, so both ends
are checked rather than assumed:

- `install-host-deps.R` asserts that `StanHeaders:::CxxFlags()` defines
  `NEW_RSTAN` **iff** the installed `rstan` is >= 2.36 — i.e. that the pair
  agrees on whether Stan's RNG type is `stan::rng_t` or `boost::ecuyer1988`.
  Without this the disagreement surfaces only much later, as a compile error
  deep in Stan Math.
- `fetch-wasm-deps.R` requires the wasm `StanHeaders` version to **equal** the
  host one. That holds by construction, since the host version is derived from
  the same index; the check catches the index having moved between the image
  build and the run.

## Running it

Via the Makefile (recommended; builds from the `r-source-release` tarball):

```sh
make r-binary-webr
```

The image (`rbest-wasm-native:<WEBR_TAG>`) is tagged and stays in the local
Docker image store/layer cache across runs and across `make clean` (which
only clears `build/`, not Docker's own state). The Makefile also tracks every
file `COPY`ed into the image as a prerequisite of the `build/webr/.image-built`
stamp, so `docker compose build` is only re-invoked when one of those files
(or `WEBR_TAG`) actually changes -- a re-run with nothing changed skips
straight to the compile step. Force a rebuild with
`rm build/webr/.image-built` or `docker compose build --no-cache` directly.

Directly, for iterating on the Docker setup itself:

```sh
cd tools/webr/docker
export RBEST_SRC="$(git rev-parse --show-toplevel)"
docker compose build
docker compose run --rm build                      # baseline variant
RBEST_DEPS=false RBEST_VARIANTS=notbb \
  docker compose run --rm build                    # the variant that works
```

Behind a corporate intercepting proxy, export `http_proxy` / `https_proxy`
first; if the network also intercepts TLS (the image build's CRAN-index
preflight fails with "Cannot reach the CRAN index" even though the proxy was
detected -- that error is TLS trust, not connectivity), point `CORP_CERTS_DIR`
at a directory of `*.pem` CA roots (see
`../../../design/wasm-build/README.md`, "Optional corporate CA interception"
-- same mechanism here). This repo ships only the empty `certs.d/`
placeholder; real corporate CA material is never committed here and must be
supplied externally. Via the Makefile:

```sh
make CORP_CERTS_DIR=/path/to/private/certs r-binary-webr
```

`CORP_CERTS_DIR` may be given as a relative path (resolved against the repo
root, where `make` is invoked from) or absolute; both work through the
Makefile. Left unset, the build has no knowledge of any specific
organization's PKI.

Outputs land in `output/` (override with `RBEST_OUT`):

- `repo/bin/emscripten/contrib/4.6/*.tgz` — the wasm package repo
- `vfs/rbest-library.data.gz` + `.js.metadata` — the VFS library image
- `build-<variant>.log`, `build-summary.csv`, `wasm-deps.csv`

### On GitHub Actions

`.github/workflows/webr-image.yaml` runs exactly this setup on an
`ubuntu-24.04` runner, on manual dispatch and on any pushed tag. It differs
from a local run in four ways only:

- `/src` is the git checkout rather than an extracted `r-source-release`
  tarball. `stage_source()` copies the subset `R CMD INSTALL` reads and strips
  native build artefacts, so no packaging step is needed first.
- Because a checkout is not a prepared working tree, the workflow runs
  `tools/make-ds.R` first: `data/*.rda` and `R/sysdata.rda` are generated, not
  tracked, and `.onAttach()` reads `pkg_sha` out of the latter, so without it
  the built package installs and then fails on `library(RBesT)`. Every
  `R CMD build` in the `Makefile` does the same thing to its staged tree.
- The image is built from scratch every time — no registry, no layer cache.
  The workflow runs a handful of times a year, so a persistent cache would cost
  storage continuously to save wall clock rarely, and a from-scratch build is
  the reproducible one. It also means CI always resolves the *current* Stan
  versions, where a local rerun keeps whatever the cached image was built with.
- The webR image tag, the emsdk version, the host R version and the strip set
  are workflow inputs, defaulting to the values in `compose.yaml`.

Afterwards the workflow verifies the image (same `verify-rbest.R` fit, in stock
webR diffed against the published npm tarball) and, for a tag, attaches the
image and the packed wasm repository to that tag's release. `wasm-deps.csv`,
`build-summary.csv` and `build-rstan.log` are uploaded as an artifact and
summarised in the job log, including on failure.

### End-to-end check

Two harnesses run against the image, both driven by
`tools/webr/run-webr-vfs.cjs`, which mounts the image into a stock webR under
node and points the package repository at a dead address so nothing can be
downloaded to cover a gap:

- **`tools/webr/verify-rbest.R`** — the CI assertion script. Checks provenance
  (library path, patch markers), that `detectCores()` is still `NA`, and that a
  small `gMAP()` fit produces finite summaries. Silent when it passes.
- **`tools/webr/demo/as-fit.R`** — a readable demonstration: the full MAP
  workflow on the `AS` data set (`gMAP` → `automixfit` → `ess` → `robustify`
  → `postmix`), printing everything as it goes.

The demo has a runner that fetches the matching `webr` npm package into
`build/webr/node` on first use and finds the newest image under `build/`:

```sh
make r-binary-webr        # once, to produce the image
make webr-demo            # fit the AS data set inside webR
make webr-demo-compare    # ... and under native R, then diff the numbers
```

`webr-demo-compare` runs the *same* script in both engines and compares the
posterior summaries it prints (`tools/webr/demo/compare-runs.R`). The seed is
fixed, but wasm32 and the native platform do not produce identical draws, so
the check is agreement to Monte Carlo error, not equality. Observed on arm64:

```
      quantity       run1       run2 rel.diff tolerance
     pred_mean  0.2551540  0.2567600  0.00625      0.05
       pred_sd  0.0754033  0.0761423  0.00971      0.05
     pred_q2.5  0.1244050  0.1221830  0.01790      0.05
    pred_q97.5  0.4303720  0.4304240  0.00012      0.05
      tau_mean  0.3296310  0.3311250  0.00451      0.05
      map_mean  0.2552500  0.2568250  0.00613      0.05
        map_sd  0.0752713  0.0762187  0.01240      0.05
      ess_elir 43.7074000 40.7718000  0.06720      0.20
     post_mean  0.2657290  0.2664470  0.00269      0.05
       post_sd  0.0453977  0.0464005  0.02160      0.05
 post_p_lt_0.3  0.7945610  0.7841490  0.01310      0.05
```

`ess_elir` gets a looser tolerance because it inherits the variability of the
EM mixture approximation on top of the Monte Carlo error; everything else
agrees to about 1%.

Cost of the demo fit (4 chains, 6000 iterations, `cores = 1`), arm64: **29 s**
sampling in webR against **2.4 s** natively, 46 s wall clock for `make
webr-demo` including webR startup and mounting the image. No modification of
`node_modules/webr` is needed -- the TBB no-op stubs live in the RBesT and
rstan wasm binaries, not in webR itself.

## The four blockers, in the order they appeared

Each was classified before being fixed, per the triage gate in the plan.

1. **webR's host-R launchers are amd64 binaries.**
   `$WEBR_ROOT/wasm/R-4.6.0/lib/R/bin/Rscript` is a native amd64 ELF that
   `execv`s `/opt/webr/host/R-4.6.0/lib/R`, and `bin/R` is a bash wrapper with
   the same hardcoded path. `webr-vars.mk` sets `R_HOME` to that wasm R home
   and RBesT's `src/Makevars` shells out to `$(R_HOME)/bin/Rscript`, so this is
   squarely on the critical path. `webr-host-shims.sh` replaces both with `sh`
   wrappers onto the native host R. (`FINDINGS` §12 only cleared
   `/opt/webr/wasm/bin`, a different directory.)

2. **Host-flag contamination — failure class (c).**
   Because `R_HOME/bin/Rscript` now resolves to native R,
   `StanHeaders:::LdFlags()` emits **native** TBB libraries into a wasm link:

   ```
   -L/usr/local/lib/R/site-library/RcppParallel/lib/ -ltbb -ltbbmalloc
   wasm-ld: error: unknown file type: .../libtbb.so
   ```

   This affects both modules. RBesT receives the native flags through
   `StanHeaders:::LdFlags()`; current rstan source sets `PKG_LIBS` from
   `RcppParallel::RcppParallelLibs()` directly. Their `Makevars.webr`
   overrides drop `PKG_LIBS` and **keep** the compile flags: Stan math's
   `init_threadpool_tbb.hpp` needs TBB's *headers* and `-DTBB_INTERFACE_NEW` to
   compile at all. Dropping the link is correct rather than a workaround — the
   prebuilt wasm `RcppParallel` ships an **empty `lib/`**, i.e. upstream links
   no TBB either.

   This is **not** arm64-specific; it is latent on amd64 too.

3. **Generated Stan C++ / Stan library mismatch — failure class (d).**

   ```
   stanExports_gMAP.h:367: no viable conversion from 'rng_t'
     (aka 'mixmax_engine<17, 36, 0>') to 'boost::ecuyer1988'
   ```

   Seen when the tracked `src/stanExports_gMAP.{cc,h}` are compiled against a
   Stan library from a different release (Stan changed its default RNG). It
   does not arise with the CRAN sources, whose stanc is the one that generated the
   tracked sources; staging re-runs `rstantools::rstan_config()` with the host
   rstan regardless. Architecture-independent.

4. **StanHeaders `std::from_chars` portability.**

   Development Stan Math passes `std::string_view` iterators to
   `std::from_chars`. libstdc++ represents those iterators as pointers, but
   Emscripten's libc++ uses wrapped iterators and rejects the call. The build
   applies an exact-match patch using `value.data()` pointer bounds to both the
   host headers used for compilation and the wasm StanHeaders package retained
   for runtime compilation. With the CRAN sources this is inert: StanHeaders
   2.32.10 uses `boost::lexical_cast`, which the patcher recognises and leaves
   alone.

5. **Missing TBB symbols at load time** — the known `FINDINGS` §3 gap, still
   present: a wasm `rstan` imports
   `_ZN3tbb...task_scheduler_observer...` symbols that no wasm build provides.
   The `e2e/` harness patches emscripten's `resolveSymbol` to return a no-op for
   `_ZN3tbb*`. (The wasm `RcppParallel` is no longer the 575-byte stub described
   in §3 — it is 117.6 kB — but its `lib/` is still empty, so the gap remains.)

## Two traps worth remembering

- **`rwasm` caches the source tarball.** `rwasm:::update_repo()` skips
  `make_remote_tarball()` when `repo/src/contrib/<pkg>_<ver>.tar.gz` already
  exists, and `repo/` lives on the persistent `/out` mount. Every rerun
  silently rebuilt the *first* run's staged sources; the build script now
  deletes the tarball before each variant. This cost several confusing
  iterations.
- **`configure` regenerates `src/Makevars`.** RBesT's rstantools `configure`
  runs `rstantools::rstan_config()` at install time, undoing any edit to
  `src/Makevars`. The fix is `src/Makevars.webr`, rwasm's supported override
  hook: it copies the file over `src/Makevars` *and* passes `--no-configure`.
- **`rwasm::add_pkg()` swallows a per-package build failure as a `warning()`,
  not an error.** `update_repo()` wraps its call to `wasm_build()` in
  `tryCatch(..., error = function(cnd) cnd)` and only warns on failure, so
  `add_pkg()` itself returns normally even when the one package you actually
  asked for (RBesT) never got built. `build_variant()` used to trust
  `add_pkg()`'s own success/failure — a `docker compose run` could therefore
  exit 0, pack a VFS image, and still leave no `RBesT_*.tgz` behind, which
  only surfaced later as a confusing `cp: ... No such type` from the
  Makefile. Fixed by verifying the artefact actually exists on disk after
  `add_pkg()` returns, regardless of what it reported, and by capturing any
  warnings via `withCallingHandlers()` so the real per-package failure
  reason ends up in `build-summary.csv` instead of being discarded.
  The patched rstan build is likewise captured in `build-rstan.log`, whose tail
  is included in the fatal missing-artifact error.

## Other findings

- The webR v0.6.0 image ships **Emscripten 5.0.7**, not 4.0.8 as the howto
  states. The emsdk pin must track the image.
- RBesT's `DESCRIPTION` used to carry `Config/Needs/wasm` entries as bare
  `url::https://…tar.gz` refs. pkgdepends treats `Config/Needs/*` as a
  dependency type for local refs and fails with *"Cannot determine package
  names… Maybe you need to add a `<packagename>=` prefix?"*. The field is gone
  now that wasm package sources are resolved by the build rather than declared
  in `DESCRIPTION`; staging still strips it defensively, should one reappear.
- The r-universe development tarballs are **not version-immutable**:
  `StanHeaders_2.39.0.9000.tar.gz` / `rstan_2.39.0.9000.tar.gz` have changed
  while retaining the same version. That is why no route uses them any more
  (see "Where the Stan packages come from"). CRAN releases are immutable; the
  build additionally validates each tarball's package/version metadata and
  records the observed checksum. The exact host `rstan` source is reused for the
  patched wasm build, and host/wasm package versions must agree.
  `DESCRIPTION`'s `Additional_repositories: https://stan-dev.r-universe.dev`
  was removed with the same change: every declared dependency is on CRAN.
- `build/webr/.image-built` and Docker's layers intentionally cache the host
  dependencies. A rerun therefore keeps whatever Stan versions the image was
  built with, even after CRAN or `repo.r-wasm.org` has moved on; to pick those
  up, remove that stamp or run the documented `docker compose build --no-cache`
  before rebuilding. (The GitHub Actions route builds the image from scratch on
  every run, so it always resolves current versions.) Within the persistent
  output repository, patched rstan is rebuilt whenever its upstream checksum or
  injected TBB stub changes.
- Host/wasm Stan version alignment is not something package metadata can
  express, which is why `fetch-wasm-deps.R` checks it explicitly: a mismatched
  Stan library breaks the generated Stan C++ (`no viable conversion from
  'rng_t' ... to 'boost::ecuyer1988'`).

## Is the amd64 copy a cheat?

Not really, but it is a shortcut.

`--platform=linux/amd64` on the `webr` stage only selects a manifest. That
stage is never `RUN`, and what it yields is wasm or text, so nothing amd64 ever
executes. The one visible cost is `webr-host-shims.sh`: webR's `R/Makefile`
copies `$R_HOST/lib/R/bin/{R,Rscript}` from its stage-1 *native* R over the
wasm tree's launchers, so in an amd64 image those two files are amd64 — and we
have to replace them.

A fully from-source, natively built toolchain is possible and needs no shims:
webR's `./configure` builds LLVM flang itself when `EMFC` is unset, and its
`R/Makefile` builds a native stage-1 R before cross-building the wasm R. Both
`ubuntu:26.04` and `emscripten/emsdk:5.0.7` are multi-arch, so such a
Dockerfile has no architecture branching at all.

`Dockerfile.source` is that recipe, complete and commented. It is **untested** —
verified only with `docker build --check` and against the upstream build files —
because it costs 3–7 hours and 25–40 GB (LLVM, then `make -C libs all`, then
two R builds) versus about ten minutes here. See
`../../../design/howto-build-rbest-webr.md` §4.6 for when that is worth it.

## amd64 portability

Reviewed statically, not executed. Nothing branches on architecture except
`check-toolchain.sh`, which accepts both `aarch64` and `x86_64`. The only
`--platform` pin is on the never-executed `webr` stage. Both other base images
are manifest lists covering amd64 and arm64.

```sh
docker buildx build --platform linux/amd64,linux/arm64 .
```

## Verification numbers (arm64, colima, 10 CPU / 32 GB)

Recorded before `include,lib` entered the strip set, so the mount figures below
are the 204 MB image; the posterior comparison is unaffected by stripping.

```
webR ready                    0.6 s
VFS mount (229.1 MB, 17 340 files)   1.0 s
gMAP(AS, chains = 2, iter = 2000, warmup = 1000, cores = 1)   5 s
```

| quantity | webR (wasm, RBesT 1.12.0) | native (RBesT 1.11.0) |
| --- | --- | --- |
| `tau[1]` mean | 0.333 | 0.332 |
| `(Intercept)` mean | −1.12 | −1.11 |
| `theta_resp` mean | 0.248 | 0.249 |
| `theta_resp_pred` mean | 0.256 | 0.252 |

Agreement is within MCMC error (different RNG streams; `Rhat` = 1).
