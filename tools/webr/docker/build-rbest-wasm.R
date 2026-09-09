#!/usr/bin/env Rscript
##
## Cross-compile RBesT to WebAssembly on the *native* architecture.
##
## Differences from design/wasm-build/build-rbest-wasm.R (repo-root relative;
## the amd64/qemu harness):
##
##   * dependencies are never cross-compiled -- their prebuilt wasm binaries are
##     fetched by fetch-wasm-deps.R, so no Fortran toolchain is needed. The one
##     exception is rstan, for which no loadable wasm binary is published
##     anywhere: it is built here from a patched source tarball (see below);
##   * only `baseline` is built by default; the notbb / nolto variants are
##     opt-in, to be enabled *after* a failure has been classified.
##
## Expects the RBesT source tree bind-mounted read-only at /src, writes to /out.
##
## Environment variables:
##   RBEST_VARIANTS  comma separated subset of: baseline,notbb,nolto,nolto-notbb
##                   (default: baseline)
##   RBEST_VFS       set to "false" to skip the VFS library image
##   RBEST_STRIP     comma separated strip set for the VFS image
##

if (file.exists("/usr/local/bin/proxy-env.R")) source("/usr/local/bin/proxy-env.R")

src_dir <- "/src"
out_dir <- "/out"
stopifnot(dir.exists(src_dir))

options(timeout = 1800)
Sys.setenv(PKG_USE_BIOCONDUCTOR = "false")
options(pkg.use_bioconductor = FALSE)

library(rwasm)

variants <- trimws(strsplit(Sys.getenv("RBEST_VARIANTS", "baseline"), ",")[[1]])

message("== rwasm ", as.character(packageVersion("rwasm")), " ==")
message("== host R ", R.version.string, " on ", R.version$platform, " ==")
message("== emcc: ", system("emcc --version | head -1", intern = TRUE), " ==")
message("== target webR R ", getOption("rwasm.webr_version"), " ==")

## ---------------------------------------------------------------------------
## Prebuilt dependency closure (no compilation).
## ---------------------------------------------------------------------------
options(rbest.fetch_wasm_deps.no_main = TRUE)
source("/usr/local/bin/fetch-wasm-deps.R")

repo_dir <- file.path(out_dir, "repo")
manifest_path <- file.path(out_dir, "wasm-deps.csv")
if (!identical(Sys.getenv("RBEST_DEPS", "true"), "false")) {
  manifest <- fetch_wasm_deps(repo_dir, src_dir = src_dir, exclude = "rstan")
  write.csv(manifest, manifest_path, row.names = FALSE)
  rwasm::write_packages(repo_dir)
} else {
  if (!file.exists(manifest_path)) {
    stop(
      "RBEST_DEPS=false requires an existing ", manifest_path,
      " from a previous dependency fetch"
    )
  }
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
}

stan_packages <- c("StanHeaders", "rstan")
resolved <- setNames(
  manifest$version[match(stan_packages, manifest$package)],
  stan_packages
)
installed <- vapply(
  stan_packages,
  function(pkg) as.character(packageVersion(pkg)),
  character(1)
)
if (anyNA(resolved) || !identical(unname(resolved), unname(installed))) {
  stop(
    "host/wasm Stan version mismatch\n",
    "  host: ", paste(names(installed), installed, collapse = ", "), "\n",
    "  wasm: ", paste(names(resolved), resolved, collapse = ", "), "\n",
    "Rebuild the image and dependency repo with `make r-binary-webr` after ",
    "removing build/webr/.image-built."
  )
}
message(
  "== host/wasm Stan versions aligned: ",
  paste(names(installed), installed, collapse = ", "), " =="
)
patch_columns <- c("patch_before", "patch_after")
stanheaders_row <- manifest[manifest$package == "StanHeaders", , drop = FALSE]
if (nrow(stanheaders_row) != 1L ||
    !all(patch_columns %in% names(stanheaders_row)) ||
    anyNA(stanheaders_row[1, patch_columns])) {
  stop(
    "wasm dependency manifest has no StanHeaders patch provenance; ",
    "refresh dependencies with RBEST_DEPS=true"
  )
}
stanheaders_patch <- c(
  before = stanheaders_row$patch_before[[1]],
  after = stanheaders_row$patch_after[[1]]
)

## ---------------------------------------------------------------------------
## rstan, cross-compiled here from a patched source tarball.
## ---------------------------------------------------------------------------
## The one dependency that is *not* taken prebuilt. Stan math's headers
## instantiate a TBB `task_scheduler_observer` as a namespace-scope global, so
## every reverse-mode AD object imports `tbb::detail::r1::observe` (and, through
## the profiling map, `deallocate_memory`). Nothing in the wasm ecosystem
## defines them, and because the observer is a *global* its constructor runs
## during module load -- so a wasm rstan without them cannot be `dlopen()`ed at
## all and `library(rstan)` fails in webR. repo.r-wasm.org publishes no wasm
## rstan at all, so it has to be built here regardless.
##
## The symbols have to be defined inside each module that imports them: R loads
## package DLLs with `dyn.load(local = TRUE)`, and webR's Emscripten loader
## gives every such module a private symbol scope, so one package's DLL can
## never satisfy another's undefined symbols. RBesT carries its own copy --
## inst/webr/tbb-stubs.cpp, which `configure` copies into its generated src/ --
## and rstan gets one here. See
## design/howto-build-rbest-webr.md section 9.
source("/usr/local/bin/tbb-patch-common.R")
source("/usr/local/bin/patch-rstan-tarball.R")
stub_file <- webr_stub_source(src_dir)
rstan_url <- Sys.getenv(
  "RBEST_RSTAN_SRC",
  "/usr/local/share/rbest/rstan-source.tar.gz"
)
rstan_patch <- patch_rstan_tarball(
  url = rstan_url,
  out_dir = file.path(out_dir, "patched"),
  stub_file = stub_file
)
## Keep the expensive wasm build only when the persistent repo was built from
## the same upstream source and TBB stub; a matching version is insufficient,
## since the tarball is patched here before it is compiled.
rstan_input_stamp <- file.path(repo_dir, ".rstan-build-inputs")
rstan_build_inputs <- c(
  paste("upstream", rstan_patch$checksum_before),
  paste("stub", webr_checksum(stub_file)),
  paste("rstan patcher", webr_checksum("/usr/local/bin/patch-rstan-tarball.R")),
  paste(
    "StanHeaders patcher",
    webr_checksum("/usr/local/bin/patch-stanheaders.R")
  )
)
cached_rstan_inputs <- if (file.exists(rstan_input_stamp)) {
  readLines(rstan_input_stamp, warn = FALSE)
} else {
  character()
}
if (!identical(cached_rstan_inputs, rstan_build_inputs)) {
  stale_rstan <- list.files(
    repo_dir,
    pattern = "^rstan_.*[.](tar[.]gz|tgz)$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(stale_rstan)) {
    unlink(stale_rstan)
    message("== removed stale cached rstan built from different source ==")
  }
  rwasm::write_packages(repo_dir)
}
rstan_log_path <- file.path(out_dir, "build-rstan.log")
rstan_log <- file(rstan_log_path, open = "wt")
sink(rstan_log, type = "output", split = TRUE)
sink(rstan_log, type = "message")
rstan_warnings <- character()
rstan_status <- withCallingHandlers(
  tryCatch(
    {
      rwasm::add_pkg(
        paste0("local::", rstan_patch$path),
        repo_dir = repo_dir,
        remotes = NULL,
        dependencies = FALSE
      )
      NULL
    },
    error = identity
  ),
  warning = function(w) {
    rstan_warnings <<- c(rstan_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
sink(type = "message")
sink(type = "output")
close(rstan_log)
if (inherits(rstan_status, "error")) {
  stop(conditionMessage(rstan_status), "\nFull rstan build log: ", rstan_log_path)
}
rstan_tgz <- list.files(
  file.path(repo_dir, "bin", "emscripten"),
  pattern = "^rstan_.*[.]tgz$", recursive = TRUE, full.names = TRUE
)
if (!length(rstan_tgz)) {
  log_tail <- tail(readLines(rstan_log_path, warn = FALSE), 40L)
  stop(
    "the patched rstan did not produce a wasm binary. Without it the image ",
    "cannot be loaded in webR at all.\n",
    if (length(rstan_warnings)) {
      paste0("rwasm warning: ", paste(rstan_warnings, collapse = "; "), "\n")
    } else {
      ""
    },
    "Tail of ", rstan_log_path, ":\n",
    paste(log_tail, collapse = "\n")
  )
}
writeLines(rstan_build_inputs, rstan_input_stamp)
message("== patched rstan: ", basename(rstan_tgz[1]), " (",
  round(file.size(rstan_tgz[1]) / 1e6, 2), " MB) ==")
rwasm::write_packages(repo_dir)

## ---------------------------------------------------------------------------
## Stage a per-variant copy of the source tree with the relevant knobs flipped.
## ---------------------------------------------------------------------------
## Whitelist: the working tree can be tens of GB (caches, revdep, docs, check
## dirs), and copying it wholesale is both slow and pointless. Only what
## `R CMD INSTALL` actually reads is staged.
stage_files <- c(
  "DESCRIPTION", "NAMESPACE", "LICENSE", "NEWS.md", ".Rbuildignore",
  "configure", "configure.win", "cleanup",
  "R", "src", "inst", "man", "data", "tools"
)

stage_source <- function(variant) {
  dest <- file.path(out_dir, "src", variant)
  unlink(dest, recursive = TRUE)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  present <- file.path(src_dir, stage_files)
  present <- present[file.exists(present)]
  file.copy(present, dest, recursive = TRUE)

  ## Drop native build artefacts that would confuse a cross-compile.
  unlink(list.files(
    file.path(dest, "src"),
    pattern = "[.](o|so|dll|dylib)$",
    full.names = TRUE
  ))
  unlink(file.path(dest, "src", "package-binary"))

  desc_path <- file.path(dest, "DESCRIPTION")
  desc <- readLines(desc_path, warn = FALSE)

  ## Drop `Config/Needs/wasm` if the package ever grows one again. pkgdepends
  ## treats Config/Needs/* as a dependency type when resolving a `local::` ref,
  ## which is how RBesT is built below, and bare `url::<tarball>` entries in
  ## such a field cannot be resolved:
  ##
  ##   Cannot determine package names for 2 packages: "url::https://.../
  ##   StanHeaders_2.39.0.9000.tar.gz" ... Maybe you need to add a
  ##   `<packagename>=` prefix?
  ##
  ## DESCRIPTION carries no such field today, and wasm package sources are not
  ## its business anyway: StanHeaders comes in as a prebuilt wasm binary and
  ## rstan is built above.
  field_start <- grep("^Config/Needs/wasm:", desc)
  if (length(field_start)) {
    idx <- field_start
    j <- field_start + 1
    while (j <= length(desc) && grepl("^[[:space:]]", desc[j])) {
      idx <- c(idx, j)
      j <- j + 1
    }
    desc <- desc[-idx]
  }

  if (variant %in% c("nolto", "nolto-notbb")) {
    ## Link-time optimisation is a common emscripten failure mode.
    desc <- desc[!grepl("^UseLTO:", desc)]
  }
  writeLines(desc, desc_path)

  ## Regenerate src/stanExports_*.{cc,h} from inst/stan with the host rstan's
  ## stanc, so the generated C++ and the Stan library it is compiled against
  ## are from the same release. With the CRAN sources resolved by
  ## install-host-deps.R (rstan 2.32.7 today) this reproduces the tracked
  ## sources, which stanc 2.32.2 generated; it is kept because nothing pins
  ## that to remain true, and because a mismatch is otherwise found only as a
  ## compile error deep in Stan Math:
  ##
  ##   stanExports_gMAP.h:367: no viable conversion from 'rng_t'
  ##   (aka 'mixmax_engine<17, 36, 0>') to 'boost::ecuyer1988'
  ##
  ## It also generates src/Makevars, which is why the TBB linkage override
  ## below must be written as src/Makevars.webr rather than edited in.
  ## This is architecture-independent -- it would bite the amd64 route too.
  rstantools::rstan_config(dest)
  message("regenerated stanExports with rstan ", packageVersion("rstan"))

  ## The same TBB no-ops that go into rstan, for RBesT's own binary: with
  ## RTLD_LOCAL every module has to define what it imports. `configure` would
  ## do this copy, but it cannot be relied on here -- src/ has just been
  ## regenerated by the rstan_config() call above, and the `notbb` variants
  ## write src/Makevars.webr, which makes rwasm pass --no-configure. So do it
  ## explicitly, for every variant.
  stub_target <- file.path(dest, "src", "webr-tbb-stubs.cpp")
  if (!file.copy(stub_file, stub_target, overwrite = TRUE)) {
    stop("failed to stage the TBB stub source at ", stub_target)
  }
  message("added src/webr-tbb-stubs.cpp from ", stub_file)

  if (variant %in% c("notbb", "nolto-notbb")) {
    ## Remove the RcppParallel/TBB linkage from the rstantools Makevars.
    ##
    ## src/Makevars asks the *host* RcppParallel for its flags, and on any
    ## architecture that yields native libraries:
    ##
    ##   -L.../RcppParallel/lib/ -ltbb -ltbbmalloc
    ##   wasm-ld: error: unknown file type: .../libtbb.so
    ##
    ## There is no wasm TBB to link against either (see FINDINGS section 3), so
    ## the linkage is dropped and RcppParallel's TBB backend compiled out.
    mk_path <- file.path(dest, "src", "Makevars")
    mk <- readLines(mk_path, warn = FALSE)
    ## Both `RcppParallel::*` **and** `StanHeaders:::CxxFlags()/LdFlags()`
    ## expand to native TBB paths, because they are evaluated by the *host* R.
    ## The prebuilt wasm RcppParallel ships an empty lib/, i.e. upstream builds
    ## the wasm stack without TBB, so drop it here as well.
    ## `PKG_CXXFLAGS` is kept: it only adds RcppParallel's *headers* and the
    ## `-DTBB_INTERFACE_NEW` / `-DSTAN_THREADS` defines, which Stan math needs
    ## to compile at all (`init_threadpool_tbb.hpp`). Only `PKG_LIBS` -- the
    ## native `-L .../RcppParallel/lib -ltbb -ltbbmalloc` -- is dropped; the
    ## prebuilt wasm RcppParallel ships an empty `lib/`, i.e. upstream links no
    ## TBB either, and side modules resolve symbols dynamically.
    mk <- mk[!grepl("^PKG_LIBS[[:space:]]*=", mk)]
    mk <- c(mk, "PKG_LIBS =")
    ## Written as `src/Makevars.webr`, which is rwasm's supported override
    ## hook: it copies the file over `src/Makevars` *and* passes
    ## `--no-configure`, so `rstantools::rstan_config()` cannot regenerate the
    ## file and put the TBB flags back.
    writeLines(mk, file.path(dest, "src", "Makevars.webr"))
    if (any(grepl("^PKG_LIBS[[:space:]]*=.+", mk))) {
      stop("failed to strip the native TBB linkage from src/Makevars")
    }
    message(
      "staged src/Makevars for '", variant, "':\n",
      paste0("  ", grep("^PKG_", mk, value = TRUE), collapse = "\n")
    )
  }

  dest
}

## ---------------------------------------------------------------------------
## Build one variant, capturing the full log.
## ---------------------------------------------------------------------------
build_variant <- function(variant) {
  message("\n############ variant: ", variant, " ############")
  pkg_dir <- stage_source(variant)

  ## `rwasm:::update_repo()` decides whether to build a package by reading the
  ## *binary* index, not by looking for the source tarball: it calls
  ## `available.packages()` on bin/emscripten/contrib/<R minor> and `next`s
  ## whenever the index already lists the package at the same version (rwasm
  ## R/repo.R). `repo/` lives on the persistent /out mount, so on any rerun
  ## after a successful build the freshly staged sources are skipped and the
  ## previous binary kept -- silently: no error, no warning, and nothing in the
  ## per-variant log beyond "Processing 1 package(s).". The check below then
  ## reports "reported success but produced no RBesT wasm binary" with no
  ## diagnostic, because the tarballs were deleted but the index still
  ## advertised them.
  ##
  ## So drop the artefacts *and* refresh the indexes, exactly as the rstan
  ## caching above does: `rwasm::write_packages()` rewrites both `PACKAGES`
  ## files from whatever tarballs are actually present. The 60-odd prebuilt
  ## dependencies stay indexed; only RBesT disappears.
  unlink(list.files(
    repo_dir,
    pattern = "^RBesT_.*[.](tar[.]gz|tgz)$",
    recursive = TRUE,
    full.names = TRUE
  ))
  rwasm::write_packages(repo_dir)

  log_path <- file.path(out_dir, paste0("build-", variant, ".log"))

  con <- file(log_path, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
  ## `rwasm::update_repo()` catches a per-package build failure internally
  ## (`tryCatch(wasm_build(...), error = function(cnd) cnd)`) and only
  ## `warning()`s -- it never propagates the error out of `add_pkg()`. So a
  ## "SUCCESS" result here does *not* guarantee RBesT itself was built; it
  ## only means `add_pkg()` didn't throw. Capture warnings so a silent
  ## per-package failure is still diagnosable, and verify the artefact
  ## below rather than trusting this status alone.
  warnings_seen <- character()
  status <- withCallingHandlers(
    tryCatch(
      {
        ## `dependencies = FALSE`: RBesT is the only thing compiled here.
        ## `remotes = NULL`: no remote preference is needed, since nothing but
        ## RBesT is resolved (rwasm's default `NA` resolves its whole
        ## webr-remotes list and fails on refs unrelated to us).
        rwasm::add_pkg(
          paste0("local::", pkg_dir),
          repo_dir = repo_dir,
          remotes = NULL,
          dependencies = FALSE
        )
        "SUCCESS"
      },
      error = function(e) paste("FAILED:", conditionMessage(e))
    ),
    warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  sink(type = "message")
  sink(type = "output")
  close(con)

  built <- list.files(
    repo_dir,
    pattern = "^RBesT_.*[.]tgz$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (identical(status, "SUCCESS") && !length(built)) {
    detail <- if (length(warnings_seen)) {
      paste(warnings_seen, collapse = "; ")
    } else {
      paste0(
        "no diagnostic captured. A log holding nothing but \"Processing 1 ",
        "package(s).\" means update_repo() skipped the build because the ",
        "binary index still listed RBesT; see ", log_path
      )
    }
    status <- paste0(
      "FAILED: rwasm::add_pkg() reported success but produced no RBesT ",
      "wasm binary -- ", detail
    )
  }

  message("variant ", variant, ": ", status)
  if (length(built)) {
    message(
      "  binary: ", basename(built[1]),
      " (", round(file.size(built[1]) / 1e6, 2), " MB)"
    )
  }
  message("  log: ", log_path)

  data.frame(
    variant = variant,
    status = status,
    binary = if (length(built)) basename(built[1]) else NA_character_,
    size_mb = if (length(built)) round(file.size(built[1]) / 1e6, 2) else NA_real_,
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, lapply(variants, build_variant))

message("\n================ SUMMARY ================")
print(results, row.names = FALSE)
write.csv(results, file.path(out_dir, "build-summary.csv"), row.names = FALSE)

if (!any(results$status == "SUCCESS")) {
  message(
    "\nNo variant succeeded. Classify the failure before changing anything:\n",
    "  (a) toolchain/architecture   -- emcc/clang errors, missing sysroot\n",
    "  (b) UseLTO / CXX_STD         -- retry with RBEST_VARIANTS=nolto\n",
    "  (c) host-flag contamination  -- native -L/-l paths in the wasm link;\n",
    "                                  retry with RBEST_VARIANTS=notbb\n",
    "  (d) stanExports_gMAP.* mismatch -- generated C++ vs the Stan library\n"
  )
  quit(status = 1L)
}

## ---------------------------------------------------------------------------
## Self-contained delivery: the whole library as one mountable image.
## ---------------------------------------------------------------------------
if (!identical(Sys.getenv("RBEST_VFS", "true"), "false")) {
  strip <- trimws(strsplit(
    Sys.getenv(
      "RBEST_STRIP",
      "demo,doc,examples,help,html,tests,vignette"
    ),
    ","
  )[[1]])
  vfs_dir <- file.path(out_dir, "vfs")
  message("\n== packing VFS library image (strip: ", paste(strip, collapse = ", "), ") ==")
  status <- tryCatch(
    {
      rwasm::make_vfs_library(
        out_dir = vfs_dir,
        out_name = "rbest-library.data",
        repo_dir = repo_dir,
        compress = TRUE,
        strip = strip
      )
      "SUCCESS"
    },
    error = function(e) paste("FAILED:", conditionMessage(e))
  )
  message("VFS image: ", status)
  ## Record both patched Stan packages next to the image, so neither artefact
  ## can be mistaken for a stock build of the version it reports.
  webr_write_patch_manifest(
    vfs_dir,
    list(
      stanheaders_patch_manifest_entry(
        resolved[["StanHeaders"]],
        stanheaders_patch
      ),
      patch_manifest_entry(rstan_patch)
    ),
    stub_file
  )
  for (f in list.files(vfs_dir, full.names = TRUE)) {
    message("  ", basename(f), " (", round(file.size(f) / 1e6, 2), " MB)")
  }
}
