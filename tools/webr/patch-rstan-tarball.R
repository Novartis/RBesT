##
## Add the TBB no-op stubs to rstan's `src/` before it is cross-compiled.
##
## Used by both build routes:
##
##   * tools/webr/build-rwasm.R (GitHub Actions) resolves rstan to a `url::`
##     ref against the stan-dev r-universe and hands that URL here;
##   * tools/webr/docker/build-rbest-wasm.R hands over the tarball the image
##     already resolved and validated at image build time.
##
## Either way the patched tarball is written next to the build output and built
## from there, so the wasm rstan that ends up in the image defines
##
##   tbb::detail::r1::observe(tbb::detail::d1::task_scheduler_observer&, bool)
##   tbb::detail::r1::deallocate_memory(void*)
##
## itself. It has to define them itself: R loads package DLLs with
## `dyn.load(local = TRUE)` and webR's loader gives each such module a private
## symbol scope, so no package's DLL can satisfy another's undefined symbols.
## RBesT carries its own guarded copy: `configure` copies
## inst/webr/tbb-stubs.cpp into its generated src/.
##
## Everything here is explicit and logged -- the checksums of the tarball before
## and after, and the file added. Nothing is edited in place with a `sed`.
##

## `src/Makevars` is expected to glob its sources:
##
##   SOURCES = $(filter-out stan_fit.cpp, $(wildcard *.cpp))
##   OBJECTS = $(SOURCES:.cpp=.o)
##
## which picks up a newly added `src/*.cpp` with no further edit. That is true
## of rstan 2.32 through 2.39, but it is an assumption about a file this
## repository does not control, so it is checked rather than trusted.
assert_makevars_globs_sources <- function(makevars) {
  if (!file.exists(makevars)) {
    stop("rstan tarball has no src/Makevars -- cannot verify how sources are selected")
  }
  mk <- readLines(makevars, warn = FALSE)
  sources <- grep("^[[:space:]]*SOURCES[[:space:]]*(:|\\+)?=", mk, value = TRUE)
  if (!length(sources)) {
    ## No SOURCES at all: R CMD INSTALL's default rule compiles every src/*.cpp.
    message("  src/Makevars sets no SOURCES; R CMD INSTALL globs src/*.cpp")
    return(invisible(TRUE))
  }
  if (!any(grepl("wildcard", sources))) {
    stop(
      "src/Makevars sets SOURCES explicitly rather than globbing:\n  ",
      paste(sources, collapse = "\n  "),
      "\nThe stub would not be compiled. Add webr-tbb-stubs.o to OBJECTS here ",
      "and re-check, rather than assuming."
    )
  }
  message("  src/Makevars globs its sources: ", trimws(sources[1]))
  invisible(TRUE)
}

write_wasm_makevars <- function(makevars) {
  mk <- readLines(makevars, warn = FALSE)
  pkg_libs <- grep(
    "^[[:space:]]*PKG_LIBS[[:space:]]*(\\+|:)?=",
    mk
  )
  if (!length(pkg_libs)) {
    message("  src/Makevars sets no PKG_LIBS; no wasm override needed")
    return(invisible(NULL))
  }
  if (any(!grepl("RcppParallel::RcppParallelLibs", mk[pkg_libs], fixed = TRUE))) {
    stop(
      "rstan src/Makevars contains unexpected PKG_LIBS entries:\n  ",
      paste(mk[pkg_libs], collapse = "\n  "),
      "\nOnly RcppParallel's native TBB linkage can be removed safely."
    )
  }
  mk <- c(mk[-pkg_libs], "", "PKG_LIBS =")
  target <- file.path(dirname(makevars), "Makevars.webr")
  writeLines(mk, target)
  message("  wrote src/Makevars.webr without native RcppParallel/TBB linkage")
  invisible(target)
}

## Download or copy `src` (a URL, or a path to an already-downloaded tarball),
## add the stub and wasm Makevars override to `src/`, repack under `out_dir`,
## and return everything the caller needs to build against it and report it.
patch_rstan_tarball <- function(url, out_dir, stub_file, work_dir = tempfile("rstan-patch")) {
  stopifnot(file.exists(stub_file))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

  message("== patching rstan for wasm ==")
  message("  source: ", url)

  tarball <- file.path(work_dir, basename(url))
  if (file.exists(url)) {
    ## The local docker route hands over the tarball the image already
    ## resolved and metadata-validated, so the rolling r-universe artifact is
    ## fetched once, at image build time.
    if (!file.copy(url, tarball, overwrite = TRUE)) {
      stop("could not copy ", url)
    }
  } else {
    utils::download.file(url, tarball, mode = "wb", quiet = TRUE)
  }
  before <- webr_checksum(tarball)
  message("  source tarball ", basename(tarball), " (",
    round(file.size(tarball) / 1e6, 2), " MB) ", before)

  unpacked <- file.path(work_dir, "src")
  dir.create(unpacked, showWarnings = FALSE)
  utils::untar(tarball, exdir = unpacked)
  pkg_dir <- file.path(unpacked, "rstan")
  if (!dir.exists(pkg_dir)) {
    stop("the rstan tarball does not unpack to a directory called 'rstan'")
  }

  version <- unname(read.dcf(file.path(pkg_dir, "DESCRIPTION"), fields = "Version")[1, 1])
  makevars <- file.path(pkg_dir, "src", "Makevars")
  assert_makevars_globs_sources(makevars)
  write_wasm_makevars(makevars)

  target <- file.path(pkg_dir, "src", "webr-tbb-stubs.cpp")
  if (!file.copy(stub_file, target, overwrite = TRUE)) {
    stop("could not add ", target)
  }
  message("  added src/webr-tbb-stubs.cpp from ", stub_file)

  ## A marker inside the installed package, so the patch is discoverable from
  ## the artefact itself and not only from the build log.
  dir.create(file.path(pkg_dir, "inst"), showWarnings = FALSE)
  writeLines(
    c(
      paste0("rstan ", version, " -- PATCHED for the RBesT webR build"),
      "",
      "src/webr-tbb-stubs.cpp was added before this package was cross-compiled",
      "to WebAssembly. It defines the two Intel TBB runtime entry points that",
      "Stan math imports through headers as no-ops:",
      "",
      "  tbb::detail::r1::observe(tbb::detail::d1::task_scheduler_observer&, bool)",
      "  tbb::detail::r1::deallocate_memory(void*)",
      "",
      "Nothing in the wasm ecosystem provides them, so without this the module",
      "cannot be dlopen()ed at all. webR is single-threaded, so the no-ops are",
      "behaviourally inert. See design/howto-build-rbest-webr.md section 9 in",
      "the RBesT sources.",
      "",
      "src/Makevars.webr removes RcppParallel's native TBB libraries from the",
      "wasm link while retaining its headers and compile-time definitions.",
      "",
      paste("upstream tarball:", url),
      paste("upstream checksum:", before)
    ),
    file.path(pkg_dir, "inst", "WEBR-PATCHES")
  )

  patched <- file.path(normalizePath(out_dir), paste0("rstan_", version, ".tar.gz"))
  ## Packed from the parent of `rstan/` so the tarball keeps the single
  ## top-level directory that pkgdepends and R CMD INSTALL expect.
  old <- setwd(unpacked)
  on.exit(setwd(old), add = TRUE)
  utils::tar(patched, "rstan", compression = "gzip", tar = "internal")
  setwd(old)

  after <- webr_checksum(patched)
  message("  repacked ", patched, " (", round(file.size(patched) / 1e6, 2), " MB) ", after)

  list(
    package = "rstan",
    version = version,
    url = url,
    path = patched,
    checksum_before = before,
    checksum_after = after
  )
}

## The ref form used to feed a local package archive back to pkgdepends.
## `local::<path>` is the documented form for a local file or directory; the
## `url::file://` form is kept as an escape hatch because this is resolved deep
## inside pkgdepends and is the one part of the patch that cannot be verified
## without running a build.
patched_ref <- function(patch, form = Sys.getenv("RBEST_PATCHED_REF_FORM", "local")) {
  switch(form,
    local = sprintf("%s=local::%s", patch$package, patch$path),
    url = sprintf("%s=url::file://%s", patch$package, patch$path),
    stop("unknown RBEST_PATCHED_REF_FORM: ", form)
  )
}

## A manifest entry describing the patch, for webr_write_patch_manifest().
patch_manifest_entry <- function(patch) {
  c(
    package = patch$package,
    version = patch$version,
    mechanism = paste(
      "source patch: TBB stubs added and native TBB linkage removed",
      "before cross-compiling"
    ),
    `upstream tarball` = patch$url,
    `upstream checksum` = patch$checksum_before,
    `patched checksum` = patch$checksum_after
  )
}
