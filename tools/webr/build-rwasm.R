#!/usr/bin/env Rscript
##
## Cross-compile RBesT and its dependency closure to WebAssembly and pack the
## result into a mountable webR library image.
##
## Runs *inside* the webR development container, which provides the emscripten
## toolchain, a wasm R build and the `rwasm` package:
##
##   docker run --rm -v "$PWD":/github/workspace -w /github/workspace \
##     --entrypoint Rscript ghcr.io/r-wasm/webr:v0.6.0 tools/webr/build-rwasm.R
##
## Environment variables (all optional):
##   RBEST_REPO_DIR   CRAN-like wasm repository output   (default _rwasm/repo)
##   RBEST_IMAGE_DIR  VFS library image output           (default _rwasm/image)
##   RBEST_IMAGE_NAME VFS image base name                (default rbest-library.data)
##   RBEST_STRIP      comma separated dirs to strip from the library
##   RBEST_REMOTES    comma separated webR-patched package refs
##   RBEST_WASM_PKGS  package list  (default tools/webr/wasm-packages.dcf)
##   CRAN_MIRROR      CRAN mirror                        (default cloud.r-project.org)
##   STAN_REPO        repo carrying rstan >= 2.39 (default stan-dev.r-universe.dev)
##
## This deliberately does *not* use `r-wasm/actions/build-rwasm`. That action
## calls `rwasm::add_pkg()` with the package defaults, and two of them are
## wrong for us -- see the comments on `dependencies` and `remotes` below.

repo_dir <- Sys.getenv("RBEST_REPO_DIR", "_rwasm/repo")
image_dir <- Sys.getenv("RBEST_IMAGE_DIR", "_rwasm/image")
image_name <- Sys.getenv("RBEST_IMAGE_NAME", "rbest-library.data")

split_env <- function(name, default) {
  value <- Sys.getenv(name, default)
  value <- strsplit(value, "[[:space:],]+")[[1]]
  value[nzchar(value)]
}

strip_dirs <- split_env(
  "RBEST_STRIP",
  ## `include` is deliberately absent: it saves ~3% of the download and
  ## rstan reads StanHeaders' installed headers on some code paths.
  "demo, doc, examples, help, html, tests, vignette"
)

## `rwasm`'s default `remotes = NA` resolves the whole `inst/webr-remotes` list
## shipped with the package -- 17 refs covering everything webR patches.
## Several of those cannot be resolved against Posit Package Manager at all
## (`Matrix`, and through it `TMB`, `rigraph` and `mmrm`) and abort pkgdepends
## with the opaque "`nrow(out)` must equal `1`", killing the resolution before
## it looks at the packages we asked for. None of them are in RBesT's closure.
## The only webR-patched package that is, is mvtnorm, so name it explicitly.
remotes <- split_env("RBEST_REMOTES", "r-wasm/mvtnorm@webr")

## Resolution runs against CRAN *only*, deliberately. Adding a second
## repository makes every package carried by both resolve twice -- pkgdepends
## returns one candidate per repository and `rwasm:::update_repo()` builds
## every row it is given, so bayesplot, loo, posterior, rstan, rstantools and
## StanHeaders would each be cross-compiled twice. The packages that must come
## from elsewhere get an explicit `url::` ref instead (see wasm-packages.dcf),
## which resolves to exactly one candidate and outranks the CRAN version that
## RBesT's `Imports` would otherwise pull in.
cran_mirror <- Sys.getenv("CRAN_MIRROR", "https://cloud.r-project.org")
stan_repo <- Sys.getenv("STAN_REPO", "https://stan-dev.r-universe.dev")

options(
  repos = c(CRAN = cran_mirror),
  timeout = 1800
)
Sys.setenv(PKG_USE_BIOCONDUCTOR = "false")
options(pkg.use_bioconductor = FALSE)
source("tools/webr/patch-stanheaders.R")

## pkgdepends auto-installs OS system requirements (e.g. `apt-get install`)
## whenever it runs as root, which it does inside the webR container. wasm
## builds cross-compile against the emscripten toolchain, not the host's
## system libraries, so this is both unnecessary and liable to fail wherever
## the container lacks apt-mirror network access. Disable it; system
## requirements are still printed, just not installed.
Sys.setenv(PKG_SYSREQS = "false")
options(pkg.sysreqs = FALSE)

stopifnot(file.exists("DESCRIPTION"))

## Look `package` up in the PACKAGES index of each repository and return a
## pinned `<pkg>=url::<artifact>` ref for the highest
## version found anywhere. This is what removes the hardcoded versions: the
## r-universe rebuilds rstan continuously and deletes the superseded tarball,
## so a literal URL 404s within days. Taking the maximum across repositories
## also means no edit is needed once CRAN ships an rstan new enough for the
## wasm build -- CRAN simply starts winning the comparison.
wasm_pkg_url <- function(package, repos) {
  best <- NULL
  for (repo in repos) {
    db <- tryCatch(
      available.packages(repos = repo, type = "source", filters = character()),
      error = function(e) NULL
    )
    if (is.null(db) || !package %in% rownames(db)) {
      next
    }
    rows <- db[rownames(db) == package, , drop = FALSE]
    for (i in seq_len(nrow(rows))) {
      ## Compare as a version, but keep the string for the URL: the tarball is
      ## named with the literal `Version` field, and `package_version()`
      ## normalises the dash in a Recommended package's version (codetools
      ## 0.2-20 becomes 0.2.20), which would 404.
      version <- rows[i, "Version"]
      if (is.null(best) || package_version(version) > package_version(best$version)) {
        file <- paste0(package, "_", version, ".tar.gz")
        repository <- unname(rows[i, "Repository"])
        if (is.na(repository) || !nzchar(repository)) {
          stop(repo, " provided no download URL for ", package, " ", version)
        }
        url <- if (identical(
          basename(sub("[?].*$", "", repository)),
          file
        )) {
          repository
        } else {
          paste0(
            sub("/+$", "", sub("[?].*$", "", repository)),
            "/",
            file
          )
        }
        best <- list(version = version, url = url)
      }
    }
  }
  if (is.null(best)) {
    stop(
      "package '", package, "' is in none of the configured repositories: ",
      paste(repos, collapse = ", ")
    )
  }
  stan_min_versions <- c(StanHeaders = "2.39", rstan = "2.39")
  if (package %in% names(stan_min_versions) &&
      package_version(best$version) <
        package_version(stan_min_versions[[package]])) {
    stop(
      "resolved ", package, " ", best$version,
      "; the webR build requires >= ", stan_min_versions[[package]]
    )
  }
  sprintf("%s=url::%s", package, best$url)
}

## Which packages need an explicit ref, and why, is declared in
## `tools/webr/wasm-packages.dcf` -- deliberately not in `DESCRIPTION`, since
## pkgdepends parses every `Config/Needs/*` field of its own accord whenever it
## resolves `local::.`, which imposes rules (a `<pkg>=` prefix on every entry)
## on a field that is not otherwise its business.
pkg_list_file <- Sys.getenv("RBEST_WASM_PKGS", "tools/webr/wasm-packages.dcf")
stopifnot(file.exists(pkg_list_file))
pkg_list <- read.dcf(pkg_list_file)

## Fields are collected across every DCF record, not just the first: the blank
## lines that separate the commented blocks in wasm-packages.dcf also start a
## new record, and keeping each field next to the comment explaining it is
## worth more than packing them into one record.
dcf_field <- function(name) {
  if (!name %in% colnames(pkg_list)) {
    return(character(0))
  }
  value <- pkg_list[, name]
  value <- value[!is.na(value)]
  value <- unlist(strsplit(value, "[[:space:],]+"))
  value[nzchar(value)]
}

lookup_repos <- c(cran_mirror, stan_repo)
pinned_ref <- function(name) wasm_pkg_url(name, lookup_repos)

## Overrides for packages already in RBesT's dependency graph go to `remotes`,
## which redirects the download source of a matching row. Adding them as direct
## refs instead would give each of them a second row -- direct *and*
## transitive -- and every row gets built.
remotes <- c(
  remotes,
  vapply(dcf_field("Remote-Refs"), pinned_ref, character(1), USE.NAMES = FALSE)
)

## rstan is not built from the upstream tarball as it stands: Stan math's
## headers instantiate a TBB `task_scheduler_observer` as a namespace-scope
## global, so every reverse-mode AD object imports two TBB runtime entry points
## that nothing in the wasm ecosystem provides, and the module cannot be
## `dlopen()`ed at all. The stubs are added to rstan's own `src/`, so that
## `rstan.so` defines them itself: R loads package DLLs with
## `dyn.load(local = TRUE)` and webR's loader gives each such module a private
## symbol scope, so no package's DLL can satisfy another's undefined symbols.
## RBesT's own binary carries the same stubs, guarded by `__EMSCRIPTEN__`:
## `configure` copies inst/webr/tbb-stubs.cpp into its generated src/. See that
## file and design/howto-build-rbest-webr.md section 9.
source("tools/webr/tbb-patch-common.R")
source("tools/webr/patch-rstan-tarball.R")
stub_file <- webr_stub_source(".")
rstan_ref <- grep("^rstan=", remotes)
if (length(rstan_ref) != 1) {
  stop(
    "expected exactly one rstan ref in `remotes`, found ", length(rstan_ref),
    ". The TBB patch below has nothing to attach to."
  )
}
rstan_patch <- patch_rstan_tarball(
  url = sub("^rstan=url::", "", remotes[rstan_ref]),
  out_dir = Sys.getenv("RBEST_PATCH_DIR", "_rwasm/patched"),
  stub_file = stub_file
)
remotes[rstan_ref] <- patched_ref(rstan_patch)

## Packages that are not in the graph at all have to be direct refs; `remotes`
## ignores names it does not match.
extra <- c(
  vapply(dcf_field("Source-Refs"), pinned_ref, character(1), USE.NAMES = FALSE),
  dcf_field("Extra-Refs")
)

## `dependencies = "hard"` rather than the `rwasm` default of `FALSE`, which
## would build RBesT alone and produce a library image with none of ggplot2,
## dplyr, posterior or bayesplot in it. Not `TRUE`, which additionally pulls
## the *Suggests* of every direct ref; somewhere in that much larger graph
## pkgdepends hits a malformed ref and aborts. Suggests are not wanted in a
## wasm library anyway.
packages <- c(extra, "local::.")

message("== rwasm ", as.character(packageVersion("rwasm")), " ==")
message("== ", R.version.string, " (host) ==")
message("== packages ==")
for (p in packages) message("  ", p)
message("== remotes ==")
for (r in remotes) message("  ", r)
message("== host pre-install: ", paste(dcf_field("Host-Refs"), collapse = ", "), " ==")
message("== strip: ", paste(strip_dirs, collapse = ", "), " ==")

dir.create(repo_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(image_dir, recursive = TRUE, showWarnings = FALSE)

## `rwasm:::wasm_build()` also installs each package's *host* dependencies, via
## `pak::pkg_install("deps::<tarball>")`, so rstan's configure and build
## scripts can find StanHeaders' headers. That call consults only
## `getOption("repos")` -- which is CRAN-only here, and CRAN's StanHeaders is
## too old for rstan 2.39. Install it from the same resolved URL the wasm build
## uses, before that happens.
host_refs <- dcf_field("Host-Refs")
for (pkg in host_refs) {
  ref <- wasm_pkg_url(pkg, lookup_repos)
  url <- sub("^[^=]*=url::", "", ref)
  installed <- tryCatch(packageVersion(pkg), error = function(e) NULL)
  archive <- basename(sub("[?].*$", "", url))
  wanted <- package_version(sub(
    "\\.tar\\.gz$", "",
    substring(archive, nchar(pkg) + 2L)
  ))
  if (!is.null(installed) && installed == wanted) {
    message("  host ", pkg, " ", installed, " matches ", wanted)
    next
  }
  message("  installing host ", pkg, " ", wanted)
  ## `pak::pkg_install()` rather than `install.packages(repos = NULL)`: the
  ## latter cannot fetch dependencies, and StanHeaders needs RcppParallel and
  ## RcppEigen built on the host first. pak resolves those from CRAN, and it is
  ## what `rwasm:::wasm_build()` uses for its own host-side installs.
  pak::pkg_install(ref, ask = FALSE)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("host install of '", pkg, "' failed")
  }
  installed <- packageVersion(pkg)
  if (installed != wanted) {
    stop(
      "host install of ", pkg, " produced ", installed,
      "; expected exactly ", wanted
    )
  }
}
patch_stanheaders_charconv(
  file.path(find.package("StanHeaders"), "include")
)

rwasm::add_pkg(
  packages,
  repo_dir = repo_dir,
  remotes = remotes,
  dependencies = "hard",
  compress = TRUE
)

stanheaders_tgz <- list.files(
  repo_dir,
  pattern = "^StanHeaders_.*[.]tgz$",
  recursive = TRUE,
  full.names = TRUE
)
if (length(stanheaders_tgz) != 1L) {
  stop(
    "expected exactly one StanHeaders wasm binary to patch; found ",
    length(stanheaders_tgz)
  )
}
stanheaders_patch <- patch_stanheaders_archive(stanheaders_tgz)
rwasm::write_packages(repo_dir)

message("\n== packing VFS library image ==")
rwasm::make_vfs_library(
  out_dir = image_dir,
  out_name = image_name,
  repo_dir = repo_dir,
  strip = strip_dirs,
  compress = TRUE
)

## `compress = TRUE` makes `rwasm` write `<name>.data.gz` and drop the
## uncompressed `.data` itself, so there is nothing further to clean up here.

## Record both patched Stan packages next to the image, so neither artefact can
## be mistaken for a stock build of the version it reports.
webr_write_patch_manifest(
  image_dir,
  list(
    stanheaders_patch_manifest_entry(
      as.character(packageVersion("StanHeaders")),
      stanheaders_patch
    ),
    patch_manifest_entry(rstan_patch)
  ),
  stub_file
)

built <- list.files(repo_dir, pattern = "\\.tgz$", recursive = TRUE)
message("\n== built ", length(built), " wasm binaries ==")
if (!any(grepl("^RBesT_", basename(built)))) {
  stop("no RBesT wasm binary was produced")
}
message("== image ==")
for (f in list.files(image_dir, full.names = TRUE)) {
  message("  ", basename(f), " (", round(file.size(f) / 1e6, 2), " MB)")
}
