#!/usr/bin/env Rscript
##
## Install RBesT's *native* dependency closure into the image.
##
## The wasm cross-compile runs `R CMD INSTALL` under host R, so:
##   * R CMD INSTALL refuses to start unless Depends/Imports/LinkingTo are
##     installed for the host;
##   * LinkingTo headers (BH, Rcpp, RcppEigen, RcppParallel, rstan,
##     StanHeaders) are taken from the *host* library;
##   * RBesT's src/Makevars shells out to Rscript for
##     StanHeaders:::CxxFlags()/LdFlags() and RcppParallel::CxxFlags(), which
##     therefore also resolve against the host library.
##
## rstan/StanHeaders come from CRAN, at versions derived rather than pinned, so
## the host headers are the same release as the wasm StanHeaders binary the
## image ships (howto-build-rbest-webr.md section 7).

source("/usr/local/bin/proxy-env.R")
source("/usr/local/bin/patch-stanheaders.R")

cran <- Sys.getenv("CRAN_MIRROR", "https://cloud.r-project.org")
options(
  repos = c(CRAN = cran),
  warn = 1,
  timeout = 1800,
  Ncpus = max(1L, parallel::detectCores())
)
Sys.setenv(MAKEFLAGS = paste0("-j", getOption("Ncpus")))

say <- function(...) message("[install-host-deps] ", ...)
say("R: ", R.version.string, " on ", R.version$platform)
say("Ncpus: ", getOption("Ncpus"))

## ---------------------------------------------------------------------------
## The Stan tarballs.
##
## CRAN is the reference for the *source*, and it is the only source: the
## stan-dev r-universe is not consulted at all. CRAN tarballs are immutable,
## which removes the rolling-development-version hazard documented in
## howto-build-rbest-webr.md section 7 (the same `2.39.0.9000` URL served
## different bytes hours apart).
##
## Which CRAN release is usable is *derived*, not pinned:
##
##   StanHeaders   the wasm binary is taken prebuilt, and the host headers must
##                 be the same release as that binary, because rstan's
##                 generated C++ and the Stan library it calls have to agree.
##                 So the version is whatever repo.r-wasm.org publishes for the
##                 target R minor. CRAN keeps only the current release under
##                 src/contrib, so a superseded one comes from
##                 src/contrib/Archive/<pkg>/.
##   rstan         no wasm rstan is published anywhere -- it is cross-compiled
##                 from this very tarball -- so nothing constrains it, and
##                 CRAN's current release is taken.
##
## Deriving rather than pinning means no edit is needed when repo.r-wasm.org
## catches up with CRAN. It does let the two versions drift apart, which
## section 3 below turns into an immediate, diagnosed failure.
## ---------------------------------------------------------------------------
webr_root <- Sys.getenv("WEBR_ROOT", "/opt/webr")
webr_version <- readLines(file.path(webr_root, "R", "R-VERSION"), warn = FALSE)
r_mm <- paste(
  R_system_version(webr_version)$major,
  R_system_version(webr_version)$minor,
  sep = "."
)
wasm_repo <- Sys.getenv("WASM_REPO", "https://repo.r-wasm.org")
cran_contrib <- paste0(sub("/+$", "", cran), "/src/contrib")
wasm_contrib <- file.path(
  sub("/+$", "", wasm_repo), "bin", "emscripten", "contrib", r_mm
)

## The version a CRAN-like repository currently advertises for one package.
## A PACKAGES index carries a single row per package, so this is unambiguous.
indexed_version <- function(name, contriburl, what) {
  db <- tryCatch(
    available.packages(contriburl = contriburl),
    error = function(e) {
      stop("could not read the ", what, " index at ", contriburl, ": ",
           conditionMessage(e))
    }
  )
  ## `available.packages()` only warns when the index is unreachable and
  ## returns an empty matrix, which would otherwise be reported as the package
  ## being absent from an index that does not exist.
  if (!nrow(db)) {
    stop("the ", what, " index at ", contriburl, " is empty or unreachable")
  }
  if (!name %in% rownames(db)) {
    stop(name, " is not in the ", what, " index at ", contriburl)
  }
  unname(db[name, "Version"])
}

stan_versions <- c(
  StanHeaders = indexed_version("StanHeaders", wasm_contrib, "wasm"),
  rstan = indexed_version("rstan", cran_contrib, "CRAN source")
)
say(
  "StanHeaders ", stan_versions[["StanHeaders"]],
  " (the wasm binary published for R ", r_mm, ")"
)
say("rstan ", stan_versions[["rstan"]], " (current on CRAN)")

## A CRAN source tarball has exactly two plausible homes: the current release
## under src/contrib, a superseded one under src/contrib/Archive/<pkg>/.
stan_urls <- function(name, version) {
  file <- paste0(name, "_", version, ".tar.gz")
  c(
    file.path(cran_contrib, file),
    file.path(cran_contrib, "Archive", name, file)
  )
}

resolve_stan <- function(name) {
  list(
    package = name,
    version = unname(stan_versions[[name]]),
    urls = stan_urls(name, stan_versions[[name]])
  )
}

tmp <- tempfile()
dir.create(tmp)

tarball_description <- function(tarball) {
  td <- tempfile()
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  entries <- untar(tarball, list = TRUE)
  description <- grep("^[^/]+/DESCRIPTION$", entries, value = TRUE)
  if (length(description) != 1L) {
    stop("expected one top-level DESCRIPTION in ", basename(tarball))
  }
  untar(tarball, files = description, exdir = td)
  read.dcf(file.path(td, description))
}

download_stan <- function(source, dir) {
  dest <- file.path(
    dir,
    paste0(source$package, "_", source$version, ".tar.gz")
  )
  ok <- FALSE
  for (url in source$urls) {
    say("downloading ", source$package, " ", source$version, ": ", url)
    status <- tryCatch(
      utils::download.file(url, dest, mode = "wb", quiet = TRUE),
      error = function(e) {
        say("  not available here: ", conditionMessage(e))
        1L
      }
    )
    if (identical(as.integer(status), 0L) && file.exists(dest) &&
        file.size(dest) > 0) {
      ok <- TRUE
      break
    }
  }
  if (!ok) {
    stop(
      "could not download ", source$package, " ", source$version, " from:\n  ",
      paste(source$urls, collapse = "\n  ")
    )
  }
  desc <- tarball_description(dest)
  if (!identical(unname(desc[1, "Package"]), source$package) ||
      !identical(unname(desc[1, "Version"]), source$version)) {
    stop(
      "downloaded tarball identity mismatch for ", source$package, "\n",
      "  expected: ", source$package, " ", source$version, "\n",
      "  got:      ", desc[1, "Package"], " ", desc[1, "Version"]
    )
  }
  say("observed sha256 for ", source$package, ": ", tools::sha256sum(dest))
  dest
}

## ---------------------------------------------------------------------------
## 1. Everything except the Stan packages, from CRAN.
## ---------------------------------------------------------------------------
cran_pkgs <- c(
  ## LinkingTo
  "BH", "Rcpp", "RcppEigen", "RcppParallel",
  ## Imports (R CMD INSTALL checks these are present)
  "rstantools", "posterior", "assertthat", "mvtnorm", "Formula", "checkmate",
  "bayesplot", "ggplot2", "dplyr", "matrixStats", "statmod", "abind", "rlang",
  "jsonlite", "lifecycle", "Rdpack",
  ## Recommended, needed in the wasm image (howto section 7)
  "codetools"
)
missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing)) {
  say("installing from CRAN: ", paste(missing, collapse = ", "))
  install.packages(missing)
}
still_missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(still_missing)) {
  stop("failed to install: ", paste(still_missing, collapse = ", "))
}

## ---------------------------------------------------------------------------
## 2. The resolved Stan packages, in dependency order.
## ---------------------------------------------------------------------------
base_pkgs <- rownames(installed.packages(priority = "base"))

## The resolved tarballs are installed with repos = NULL, which does not
## resolve dependencies, so their own hard dependencies must be present first.
install_tarball_deps <- function(tarball) {
  desc <- tarball_description(tarball)
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))
  deps <- unlist(strsplit(paste(desc[1, fields], collapse = ","), ","))
  deps <- trimws(sub("\\(.*", "", deps))
  deps <- setdiff(deps[nzchar(deps)], c("R", base_pkgs, rownames(installed.packages())))
  if (length(deps)) {
    say("dependencies of ", basename(tarball), ": ", paste(deps, collapse = ", "))
    install.packages(deps)
    still <- setdiff(deps, rownames(installed.packages()))
    if (length(still)) stop("failed to install: ", paste(still, collapse = ", "))
  }
  invisible(NULL)
}

for (name in c("StanHeaders", "rstan")) {
  source <- resolve_stan(name)
  tarball <- download_stan(source, tmp)
  install_tarball_deps(tarball)
  say("installing ", name, " from ", basename(tarball))
  install.packages(tarball, repos = NULL, type = "source")
  if (!name %in% rownames(installed.packages())) {
    stop(
      "failed to install ", name, " ", source$version, ".\n",
      "rstan is compiled against the StanHeaders installed just before it, so ",
      "the likely cause is that the two derived versions do not agree on ",
      "Stan's RNG API -- look for 'ecuyer1988' or 'rng_t' in the compile ",
      "output above. See howto-build-rbest-webr.md section 7."
    )
  }
  installed <- packageVersion(name)
  if (installed != package_version(source$version)) {
    stop(
      "installed ", name, " ", installed,
      " does not match the resolved version ", source$version
    )
  }
  say(name, " ", installed, " installed")
  if (identical(name, "StanHeaders")) {
    patch_stanheaders_charconv(
      file.path(find.package("StanHeaders"), "include")
    )
  }
  if (identical(name, "rstan")) {
    ## Kept: build-rbest-wasm.R cross-compiles rstan from this very tarball
    ## (with the TBB stubs added), rather than taking a published wasm binary
    ## -- no wasm rstan is published on repo.r-wasm.org at all, and stan-dev's
    ## cannot be loaded in webR. Reusing the file that was just validated is
    ## what guarantees the wasm rstan and the host rstan that
    ## `rstantools::rstan_config()` generates against are the same content.
    dir.create("/usr/local/share/rbest", recursive = TRUE, showWarnings = FALSE)
    file.copy(tarball, "/usr/local/share/rbest/rstan-source.tar.gz", overwrite = TRUE)
    say("kept the rstan source at /usr/local/share/rbest/rstan-source.tar.gz")
  }
}

unlink(tmp, recursive = TRUE)

## ---------------------------------------------------------------------------
## 3. Prove the flags RBesT's src/Makevars depends on are obtainable.
## ---------------------------------------------------------------------------
stan_cxxflags <- StanHeaders:::CxxFlags(as_character = TRUE)
say("StanHeaders:::CxxFlags(): ", stan_cxxflags)
say("StanHeaders:::LdFlags():  ", StanHeaders:::LdFlags(as_character = TRUE))
say("RcppParallel::CxxFlags(): ", RcppParallel::CxxFlags())

## ---------------------------------------------------------------------------
## 4. The two derived versions must agree on Stan's RNG API.
##
## rstan 2.36 replaced `boost::ecuyer1988` with `stan::rng_t`, and StanHeaders
## selects between them by defining `NEW_RSTAN` in its CxxFlags() -- keyed on
## the version of rstan it finds installed. Since the pair is derived rather
## than pinned, the two can move independently, and a StanHeaders whose
## selection logic does not match the installed rstan surfaces only much later,
## as a compile error deep in Stan Math:
##
##   stanExports_gMAP.h: no viable conversion from 'rng_t'
##   (aka 'mixmax_engine<17, 36, 0>') to 'boost::ecuyer1988'
##
## Check the selection here, where both versions are still in hand, and say
## which type the pair has actually settled on.
## ---------------------------------------------------------------------------
new_rstan <- grepl("-DNEW_RSTAN([[:space:]=]|$)", stan_cxxflags)
expected_new_rstan <- packageVersion("rstan") >= package_version("2.36")
say(
  "Stan RNG type: ",
  if (new_rstan) "stan::rng_t (NEW_RSTAN defined)" else "boost::ecuyer1988"
)
if (!identical(new_rstan, expected_new_rstan)) {
  stop(
    "StanHeaders ", packageVersion("StanHeaders"), " ",
    if (new_rstan) "defines" else "does not define", " NEW_RSTAN for rstan ",
    packageVersion("rstan"), ", which expects the opposite.\n",
    "The two derived versions do not agree on Stan's RNG API; RBesT's ",
    "generated stanExports_gMAP.* would not compile. See ",
    "howto-build-rbest-webr.md section 7."
  )
}
say("done")
