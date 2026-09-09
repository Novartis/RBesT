#!/usr/bin/env Rscript
##
## Populate a CRAN-like wasm repo with *prebuilt* binaries for RBesT's runtime
## dependency closure.
##
## Nothing here is compiled: every dependency taken from this repo already has
## a published wasm binary, so the local toolchain only has to build RBesT and
## rstan (the one package with no usable published binary -- see `exclude`).
## That is what removes the need for `flang` (which is amd64-only) and keeps
## this route viable on arm64.
##
## Sourced by build-rbest-wasm.R; can also be run standalone:
##   R -q -f /usr/local/bin/fetch-wasm-deps.R
##

if (file.exists("/usr/local/bin/proxy-env.R")) {
  source("/usr/local/bin/proxy-env.R")
}
source("/usr/local/bin/patch-stanheaders.R")

## One repository. Every wasm binary the closure needs is published here: the
## full closure was walked against this index and rstan is the only miss, and
## rstan is cross-compiled locally regardless (see `exclude` below), because no
## wasm rstan that can be loaded in webR is published anywhere.
wasm_repos <- c(
  r_wasm = "https://repo.r-wasm.org"
)

webr_root <- Sys.getenv("WEBR_ROOT", "/opt/webr")
webr_version <- readLines(file.path(webr_root, "R", "R-VERSION"), warn = FALSE)
r_mm <- paste(
  R_system_version(webr_version)$major,
  R_system_version(webr_version)$minor,
  sep = "."
)

contrib_path <- function(repo) {
  file.path(repo, "bin", "emscripten", "contrib", r_mm)
}

say <- function(...) message("[fetch-wasm-deps] ", ...)

## StanHeaders is the one package whose wasm binary must be the *same release*
## as the host headers RBesT and rstan are compiled against, because rstan's
## generated C++ and the Stan library it calls have to agree. install-host-deps.R
## derives the host version *from* this index, so equality holds by
## construction; the check stays because it is cheap and because it is the one
## thing that would catch the index having moved on between the image build and
## this run.
host_stanheaders_version <- function() {
  as.character(packageVersion("StanHeaders"))
}

## Packages already present in webR's wasm R installation; never fetch these.
wasm_base_pkgs <- function() {
  lib <- file.path(webr_root, "wasm", paste0("R-", webr_version), "lib", "R", "library")
  c(list.dirs(lib, full.names = FALSE, recursive = FALSE), "R")
}

## Hard dependencies declared by a package DESCRIPTION or an available.packages
## row.
hard_deps <- function(fields) {
  deps <- unlist(strsplit(paste(fields[!is.na(fields)], collapse = ","), ","))
  deps <- trimws(sub("\\(.*", "", deps))
  deps[nzchar(deps)]
}

fetch_wasm_deps <- function(repo_dir, src_dir = "/src", extra = "codetools",
                            exclude = character()) {
  dest <- file.path(repo_dir, "bin", "emscripten", "contrib", r_mm)
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(repo_dir, "src", "contrib"), recursive = TRUE, showWarnings = FALSE)

  ## One available.packages() database per repo, kept separate so that the
  ## priority order above is honoured per package.
  dbs <- lapply(wasm_repos, function(repo) {
    db <- tryCatch(
      available.packages(contriburl = contrib_path(repo)),
      error = function(e) NULL
    )
    if (is.null(db) || nrow(db) == 0) {
      say("WARNING: no wasm PACKAGES index at ", contrib_path(repo))
      NULL
    } else {
      say(repo, ": ", nrow(db), " wasm binaries indexed")
      db
    }
  })
  dbs <- dbs[!vapply(dbs, is.null, logical(1))]
  if (!length(dbs)) stop("no wasm package index could be read")

  lookup <- function(pkg) {
    for (name in names(dbs)) {
      db <- dbs[[name]]
      if (pkg %in% rownames(db)) {
        row <- db[pkg, ]
        if (identical(pkg, "StanHeaders")) {
          host <- host_stanheaders_version()
          if (!identical(unname(row[["Version"]]), host)) {
            stop(
              name, " offers StanHeaders ", row[["Version"]],
              ", but the host has ", host, ".\n",
              "The wasm and host Stan libraries must be the same release. ",
              "The host version is derived from this same index at image ",
              "build time, so the index has moved since: rebuild the image ",
              "(remove build/webr/.image-built)."
            )
          }
        }
        return(list(repo = name, row = row))
      }
    }
    NULL
  }

  ## A package that is built here rather than downloaded still belongs in the
  ## closure walk -- its own dependencies do have to be fetched. It need not be
  ## in any wasm index, and rstan is not in one: no loadable wasm rstan is
  ## published. Its metadata therefore comes from the host installation, which
  ## install-host-deps.R has pinned to the very tarball that will be
  ## cross-compiled, so the recorded version is by construction the one built.
  host_built <- function(pkg) {
    desc <- read.dcf(file.path(find.package(pkg), "DESCRIPTION"))
    fields <- c("Version", "Depends", "Imports", "LinkingTo")
    row <- vapply(
      fields,
      function(f) if (f %in% colnames(desc)) unname(desc[1, f]) else NA_character_,
      character(1)
    )
    list(repo = "host", row = row)
  }

  ## Seed the closure from RBesT's own DESCRIPTION.
  desc <- read.dcf(file.path(src_dir, "DESCRIPTION"))
  seeds <- hard_deps(desc[1, intersect(c("Depends", "Imports"), colnames(desc))])
  seeds <- unique(c(seeds, extra))

  base_pkgs <- wasm_base_pkgs()
  todo <- setdiff(seeds, base_pkgs)
  seen <- character()
  resolved <- list()
  missing <- character()

  while (length(todo)) {
    pkg <- todo[[1]]
    todo <- todo[-1]
    if (pkg %in% seen || pkg %in% base_pkgs) next
    seen <- c(seen, pkg)

    hit <- if (pkg %in% exclude) host_built(pkg) else lookup(pkg)
    if (is.null(hit)) {
      missing <- c(missing, pkg)
      next
    }
    resolved[[pkg]] <- hit
    todo <- unique(c(
      todo,
      setdiff(hard_deps(hit$row[c("Depends", "Imports", "LinkingTo")]), c(seen, base_pkgs))
    ))
  }

  if (length(missing)) {
    stop(
      "no wasm binary published for: ", paste(sort(missing), collapse = ", "),
      "\nEither the package must be cross-compiled locally too, or the repo ",
      "list in fetch-wasm-deps.R needs extending."
    )
  }

  say("closure: ", length(resolved), " packages")
  manifest <- data.frame(
    package = character(), version = character(), repo = character(),
    size_kb = numeric(), patch_before = character(), patch_after = character(),
    stringsAsFactors = FALSE
  )

  for (pkg in sort(names(resolved))) {
    hit <- resolved[[pkg]]
    ver <- unname(hit$row[["Version"]])
    ## Excluded packages are still walked -- their own dependencies belong in
    ## the closure -- but no binary is downloaded, because it is built here
    ## instead. rstan is excluded for exactly one reason: no wasm rstan that
    ## can be loaded in webR is published anywhere (see the TBB stubs in
    ## build-rbest-wasm.R and design/howto-build-rbest-webr.md section 9), so
    ## it has no index row either and its metadata comes from the host.
    if (pkg %in% exclude) {
      say(sprintf("  %-14s %-16s %-8s %s", pkg, ver, hit$repo, "SKIPPED, built from patched source"))
      manifest <- rbind(manifest, data.frame(
        package = pkg, version = ver, repo = hit$repo,
        size_kb = NA_real_, patch_before = NA_character_,
        patch_after = NA_character_, stringsAsFactors = FALSE
      ))
      next
    }
    file <- paste0(pkg, "_", ver, ".tgz")
    target <- file.path(dest, file)
    repository <- unname(hit$row[["Repository"]])
    if (is.na(repository) || !nzchar(repository)) {
      stop(hit$repo, " provided no download URL for ", pkg, " ", ver)
    }
    ## repo.r-wasm.org exposes a contrib directory rather than a
    ## content-addressed artifact URL, but the check costs nothing and keeps
    ## the function usable against either shape.
    url <- if (identical(basename(sub("[?].*$", "", repository)), file)) {
      repository
    } else {
      paste0(sub("/+$", "", sub("[?].*$", "", repository)), "/", file)
    }
    ## Every source is now a versioned CRAN-like artifact, so a file already
    ## present under the persistent /out mount is the right one.
    if (!file.exists(target)) {
      utils::download.file(url, target, mode = "wb", quiet = TRUE)
    }
    stanheaders_patch <- if (identical(pkg, "StanHeaders")) {
      patch_stanheaders_archive(target)
    } else {
      c(before = NA_character_, after = NA_character_)
    }
    size_kb <- file.size(target) / 1024
    say(sprintf("  %-14s %-16s %-8s %8.1f kB", pkg, ver, hit$repo, size_kb))
    ## A binary this small is a stub, not a real package (see FINDINGS section 3
    ## on repo.r-wasm.org's RcppParallel).
    if (size_kb < 5) {
      say("    WARNING: suspiciously small -- likely a stub binary")
    }
    manifest <- rbind(manifest, data.frame(
      package = pkg, version = ver, repo = hit$repo,
      size_kb = round(size_kb, 1),
      patch_before = unname(stanheaders_patch[["before"]]),
      patch_after = unname(stanheaders_patch[["after"]]),
      stringsAsFactors = FALSE
    ))
  }

  manifest
}

if (!isTRUE(getOption("rbest.fetch_wasm_deps.no_main"))) {
  options(timeout = 1800)
  repo_dir <- Sys.getenv("RBEST_REPO_DIR", "/out/repo")
  ## Same exclusion as build-rbest-wasm.R: rstan has no wasm binary to fetch,
  ## it is cross-compiled from the patched host source.
  manifest <- fetch_wasm_deps(repo_dir, exclude = "rstan")
  write.csv(
    manifest,
    file.path(dirname(repo_dir), "wasm-deps.csv"),
    row.names = FALSE
  )
  library(rwasm)
  rwasm::write_packages(repo_dir)
}
