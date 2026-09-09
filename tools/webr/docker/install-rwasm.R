#!/usr/bin/env Rscript
##
## Install rwasm (and its dependency closure) into the webR build container.
##
## Split out of the Dockerfile because doing this robustly behind a corporate
## proxy needs more than a one-liner:
##
##   * The container inherits no proxy settings unless they are passed as
##     build args, and the resulting failure ("package 'pak' is not available
##     for this version of R") points at the wrong thing entirely -- the real
##     cause is that the CRAN index could not be fetched at all.
##   * pak/pkgdepends contacts bioconductor.org on every run, even when
##     nothing in the dependency graph is a Bioconductor package. On a
##     locked-down network that is a gratuitous failure point.
##   * rwasm is not on CRAN and is not in any r-universe, so it has to come
##     from GitHub. Different firewalls block different GitHub hosts, so we
##     try more than one route.
##
## Environment variables (set from Dockerfile build args):
##   CRAN_MIRROR   CRAN mirror to use          (default https://cloud.r-project.org)
##   RWASM_REF     git ref of r-wasm/rwasm     (default HEAD)

source("/usr/local/bin/proxy-env.R")

cran <- Sys.getenv("CRAN_MIRROR", "https://cloud.r-project.org")
ref <- Sys.getenv("RWASM_REF", "HEAD")
options(repos = c(CRAN = cran), warn = 1, timeout = 600)

## pkgdepends' env prefix is "pkg.", so this is the documented way to switch
## the Bioconductor lookup off.
Sys.setenv(PKG_USE_BIOCONDUCTOR = "false")
options(pkg.use_bioconductor = FALSE)

say <- function(...) message("[install-rwasm] ", ...)

say("R:    ", R.version.string)
say("CRAN: ", cran)
say("ref:  ", ref)
proxy <- Sys.getenv("https_proxy", Sys.getenv("http_proxy", ""))
say("proxy: ", if (nzchar(proxy)) proxy else "(none)")

## ---------------------------------------------------------------------------
## Preflight: fail loudly and usefully if there is simply no network.
## ---------------------------------------------------------------------------
reachable <- function(url) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch(
    {
      suppressWarnings(utils::download.file(url, tmp, quiet = TRUE))
      file.exists(tmp) && file.size(tmp) > 0
    },
    error = function(e) FALSE
  )
  isTRUE(ok)
}

if (!reachable(paste0(cran, "/src/contrib/PACKAGES.gz"))) {
  stop(
    "\n\n",
    "Cannot reach the CRAN index at ", cran, "\n\n",
    "This is almost always a proxy problem rather than an R problem: if\n",
    "the index cannot be fetched, install.packages() reports the\n",
    "misleading 'package is not available for this version of R'.\n\n",
    "If you are behind a corporate proxy, pass it into the *build*:\n\n",
    "    export HTTP_PROXY=http://proxy.example.com:8080\n",
    "    export HTTPS_PROXY=$HTTP_PROXY\n",
    "    docker compose build \\\n",
    "        --build-arg http_proxy=$HTTP_PROXY \\\n",
    "        --build-arg https_proxy=$HTTPS_PROXY\n\n",
    "(compose.yaml already forwards these from the environment, so\n",
    " exporting them and re-running `docker compose run --rm build`\n",
    " is normally enough.)\n\n",
    "If you have an internal CRAN mirror instead, use:\n\n",
    "    docker compose build --build-arg CRAN_MIRROR=https://internal/cran\n",
    call. = FALSE
  )
}
say("CRAN index reachable")

## ---------------------------------------------------------------------------
## rwasm is GitHub-only. Try the routes in order of how likely they are to be
## permitted by a restrictive firewall.
## ---------------------------------------------------------------------------
install_via_remotes <- function() {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    utils::install.packages("remotes")
  }
  remotes::install_github(
    paste0("r-wasm/rwasm", if (!identical(ref, "HEAD")) paste0("@", ref) else ""),
    upgrade = "never"
  )
}

install_via_codeload <- function() {
  ## Some networks allow codeload.github.com but block api.github.com.
  tar <- tempfile(fileext = ".tar.gz")
  br <- if (identical(ref, "HEAD")) "main" else ref
  utils::download.file(
    sprintf("https://codeload.github.com/r-wasm/rwasm/tar.gz/%s", br),
    tar, quiet = TRUE
  )
  ex <- tempfile()
  dir.create(ex)
  utils::untar(tar, exdir = ex)
  pkg <- list.files(ex, full.names = TRUE)[1]
  ## Dependencies from CRAN first, then the package itself.
  desc <- read.dcf(file.path(pkg, "DESCRIPTION"))
  fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(desc))
  deps <- unlist(strsplit(paste(desc[1, fields], collapse = ","), ","))
  deps <- trimws(sub("\\(.*", "", deps))
  deps <- setdiff(deps, c("R", "", NA))
  deps <- deps[!deps %in% rownames(utils::installed.packages())]
  if (length(deps)) utils::install.packages(deps)
  utils::install.packages(pkg, repos = NULL, type = "source")
}

install_via_pak <- function() {
  if (!requireNamespace("pak", quietly = TRUE)) {
    utils::install.packages("pak")
  }
  pak::pak(paste0("r-wasm/rwasm", if (!identical(ref, "HEAD")) paste0("@", ref) else ""))
}

routes <- list(
  remotes = install_via_remotes,
  codeload = install_via_codeload,
  pak = install_via_pak
)

errs <- character()
for (nm in names(routes)) {
  if (requireNamespace("rwasm", quietly = TRUE)) break
  say("trying route: ", nm)
  ok <- tryCatch(
    {
      routes[[nm]]()
      TRUE
    },
    error = function(e) {
      errs[[nm]] <<- conditionMessage(e)
      say("route '", nm, "' failed: ", conditionMessage(e))
      FALSE
    }
  )
  if (ok && requireNamespace("rwasm", quietly = TRUE)) {
    say("installed rwasm via '", nm, "'")
    break
  }
}

if (!requireNamespace("rwasm", quietly = TRUE)) {
  stop(
    "\n\nCould not install rwasm by any route.\n\n",
    paste(sprintf("  %-9s %s", names(errs), errs), collapse = "\n"),
    "\n\nIf GitHub is unreachable from this container, clone rwasm on the\n",
    "host and mount it, then install with:\n",
    "    R -e 'install.packages(\"/path/to/rwasm\", repos = NULL, type = \"source\")'\n",
    call. = FALSE
  )
}

say("rwasm ", as.character(packageVersion("rwasm")), " OK")
