# Verification run for the wasm/webR build of RBesT.
#
# Executed inside a webR session by tools/webr/run-webr-vfs.cjs, against a
# freshly built VFS library image and with the package repository pointed at a
# dead address, so nothing can be silently downloaded to cover a gap in the
# image.
#
# Any stop() here fails the CI job. Keep the assertions meaningful: a build can
# succeed and still produce an image that cannot be loaded (see the TBB symbol
# discussion in the troubleshooting section of
# design/howto-build-rbest-webr.md).

report <- function(...) cat("[verify]", ..., "\n")

report("R", as.character(getRversion()), "|", R.version$platform)

## 1. Loading. This is where an unresolved-symbol problem surfaces -- the
##    wasm binary links fine and only fails at dlopen() time.
report("loading rstan ...")
library(rstan)
report("rstan", as.character(packageVersion("rstan")), "loaded")

report("loading RBesT ...")
library(RBesT)
report("RBesT", as.character(packageVersion("RBesT")), "loaded")

## 2. The library must have come from the mounted image, not from anywhere
##    webR might have fetched it.
report("RBesT lib:", dirname(system.file(package = "RBesT")))

## 2b. The wasm rstan is a *patched* build: the two Intel TBB runtime entry
##     points Stan math imports through headers are supplied as no-ops, without
##     which this script cannot get past `library(rstan)` above. Assert the
##     marker is there, so a green verify can never be read as evidence that
##     stock rstan works in webR.
patch_marker <- system.file("WEBR-PATCHES", package = "rstan")
if (!nzchar(patch_marker)) {
  stop("rstan loaded but carries no WEBR-PATCHES marker: which rstan is this?")
}
report("rstan patch marker:", readLines(patch_marker, warn = FALSE)[1])

stanheaders_marker <- system.file("WEBR-PATCHES", package = "StanHeaders")
if (!nzchar(stanheaders_marker)) {
  stop(
    "StanHeaders loaded but carries no WEBR-PATCHES marker: ",
    "the libc++ header fix is not proven"
  )
}
stanheaders_header <- system.file(
  "include", "stan", "math", "prim", "core", "init_threadpool_tbb.hpp",
  package = "StanHeaders"
)
stanheaders_lines <- readLines(stanheaders_header, warn = FALSE)
stan_charconv_fixed <- any(grepl(
  "std::from_chars(value.data(), value.data() + value.size(), num_threads)",
  stanheaders_lines,
  fixed = TRUE
)) && any(grepl(
  "end != value.data() + value.size()",
  stanheaders_lines,
  fixed = TRUE
))
stan_boost_parser <- any(
  stanheaders_lines ==
    "          = boost::lexical_cast<int>(env_stan_num_threads);"
) && any(
  stanheaders_lines == "    } catch (const boost::bad_lexical_cast&) {"
)
if (!xor(stan_charconv_fixed, stan_boost_parser)) {
  stop("StanHeaders does not contain a recognized libc++-compatible parser")
}
report(
  "StanHeaders compatibility marker:",
  readLines(stanheaders_marker, warn = FALSE)[1]
)

## 3. There is no forking and detectCores() is NA in webR, so chains run
##    sequentially. Assert it rather than letting a future webR silently
##    change the answer.
if (!identical(parallel::detectCores(), NA_integer_)) {
  report("note: detectCores() is no longer NA:", parallel::detectCores())
}

## 4. A real fit. Precompiled Stan model, actual sampling.
report("fitting gMAP ...")
elapsed <- system.time(
  map_mc <- gMAP(cbind(r, n - r) ~ 1 | study,
    data = AS,
    family = binomial,
    tau.dist = "HalfNormal",
    tau.prior = 0.5,
    beta.prior = 2,
    chains = 2,
    iter = 2000,
    warmup = 500,
    cores = 1
  )
)
report("gMAP took", sprintf("%.1f s", elapsed[["elapsed"]]))

print(map_mc)

## `summary()` returns a data.frame per quantity, and `is.finite()` has no
## data.frame method, so reduce to a numeric matrix before testing.
fit_summary <- as.matrix(summary(map_mc)$theta.pred)
if (!all(is.finite(fit_summary))) {
  stop("gMAP posterior summary contains non-finite values")
}

## 5. The non-Stan half of the package, on the fit's predictive draws.
report("running automixfit ...")
map <- automixfit(map_mc)
print(map)

if (!all(is.finite(summary(map)))) {
  stop("automixfit summary contains non-finite values")
}

report("OK")
