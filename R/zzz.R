.onLoad <- function(libname, pkgname) {
  if (!("methods" %in% .packages())) attachNamespace("methods")
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}

.onAttach <- function(...) {
  ver <- utils::packageVersion("RBesT")
  packageStartupMessage(
    "This is RBesT version ",
    ver,
    " (released ",
    format(pkg_create_date, "%F"),
    ", git-sha ",
    pkg_sha,
    ")"
  )
}
