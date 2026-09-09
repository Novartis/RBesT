## Normalise proxy environment variables.
##
## compose.yaml forwards both the lower- and upper-case spellings, so one of
## each pair is typically the empty string. An empty `http_proxy` is not the
## same as an unset one -- libcurl treats it as "no proxy" and it therefore
## masks a perfectly good `HTTP_PROXY`. Collapse each pair to the non-empty
## value and unset it entirely if neither is set.
##
## Sourced by install-rwasm.R and build-rbest-wasm.R.

local({
  pairs <- list(
    c("http_proxy", "HTTP_PROXY"),
    c("https_proxy", "HTTPS_PROXY"),
    c("no_proxy", "NO_PROXY")
  )
  for (p in pairs) {
    vals <- Sys.getenv(p)
    vals <- vals[nzchar(vals)]
    if (length(vals)) {
      args <- as.list(rep(vals[[1]], 2))
      names(args) <- p
      do.call(Sys.setenv, args)
    } else {
      Sys.unsetenv(p)
    }
  }
})
