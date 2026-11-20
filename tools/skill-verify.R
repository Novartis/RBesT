#!/usr/bin/env Rscript
## Executes every fenced R block in the skill's reference files so that the
## documented API cannot drift from the package. Fails loudly: no tryCatch, the
## first broken block aborts the run.

args <- commandArgs(trailingOnly = TRUE)
skill_dir <- if (length(args) >= 1) args[1] else "."
ref_dir <- file.path(skill_dir, "references")
skill_md <- file.path(skill_dir, "SKILL.md")

if (!dir.exists(ref_dir)) {
  stop("no references/ directory under '", skill_dir, "'", call. = FALSE)
}

## Snippets are checked for correctness, not for statistical precision.
options(
  RBesT.MC.chains = 2,
  RBesT.MC.iter = 500,
  RBesT.MC.warmup = 250,
  RBesT.MC.thin = 1
)
message(
  "Reduced sampling for verification: ",
  paste(
    sprintf(
      "%s=%s",
      c("chains", "iter", "warmup", "thin"),
      c(
        getOption("RBesT.MC.chains"),
        getOption("RBesT.MC.iter"),
        getOption("RBesT.MC.warmup"),
        getOption("RBesT.MC.thin")
      )
    ),
    collapse = " "
  )
)

## --- target version --------------------------------------------------------

skill_lines <- readLines(skill_md, warn = FALSE)
declared <- regmatches(
  skill_lines,
  regexpr("RBesT [0-9]+[.-][0-9]+[.-][0-9]+", skill_lines)
)
declared <- unique(sub("^RBesT ", "", unlist(declared)))
if (length(declared) != 1) {
  stop(
    "SKILL.md must declare exactly one target version, found: ",
    paste(declared, collapse = ", "),
    call. = FALSE
  )
}
installed <- as.character(utils::packageVersion("RBesT"))
if (package_version(gsub("-", ".", declared)) != package_version(installed)) {
  stop(
    "skill targets RBesT ",
    declared,
    " but RBesT ",
    installed,
    " is installed; run `make skill-verify`, which verifies against a dev ",
    "install of the working tree",
    call. = FALSE
  )
}
message("Target version ", declared, " matches the installed package.")

## --- link integrity --------------------------------------------------------

ref_files <- sort(list.files(ref_dir, pattern = "[.]md$", full.names = TRUE))
if (length(ref_files) == 0) {
  stop("no reference files found in ", ref_dir, call. = FALSE)
}

links <- unlist(regmatches(
  skill_lines,
  gregexpr("[.]/references/[A-Za-z0-9._-]+[.]md", skill_lines)
))
linked <- unique(basename(links))
missing_target <- setdiff(linked, basename(ref_files))
if (length(missing_target) > 0) {
  stop(
    "SKILL.md links to missing reference file(s): ",
    paste(missing_target, collapse = ", "),
    call. = FALSE
  )
}
unlinked <- setdiff(basename(ref_files), linked)
if (length(unlinked) > 0) {
  stop(
    "reference file(s) not reachable from SKILL.md: ",
    paste(unlinked, collapse = ", "),
    call. = FALSE
  )
}
message("All ", length(ref_files), " reference files are linked from SKILL.md.")

## --- block extraction ------------------------------------------------------

extract_blocks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fences <- grep("^```", lines)
  blocks <- list()
  i <- 1
  while (i < length(fences)) {
    open <- fences[i]
    close <- fences[i + 1]
    if (grepl("^```r\\s*$", lines[open])) {
      body <- if (close > open + 1) lines[(open + 1):(close - 1)] else character()
      blocks[[length(blocks) + 1]] <- list(line = open, code = body)
    }
    i <- i + 2
  }
  blocks
}

## Function names used in call position, so an invented API fails here rather
## than in whatever the snippet happens to do with its result.
is_empty_symbol <- function(x) is.symbol(x) && !nzchar(as.character(x))

called_names <- function(expr, acc = character()) {
  parts <- NULL
  if (is.call(expr)) {
    head <- expr[[1]]
    if (is.name(head)) {
      acc <- c(acc, as.character(head))
    }
    parts <- as.list(expr)[-1]
  } else if (is.pairlist(expr) || is.expression(expr)) {
    parts <- as.list(expr)
  }
  for (i in seq_along(parts)) {
    ## empty symbols are missing formal defaults and empty index arguments
    if (is_empty_symbol(parts[[i]])) next
    acc <- called_names(parts[[i]], acc)
  }
  acc
}

allowed_namespaces <- c(
  "RBesT",
  "base",
  "stats",
  "utils",
  "graphics",
  "grDevices",
  "methods",
  "tools",
  "R_GlobalEnv",
  ""
)

check_symbols <- function(exprs, env, where) {
  names <- setdiff(unique(unlist(lapply(exprs, called_names))), c("::", ":::"))
  reexported <- getNamespaceExports("RBesT")
  for (nm in names) {
    if (!exists(nm, envir = env, mode = "function")) {
      stop(where, ": no such function `", nm, "()`", call. = FALSE)
    }
    if (nm %in% reexported) next
    fn <- get(nm, envir = env, mode = "function")
    ns <- if (is.null(environment(fn))) "base" else environmentName(environment(fn))
    if (!ns %in% allowed_namespaces) {
      stop(
        where,
        ": `",
        nm,
        "()` comes from '",
        ns,
        "', which snippets must reference explicitly with `::`",
        call. = FALSE
      )
    }
  }
}

## --- run -------------------------------------------------------------------

## Some files document behaviour-changing options, so each file starts from the
## same global state regardless of the order they are processed in.
restore_options <- function(snapshot) {
  added <- setdiff(names(options()), names(snapshot))
  if (length(added) > 0) {
    options(stats::setNames(vector("list", length(added)), added))
  }
  options(snapshot)
}

baseline_options <- options()
n_run <- 0L
n_skipped <- 0L

for (path in ref_files) {
  blocks <- extract_blocks(path)
  env <- new.env(parent = globalenv())
  message("\n== ", basename(path), ": ", length(blocks), " block(s)")
  for (k in seq_along(blocks)) {
    block <- blocks[[k]]
    where <- sprintf("%s block %d (line %d)", basename(path), k, block$line)
    exprs <- parse(text = paste(block$code, collapse = "\n"))
    if (length(block$code) > 0 && grepl("^\\s*#\\s*verify:\\s*skip", block$code[1])) {
      message("   skip  ", where, " (marked `# verify: skip`)")
      n_skipped <- n_skipped + 1L
      next
    }
    message("   run   ", where)
    for (e in exprs) {
      eval(e, envir = env)
    }
    ## after evaluation, so functions the block defines itself are in scope
    check_symbols(exprs, env, where)
    n_run <- n_run + 1L
  }
  restore_options(baseline_options)
}

message(
  "\nOK: ",
  n_run,
  " block(s) executed, ",
  n_skipped,
  " skipped, across ",
  length(ref_files),
  " reference file(s)."
)
