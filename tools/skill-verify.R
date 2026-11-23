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

## The quasi-Monte Carlo assurance recipe in references/design.md offers three
## point generators; all of them are executed here.
qmc_pkgs <- c("spacefillr", "qrng", "randtoolbox")
qmc_missing <- qmc_pkgs[!vapply(qmc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(qmc_missing) > 0) {
  stop(
    "skill verification needs the package(s): ",
    paste(qmc_missing, collapse = ", "),
    call. = FALSE
  )
}

## Snippets print plots; keep them out of the working directory.
grDevices::pdf(NULL)

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

## --- permitted exceptions --------------------------------------------------

## The table under "## Permitted exceptions" in SKILL.md: one row per case in
## which the skill allows numerical code outside RBesT.
section_lines <- function(lines, heading) {
  start <- which(lines == heading)
  if (length(start) != 1) {
    stop("SKILL.md must have exactly one '", heading, "' section", call. = FALSE)
  }
  later <- grep("^## ", lines)
  later <- later[later > start]
  end <- if (length(later) > 0) later[1] - 1 else length(lines)
  lines[start:end]
}
table_rows <- grep(
  "^[|]",
  section_lines(skill_lines, "## Permitted exceptions"),
  value = TRUE
)
split_row <- function(row) trimws(strsplit(row, "|", fixed = TRUE)[[1]][-1])
exception_cols <- c("Exception", "When", "Recipe", "Obsolete when", "Confirmed for")
if (length(table_rows) < 2 || !identical(split_row(table_rows[1]), exception_cols)) {
  stop(
    "the 'Permitted exceptions' table in SKILL.md must have the columns: ",
    paste(exception_cols, collapse = " | "),
    call. = FALSE
  )
}
exceptions <- lapply(table_rows[-(1:2)], function(row) {
  cells <- split_row(row)
  if (length(cells) != length(exception_cols)) {
    stop("malformed 'Permitted exceptions' row: ", row, call. = FALSE)
  }
  stats::setNames(as.list(cells), exception_cols)
})
exception_summary <- vapply(
  exceptions,
  function(e) {
    paste0(
      "  * ", e[["Exception"]],
      "\n    obsolete when: ", e[["Obsolete when"]],
      "\n    confirmed for: ", e[["Confirmed for"]]
    )
  },
  character(1)
)
review_hint <- paste0(
  "\nThe skill has to be reviewed against the package for each release; see ",
  "'Agent Skill Release Review' in AGENTS.md. The target version in SKILL.md ",
  "and the 'Confirmed for' entry of every permitted exception are the ",
  "sign-off of that review: change them only after doing it.\n",
  "Permitted exceptions to check",
  if (length(exceptions) == 0) ": none" else paste0(":\n", paste(exception_summary, collapse = "\n"))
)

installed <- as.character(utils::packageVersion("RBesT"))
same_version <- function(a, b) {
  package_version(gsub("-", ".", a)) == package_version(gsub("-", ".", b))
}
if (!same_version(declared, installed)) {
  stop(
    "the skill has not been reviewed for this version of the package: it ",
    "targets RBesT ", declared, ", the package under verification is RBesT ",
    installed, ". (If this is not the dev install of the working tree, run ",
    "`make skill-verify`.)",
    review_hint,
    call. = FALSE
  )
}
message("Target version ", declared, " matches the installed package.")

for (e in exceptions) {
  confirmed <- e[["Confirmed for"]]
  if (!grepl("^[0-9]+[.-][0-9]+[.-][0-9]+$", confirmed) || !same_version(confirmed, declared)) {
    stop(
      "permitted exception '", e[["Exception"]], "' is confirmed for '",
      confirmed, "', not for the target version ", declared, ".",
      review_hint,
      call. = FALSE
    )
  }
  ## the recipe cell names a reference file and a section heading in it
  ref <- regmatches(e[["Recipe"]], regexpr("[.]/references/[A-Za-z0-9._-]+[.]md", e[["Recipe"]]))
  heading <- regmatches(e[["Recipe"]], regexpr("section \"[^\"]+\"", e[["Recipe"]]))
  if (length(ref) != 1 || length(heading) != 1) {
    stop(
      "permitted exception '", e[["Exception"]], "': the Recipe cell must link ",
      "a reference file and name its section as: section \"...\"",
      call. = FALSE
    )
  }
  heading <- sub("^section \"(.*)\"$", "\\1", heading)
  ref_path <- file.path(skill_dir, ref)
  if (!file.exists(ref_path) || !any(readLines(ref_path, warn = FALSE) == paste("##", heading))) {
    stop(
      "permitted exception '", e[["Exception"]], "': no section '## ", heading,
      "' in ", ref,
      call. = FALSE
    )
  }
}
message(
  length(exceptions),
  " permitted exception(s) confirmed for ", declared, ", recipes found."
)

## --- in-file anchors -------------------------------------------------------

anchors <- unique(unlist(regmatches(
  skill_lines,
  gregexpr("[]][(]#[a-z0-9-]+[)]", skill_lines)
)))
anchors <- gsub("^[]][(]#|[)]$", "", anchors)
headings <- sub("^#+ +", "", grep("^#+ ", skill_lines, value = TRUE))
slugs <- gsub(" ", "-", gsub("[^a-z0-9 -]", "", tolower(headings)))
dangling <- setdiff(anchors, slugs)
if (length(dangling) > 0) {
  stop(
    "SKILL.md links to missing section anchor(s): ",
    paste0("#", dangling, collapse = ", "),
    call. = FALSE
  )
}

## --- script template version -----------------------------------------------

## The script template repeats the target version; keep it in step.
template_md <- file.path(ref_dir, "running-code.md")
template_lines <- readLines(template_md, warn = FALSE)
template_versions <- c(
  sub("^RBesT ", "", unlist(regmatches(
    template_lines,
    gregexpr("RBesT [0-9]+[.-][0-9]+[.-][0-9]+", template_lines)
  ))),
  gsub("^.*< \"|\"$", "", unlist(regmatches(
    template_lines,
    gregexpr("packageVersion[(]\"RBesT\"[)] < \"[0-9.-]+\"", template_lines)
  )))
)
if (length(template_versions) == 0) {
  stop("running-code.md: the script template states no RBesT version", call. = FALSE)
}
stale <- unique(template_versions[!vapply(template_versions, same_version, logical(1), b = declared)])
if (length(stale) > 0) {
  stop(
    "running-code.md states RBesT version(s) ",
    paste(stale, collapse = ", "),
    " but the skill targets ", declared,
    call. = FALSE
  )
}

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
