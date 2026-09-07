##
## Patch Stan Math's STAN_NUM_THREADS parser for Emscripten's libc++.
##
## std::from_chars() takes raw character pointers. Development Stan Math passes
## std::string_view iterators, which happen to work with libstdc++ but fail
## with libc++ because its iterators are wrapped types. Use string_view::data()
## for the equivalent, standard-conforming pointer range. Released StanHeaders
## versions may instead use boost::lexical_cast(), which needs no patch.
##

patch_stanheaders_charconv <- function(root) {
  suffix <- file.path(
    "stan", "math", "prim", "core", "init_threadpool_tbb.hpp"
  )
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE)
  header <- candidates[endsWith(candidates, suffix)]
  if (length(header) != 1L) {
    stop(
      "expected exactly one Stan Math header ending in ", suffix,
      " under ", root, "; found ", length(header)
    )
  }

  lines <- readLines(header, warn = FALSE)
  old_call <- "      = std::from_chars(value.begin(), value.end(), num_threads);"
  new_call <- paste0(
    "      = std::from_chars(value.data(), value.data() + value.size(), ",
    "num_threads);"
  )
  old_end <- "  if (error != std::errc() || end != value.end()"
  new_end <- "  if (error != std::errc() || end != value.data() + value.size()"
  boost_call <- "          = boost::lexical_cast<int>(env_stan_num_threads);"
  boost_catch <- "    } catch (const boost::bad_lexical_cast&) {"

  if (sum(lines == new_call) == 1L && sum(lines == new_end) == 1L &&
      !any(lines %in% c(old_call, old_end))) {
    message("  StanHeaders charconv patch already present: ", header)
    return(invisible(FALSE))
  }
  if (sum(lines == "#include <boost/lexical_cast.hpp>") == 1L &&
      sum(lines == boost_call) == 1L &&
      sum(lines == boost_catch) == 1L &&
      !any(grepl("from_chars", lines, fixed = TRUE))) {
    message("  StanHeaders uses the libc++-compatible Boost parser: ", header)
    return(invisible(FALSE))
  }
  if (sum(lines == old_call) != 1L || sum(lines == old_end) != 1L ||
      any(lines %in% c(new_call, new_end))) {
    stop(
      "unexpected init_threadpool_tbb.hpp charconv implementation in ", header,
      "\nRefusing to patch code that does not match the known upstream form."
    )
  }

  lines[lines == old_call] <- new_call
  lines[lines == old_end] <- new_end
  writeLines(lines, header)
  message("  patched StanHeaders charconv pointer range: ", header)
  invisible(TRUE)
}

patch_stanheaders_archive <- function(archive) {
  archive <- normalizePath(archive)
  before <- if (exists("webr_checksum")) webr_checksum(archive) else {
    paste0("sha256:", unname(tools::sha256sum(archive)))
  }
  td <- tempfile("stanheaders-patch")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::untar(archive, exdir = td)
  source_changed <- patch_stanheaders_charconv(td)

  entries <- list.files(td, all.files = TRUE, no.. = TRUE)
  if (!length(entries)) {
    stop("StanHeaders archive is empty: ", archive)
  }
  package_root <- file.path(td, "StanHeaders")
  if (!dir.exists(package_root)) {
    stop("StanHeaders archive does not contain a StanHeaders/ package root")
  }
  marker <- c(
    "StanHeaders -- CHECKED/PATCHED for the RBesT webR build",
    "",
    "stan/math/prim/core/init_threadpool_tbb.hpp is libc++ compatible. The",
    "build accepts upstream's boost::lexical_cast parser or ensures that",
    "std::from_chars receives raw pointers from std::string_view::data()."
  )
  marker_path <- file.path(package_root, "WEBR-PATCHES")
  marker_changed <- !file.exists(marker_path) ||
    !identical(readLines(marker_path, warn = FALSE), marker)
  if (!source_changed && !marker_changed) {
    message("  StanHeaders archive patch already present: ", archive)
    return(invisible(c(before = before, after = before)))
  }
  writeLines(marker, marker_path)
  old <- setwd(td)
  on.exit(setwd(old), add = TRUE)
  utils::tar(
    archive,
    entries,
    compression = "gzip",
    tar = "internal"
  )
  setwd(old)

  after <- if (exists("webr_checksum")) webr_checksum(archive) else {
    paste0("sha256:", unname(tools::sha256sum(archive)))
  }
  message("  patched StanHeaders archive: ", before, " -> ", after)
  invisible(c(before = before, after = after))
}

stanheaders_patch_manifest_entry <- function(version, patch) {
  c(
    package = "StanHeaders",
    version = version,
    mechanism = paste(
      "header check/patch: accept upstream's Boost parser or ensure",
      "std::from_chars receives string_view::data() pointers for libc++",
      "compatibility"
    ),
    `wasm package before` = unname(patch[["before"]]),
    `wasm package after` = unname(patch[["after"]])
  )
}
