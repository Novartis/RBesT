#! /usr/bin/env Rscript

## rm(list=ls())

## generate example data sets
make_ds <- function() {
  colitis <<- data.frame(
    study = c("Van_assche", "Feagan", "Rutgeerts-1", "Rutgeerts-2"),
    n = c(56, 63, 121, 123),
    r = c(6, 9, 18, 7)
  )

  AS <<- data.frame(
    study = paste("Study", 1:8),
    n = c(107, 44, 51, 39, 139, 20, 78, 35),
    r = c(23, 12, 19, 9, 39, 6, 9, 10)
  )

  transplant <<- data.frame(
    study = paste("Study", 1:11),
    n = c(33, 45, 74, 103, 140, 49, 83, 59, 22, 109, 213),
    r = c(6, 8, 17, 28, 26, 8, 22, 8, 6, 16, 53)
  )

  crohn <<- dat <- data.frame(
    study = c(
      "Gastr06",
      "AIMed07",
      "NEJM07",
      "Gastr01a",
      "APhTh04",
      "Gastr01b"
    ),
    n = c(74, 166, 328, 20, 25, 58),
    y = c(-51, -49, -36, -47, -90, -54)
  )

  asthma <<- data.frame(
    study = paste0("Study_", 1:10),
    NCT = c(
      "ISRCTN75169762",
      "NCT01000506",
      "NCT01691521",
      "NCT01287039",
      "NCT01285323",
      "NCT01238861",
      "NCT01545440; NCT01545453",
      "NCT01691508",
      "NCT00292877",
      "NCT01312961"
    ),
    d = c(0.96, 1.00, 0.62, 1.00, 1.00, 1.00, 1.00, 0.46, 0.50, 0.23),
    n = c(32, 155, 191, 244, 232, 80, 66, 66, 11, 52),
    mu_hat = c(3.73, 2.40, 1.75, 1.80, 2.11, 0.57, 0.88, 2.12, 1.99, 4.64),
    log_mu_hat = c(1.32, 0.88, 0.56, 0.59, 0.75, -0.56, -0.13, 0.75, 0.69, 1.53),
    se_log_mu_hat = c(0.18, 0.11, 0.11, 0.10, 0.13, 0.16, 0.28, 0.16, 0.31, 5.28),
    kappa_hat = c(0.73, 0.80, 1.43, 1.97, 3.38, 1.61, 2.56, 0.67, 0.00, 1.13),
    log_kappa_hat = c(-0.32, -0.22, 0.36, 0.68, 1.22, 0.48, 0.94, -0.39, -70.67, 0.12),
    se_log_kappa_hat = c(0.38, 0.17, 0.22, 0.15, 0.14, 0.48, 0.51, 0.53, 41.44, 6.67),
    phase = c(
      "not phase III / early proof-of-concept",
      "phase II",
      "phase III",
      "phase III",
      "phase III",
      "phase IIb",
      "converted from phase III to phase IIb",
      "phase III",
      "not phase III / early proof-of-concept",
      "phase II"
    ),
    stringsAsFactors = FALSE
  )

  ## use_data expects it's data sets in the global env (which is why
  ## we do <<-)
  use_data(AS, transplant, colitis, crohn, asthma, overwrite = TRUE)
}


make_internal_ds <- function() {
  if (!file.exists("inst/sbc/calibration.rds")) {
    stop("Please create calibration run data first!")
  }

  calibration <- readRDS("inst/sbc/calibration.rds")

  calibration_meta <- calibration[c("S", "B", "git_hash", "created")]
  calibration_data <- calibration$data

  calibration_md5 <- strsplit(readLines("inst/sbc/calibration.md5"), ": +")
  vals <- sapply(calibration_md5, function(x) {
    x[[2]]
  })
  keys <- sapply(calibration_md5, function(x) {
    x[[1]]
  })
  names(vals) <- keys
  calibration_meta["MD5"] <- vals["MD5"]

  pkg_create_date <- Sys.time()
  pkg_sha <- "$Format:%h$"

  if (gsub("\\$", "", pkg_sha) == "Format:%h") {
    pkg_sha <- system("git rev-parse --short HEAD", intern = TRUE)
  }

  use_data(
    calibration_data,
    calibration_meta,
    pkg_create_date,
    pkg_sha,
    internal = TRUE,
    overwrite = TRUE
  )
}

## `use_data()` is a usethis function that devtools only re-exports, and it is
## the only thing needed here. Depending on usethis directly keeps this script
## runnable in the webR build's CI job without pulling in devtools.
library(usethis)

## cleanup first
if (file.exists("R/sysdata.rda")) {
  file.remove("R/sysdata.rda")
}

for (rda in dir("data", pattern = "*rda", full.names = TRUE)) {
  file.remove(rda)
}

##load_all()

make_ds()
make_internal_ds()
