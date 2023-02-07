#! /usr/bin/env Rscript

rm(list=ls())

## generate example data sets
make_ds <- function() {
    colitis <<- data.frame(study=c("Van_assche", "Feagan", "Rutgeerts-1", "Rutgeerts-2"),
                           n = c(56,63,121,123),
                           r = c(6,9,18,7)
                           )

    AS <<- data.frame(study=paste("Study", 1:8),
                      n=c(107,44,51,39,139,20,78,35),
                      r=c(23,12,19,9,39,6,9,10))

    transplant <<- data.frame(study=paste("Study", 1:11),
                              n=c( 33, 45, 74,103,140, 49, 83, 59, 22,109,213),
                              r=c(  6,  8, 17, 28, 26,  8, 22,  8,  6, 16, 53))

    crohn <<- dat <- data.frame(study=c("Gastr06","AIMed07","NEJM07","Gastr01a","APhTh04","Gastr01b"),
                                n=c(74, 166, 328, 20, 25, 58),
                                y=c(-51, -49, -36, -47, -90, -54))

    ## use_data expects it's data sets in the global env (which is why
    ## we do <<-)
    use_data(AS, transplant, colitis, crohn, overwrite=TRUE)
}


make_internal_ds <- function() {
    if(!file.exists("inst/sbc/calibration.rds"))
        stop("Please create calibration run data first!")

    calibration  <- readRDS("inst/sbc/calibration.rds")

    calibration_meta <- calibration[c("S", "B", "git_hash", "created")]
    calibration_data <- calibration$data

    calibration_md5  <- strsplit(readLines("inst/sbc/calibration.md5"), ": +")
    vals  <- sapply(calibration_md5, function(x) { x[[2]] } )
    keys  <- sapply(calibration_md5, function(x) { x[[1]] } )
    names(vals)  <- keys
    calibration_meta["MD5"]  <- vals["MD5"]

    pkg_create_date  <- Sys.time()

    use_data(
        calibration_data,
        calibration_meta,
        pkg_create_date,
        internal=TRUE, overwrite=TRUE)
}

library(devtools)

## cleanup first
if(file.exists("R/sysdata.rda"))
    file.remove("R/sysdata.rda")

for(rda in dir("data", pattern="*rda", full.names=TRUE))
    file.remove(rda)

##load_all()

make_ds()
make_internal_ds()
