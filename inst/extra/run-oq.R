library(testthat)
library(RBesT)

cat("TEST RUN DATE:", date(), "\n")

cat("TESTING PACKAGE:\n")
print(packageDescription("RBesT"))

## enforce that all tests are run
Sys.setenv(NOT_CRAN="true")

cat("RUNNING PACKAGE TESTS:\n")
## run each section separatley to get subsequent numbering per section
## of the TAP reporter; execution order is in line with vignette steps
for(test in c("gMAP", "EM", "oc1S", "oc2S", "mixdist", "mixdiff", "preddist", "postmix", "utils", "pos1S", "pos2S")) {
    test_package("RBesT", filter=test, reporter="tap")
}

## finally run all tests once more, but with the stop reporter. This
## ensures that the last line of this script is only displayed if and
## only if all tests run successful
test_package("RBesT", reporter="stop")

cat("\n\nR SESSION INFO:\n")

print(sessionInfo())

cat("\nTEST FINISH DATE:", date(), "\n")
cat("\n\nALL TESTS SUCCESSFUL\n")

