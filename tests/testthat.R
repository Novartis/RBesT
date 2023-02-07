if(Sys.getenv("NOT_CRAN") == "true") {
    library(testthat)
    library(RBesT)

    test_check("RBesT")
}
