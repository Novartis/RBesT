
source_example <- function(example, env=parent.frame(), disable_plots=TRUE) {
    ex_source <- readLines(system.file("examples", example, package="RBesT", mustWork=TRUE))
    if(disable_plots) {
        ex_source <- grep("plot\\(", ex_source, value=TRUE, invert=TRUE)
    }
    suppressMessages(ex <- source(textConnection(ex_source),
                                  local=env, echo=FALSE, verbose=FALSE))
    invisible(ex)
}
