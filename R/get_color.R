# internal utilities to work with bayesplot -------------------------------
#' @keywords internal
get_color <- function(levels) {
    color_code <- sapply(levels, function(lev)
        switch(lev, l=1, lh=2, m=3, mh=4, d=5, dh=6, lev)
                         )
    unname(unlist(bayesplot::color_scheme_get())[color_code])
}
