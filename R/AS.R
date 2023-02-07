#' Ankylosing Spondylitis.
#'
#' Data set containing historical information for placebo for a phase
#' II trial of ankylosing spondylitis patients. The primary efficacy
#' endpoint was the percentage of patients with a 20% response
#' according to the Assessment of SpondyloArthritis international
#' Society criteria for improvement (ASAS20) at week 6.
#'
#' @format A data frame with 8 rows and 3 variables:
#' \describe{
#'   \item{study}{study}
#'   \item{n}{study size}
#'   \item{r}{number of events}
#' }
#' @references Baeten D. et. al, \emph{The Lancet}, 2013, (382), 9906, p 1705
#'
#' @template example-start
#' @examples
#' set.seed(34563)
#' map_AS <- gMAP(cbind(r, n-r) ~ 1 | study,
#'                family=binomial,
#'                data=AS,
#'                tau.dist="HalfNormal", tau.prior=1,
#'                beta.prior=2)
#' @template example-stop
"AS"
