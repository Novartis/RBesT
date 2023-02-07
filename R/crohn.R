#' Crohn's disease.
#'
#' Data set containing historical information for placebo arm of
#' relevant studies for the treatment of Crohn's disease. The primary
#' outcome is change from baseline in Crohn's Disease Activity Index
#' (CDAI) over a duration of 6 weeks. Standard deviation of change
#' from baseline endpoint is approximately 88.
#'
#' @format A data frame with 4 rows and 3 variables:
#' \describe{
#'   \item{study}{study}
#'   \item{n}{study size}
#'   \item{y}{mean CDAI change}
#' }
#'
#' @references Hueber W. et. al, \emph{Gut}, 2012, 61(12):1693-1700
#'
#' @template example-start
#' @examples
#' set.seed(546346)
#' map_crohn <- gMAP(cbind(y, y.se) ~ 1 | study,
#'                   family=gaussian,
#'                   data=transform(crohn, y.se=88/sqrt(n)),
#'                   weights=n,
#'                   tau.dist="HalfNormal", tau.prior=44,
#'                   beta.prior=cbind(0,88))
#' @template example-stop
"crohn"
