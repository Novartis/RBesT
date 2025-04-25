#' Predictions from gMAP analyses
#'
#' @name predict.gMAP
#'
#' @description
#' Produces a sample of the predictive distribution.
#'
#' @param x,object gMAP analysis object for which predictions are performed
#' @param newdata data.frame which must contain the same columns as
#' input into the gMAP analysis. If left out, then a posterior prediction for
#' the fitted data entries from the gMAP object is performed (shrinkage estimates).
#' @param probs defines quantiles to be reported.
#' @param type sets reported scale (`response` (default) or `link`).
#' @param na.action how to handle missings.
#' @param thin thinning applied is derived from the `gMAP` object.
#' @param digits number of displayed significant digits.
#' @param ... ignored.
#'
#' @details Predictions are made using the \eqn{\tau} prediction
#' stratum of the gMAP object. For details on the syntax, please refer
#' to [predict.glm()] and the example below.
#'
#' @seealso [gMAP()], [predict.glm()]
#'
#' @template example-start
#' @example inst/examples/predict_gMAP.R
#' @template example-stop
#' @rdname predict.gMAP
#' @method predict gMAP
#' @export
predict.gMAP <- function(
  object,
  newdata,
  type = c("response", "link"),
  probs = c(0.025, 0.5, 0.975),
  na.action = na.pass,
  thin,
  ...
) {
  f <- object$formula
  mf <- object$model
  tt <- terms(f, data = mf, lhs = 1, rhs = 1)
  type <- match.arg(type)
  if (missing(newdata)) {
    posterior_predict <- TRUE
    X <- model.matrix(f, mf, rhs = 1)
    log_offset <- object$log_offset
    group.factor <- model.part(f, data = mf, rhs = 2)
  } else {
    posterior_predict <- FALSE
    Terms <- delete.response(tt)
    ## replace model frame with newdata context
    m <- model.frame(
      Terms,
      newdata,
      na.action = na.action,
      xlev = .getXlevels(tt, mf)
    )
    if (!is.null(cl <- attr(Terms, "dataClasses"))) {
      .checkMFClasses(cl, m)
    }
    X <- model.matrix(
      Terms,
      m,
      contrasts.arg = attr(model.matrix(f, mf, rhs = 1), "contrasts")
    )
    log_offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) {
      for (i in off.num) {
        log_offset <- log_offset +
          eval(
            attr(
              tt,
              "variables"
            )[[i + 1]],
            newdata
          )
      }
    }
    if (!is.null(object$call$offset)) {
      log_offset <- log_offset + eval(object$call$offset, newdata)
    }

    group.factor <- model.part(f, data = newdata, rhs = 2)
  }

  if (ncol(group.factor) != 1) {
    stop("Grouping factor must be a single term (study).")
  }
  group.factor <- group.factor[, 1]

  if (!is.factor(group.factor)) {
    group.factor <- factor(group.factor)
  }
  labels <- as.character(group.factor)
  group.index <- array(as.integer(group.factor))

  ## nubmer of groups coded by the factor
  n.groups <- nlevels(group.factor)
  ## number of groups actually observed in the data
  n.groups.obs <- length(unique(group.index))

  if (missing(thin)) {
    thin <- object$thin
  }

  beta <- rstan::extract(
    object$fit,
    inc_warmup = FALSE,
    permuted = FALSE,
    pars = "beta"
  )
  n.pred <- nrow(X)
  n.iter <- dim(beta)[1]
  n.chains <- dim(beta)[2]

  if (posterior_predict) {
    pred <- aperm(
      rstan::extract(
        object$fit,
        inc_warmup = FALSE,
        permuted = FALSE,
        pars = "theta"
      ),
      c(3, 1, 2)
    )
  } else {
    pred <- apply(beta, c(1, 2), function(x) X %*% x)
    if (n.pred == 1) {
      pred <- array(pred, dim = c(1, dim(pred)))
    }
  }

  sub_ind <- seq(1, n.iter, by = thin)

  pred <- t(matrix(pred[, sub_ind, ], nrow = n.pred))

  if (!posterior_predict) {
    ## in case we make a prediction unconditional on the fitted
    ## data, we draw random effects here (one for each study per
    ## iteration)

    S <- nrow(pred)
    ## sample random effects for as many groups defined, which can
    ## be more than the ones in the data set, since we sample for
    ## all defined factor levels
    tau <- as.vector(rstan::extract(
      object$fit,
      inc_warmup = FALSE,
      permuted = FALSE,
      pars = paste0("tau[", object$tau.strata.pred, "]")
    )[sub_ind, , ])
    if (object$REdist == "normal") {
      re <- tau * matrix(rnorm(n.groups * S, 0, 1), nrow = S)
    }
    if (object$REdist == "t") {
      re <- tau * matrix(rt(n.groups * S, df = object$t.df), nrow = S)
    }

    ## ... and add it to predictions
    pred <- pred + re[, group.index]
  }

  if (type == "response") {
    pred <- object$family$linkinv(pred)
  }

  predNames <- NULL
  if (!is.null(rownames(X))) {
    predNames <- rownames(X)
  }
  dimnames(pred) <- list(NULL, predNames)

  stat <- SimSum(pred, probs = probs, margin = 2)
  attr(pred, "summary") <- stat
  attr(pred, "type") <- type
  attr(pred, "family") <- object$family
  attr(pred, "sigma_ref") <- object$sigma_ref
  invisible(structure(pred, class = c("gMAPpred")))
}

#' @rdname predict.gMAP
#' @method print gMAPpred
#' @export
print.gMAPpred <- function(x, digits = 3, ...) {
  cat("Meta-Analytic-Predictive Prior Predictions\n")
  cat("Scale:", attr(x, "type"), "\n")
  cat("\n")
  cat("Summary:\n")
  print(signif(attr(x, "summary"), digits = digits))
}

#' @rdname predict.gMAP
#' @method summary gMAPpred
#' @export
summary.gMAPpred <- function(object, ...) {
  attr(object, "summary")
}

#' @rdname predict.gMAP
#' @method as.matrix gMAPpred
#' @export
as.matrix.gMAPpred <- function(x, ...) {
  class(x) <- "matrix"
  x
}
