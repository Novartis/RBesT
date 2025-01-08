#' Diagnostic plots for EM fits
#'
#' Produce diagnostic plots of EM fits returned from \code{\link{mixfit}}.
#'
#' @param x EM fit
#' @param size Optional argument passed to \code{ggplot2} routines
#' which control line thickness.
#' @param link Choice of an applied link function. Can take one of the
#' values \code{identity} (default), \code{logit} or \code{log}.
#' @param ... Ignored.
#'
#' Overlays the fitted mixture density with a histogram and a density
#' plot of the raw sample fitted. Applying a link function can be
#' beneficial, for example a \code{logit} (\code{log}) link for beta
#' (gamma) mixtures obtained from a Binomial (Poisson)
#' \code{\link{gMAP}} analysis.
#'
#' @template plot-help
#'
#' @return A list of \code{\link[ggplot2]{ggplot}} plots for
#' diagnostics of the EM run. Detailed EM diagnostic plots are
#' included only if the global option \code{RBesT.verbose} is set to
#' \code{TRUE}. These include plots of the parameters of each
#' component vs the iteration. The plot of the mixture density with a
#' histogram and a density of the fitted sample is always returned.
#'
#' @family EM
#'
#' @examples
#'
#' bmix <- mixbeta(rob = c(0.2, 1, 1), inf = c(0.8, 10, 2))
#' bsamp <- rmix(bmix, 1000)
#' bfit <- mixfit(bsamp, type = "beta", Nc = 2)
#' pl <- plot(bfit)
#'
#' print(pl$mixdens)
#' print(pl$mix)
#'
#' \donttest{
#' # a number of additional plots are generated in verbose mode
#' .user_option <- options(RBesT.verbose = TRUE)
#' pl_all <- plot(bfit)
#'
#' # recover previous user options
#' options(.user_option)
#'
#' names(pl_all)
#' # [1] "mixdist" "a"   "b"   "w"   "m"   "N"   "Lm"  "lN"  "Lw"  "lli" "mixdens" "mixecdf" "mix"
#' }
#'
#' @method plot EM
#' @export
plot.EM <- function(x, size = 1.25, link = c("identity", "logit", "log"), ...) {
  pl <- list()
  if (!is_mixmv(x)) {
    pl$mixdist <- plot.mix(x, size = size, ...)
  }
  ## in verbose mode we output EM fit diagnostics
  if (getOption("RBesT.verbose", FALSE)) {
    ## these NULL assignments make R check happy
    a <- b <- w <- s <- comp <- iteration <- NULL
    Nc <- ncol(x)
    pseq <- lapply(
      attr(x, "traceMix"),
      function(m) {
        class(m) <- "matrix"
        m <- as.data.frame(t(m))
        m$comp <- 1:Nc
        m
      }
    )
    names(pseq) <- 1:length(pseq) - 1
    Mw <- dplyr::bind_rows(pseq, .id = "iteration")
    Mw <- Mw[c(1, 5, 2, 3, 4)]
    Mw$iteration <- as.numeric(Mw$iteration)
    if ("EMbmm" %in% class(x)) {
      Mw <- within(Mw, {
        m <- a / (a + b)
        N <- a + b
        Lm <- logit(m)
        lN <- log(N)
      })
      if (Nc != 1) {
        Mw <- within(Mw, {
          Lw <- logit(w)
        })
      }
    }
    if ("EMnmm" %in% class(x)) {
      Mw <- within(Mw, {
        ls <- log(s)
      })
      if (Nc != 1) {
        Mw <- within(Mw, {
          Lw <- logit(w)
        })
      }
    }
    if ("EMgmm" %in% class(x)) {
      Mw <- within(Mw, {
        la <- log(a)
        lb <- log(b)
      })
      if (Nc != 1) {
        Mw <- within(Mw, {
          Lw <- logit(w)
        })
      }
    }
    pars <- names(Mw)[-c(1, 2)]
    Mw <- within(Mw, {
      Comp <- factor(comp)
    })
    LL <- data.frame(iteration = 0:max(Mw$iteration), lli = attr(x, "traceLli"))
    basePl <- ggplot(Mw, aes(x = .data$iteration, colour = .data$Comp)) +
      geom_line(size = size)
    for (p in pars) {
      pl[[p]] <- basePl + aes(y = .data$p)
    }
    pl$lli <- ggplot(subset(LL, iteration > 0), aes(x = .data$iteration, y = .data$lli)) +
      geom_line(size = size) +
      ylab("log-likelihood")
  }
  ## pl$mix <- plot.mix(x, comp=TRUE, samp=attr(x, "x"), ...)
  link <- match.arg(link)
  dlink(x) <- link_map[[link]]

  cols <- bayesplot::color_scheme_get(i = 1:6)

  if (!is_mixmv(x)) {
    ## univariate case
    samp <- data.frame(Sample = mixlink(x, as.vector(attr(x, "x"))))
    ## workaround a weird bug in ggplot which enlarges the interval
    interval <- quantile(samp$Sample, c(0.025, 0.975))
    n_fun <- 501
    max_span <- diff(range(samp))
    interval_span <- diff(interval)
    n_fun <- min(5E3, round(n_fun * max_span / interval_span))

    if (!is.dlink_identity(dlink(x))) {
      subtitle <- paste("Link:", dlink(x)$name)
    } else {
      subtitle <- NULL
    }

    pl$mixdens <- bayesplot::mcmc_dens(samp) + bayesplot::facet_text(FALSE) +
      stat_function(inherit.aes = FALSE, fun = dmix, args = list(mix = x), size = size, n = n_fun) +
      ggtitle("Parametric Mixture (black line) and Kernel Estimate of Sample Density", subtitle = subtitle) +
      bayesplot::xaxis_title(FALSE)

    pl$mixecdf <- ggplot(samp, aes(x = .data$Sample)) +
      stat_ecdf(geom = "area", size = 0, fill = cols$light) +
      stat_ecdf(geom = "step", size = size, colour = cols$mid) +
      stat_function(fun = pmix, args = list(mix = x), size = size, n = n_fun) +
      ggtitle("Estimated Cumulative Density from Parametric Mixture (black line) and Sample", subtitle = subtitle) +
      bayesplot::bayesplot_theme_get() +
      bayesplot::yaxis_title(FALSE) +
      bayesplot::xaxis_title(FALSE) +
      bayesplot::facet_text(FALSE)

    pl$mix <- bayesplot::mcmc_hist(samp, binwidth = diff(interval) / 50, freq = FALSE) + bayesplot::facet_text(FALSE) +
      stat_function(inherit.aes = FALSE, fun = dmix, args = list(mix = x), size = size, n = n_fun) +
      ggtitle("Parametric Mixture Density (black line) and Histogram of Sample", subtitle = subtitle) +
      bayesplot::xaxis_title(FALSE)
  } else if (inherits(x, "mvnormMix")) {
    var1 <- var2 <- NULL
    ## multivariate case: only support mvnorm for now (the only
    ## one supported as of Aug 2023). Plot the pair-wise marginal
    ## densities, which are the pair-wise marginal mixtures.
    message("Diagnostic plots for mixture multivariate normal densities are experimental.\nPlease note that these are subject to changes in future releases.")
    samp <- attr(x, "x")
    p <- ncol(samp)
    dim_labels <- mvnorm_dim_labels(x[-1, 1])
    if (is.null(dim_labels)) {
      dim_labels <- paste0("Dimension ", 1:p)
    }
    breaks <- 60
    nbins <- 20
    pairs <- subset(expand.grid(var1 = 1:p, var2 = 1:p), var1 >= var2)
    layout <- matrix(NA, nrow = p, ncol = p)
    pl_pairs <- list()
    pl_pairs_compact <- list()
    for (i in seq_len(nrow(pairs))) {
      v1 <- pairs$var1[i]
      v2 <- pairs$var2[i]
      label <- paste0("mixpair[", v2, ",", v1, "]")
      mix_sub <- mvnorm_extract_dim(x, unique(c(v2, v1)))
      layout[v1, v2] <- i
      if (v1 == v2) {
        interval <- quantile(samp[, v1], c(0.025, 0.975))
        pl_pairs[[label]] <- bayesplot::mcmc_hist(samp[, v1, drop = FALSE], binwidth = diff(interval) / 50, freq = FALSE) +
          bayesplot::facet_text(FALSE) +
          stat_function(inherit.aes = FALSE, fun = function(mix, x) dmix(mix, matrix(x, nrow = length(x))), args = list(mix = mix_sub)) +
          ylab(dim_labels[v2]) +
          xlab(dim_labels[v2])
      } else {
        data_ranges <- apply(samp[, c(v2, v1), drop = FALSE], 2, range)
        colnames(data_ranges) <- c("x", "y")
        data_grid <- expand.grid(apply(data_ranges, 2, function(r) seq(r[1], r[2], length = breaks), simplify = FALSE))
        data_grid$z <- dmix(mix_sub, as.matrix(data_grid), log = TRUE)
        pl_pairs[[label]] <- bayesplot::mcmc_scatter(samp[, c(v2, v1), drop = FALSE], alpha = 0.1) +
          bayesplot::facet_text(FALSE) +
          geom_contour(aes(z = .data$z), data = data_grid, bins = nbins, colour = "black") +
          xlab(dim_labels[v2]) +
          ylab(dim_labels[v1])
      }
      pl_pairs_compact[[label]] <- pl_pairs[[label]]
      if (v1 != p) {
        pl_pairs_compact[[label]] <- pl_pairs_compact[[label]] + bayesplot::xaxis_title(FALSE) + bayesplot::xaxis_ticks(FALSE) + bayesplot::xaxis_text(FALSE)
      }
      if (v2 != 1) {
        pl_pairs_compact[[label]] <- pl_pairs_compact[[label]] + bayesplot::yaxis_title(FALSE) + bayesplot::yaxis_ticks(FALSE) + bayesplot::yaxis_text(FALSE)
      }
    }
    pl$mixpairs <- bayesplot::bayesplot_grid(plots = pl_pairs_compact, grid_args = list(nrow = p, ncol = p, layout_matrix = layout))
    pl$mixpairs$bayesplots <- pl_pairs
  }

  pl
}
