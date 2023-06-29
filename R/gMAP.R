#' Meta-Analytic-Predictive Analysis for Generalized Linear Models
#'
#' Meta-Analytic-Predictive (MAP) analysis for generalized linear
#' models suitable for normal, binary, or Poisson data. Model
#' specification and overall syntax follows mainly
#' \code{\link[stats:glm]{glm}} conventions.
#'
#' @param formula the model formula describing the linear predictor
#' and encoding the grouping; see details
#' @param family the family of distributions defining the statistical
#' model (\code{binomial}, \code{gaussian}, or \code{poisson})
#' @param data optional data frame containing the variables of the
#' model. If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}.
#' @param weights optional weight vector; see details below.
#' @param offset offset term in statistical model used for Poisson
#' data
#' @param tau.strata sets the exchangability stratum per study. That
#' is, it is expected that each study belongs to a single
#' stratum. Default is to assign all studies to stratum 1. See section
#' differential heterogeniety below.
#' @param tau.strata.pred the index for the prediction stratum; default is 1.
#' @param tau.dist type of prior distribution for \code{tau};
#' supported priors are \code{HalfNormal} (default),
#' \code{TruncNormal}, \code{Uniform}, \code{Gamma}, \code{InvGamma},
#' \code{LogNormal}, \code{TruncCauchy}, \code{Exp} and \code{Fixed}.
#' @param tau.prior parameters of prior distribution for \code{tau};
#' see section prior specification below.
#' @param beta.prior mean and standard deviation for normal priors of
#' regression coefficients, see section prior specification below.
#' @param prior_PD logical to indicate if the prior predictive distribution should be sampled (no conditioning on the data). Defaults to \code{FALSE}.
#' @param REdist type of random effects distribution. \code{Normal} (default) or \code{t}.
#' @param t.df degrees of freedom if random-effects distribution is \code{t}.
#' @param contrasts an optional list; See \code{contrasts.arg} from
#' \code{\link[stats:model.matrix.default]{model.matrix.default}}.
#' @template args-sampling
#' @param digits number of displayed significant digits.
#' @param probs defines quantiles to be reported.
#' @param type sets reported scale (\code{response} (default) or \code{link}).
#' @param x,object \code{gMAP} analysis object created by \code{gMAP} function
#' @param ... optional arguments are ignored
#'
#' @details
#'
#' The meta-analytic-predictive (MAP) approach derives a prior from
#' historical data using a hierarchical model.  The statistical model is
#' formulated as a generalized linear mixed model for binary, normal
#' (with fixed \eqn{\sigma}) and Poisson endpoints:
#'
#' \deqn{y_{ih}|\theta_{ih} \sim f(y_{ih} | \theta_{ih})}{y_ih|\theta_ih ~ f(y_ih | \theta_ih)}
#'
#' Here, \eqn{i=1,\ldots,N} is the index for observations, and
#' \eqn{h=1,\ldots,H} is the index for the grouping (usually studies).
#' The model assumes the linear predictor for a transformed mean as
#'
#' \deqn{g(\theta_{ih}; x_{ih},\beta) = x_{ih} \, \beta + \epsilon_h}{g(\theta_ih; x_ih,\beta) = x_ih \beta + \epsilon_h}
#'
#' with \eqn{x_{ih}}{x_ih} being the row vector of \eqn{k} covariates for
#' observation \eqn{i}.  The variance component is assumed by default
#' normal
#'
#' \deqn{\epsilon_h \sim N(0,\tau^2), \qquad h=1,\ldots,H}{\epsilon_h ~ N(0,\tau^2), h=1,...,H}
#'
#' Lastly, the Bayesian implementation assumes independent normal
#' priors for the \eqn{k} regression coefficients and a prior for the
#' between-group standard deviation \eqn{\tau} (see \code{taud.dist}
#' for available distributions).
#'
#' The MAP prior will then be derived from the above model as the
#' conditional distribution of \eqn{\theta_{\star}}{\theta_*} given the
#' available data and the vector of covariates \eqn{x_{\star}}{x_*}
#' defining the overall intercept
#'
#' \deqn{\theta_{\star}| x_{\star},y .}{\theta_*| x_*,y .}
#'
#' A simple and common case arises for one observation (summary
#' statistic) per trial. For a normal endpoint, the model then simplifies
#' to the standard normal-normal hierarchical model. In the above
#' notation, \eqn{i=h=1,\ldots,H} and
#'
#' \deqn{y_h|\theta_h \sim N(\theta_h,s_h^2)}{y_h|\theta_h ~ N(\theta_h,s_h^2)}
#' \deqn{\theta_h = \mu + \epsilon_h}{\theta_h = \mu + \epsilon_h}
#' \deqn{\epsilon_h \sim N(0,\tau^2),}{\epsilon_h ~ N(0,\tau^2),}
#'
#' where the more common \eqn{\mu} is used for the only (intercept)
#' parameter \eqn{\beta_1}. Since there are no covariates, the MAP
#' prior is simply \eqn{Pr(\theta_{\star} |
#' y_1,\ldots,y_H)}{Pr(\theta_* | y_1,\ldots,y_H)}.
#'
#' The hierarchical model is a compromise between the two extreme
#' cases of full pooling (\eqn{\tau=0}, full borrowing, no
#' discounting) and no pooling (\eqn{\tau=\infty}, no borrowing,
#' stratification). The information content of the
#' historical data grows with H (number of historical data items)
#' indefinitely for full pooling whereas no information is
#' gained in a stratified analysis. For a fixed
#' \eqn{\tau}, the maximum effective sample
#' size of the MAP prior is \eqn{n_\infty} (\eqn{H\rightarrow
#' \infty}{H->\infty}), which for a normal endpoint with fixed
#' \eqn{\sigma} is
#'
#' \deqn{n_\infty = \left(\frac{\tau^2}{\sigma^2}\right)^{-1},}{n_\infty = (\tau^2/\sigma^2)^-1}
#'
#' (\emph{Neuenschwander et al., 2010}). Hence, the ratio
#' \eqn{\tau/\sigma} limits the amount of information a MAP prior is
#' equivalent to. This allows for a classification of \eqn{\tau}
#' values in relation to \eqn{\sigma}, which is crucial to define a
#' prior \eqn{P_\tau}. The following classification is useful in a
#' clinical trial setting:
#'
#' \tabular{lcc}{
#' Heterogeneity \tab \eqn{\tau/\sigma} \tab \eqn{n_\infty} \cr
#' small \tab 0.0625 \tab 256 \cr
#' moderate \tab 0.125 \tab 64 \cr
#' substantial \tab 0.25 \tab 16 \cr
#' large \tab 0.5 \tab 4 \cr
#' very large \tab 1.0 \tab 1
#' }
#'
#' The above formula for \eqn{n_\infty} assumes a known
#' \eqn{\tau}. This is unrealistic as the between-trial heterogeneity
#' parameter is often not well estimable, in particular if the number
#' of trials is small (H small). The above table helps to specify a
#' prior distribution for \eqn{\tau} appropriate for the given context
#' which defines the crucial parameter \eqn{\sigma}. For binary and
#' Poisson endpoints, normal approximations can be used to determine
#' \eqn{\sigma}. See examples below for concrete cases.
#'
#' The design matrix \eqn{X} is defined by the formula for the linear
#' predictor and is always of the form \code{response ~ predictor |
#' grouping}, which follows \code{\link[stats:glm]{glm}}
#' conventions. The syntax has been extended to include a
#' specification of the grouping (for example study) factor of the
#' data with a horizontal bar, \code{|}. The bar separates the
#' optionally specified grouping level, i.e. in the binary endpoint
#' case \code{cbind(r, n-r) ~ 1 | study}. By default it is assumed
#' that each row corresponds to an individual group (for which an
#' individual parameter is estimated). Specifics for the different
#' endpoints are:
#'
#' \describe{
#'
#' \item{normal}{\code{family=gaussian} assumes an identity link
#' function. The \code{response} should be given as matrix with two
#' columns with the first column being the observed mean value
#' \eqn{y_{ih}}{y_ih} and the second column the standard error
#' \eqn{se_{ih}}{se_ih} (of the mean). Additionally, it is recommended
#' to specify with the \code{weight} argument the number of units
#' which contributed to the (mean) measurement
#' \eqn{y_{ih}}{y_ih}. This information is used to estimate
#' \eqn{\sigma}.}
#'
#' \item{binary}{\code{family=binomial} assumes a logit link
#' function. The \code{response} must be given as two-column matrix
#' with number of responders \eqn{r} (first column) and non-responders
#' \eqn{n-r} (second column).}
#'
#' \item{Poisson}{\code{family=poisson} assumes a log link
#' function. The \code{response} is a vector of counts. The total
#' exposure times can be specified by an \code{offset}, which will be
#' linearly added to the linear predictor. The \code{offset} can be
#' given as part of the formula, \code{y ~ 1 + offset(log(exposure))}
#' or as the \code{offset} argument to \code{gMAP}. Note that the
#' exposure unit must be given as log-offset.}
#'
#' }
#'
#' @section Differential Discounting:
#'
#' The above model assumes the same between-group standard deviation
#' \eqn{\tau}, which implies that the data are equally relevant. This
#' assumption can be relaxed to more than one \eqn{\tau}. That is,
#'
#' \deqn{\epsilon_h \sim N(0,\tau_{s(h)}^2)}{\epsilon_h ~ N(0,\tau_s(h)^2)}
#'
#' where \eqn{s(h)} assignes group \eqn{h} to one of \eqn{S}
#' between-group heterogeneity strata.
#'
#' For example, in a situation with two randomized and four
#' observational studies, one may want to assume \eqn{\tau_1} (for
#' trials 1 and 2) and \eqn{\tau_2} (for trials 3-6) for the
#' between-trial standard deviations of the control means. More
#' heterogeneity (less relevance) for the observational studies can
#' then be expressed by appropriate priors for \eqn{\tau_1} and
#' \eqn{\tau_2}. In this case, \eqn{S=2} and the strata assignments
#' (see \code{tau.strata} argument) would be \eqn{s(1)=s(2)=1,
#' s(3)=\ldots=s(6)=2}.
#'
#' @section Prior Specification:
#'
#' The prior distribution for the regression coefficients \eqn{\beta}
#' is normal.
#'
#' \itemize{
#' \item If a single number is given, then this is used as the standard
#' deviation and the default mean of 0 is used.
#'
#' \item If a vector is given, it must be of the same length
#' as number of covariates defined and is used as standard
#' deviation.
#'
#' \item If a matrix with a single row is given, its first row will be
#' used as mean and the second row will be used as standard deviation
#' for all regression coefficients.
#'
#' \item Lastly, a two-column matrix (mean and standard deviation columns)
#' with as many columns as regression coefficients can be given.
#' }
#'
#' It is recommended to always specify a \code{beta.prior}. Per
#' default a mean of 0 is set. The standard deviation is set to 2 for
#' the binary case, to 100 * \code{sd(y)} for the normal case and to
#' \code{sd(log(y + 0.5 + offset))} for the Poisson case.
#'
#' For the between-trial heterogeniety \eqn{\tau} prior, a dispersion
#' parameter must always be given for each exchangeability
#' stratum. For the different \code{tau.prior} distributions, two
#' parameters are needed out of which one is set to a default value if
#' applicable:
#'
#' \tabular{lccl}{
#' Prior \tab \eqn{a} \tab \eqn{b} \tab default \cr
#' \code{HalfNormal}  \tab \eqn{\mu = 0} \tab  \eqn{\sigma} \tab \cr
#' \code{TruncNormal} \tab \eqn{\mu} \tab  \eqn{\sigma} \tab \eqn{\mu = 0} \cr
#' \code{Uniform}     \tab a \tab b \tab a = 0 \cr
#' \code{Gamma}       \tab \eqn{\alpha} \tab \eqn{\beta} \tab \cr
#' \code{InvGamma}    \tab \eqn{\alpha} \tab \eqn{\beta} \tab \cr
#' \code{LogNormal}   \tab \eqn{\mu_{\log}}{\mu_log} \tab \eqn{\sigma_{\log}}{\sigma_log} \tab \cr
#' \code{TruncCauchy} \tab \eqn{\mu} \tab \eqn{\sigma} \tab \eqn{\mu = 0} \cr
#' \code{Exp}         \tab \eqn{\beta} \tab 0 \tab \cr
#' \code{Fixed}       \tab a \tab 0 \tab \cr
#' }
#'
#' For a prior distribution with a default location parameter, a
#' vector of length equal to the number of exchangability strata can
#' be given. Otherwise, a two-column matrix with as many rows as
#' exchangability strata must be given, except for a single \eqn{\tau}
#' stratum, for which a vector of length two defines the parameters a
#' and b.
#'
#' @section Random seed: The MAP analysis is performed using
#' Markov-Chain-Monte-Carlo (MCMC) in \code{\link[rstan]{rstan}}. MCMC
#' is a stochastic algorithm. To obtain exactly reproducible results
#' you must use the \code{\link[base:set.seed]{set.seed}} function
#' before calling \code{gMAP}. See \code{\link[=RBesT-package]{RBesT}}
#' overview page for global options on setting further MCMC simulation
#' parameters.
#'
#' @return The function returns a S3 object of type \code{gMAP}. See
#' the methods section below for applicable functions to query the
#' object.
#'
#' @references Neuenschwander B, Capkun-Niggli G, Branson M,
#' Spiegelhalter DJ. Summarizing historical information on controls in
#' clinical trials. \emph{Clin Trials}. 2010; 7(1):5-18
#'
#' Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D,
#' Neuenschwander B.  Robust meta-analytic-predictive priors in
#' clinical trials with historical control information.
#' \emph{Biometrics} 2014;70(4):1023-1032.
#'
#' Weber S, Li Y, Seaman III J.W., Kakizume T, Schmidli H. Applying
#' Meta-Analytic Predictive Priors with the {R} {B}ayesian evidence
#' synthesis tools. \emph{JSS} 2021; 100(19):1-32
#'
#' @seealso \code{\link{plot.gMAP}}, \code{\link{forest_plot}}, \code{\link{automixfit}}, \code{\link{predict.gMAP}}
#'
#' @template example-start
#' @examples
#' # Binary data example 1
#'
#' # Mean response rate is ~0.25. For binary endpoints
#' # a conservative choice for tau is a HalfNormal(0,1) as long as
#' # the mean response rate is in the range of 0.2 to 0.8. For
#' # very small or large rates consider the n_infinity approach
#' # illustrated below.
#' # for exact reproducible results, the seed must be set
#' set.seed(34563)
#' map_AS <- gMAP(cbind(r, n-r) ~ 1 | study,
#'                family=binomial,
#'                data=AS,
#'                tau.dist="HalfNormal", tau.prior=1,
#'                beta.prior=2)
#' print(map_AS)
#'
#' # obtain numerical summaries
#' map_sum <- summary(map_AS)
#' print(map_sum)
#' names(map_sum)
#' # [1] "tau"        "beta"       "theta.pred" "theta"
#' map_sum$theta.pred
#'
#' \donttest{
#' # graphical model checks (returns list of ggplot2 plots)
#' map_checks <- plot(map_AS)
#' # forest plot with shrinkage estimates
#' map_checks$forest_model
#' # density of MAP prior on response scale
#' map_checks$densityThetaStar
#' # density of MAP prior on link scale
#' map_checks$densityThetaStarLink
#' }
#'
#' # obtain shrinkage estimates
#' fitted(map_AS)
#'
#' # regression coefficients
#' coef(map_AS)
#'
#' # finally fit MAP prior with parametric mixture
#' map_mix <- mixfit(map_AS, Nc=2)
#' plot(map_mix)$mix
#'
#' \donttest{
#' # optionally select number of components automatically via AIC
#' map_automix <- automixfit(map_AS)
#' plot(map_automix)$mix
#' }
#'
#' # Normal example 2, see normal vignette
#'
#' # Prior considerations
#'
#' # The general principle to derive a prior for tau can be based on the
#' # n_infinity concept as discussed in Neuenschwander et al., 2010.
#' # This assumes a normal approximation which applies for the colitis
#' # data set as:
#' p_bar <- mean(with(colitis, r/n))
#' s <- round(1/sqrt(p_bar * (1-p_bar)), 1)
#' # s is the approximate sampling standard deviation and a
#' # conservative prior is tau ~ HalfNormal(0,s/2)
#' tau_prior_sd <- s/2
#'
#' # Evaluate HalfNormal prior for tau
#' tau_cat <- c(pooling=0
#'             ,small=0.0625
#'             ,moderate=0.125
#'             ,substantial=0.25
#'             ,large=0.5
#'             ,veryLarge=1
#'             ,stratified=Inf)
#' # Interval probabilites (basically saying we are assuming
#' # heterogeniety to be smaller than very large)
#' diff(2*pnorm(tau_cat * s, 0, tau_prior_sd))
#' # Cumulative probabilities as 1-F
#' 1 - 2*(pnorm(tau_cat * s, 0, tau_prior_sd) - 0.5)
#'
#' @template example-stop
#' @export
gMAP <- function (formula,
                  family = gaussian,
                  data,
                  weights,
                  offset,
                  tau.strata,
                  tau.dist=c("HalfNormal","TruncNormal","Uniform","Gamma","InvGamma","LogNormal","TruncCauchy","Exp", "Fixed"),
                  tau.prior,
                  tau.strata.pred=1,
                  beta.prior,
                  prior_PD=FALSE,
                  REdist=c("normal","t"),
                  t.df=5,
                  contrasts=NULL,
                  iter=getOption("RBesT.MC.iter" , 6000),
                  warmup=getOption("RBesT.MC.warmup", 2000),
                  thin=getOption("RBesT.MC.thin", 4),
                  init=getOption("RBesT.MC.init", 1),
                  chains=getOption("RBesT.MC.chains", 4),
                  cores=getOption("mc.cores", 1L)
                  ) {
    call <- match.call()

    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }

    if (missing(data))
        data <- environment(formula)

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "offset", "weights", "tau.strata", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]

    f <- Formula::Formula(formula)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- f
    ##mf$weights <- weights
    mf <- eval(mf, parent.frame())
    mt <- terms(f, rhs=1)

    ## we only support a single LHS
    assert_that(length(f)[1] == 1)
    ## check that we have an overall intercept, otherwise the
    ## RBesT approach is not straightforward
    ##assert_that(attr(terms(mf, rhs=1), "intercept") == 1)
    has_intercept <- attr(terms(f, rhs=1), "intercept") == 1

    ## if the model response is a matrix, then the first column is
    ## interpreted as y and the second as n
    y <- model.response(mf)
    if(is.matrix(y)) {
        ## if y is a matrix, then we expect it to have 2 columns
        assert_that(ncol(y) == 2)
        ## y.aux is the standard deviation for the case of a normal
        ## y.aux is the number of non-responders for the binary case
        y.aux <- array(y[,2])
        y     <- array(y[,1])
    } else {
        y <- array(y)
        y.aux <- NULL
    }

    weights <- model.weights(mf)
    H <- NROW(y)

    ## todo: offset right now only taken care of for poisson
    ## regression, should be handled for all cases
    if(missing(offset))
        log_offset <- model.offset(model.part(f, data = mf, rhs = 1, terms = TRUE))
    else
        log_offset <- model.offset(mf)
    if(is.null(log_offset)) log_offset <- rep(0, H)
    log_offset <- array(log_offset)

    ## first define dummy data for all cases, which get overwritten for
    ## the case calculated below

    y_se <- array(rep(0,H))

    r   <- array(rep(0,H))
    r_n <- array(rep(1,H))

    count <- array(rep(0, H))

    if(length(f)[2] == 1) {
        ## no grouping has been given, then treat each data row as
        ## individual study
        group.factor <- 1:H
    } else {
        group.factor <- model.part(f, data = mf, rhs = 2)
        if(ncol(group.factor) != 1)
            stop("Grouping factor must be a single term (study).")
        group.factor <- group.factor[,1]
    }
    if (!is.factor(group.factor)) {
        group.factor <- factor(group.factor)
    }
    labels <- as.character(group.factor)
    group.index <- array(as.integer(group.factor))

    ## nubmer of groups coded by the factor
    n.groups <- nlevels(group.factor)
    ## number of groups actually observed in the data
    n.groups.obs <- length(unique(group.index))

    ## estimate of the reference scale, used in the normal case
    sigma_ref <- 0

    ## guess of the sd on the link scale, used to scale variables
    sigma_guess <- 1

    ## guessed tau
    tau_guess <- 1

    ## from here on everything should be defined which is needed for
    ## the given cases. Hence check if all inputs are given for the
    ## particular cases.
    if(family$family == "gaussian") {
        assert_that(family$link == "identity")
        if(is.null(y.aux)) {
            message("No standard error specified for normal data. Assuming standard error of 1 for all data items.")
            y.aux <- rep(1, H)
        }
        y_se <- y.aux
        y_n <- weights
        ## reference scale is the pooled variance estimate scaled by
        ## the total sample size
        if(!is.null(y_se) & !is.null(y_n)) {
            sigma_ref <- sqrt( sum(y_n) * 1/sum(1/y_se^2) )
            sigma_guess <- sigma_ref
        } else {
            sigma_ref <- NULL
            sigma_guess <- tau_guess
        }
        if(n.groups.obs > 1) {
            tau_guess <- max(sigma_guess/10, sd( tapply(y, group.index, mean) ) )
        }
    }
    if(family$family == "binomial") {
        assert_that(family$link == "logit")
        r   <- y
        nr  <- y.aux
        r_n <- y + y.aux
        lodds <- log((y + 0.5)/(nr + 0.5))
        ## p_bar would be 1 of r=n
        ##p_bar <- (sum(r) + 0.5)/ (sum(r_n) + 0.5)
        p_bar <- mean(inv_logit(lodds))
        sigma_guess <- 1/sqrt(p_bar * (1-p_bar))
        if(n.groups.obs > 1) {
            tau_guess <- max(sigma_guess/10, sd( tapply(lodds, group.index, mean) ) )
        }
    }
    if(family$family == "poisson") {
        assert_that(family$link == "log")
        count <- y
        sigma_guess <- 1/exp(mean(log(y + 0.5) - log_offset))
        if(n.groups.obs > 1) {
            tau_guess <- max(sigma_guess/10, sd( tapply(log(y + 0.5) - log_offset, group.index, mean) ) )
        } else {
            tau_guess <- sigma_guess
        }
    }

    ## create a unique label vector
    ulabels <- labels
    if(length(unique(ulabels)) != length(ulabels)) {
        group_column <- names(model.part(f, data = mf, rhs = 2))[1]
        data_factors <- setdiff(names(.getXlevels(mt, mf)), group_column)
        if(!is.null(model.extract(mf, "tau.strata"))) {
            group_column <- c(group_column, "(tau.strata)")
        }
        label_columns <- c(group_column, data_factors)
        if(length(label_columns) > 1)
            ulabels <- do.call(paste, c(mf[,label_columns], list(sep="/")))
        ## if now labels are still not unique, we label them sequentially
        if(length(unique(ulabels)) != length(ulabels)) {
            for(l in unique(labels)) {
                ind <- labels == l
                if(sum(ind) > 1) {
                    ulabels[ind] <- paste(ulabels[ind], seq(1,sum(ind)), sep="/")
                }
            }
        }
    }

    ## per group we must have an assignment to a tau stratum
    tau.strata.factor <- model.extract(mf, "tau.strata")
    if(is.null(tau.strata.factor)) {
        tau.strata.factor <- rep(1, H)
    } else {
        ## check that per group the tau stratum is unique
        for(g in levels(group.factor)) {
            gind <- group.factor == g
            if(length(unique(tau.strata.factor[gind])) != 1)
                stop("Found multiple tau strata defined for group", g, "!\nEach tau stratum must correspond to a unique group.")
        }
    }
    if(!is.factor(tau.strata.factor))
        tau.strata.factor <- factor(tau.strata.factor)
    tau.strata.index <- as.integer(tau.strata.factor)
    n.tau.strata <- max(nlevels(tau.strata.factor), tau.strata.pred)
    tau.strata.index <- array(as.integer(tau.strata.factor))

    ## setup design matrix
    X <- model.matrix(f, mf, rhs=1, contrasts.arg=contrasts)

    ## stratified estimates
    est_strat <- function(alpha) {
        z <- qnorm(1-alpha/2)
        theta_resp.strat <- switch(family$family,
                                   gaussian = cbind(y, y_se ,y - z * y_se,y + z * y_se),
                                   binomial = cbind(r/r_n, sqrt(r/(r_n) * (1-r/r_n) / r_n), BinaryExactCI(r, r_n, alpha, drop=FALSE) ),
                                   poisson = cbind(count/exp(log_offset), sqrt(count/exp(2*log_offset)), do.call(cbind, lapply(c(low=alpha/2, high=1-alpha/2), qgamma, shape=count + 0.5 * (count == 0), rate=exp(log_offset))))
                                   )
        dimnames(theta_resp.strat) <- list(ulabels,c("mean","se",paste0(c(100 * alpha/2, 100 *(1-alpha/2)), "%" )))
        theta_resp.strat
    }
    theta_resp.strat <- est_strat(0.05)
    theta.strat <- family$linkfun(theta_resp.strat)

    ## pooled estimates via glm fit
    fit.pooled <- if(family$family == "gaussian") {
        glm.fit(X, y, weights=as.vector(1/y_se^2), offset=log_offset, family=family)
    } else {
        glm.fit(X, model.response(mf), weights=as.vector(weights), offset=as.vector(log_offset), family=family)
    }
    theta.pooled <- fit.pooled$fitted.values
    theta_resp.pooled <- family$linkinv(fit.pooled$fitted.values)

    ## pooled fit should be replaced in the future with call to Stan
    ## optimizer:
    ##dataFixedL <- modifyList(dataL, list(tau_prior_dist=-1, tau_prior=matrix(1e-5, nrow=n.tau.strata, ncol=2)))
    ##fit_pooled <- optimizing(gMAP_model@stanmodel, data=dataFixedL)

    mX <- NCOL(X)

    if (missing(beta.prior)) {
        if(family$family == "gaussian") {
            beta.prior <- c(1e2*tau_guess)
        }
        if(family$family == "poisson") {
            beta.prior <- log(1e2) + tau_guess
        }
        if(family$family == "binomial") {
            beta.prior <- c(2)
        }
        message(paste("Assuming default prior dispersion for beta:", paste(beta.prior, collapse=", ")))
    }
    if(NCOL(beta.prior) == 1) {
        if(length(beta.prior) != 1) {
            assert_that(length(beta.prior) == mX)
        } else {
            beta.prior <- rep(beta.prior, mX)
        }
        beta.prior.location <- rep(0, mX)
        message("Assuming default prior location   for beta: ", paste(beta.prior.location, collapse=", "))
        if(mX > 1)
            warning("Check default prior location for intercept and regression coefficients!")
        beta.prior <- cbind(mean=beta.prior.location, sd=beta.prior)
    }
    if(!is.matrix(beta.prior))
        beta.prior <- matrix(beta.prior, mX, 2, byrow=TRUE, list(NULL, c("mean", "sd")))

    tau.dist <- match.arg(tau.dist)

    if(missing(tau.prior)) {
        ## abort execution if tau.prior not given
        stop("tau.prior must be set. This parameter is problem specific. Please consult documentation for details.")
        tau.prior <- switch(tau.dist,
                            Fixed = c(1, 0),
                            HalfNormal = c(0, 1),
                            TruncNormal = c(0, 1),
                            Uniform = c(0, 1),
                            Gamma = c(1, 1),
                            InvGamma = c(2, 1),
                            LogNormal = c(0, 1),
                            TruncCauchy = c(0, 1),
                            Exp = c(1, 0))
        tau.prior <- matrix(tau.prior, nrow=1, ncol=2)
    }

    assert_that(is.numeric(tau.prior))

    ## in case the user did not provide a matrix as prior.tau, try to
    ## guess if possible
    if(NCOL(tau.prior) == 1) {
        if(n.tau.strata == 1 & length(tau.prior) == 2)
            tau.prior <- matrix(tau.prior, nrow=1, ncol=2)
        if(n.tau.strata > 1 & !is.matrix(tau.prior) & tau.dist %in% c("LogNormal", "Gamma", "InvGamma")) {
            stop("Random effects dispersion distribution LogNormal, Gamma and InvGamma require matrix for tau.prior.")
        }
        if(!is.matrix(tau.prior)) {
            if(tau.dist %in% c("Fixed", "Exp")) {
                tau.prior <- cbind(tau.prior, 0)
            } else {
                tau.prior <- cbind(0, tau.prior)
            }
        }
    }

    if(NROW(tau.prior) < n.tau.strata) {
        stop("Multiple tau.strata defined, but tau.prior parameter not set for all strata.")
    }
    if(NROW(tau.prior) > n.tau.strata) {
        stop("More tau priors defined than tau.strata defined.")
    }

    ## code prior distribution
    tau_prior_dist <- switch(tau.dist
                            ,Fixed=-1
                            ,HalfNormal=0
                            ,TruncNormal=1
                            ,Uniform=2
                            ,Gamma=3
                            ,InvGamma=4
                            ,LogNormal=5
                            ,TruncCauchy=6
                            ,Exp=7
                             )

    if(tau.dist == "HalfNormal")
        assert_that(all(tau.prior[,1] == 0))

    REdist <- match.arg(REdist)

    re_dist <- ifelse(REdist == "normal", 0, 1)
    re_dist_t_df <- t.df

    link <- switch(family$family,
                   gaussian = 1,
                   binomial = 2,
                   poisson  = 3
                   )

    ## Model parametrization
    ## 0 = Use CP
    ## 1 = Use NCP
    ## 2 = Automatically detect which param to take
    ncp <- getOption("RBesT.MC.ncp",      1)

    assert_number(ncp, lower=0, upper=2)

    ## automatically detect if we have a sparse or rich data situation
    ## (very experimental detection, default is to use NCP)
    if(ncp == 2) {
        ncp <- 1
        ## we only have a tau_guess for H>1, then we set the CP
        ## parametrization whenever on average of the standard error is
        ## much smaller than the guessed tau in which case each group
        ## is estimated with high precision from the data in
        ## comparison to the between-group variation
        if( H>1 & sqrt(tau_guess^2 / max(theta_resp.strat[,"se"]^2)) > 20)
            ncp <- 0
    }

    ## calculate very roughly the scale of tau and mu; tau is
    ## calculated on the log-scale

    ## approximate maximal sample size we may get
    nInf <- 0.9 * (sigma_guess / tau_guess)^2

    ## consider here that the ss is chisq distributed; however, we
    ## need the square root transformed distribution; now this
    ## estimate will be over-confident since the between-group
    ## variation decreases the information we have. Hence we inflate
    ## the resulting sd according to the ratio of
    ## n_inf/n_(n.groups-1)/2 (eq. 11, Neuenschwander 2010)
    if(n.groups > 1) {
        ms <- square_root_gamma_stats((n.groups-1)/2, 2 * tau_guess^2 /(n.groups-1))
        ms[2] <- sqrt(1 + 2/(n.groups-1)) * ms[2]
        tau_raw_guess <- c(log(ms[1]) - log( sqrt( (1 + ms[2]^2/ms[1]^2) ) ),
                           sqrt(log(1 + ms[2]^2/ms[1]^2)))
    } else {
        tau_raw_guess <- c(log(tau_guess), 1)
    }

    beta_raw_guess <- rbind(mean=fit.pooled$coefficients,
                            sd=rep(sigma_guess/sqrt(nInf), mX))

    assert_logical(prior_PD, FALSE, len=1)

    fitData <- list("H", "X", "mX", "link",
                    "y", "y_se",
                    "r", "r_n",
                    "count", "log_offset",
                    "tau_prior_dist",
                    "re_dist", "re_dist_t_df",
                    "group.index", "n.groups",
                    "tau.strata.index", "n.tau.strata", "tau.strata.pred",
                    "beta.prior", "tau.prior",
                    "ncp", "tau_raw_guess", "beta_raw_guess",
                    "prior_PD")

    ## MODEL SETUP

    ##para <- c("beta", "tau", "theta", "theta_pred", "theta_resp_pred", "beta_raw", "tau_raw", "lp__")

    ## place data in a named list for safer passing it around in R
    dataL <- mget(unlist(fitData), envir=as.environment(-1))

    ## convert to Stan's 0/1 convention
    dataL$prior_PD <- as.integer(dataL$prior_PD)

    ## change variable naming conventions, replace forbidden "." to
    ## "_"
    names(dataL) <- gsub("\\.", "_", names(dataL))

    ## run model with Stan

    rescale  <- getOption("RBesT.MC.rescale",      TRUE)
    control_user <- getOption("RBesT.MC.control", list())
    control <- modifyList(list(adapt_delta=0.99, stepsize=0.01, max_treedepth=20), control_user)
    verbose  <- getOption("RBesT.verbose", FALSE)

    assert_flag(rescale)
    assert_number(init, lower=0, finite=TRUE)

    if(!rescale) {
        dataL$tau_raw_guess[2] <- 1
        dataL$beta_raw_guess[2,] <- 1
    }

    exclude_pars <- c("beta_raw", "tau_raw", "xi_eta")
    ## in absence of an overall intercept we drop the MAP posterior
    if(!has_intercept)
        exclude_pars <- c(exclude_pars, "theta_pred", "theta_resp_pred")

    if(verbose)
        exclude_pars <- c()

    ## MODEL RUN
    stan_msg <- capture.output(fit <- rstan::sampling(stanmodels$gMAP,
                                                      data=dataL,
                                                      ##pars=para,
                                                      warmup=warmup,
                                                      iter=iter,
                                                      chains=chains,
                                                      cores=cores,
                                                      thin=thin,
                                                      init=init,
                                                      refresh=0,
                                                      control=control,
                                                      algorithm = "NUTS",
                                                      open_progress=FALSE,
                                                      pars=exclude_pars,
                                                      include=FALSE,
                                                      save_warmup=TRUE
                                                      ))

    if(attributes(fit)$mode != 0)
        stop("Stan sampler did not run successfully!")

    ## only display Stan messages in verbose mode
    if(verbose) {
        cat(paste(c(stan_msg, ""), collapse="\n"))
    }

    ## MODEL FINISHED
    fit_sum <- rstan::summary(fit)$summary

    vars <- rownames(fit_sum)

    beta_ind <- grep("^beta\\[", vars)
    tau_ind  <- grep("^tau\\[", vars)
    lp_ind <- grep("^lp__", vars)

    beta <- fit_sum[beta_ind, "mean"]
    tau  <- fit_sum[tau_ind , "mean"]

    names(beta) <- colnames(X)
    names(tau)  <- paste0("tau", seq(n.tau.strata))

    Rhat.max <- max(fit_sum[,"Rhat"], na.rm=TRUE)

    if(Rhat.max > 1.1)
        warning("Maximal Rhat > 1.1. Consider increasing RBesT.MC.warmup MCMC parameter.")

    Neff.min <- min(fit_sum[c(beta_ind, tau_ind, lp_ind),"n_eff"], na.rm=TRUE)

    if(Neff.min < 1e3)
        message("Final MCMC sample equivalent to less than 1000 independent draws.\nPlease consider increasing the MCMC simulation size.")

    ## set internal RBesT thinning to 1
    thin <- 1

    ## finally include a check if the Stan NuTS sample had any
    ## divergence in the sampling phase, these are not supposed to
    ## happen and can often be avoided by increasing adapt_delta
    sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
    n_divergent <- sum(sapply(sampler_params, function(x) sum(x[,'divergent__'])) )
    if(n_divergent > 0) {
        warning(paste("In total", n_divergent, "divergent transitions occured during the sampling phase.\nPlease consider increasing adapt_delta closer to 1 with the following command prior to gMAP:\noptions(RBesT.MC.control=list(adapt_delta=0.999))"))
    }

    Out <- list(theta.strat=theta.strat,
                theta_resp.strat=theta_resp.strat,
                theta.pooled=theta.pooled,
                theta_resp.pooled=theta_resp.pooled,
                n.tau.strata=n.tau.strata,
                sigma_ref=sigma_ref,
                tau.strata.pred=tau.strata.pred,
                has_intercept=has_intercept,
                tau = tau,
                beta = beta,
                REdist = REdist,
                t.df=t.df,
                X = X,
                Rhat.max = Rhat.max,
                thin = thin,
                call = call,
                family = family,
                formula = f,
                model = mf,
                terms = mt,
                xlevels = .getXlevels(mt, mf),
                group.factor=group.factor,
                tau.strata.factor=tau.strata.factor,
                data = data,
                log_offset = log_offset,
                est_strat=est_strat,
                fit=fit,
                fit.data=dataL
                )

    structure(Out, class=c("gMAP"))
}


#' @describeIn gMAP displays a summary of the gMAP analysis.
#' @export
print.gMAP <- function(x, digits=3, probs=c(0.025, 0.5, 0.975), ...) {
    cat("Generalized Meta Analytic Predictive Prior Analysis\n")
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    cat("Exchangeability tau strata:", x$n.tau.strata,"\n")
    cat("Prediction tau stratum    :", x$tau.strata.pred,"\n")
    cat("Maximal Rhat              :", signif(x$Rhat.max, digits=digits),"\n")

    if(x$family$family == "gaussian" & !is.null(x$sigma_ref))
        cat("Estimated reference scale :", signif(x$sigma_ref, digits=digits), "\n")

    csum_tau <- rstan::summary(x$fit, probs=probs, pars=paste0("tau[", x$tau.strata.pred , "]"))$summary

    f <- colnames(csum_tau)[-match(c("se_mean", "n_eff", "Rhat"), colnames(csum_tau))]
    csum_tau <- csum_tau[,f]

    cat("\nBetween-trial heterogeneity of tau prediction stratum\n")
    print(signif(csum_tau, digits=digits))

    if(x$has_intercept) {
        csum_map <- rstan::summary(x$fit, probs=probs, pars="theta_resp_pred")$summary
        csum_map <- csum_map[,f]
        cat("\nMAP Prior MCMC sample\n")
        print(signif(csum_map, digits=digits))
    }

    div_trans <- sum(rstan::get_divergent_iterations(x$fit))
    num_sim <- length(rstan::get_divergent_iterations(x$fit))
    if (div_trans > 0) {
      warning(
        "The sampler detected ", div_trans, " out of ", num_sim, " transitions ending in a divergence after warmup.\n",
        "Increasing 'adapt_delta' closer to 1 may help to avoid these. Use for example: \n",
        paste0("options(RBesT.MC.control=list(adapt_delta=0.999))"),
        call.=FALSE
      )
    }
    Rhats <- bayesplot::rhat(x$fit)
    if (any(Rhats > 1.1, na.rm = TRUE)) {
      warning(
        "Parts of the model have not converged (some Rhats are > 1.1).\n",
        "Be careful when analysing the results! It is recommend to run\n",
        "more iterations and/or setting stronger priors.", call.=FALSE
      )
    }

    invisible(x)
}

#' @describeIn gMAP returns the quantiles of the posterior shrinkage
#' estimates for each data item used during the analysis of the given
#' \code{gMAP} object.
#' @export
fitted.gMAP <- function(object, type=c("response", "link"), probs = c(0.025, 0.5, 0.975), ...) {
    type <- match.arg(type)
    trans <- if(type == "response") object$family$linkinv else identity
    sim <- rstan::extract(object$fit, pars="theta")$theta
    res <- SimSum(trans(sim), probs=probs, margin=2)
    dimnames(res) <- list(rownames(object$theta_resp.strat), colnames(res))
    res
}

#' @describeIn gMAP returns the quantiles of the predictive
#' distribution. User can choose with \code{type} if the result is on
#' the response or the link scale.
#' @export
coef.gMAP <- function(object, probs = c(0.025, 0.5, 0.975), ...) {
    csum <- rstan::summary(object$fit, probs=probs, pars="beta")$summary
    f <- colnames(csum)[-match(c("se_mean", "n_eff", "Rhat"), colnames(csum))]
    csum <- subset(csum, select=f)
    rownames(csum) <- colnames(object$X)
    csum
}

#' @describeIn gMAP extracts the posterior sample of the model.
#' @method as.matrix gMAP
#' @export
as.matrix.gMAP <- function(x, ...) {
    as.matrix(x$fit, pars=c("lp__"), include=FALSE)
}

#' @method model.matrix gMAP
#' @export
model.matrix.gMAP <- function(object, ...) {
  return(model.matrix.default(object, object$data, contrasts.arg=object$contrast))
}


#' @describeIn gMAP returns the summaries of a gMAP.
#' analysis. Output is a \code{gMAPsummary} object, which is a list containing
#' \describe{
#' \item{\code{tau}}{posterior summary of the heterogeneity standard deviation}
#' \item{\code{beta}}{posterior summary of the regression coefficients}
#' \item{\code{theta.pred}}{summary of the predictive distribution (given in dependence on the \code{type} argument either on \code{response} or \code{link} scale)}
#' \item{\code{theta}}{posterior summary of the mean estimate (also depends on the \code{type} argument)}
#' }
#' @method summary gMAP
#' @export
summary.gMAP <- function(object, type=c("response", "link"), probs = c(0.025, 0.5, 0.975), ...) {
    call <- match.call()
    type <- match.arg(type)
    csum_beta <- rstan::summary(object$fit, probs=probs, pars=c("beta"))$summary
    csum_tau <- rstan::summary(object$fit, probs=probs,  pars=c("tau"))$summary
    if(object$has_intercept) {
        if(type == "response") {
            csum_pred <- rstan::summary(object$fit, probs=probs, pars=c("theta_resp_pred"))$summary
            csum_mean <- SimSum( object$family$linkinv( rstan::extract(object$fit, pars=c("beta[1]"))[[1]] ), probs=probs )
            rownames(csum_mean) <- "theta_resp"
        } else {
            csum_pred <- rstan::summary(object$fit, probs=probs, pars=c("theta_pred"))$summary
            csum_mean <- rstan::summary(object$fit, probs=probs, pars=c("beta[1]"))$summary
            rownames(csum_mean) <- "theta"
        }
    } else {
        csum_pred <- NULL
        csum_mean <- NULL
    }
    f <- colnames(csum_beta)[-match(c("se_mean", "n_eff", "Rhat"), colnames(csum_beta))]
    rownames(csum_beta) <- colnames(object$X)
    Out <- list(tau=subset(csum_tau, select=f), beta=subset(csum_beta, select=f))
    if(object$has_intercept) {
        Out <- c(Out, list(theta.pred=subset(csum_pred,select=f), theta=subset(csum_mean,select=f)))
    }
    structure(Out, class=c("gMAPsummary"), call=call)
}


#' @export
print.gMAPsummary <- function(x, digits=3, ...) {
    cat("Heterogeneity parameter tau per stratum:\n")
    print(signif(x$tau, digits=digits))
    cat("\nRegression coefficients:\n")
    print(signif(x$beta, digits=digits))
    if("theta" %in% names(x)) {
        cat("\nMean estimate MCMC sample:\n")
        print(signif(x$theta, digits=digits))
    }
    if("theta.pred" %in% names(x)) {
        cat("\nMAP Prior MCMC sample:\n")
        print(signif(x$theta.pred, digits=digits))
    }
    invisible(x)
}

## calculate the for a gamma distribution the mean and standard
## deviation of the square root transformed variable
## see square-root-of-gamma mathematica file
square_root_gamma_stats <- function(a, b) {
    m <- sqrt(b) * exp(lgamma(0.5 + a) - lgamma(a))
    v <- b*a - m^2
    c(mean=m, sd=sqrt(v))
}
