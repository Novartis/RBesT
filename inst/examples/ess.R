# Conjugate Beta example
a <- 5
b <- 15
prior <- mixbeta(c(1, a, b))

ess(prior)
(a+b)

# Beta mixture example
bmix <- mixbeta(rob=c(0.2, 1, 1), inf=c(0.8, 10, 2))

ess(bmix, "elir")

ess(bmix, "moment")
# moments method is equivalent to
# first calculate moments
bmix_sum <- summary(bmix)
# then calculate a and b of a matching beta
ab_matched <- ms2beta(bmix_sum["mean"], bmix_sum["sd"])
# finally take the sum of a and b which are equivalent
# to number of responders/non-responders respectivley
round(sum(ab_matched))

ess(bmix, method="morita")

# Predictive consistency of elir
n_forward <- 1E1
bmixPred <- preddist(bmix, n=n_forward)
pred_samp <- rmix(bmixPred, 1E2)
# use more samples here for greater accuracy
# pred_samp <- rmix(bmixPred, 1E3) 
pred_ess <- sapply(pred_samp, function(r) ess(postmix(bmix, r=r, n=n_forward), "elir") )
ess(bmix, "elir")
mean(pred_ess) - n_forward

# Normal mixture example
nmix <- mixnorm(rob=c(0.5, 0, 2), inf=c(0.5, 3, 4), sigma=10)

ess(nmix, "elir")

ess(nmix, "moment")

## the reference scale determines the ESS
sigma(nmix) <- 20
ess(nmix)

# Gamma mixture example
gmix <- mixgamma(rob=c(0.3, 20, 4), inf=c(0.7, 50, 10))

ess(gmix) ## interpreted as appropriate for a Poisson likelihood (default)

likelihood(gmix) <- "exp"
ess(gmix) ## interpreted as appropriate for an exponential likelihood


