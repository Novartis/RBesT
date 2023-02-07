context("pos2S: 2-Sample Probability of Success")

## test the analytical OC function via brute force simulation
set.seed(12354)

## expect results to be 1% exact
eps <- 1e-2

prior1 <- mixnorm(c(0.3, -0.2, 2), c(0.7, 0, 50), sigma=1)
prior2 <- mixnorm(c(1.0, 0, 50), sigma=1)

##prior1 <- mixnorm(c(0.3, -0.2, 20), c(0.7, 0, 50), sigma=1)
##prior2 <- mixnorm(c(1.0, 0, 50), sigma=1)

N1 <- 30
N2 <- 40

## we test here that the PoS indeed averages over the predictive of an
## informative prior weighting the conditional power respectively.

s <- 1

pcrit <- 0.80
qcrit <- 0

N_samp <- 1E4

dec <- decision2S(pcrit, qcrit)
decU <- decision2S(pcrit, qcrit, lower.tail=FALSE)

test_pos2S <- function(prior1, prior2, ia_dist1, ia_dist2, n1, n2, dec, decU) {
    skip_on_cran()

    ## the PoS is the expected value of the condition power integrated
    ## over the interim density which is what we check here
    cpo_analytic <- oc2S(prior1, prior2, n1, n2, dec)
    pos_analytic <- pos2S(prior1, prior2, n1, n2, dec)
    samp1 <- rmix(ia_dist1, N_samp)
    samp2 <- rmix(ia_dist2, N_samp)
    pos_mc <- mean(cpo_analytic(samp1, samp2))
    ##print(pos_mc)
    ##print(pos_analytic(ia_dist1, ia_dist2))
    expect_true(all(abs(pos_mc - pos_analytic(ia_dist1, ia_dist2)) < eps))
    lower.tail <- attr(dec,"lower.tail")
    if(lower.tail) {
        ##cat("Also testing lower.tail=FALSE\n")
        test_pos2S(prior1, prior2, ia_dist1, ia_dist2, n1, n2, decU)
    }
}

ia1 <- postmix(prior1, m=0.2, se=s/sqrt(15))
ia2 <- postmix(prior2, m=0, se=s/sqrt(15))

test_that("Normal PoS 2 sample function matches MC integration of CPO",
          test_pos2S(prior1, prior2,
                     ia1, ia2,           
                     N1, N2,
                     dec, decU))

## also run a MC comparison
pos2S_normal_MC <- function(prior1, prior2, N1, N2, dtheta1, dtheta2, pcrit=0.975, qcrit=0) {
    skip_on_cran()

    mean_sd1 <- sigma(prior1) / sqrt(N1)
    mean_sd2 <- sigma(prior2) / sqrt(N2)

    mean_prior1 <- prior1
    sigma(mean_prior1) <- mean_sd1
    mean_prior2 <- prior2
    sigma(mean_prior2) <- mean_sd2

    pred_dtheta1 <- preddist(dtheta1, n=N1)##, sigma=mean_sd1)
    pred_dtheta2 <- preddist(dtheta2, n=N2)##, sigma=mean_sd1)
    
    ##mean_samp1 <- rnorm(Nsim, theta1, mean_sd1)
    ##mean_samp2 <- rnorm(Nsim, theta2, mean_sd2)
    mean_samp1 <- rmix(pred_dtheta1, N_samp)
    mean_samp2 <- rmix(pred_dtheta2, N_samp)

    dec <- rep(NA, N_samp)

    for(i in 1:N_samp) {
        post1 <- postmix(mean_prior1, m=mean_samp1[i], se=mean_sd1)
        post2 <- postmix(mean_prior2, m=mean_samp2[i], se=mean_sd2)
        dec[i] <- as.numeric(pmix(RBesT:::mixnormdiff(post1, post2), qcrit) > pcrit)
    }

    mean(dec)
}

test_that("Normal PoS 2 sample function matches MC integration",
          {
              pos_mc <- pos2S_normal_MC(prior1, prior2, N1, N2, ia1, ia2, pcrit=0.8, qcrit=0)
              pos_analytic <- pos2S(prior1, prior2, N1, N2, dec)
              expect_true(all(abs(pos_mc - pos_analytic(ia1, ia2)) < eps))
          })

beta_prior <- mixbeta(c(1, 1, 1))
beta_ia1 <- postmix(beta_prior, r=20, n=50)
beta_ia2 <- postmix(beta_prior, r=30, n=50)
test_that("Binomial PoS 2 sample function matches MC integration of CPO",
          test_pos2S(beta_prior, beta_prior,
                     beta_ia1, beta_ia2,           
                     N1, N2,
                     dec, decU))


gamma_prior <- mixgamma(c(1, 1, 1), param="mn")
alpha <- 0.05
dec_count  <- decision2S(1-alpha, 0, lower.tail=TRUE)
dec_countU  <- decision2S(1-alpha, 0, lower.tail=FALSE)
gamma_ia1 <- postmix(gamma_prior, m=0.7, n=60)
gamma_ia2 <- postmix(gamma_prior, m=1.2, n=60)
test_that("Poisson PoS 2 sample function matches MC integration of CPO",
          test_pos2S(gamma_prior, gamma_prior,
                     gamma_ia1, gamma_ia2,           
                     N1, N2,
                     dec_count, dec_countU))

