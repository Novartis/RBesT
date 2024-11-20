# create a fake data set with a covariate
trans_cov <- transform(transplant, country=cut(1:11, c(0,5,8,Inf), c("CH", "US", "DE")))
set.seed(34246)
map <- gMAP(cbind(r, n-r) ~ 1 + country | study,
            data=trans_cov,
            tau.dist="HalfNormal",
            tau.prior=1,
            # Note on priors: we make the overall intercept weakly-informative
            # and the regression coefficients must have tighter sd as these are
            # deviations in the default contrast parametrization
            beta.prior=rbind(c(0,2), c(0,1), c(0,1)),
            family=binomial,
            ## ensure fast example runtime
            thin=1, chains=1)

# posterior predictive distribution for each input data item (shrinkage estimates)
pred_cov <- predict(map)
pred_cov

# extract sample as matrix
samp <- as.matrix(pred_cov)

# predictive distribution for each input data item (if the input studies were new ones)
pred_cov_pred <- predict(map, trans_cov)
pred_cov_pred


# a summary function returns the results as matrix
summary(pred_cov)

# obtain a prediction for new data with specific covariates
pred_new <- predict(map, data.frame(country="CH", study=12))
pred_new
