#include /include/license.stan
#include /include/copyright_novartis.stan

// gMAP Stan Analysis
data {
  // number of input historical trials
  int<lower=1> H;

  // link function (1=normal, 2=binary, 3=poisson)
  int<lower=1,upper=3> link;

  // normal data, link=identity=1
  vector[H] y;
  vector[H] y_se;
  
  // binomial data, link=logit=2
  int<lower=0>    r[H];
  int<lower=1>  r_n[H];

  // count data, link=log=3
  int<lower=0> count[H];
  vector[H]    log_offset;

  // exchangeability cluster mapping
  int<lower=1> n_groups;
  int<lower=1,upper=n_groups> group_index[H];

  // tau prediction stratum
  int<lower=1,upper=n_groups> n_tau_strata;
  int<lower=1,upper=n_tau_strata> tau_strata_pred;
  // data item to tau stratum mapping
  int<lower=1,upper=n_tau_strata> tau_strata_index[H];

  // number of predictors
  int<lower=1> mX;
  // design matrix
  matrix[H,mX] X;
  
  // design matrix prediction (not used, only intercept prediction)
  //matrix[H,mX] Xpred;

  // priors
  matrix[mX,2] beta_prior;
  matrix[n_tau_strata,2] tau_prior;

  // model user choices
  int<lower=-1,upper=7> tau_prior_dist;
  int<lower= 0,upper=1> re_dist;
  real<lower=0> re_dist_t_df;

  // ncp parametrization?
  int<lower=0,upper=1> ncp;

  // guesses on the parameter location and scales
  vector[mX] beta_raw_guess[2];
  real       tau_raw_guess[2];

  // sample from prior predictive (do not add data to likelihood)
  int<lower=0,upper=1> prior_PD;
}
transformed data {
  vector[mX] beta_prior_stan[2];
  vector[n_tau_strata] tau_prior_stan[2];
  //matrix[n_groups, n_tau_strata] S;
  //matrix[H, n_groups] Z;
  matrix[H, mX] X_param;
  // group index to tau stratum mapping
  int<lower=1,upper=n_tau_strata> tau_strata_gindex[n_groups] = rep_array(tau_strata_pred, n_groups);

  for (i in 1:mX) {
    beta_prior_stan[1,i] = beta_prior[i,1];
    beta_prior_stan[2,i] = beta_prior[i,2];
  }
  
  for (i in 1:n_tau_strata) {
    tau_prior_stan[1,i] = tau_prior[i,1];
    tau_prior_stan[2,i] = tau_prior[i,2];
  }
  
  for (i in 1:H) {
    tau_strata_gindex[group_index[i]] = tau_strata_index[i];
  }
  
  /*
  // strata to group mapping
  S = rep_matrix(0, n_groups, n_tau_strata);
  for (i in 1:n_groups)
    S[i,tau_strata_index[i]] = 1.0;
  
  // groups to trial mapping
  Z = rep_matrix(0, H, n_groups);
  for (i in 1:H)
    Z[i,group_index[i]] = 1.0;
  */
  
  print("Stan gMAP analysis");

  if(link == 1)            print("likelihood:      Normal (identity link)");
  if(link == 2)            print("likelihood:      Binomial (logit link)");
  if(link == 3)            print("likelihood:      Poisson (log link)");
  
  if(tau_prior_dist == -1) print("tau distrib.:    Fixed");
  if(tau_prior_dist ==  0) print("tau distrib.:    HalfNormal");
  if(tau_prior_dist ==  1) print("tau distrib.:    TruncNormal");
  if(tau_prior_dist ==  2) print("tau distrib.:    Uniform");
  if(tau_prior_dist ==  3) print("tau distrib.:    Gamma");
  if(tau_prior_dist ==  4) print("tau distrib.:    InvGamma");
  if(tau_prior_dist ==  5) print("tau distrib.:    LogNormal");
  if(tau_prior_dist ==  6) print("tau distrib.:    TruncCauchy");
  if(tau_prior_dist ==  7) print("tau distrib.:    Exponential");

  if(re_dist == 0)         print("random effects:  Normal");
  if(re_dist == 1)         print("random effects:  Student-t, df = ", re_dist_t_df);

  if(ncp) {
    X_param = X;
    print("parametrization: Non-Centered");
  } else {
    print("parametrization: Centered");
    X_param = X;
    for (i in 1:H) {
      if(X_param[i,1] != 1)
        reject("Centered parametrization requires treatment contrast parametrization!");
      X_param[i,1] = 0;
    }
  }
  if(prior_PD)
    print("Info: Sampling from prior predictive distribution.");
}
parameters {
  vector[mX]           beta_raw;
  vector[n_tau_strata] tau_raw;
  vector[n_groups]     xi_eta;
}
transformed parameters {
  vector[H] theta;
  vector[n_groups] eta;
  vector[mX] beta;
  vector[n_tau_strata] tau;
  vector[n_groups] tau_group;

  beta = beta_raw_guess[1] + beta_raw_guess[2] .* beta_raw;

  // fixed tau distribution ignores raw_tau
  if(tau_prior_dist == -1)
    tau = tau_prior_stan[1];
  else
    tau = exp(tau_raw_guess[1] + tau_raw_guess[2] * tau_raw);

  tau_group = tau[tau_strata_gindex];
  
  if (ncp)  // NCP
    eta = xi_eta .* tau_group;
  else // CP places overall intercept into random effect
    eta = beta_raw_guess[1,1] + beta_raw_guess[2,1] * xi_eta;

  theta = X_param * beta + eta[group_index];
}
model {
  if (ncp) {
    // standardized random effect distribution (aka Matt trick)
    if(re_dist == 0) xi_eta ~ normal(                 0, 1);
    if(re_dist == 1) xi_eta ~ student_t(re_dist_t_df, 0, 1);
  } else {
    // random effect distribution
    if(re_dist == 0) xi_eta ~ normal( (beta[1] - beta_raw_guess[1,1])/beta_raw_guess[2,1], tau_group / beta_raw_guess[2,1]);
    if(re_dist == 1) xi_eta ~ student_t(re_dist_t_df, (beta[1] - beta_raw_guess[1,1])/beta_raw_guess[2,1], tau_group / beta_raw_guess[2,1]);
  }

  // assign priors to coefficients
  beta ~ normal(beta_prior_stan[1], beta_prior_stan[2]);

  // fixed (needs fake assignment)
  if(tau_prior_dist == -1) tau_raw ~ normal(                 0, 1);
  // half-normal
  if(tau_prior_dist ==  0) tau ~ normal(                     0, tau_prior_stan[2]);
  // truncated normal
  if(tau_prior_dist ==  1) tau ~ normal(     tau_prior_stan[1], tau_prior_stan[2]);
  if(tau_prior_dist ==  2) tau ~ uniform(    tau_prior_stan[1], tau_prior_stan[2]);
  if(tau_prior_dist ==  3) tau ~ gamma(      tau_prior_stan[1], tau_prior_stan[2]);
  if(tau_prior_dist ==  4) tau ~ inv_gamma(  tau_prior_stan[1], tau_prior_stan[2]);
  if(tau_prior_dist ==  5) tau ~ lognormal(  tau_prior_stan[1], tau_prior_stan[2]);
  if(tau_prior_dist ==  6) tau ~ cauchy(     tau_prior_stan[1], tau_prior_stan[2]);
  if(tau_prior_dist ==  7) tau ~ exponential(tau_prior_stan[1]);

  // add Jacobian adjustement due to shifting and transforming tau_raw
  if(tau_prior_dist != -1) target += tau_raw_guess[2] * tau_raw;
  
  // finally compute data-likelihood
  if(!prior_PD) {
    if(link == 1) y     ~ normal(              theta, y_se);
    if(link == 2) r     ~ binomial_logit(r_n,  theta);
    if(link == 3) count ~ poisson_log(log_offset + theta);
  }
}
generated quantities {
  real theta_pred;
  real theta_resp_pred;

  // make intercept only prediction
  if(re_dist == 0) theta_pred = normal_rng(                 beta[1], tau[tau_strata_pred]);
  if(re_dist == 1) theta_pred = student_t_rng(re_dist_t_df, beta[1], tau[tau_strata_pred]);

  if(link == 1) theta_resp_pred =           theta_pred;
  if(link == 2) theta_resp_pred = inv_logit(theta_pred);
  if(link == 3) theta_resp_pred = exp(      theta_pred);
}
