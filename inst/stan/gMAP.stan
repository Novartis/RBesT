#include /include/license.stan
#include /include/copyright_novartis.stan

// gMAP Stan Analysis
functions {
  /*
   * Orthonormal (Helmert) basis Q of the zero-sum subspace of R^J:
   *
   *   Q'Q = I_{J-1},   1'Q = 0,   QQ' = I_J - (1/J) 1 1'
   *
   * Used for the sum-to-zero (s2z) reparametrization which marginalizes the
   * data-free common shift of the group random effects. The dedicated
   * sum_to_zero_vector type requires Stan >= 2.36 while RBesT targets 2.32,
   * hence the explicit basis.
   *
   * Column k has k entries 1/sqrt(k(k+1)), one entry -k/sqrt(k(k+1)) and is
   * zero afterwards. For J = 1 the result is the 1 x 0 matrix.
   *
   * Building Q costs O(J^2) and applying it another O(J^2). J is the number
   * of historical trials (typically < 30), so the O(J) cumulative-sum Helmert
   * recursion is not worth the added complexity.
   */
  matrix zero_sum_basis(int J) {
    matrix[J, J - 1] Q = rep_matrix(0.0, J, J - 1);
    for (k in 1 : (J - 1)) {
      real s = inv_sqrt(k * (k + 1.0));
      for (i in 1 : k) {
        Q[i, k] = s;
      }
      Q[k + 1, k] = -k * s;
    }
    return Q;
  }
}
data {
  // number of input historical trials
  int<lower=1> H;
  
  // link function (1=normal, 2=binary, 3=poisson)
  int<lower=1, upper=3> link;
  
  // normal data, link=identity=1
  vector[H] y;
  vector[H] y_se;
  
  // binomial data, link=logit=2
  array[H] int<lower=0> r;
  array[H] int<lower=1> r_n;
  
  // count data, link=log=3
  array[H] int<lower=0> count;
  vector[H] log_offset;
  
  // exchangeability cluster mapping
  int<lower=1> n_groups;
  array[H] int<lower=1, upper=n_groups> group_index;
  
  // tau prediction stratum
  int<lower=1, upper=n_groups> n_tau_strata;
  int<lower=1, upper=n_tau_strata> tau_strata_pred;
  // data item to tau stratum mapping
  array[H] int<lower=1, upper=n_tau_strata> tau_strata_index;
  
  // number of predictors
  int<lower=1> mX;
  // design matrix
  matrix[H, mX] X;
  
  // does the model include an overall intercept? (first column of X is then
  // identically 1); required for the sum-to-zero reparametrization
  int<lower=0, upper=1> has_intercept;
  
  // user switch enabling the sum-to-zero reparametrization (option
  // RBesT.MC.s2z, default on); set to 0 to recover the legacy sampling scheme
  int<lower=0, upper=1> use_s2z;
  
  // design matrix prediction (not used, only intercept prediction)
  //matrix[H,mX] Xpred;
  
  // priors
  matrix[mX, 2] beta_prior;
  matrix[n_tau_strata, 2] tau_prior;
  
  // model user choices
  int<lower=-1, upper=7> tau_prior_dist;
  int<lower=0, upper=1> re_dist;
  real<lower=0> re_dist_t_df;
  
  // ncp parametrization?
  int<lower=0, upper=1> ncp;
  
  // guesses on the parameter location and scales
  array[2] vector[mX] beta_raw_guess;
  array[2] real tau_raw_guess;
  
  // sample from prior predictive (do not add data to likelihood)
  int<lower=0, upper=1> prior_PD;
}
transformed data {
  array[2] vector[mX] beta_prior_stan;
  array[2] vector[n_tau_strata] tau_prior_stan;
  //matrix[n_groups, n_tau_strata] S;
  //matrix[H, n_groups] Z;
  matrix[H, mX] X_param;
  // group index to tau stratum mapping
  array[n_groups] int<lower=1, upper=n_tau_strata> tau_strata_gindex = rep_array(tau_strata_pred,
                                                                    n_groups);
  
  /*
   * Sum-to-zero (s2z) reparametrization.
   *
   * The data see the group effects only through beta[1] + eps_j, so the common
   * shift mean(eps) is conditionally data-free and is marginalized out here.
   * This removes the beta[1]/mean(eps) ridge; the super-population beta[1] is
   * recovered exactly (not approximately) in transformed parameters.
   *
   * Only applicable when
   *   - the random effects are normal (a vector of iid Student-t values is not
   *     spherically symmetric, so the shift is not data-free),
   *   - there is a single tau stratum (with heterogeneous scales the exact
   *     decomposition is the precision-weighted one, which needs a basis
   *     rebuilt per gradient evaluation; see
   *     design/issue-s2z-multi-tau-strata.md), and
   *   - the model has an overall intercept which can absorb the shift.
   * Otherwise the legacy parametrization is used unchanged.
   *
   * use_s2z is a user escape hatch rather than a correctness condition: the
   * two parametrizations describe the same model, so switching it off changes
   * only the sampling geometry. It exists so that a user who hits pathological
   * behaviour in the new geometry can fall back without downgrading.
   */
  int s2z = use_s2z && (re_dist == 0) && (n_tau_strata == 1) && has_intercept;
  // number of sampled group coordinates: J-1 free subspace coordinates under
  // s2z, J group effects otherwise
  int n_re = s2z ? n_groups - 1 : n_groups;
  // dimension of the data-free recovery direction zeta (0 or 1)
  int n_rec = s2z ? 1 : 0;
  real inv_sqrt_J = inv_sqrt(n_groups);
  matrix[s2z ? n_groups : 0, s2z ? n_groups - 1 : 0] Q;
  
  for (i in 1 : mX) {
    beta_prior_stan[1, i] = beta_prior[i, 1];
    beta_prior_stan[2, i] = beta_prior[i, 2];
  }
  
  for (i in 1 : n_tau_strata) {
    tau_prior_stan[1, i] = tau_prior[i, 1];
    tau_prior_stan[2, i] = tau_prior[i, 2];
  }
  
  for (i in 1 : H) {
    tau_strata_gindex[group_index[i]] = tau_strata_index[i];
  }
  
  if (s2z) {
    /*
     * Absorbing the common shift into the intercept is only valid if the first
     * design column is *identically* one: a shift delta moves theta[h] by
     * delta * (1 - X[h,1]) otherwise, which is a different model. This holds
     * whenever has_intercept is set, since model.matrix() then emits an
     * all-ones (Intercept) column, but assert it so that a future change to
     * how X is built cannot break shift absorption silently.
     *
     * NOTE: this is a *different* invariant from the treatment-contrast guard
     * of the legacy centered parametrization below; do not conflate them.
     */
    for (i in 1 : H) {
      if (X[i, 1] != 1) {
        reject("s2z requires an all-ones intercept column!");
      }
    }
    // sd_alpha = hypot(s1, tau/sqrt(J)) must be strictly positive, which can
    // fail for s1 = 0 combined with a fixed tau = 0
    if (beta_prior_stan[2, 1] <= 0) {
      reject("s2z requires a strictly positive intercept prior sd!");
    }
    Q = zero_sum_basis(n_groups);
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
  
  if (link == 1) 
    print("likelihood:      Normal (identity link)");
  if (link == 2) 
    print("likelihood:      Binomial (logit link)");
  if (link == 3) 
    print("likelihood:      Poisson (log link)");
  
  if (tau_prior_dist == -1) 
    print("tau distrib.:    Fixed");
  if (tau_prior_dist == 0) 
    print("tau distrib.:    HalfNormal");
  if (tau_prior_dist == 1) 
    print("tau distrib.:    TruncNormal");
  if (tau_prior_dist == 2) 
    print("tau distrib.:    Uniform");
  if (tau_prior_dist == 3) 
    print("tau distrib.:    Gamma");
  if (tau_prior_dist == 4) 
    print("tau distrib.:    InvGamma");
  if (tau_prior_dist == 5) 
    print("tau distrib.:    LogNormal");
  if (tau_prior_dist == 6) 
    print("tau distrib.:    TruncCauchy");
  if (tau_prior_dist == 7) 
    print("tau distrib.:    Exponential");
  
  if (re_dist == 0) 
    print("random effects:  Normal");
  if (re_dist == 1) 
    print("random effects:  Student-t, df = ", re_dist_t_df);
  
  /*
   * X_param is LEGACY-ONLY: its centered branch zeroes the intercept column so
   * that the intercept can be folded into the group effects. The s2z path
   * builds theta from X directly and must NOT wire the intercept through
   * X_param.
   *
   * The treatment-contrast guard below still executes for s2z fits, but can
   * never trigger for them: the s2z precondition already rules out
   * X[i,1] != 1.
   * That is not a user-visible behaviour change: it fires only when
   * X[i,1] != 1, which the s2z precondition above already rules out.
   */
  if (ncp) {
    X_param = X;
    print("parametrization: Non-Centered");
  } else {
    print("parametrization: Centered");
    X_param = X;
    for (i in 1 : H) {
      if (X_param[i, 1] != 1) 
        reject("Centered parametrization requires treatment contrast parametrization!");
      X_param[i, 1] = 0;
    }
  }
  
  if (prior_PD) 
    print("Info: Sampling from prior predictive distribution.");
}
parameters {
  vector[mX] beta_raw;
  vector[n_tau_strata] tau_raw;
  // under s2z these are the J-1 free coordinates inside the zero-sum subspace,
  // otherwise the J group effects
  vector[n_re] xi_eta;
  // data-free recovery direction zeta ~ N(0,1); length 0 unless s2z
  vector[n_rec] xi_abar;
}
transformed parameters {
  vector[H] theta;
  vector[mX] beta;
  vector[n_tau_strata] tau;
  
  beta = beta_raw_guess[1] + beta_raw_guess[2] .* beta_raw;
  
  // fixed tau distribution ignores raw_tau
  if (tau_prior_dist == -1) 
    tau = tau_prior_stan[1];
  else 
    tau = exp(tau_raw_guess[1] + tau_raw_guess[2] * tau_raw);
  
  // expand random effect to groups in loop for performance reasons
  if (s2z) {
    /*
     * Sum-to-zero path. beta[1] holds the *sampled* intercept
     * alpha = beta[1] + mean(eps) until the very end of this block, where it
     * is overwritten with the recovered super-population intercept. This
     * sequencing is load-bearing: theta must be formed while beta[1] still
     * holds alpha.
     */
    real alpha = beta[1]; // local, NOT an output
    real m1 = beta_prior_stan[1, 1];
    real s1 = beta_prior_stan[2, 1];
    real sd_a = tau[1] * inv_sqrt_J; // sd of the common shift
    real sd_alpha = hypot(s1, sd_a); // widened intercept prior sd
    // r of the design doc, renamed since `r` is the binomial response data
    real r_shift = sd_a / sd_alpha; // in (0, 1)
    /*
     * Recovery of the common shift. With zeta = xi_abar[1] ~ N(0,1) entering
     * no other density statement, this reproduces
     *   mean(eps) | alpha, tau ~ N(r^2 (alpha - m1), (s1 r)^2)
     * exactly, so beta[1] is an exact draw from its conditional. Written via
     * hypot and the ratio r so that neither s1^2 nor tau^2/J materializes.
     */
    real abar = r_shift * (r_shift * (alpha - m1) + s1 * xi_abar[1]);
    // group effects inside the zero-sum subspace; marginal sd is
    // tau * sqrt(1 - 1/J), which is exactly the distribution of
    // eps - mean(eps). Do NOT rescale by (1 - 1/J)^-1/2.
    vector[n_groups] re = Q
                          * (ncp ? tau[1] * xi_eta
                             : beta_raw_guess[2, 1] * xi_eta);
    
    // note: X, not X_param -- the intercept column must be kept here
    for (h in 1 : H) {
      theta[h] = X[h] * beta + re[group_index[h]];
    }
    
    // super-population intercept; must never receive a `~` statement, see
    // the model block
    beta[1] = alpha - abar;
  } else if (ncp) {
    if (n_tau_strata == 1) {
      // most common case of just one stratum which simplifies things
      // and in ncp mode
      for (h in 1 : H) {
        theta[h] = X_param[h] * beta + xi_eta[group_index[h]] * tau[1];
      }
    } else {
      for (h in 1 : H) {
        theta[h] = X_param[h] * beta
                   + xi_eta[group_index[h]]
                     * tau[tau_strata_gindex[group_index[h]]];
      }
    }
  } else {
    for (h in 1 : H) {
      theta[h] = X_param[h] * beta + beta_raw_guess[1, 1]
                 + beta_raw_guess[2, 1] * xi_eta[group_index[h]];
    }
  }
}
model {
  if (s2z) {
    /*
     * alpha is recomputed from the raw parameter: beta[1] now holds the
     * RECOVERED super-population intercept, which must never appear in any
     * density statement. That is exactly what makes zeta = xi_abar[1] an
     * exact N(0,1) draw and beta[1] an exact draw from its conditional.
     * A test asserts the N(0,1) marginal of xi_abar to catch any regression
     * here (e.g. re-adding a vectorized `beta ~ normal(...)`).
     */
    real alpha = beta_raw_guess[1, 1] + beta_raw_guess[2, 1] * beta_raw[1];
    
    xi_abar ~ std_normal();
    
    // free subspace coordinates; centered and non-centered are the same
    // target, no Jacobian is needed since the density is stated on the
    // sampled variable itself
    if (ncp) {
      xi_eta ~ std_normal();
    } else {
      xi_eta ~ normal(0, tau[1] / beta_raw_guess[2, 1]);
    }
    
    // widened intercept prior, alpha | tau ~ N(m1, s1^2 + tau^2/J).
    // The `~` form drops only genuinely constant terms; the -log(sd) here
    // depends on tau and is therefore retained, which is what the tau
    // marginal needs. A log_prob test over a tau grid asserts this.
    alpha ~ normal(beta_prior_stan[1, 1],
                   hypot(beta_prior_stan[2, 1], tau[1] * inv_sqrt_J));
    
    // remaining coefficients keep their original priors
    if (mX > 1) {
      beta[2 : mX] ~ normal(beta_prior_stan[1][2 : mX],
                            beta_prior_stan[2][2 : mX]);
    }
  } else {
    if (ncp) {
      // standardized random effect distribution (aka Matt trick)
      if (re_dist == 0) 
        xi_eta ~ normal(0, 1);
      if (re_dist == 1) 
        xi_eta ~ student_t(re_dist_t_df, 0, 1);
    } else {
      // random effect distribution
      if (re_dist == 0) 
        xi_eta ~ normal((beta[1] - beta_raw_guess[1, 1])
                        / beta_raw_guess[2, 1],
                        tau[tau_strata_gindex] / beta_raw_guess[2, 1]);
      if (re_dist == 1) 
        xi_eta ~ student_t(re_dist_t_df,
                           (beta[1] - beta_raw_guess[1, 1])
                           / beta_raw_guess[2, 1],
                           tau[tau_strata_gindex] / beta_raw_guess[2, 1]);
    }
    
    // assign priors to coefficients
    beta ~ normal(beta_prior_stan[1], beta_prior_stan[2]);
  }
  
  // fixed (needs fake assignment)
  if (tau_prior_dist == -1) 
    tau_raw ~ normal(0, 1);
  // half-normal
  if (tau_prior_dist == 0) 
    tau ~ normal(0, tau_prior_stan[2]);
  // truncated normal
  if (tau_prior_dist == 1) 
    tau ~ normal(tau_prior_stan[1], tau_prior_stan[2]);
  if (tau_prior_dist == 2) 
    tau ~ uniform(tau_prior_stan[1], tau_prior_stan[2]);
  if (tau_prior_dist == 3) 
    tau ~ gamma(tau_prior_stan[1], tau_prior_stan[2]);
  if (tau_prior_dist == 4) 
    tau ~ inv_gamma(tau_prior_stan[1], tau_prior_stan[2]);
  if (tau_prior_dist == 5) 
    tau ~ lognormal(tau_prior_stan[1], tau_prior_stan[2]);
  if (tau_prior_dist == 6) 
    tau ~ cauchy(tau_prior_stan[1], tau_prior_stan[2]);
  if (tau_prior_dist == 7) 
    tau ~ exponential(tau_prior_stan[1]);
  
  // add Jacobian adjustement due to shifting and transforming tau_raw
  if (tau_prior_dist != -1) 
    target += tau_raw_guess[2] * tau_raw;
  
  // finally compute data-likelihood
  if (!prior_PD) {
    if (link == 1) 
      y ~ normal(theta, y_se);
    if (link == 2) 
      r ~ binomial_logit(r_n, theta);
    if (link == 3) 
      count ~ poisson_log(log_offset + theta);
  }
}
generated quantities {
  real theta_pred;
  real theta_resp_pred;
  
  // make intercept only prediction
  if (re_dist == 0) 
    theta_pred = normal_rng(beta[1], tau[tau_strata_pred]);
  if (re_dist == 1) 
    theta_pred = student_t_rng(re_dist_t_df, beta[1], tau[tau_strata_pred]);
  
  if (link == 1) 
    theta_resp_pred = theta_pred;
  if (link == 2) 
    theta_resp_pred = inv_logit(theta_pred);
  if (link == 3) 
    theta_resp_pred = exp(theta_pred);
}


