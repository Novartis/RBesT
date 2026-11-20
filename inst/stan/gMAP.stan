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

  /*
   * Partial-centering model parametrization for the group random effects,
   * conventional
   * (non-s2z) parametrization.
   *
   * For group j with fraction c_j the sampled coordinate xi_j carries the
   * location m_j and the scale sc_j returned here,
   *
   *   sc_j = (1 - c_j) + c_j * tau_j / g_j,
   *   m_j  = c_j * (beta1 - anchor - loc_j) / g_j,
   *
   * with g_j = group_scale_guess[j] and anchor = beta_raw_guess[1,1] the
   * per-group physical-scale numeraire and the location guess the centered
   * path already uses, and loc_j = group_location_guess[j] a per-group
   * posterior location guess (0 recovers the previous behavior exactly).
   * The physical group effect is
   * eps_j = tau_j * (xi_j - m_j) / sc_j, so stating xi_j ~ D(m_j, sc_j) is
   * *exactly* eps_j ~ D(0, tau_j) pushed through the parametrization: the
   * location-scale density already contains the Jacobian tau_j / sc_j, and
   * eps_j's marginal cancels m_j exactly, hence no separate adjustment and
   * no change to eps_j's D(0, tau_j) prior is introduced by either
   * numeraire. This holds for *any* m_j, sc_j (shared or per-group, zero or
   * not): loc_j only relocates where xi_j's probability mass sits, chosen
   * so the raw coordinate's natural range (where the numeric `init`
   * argument's random radius is centered) corresponds to a good starting
   * point relative to the quadrature-estimated posterior mean of eps_j.
   *
   * c = 0 gives m = 0, sc = 1 and eps = tau * xi -- the non-centered
   * model parametrization, *unconditionally* on loc_j (loc_j is itself
   * scaled by c_j, so it vanishes exactly at c = 0, preserving the shipped
   * non-centered identity even when group_location_guess is nonzero);
   * c = 1 gives m = (beta1 - anchor - loc) / g, sc = tau / g and
   * eps = g * xi - (beta1 - anchor) + loc -- the centered model
   * parametrization with an additional, exact group-specific shift.
   *
   * Interpolation is linear in the scale, matching brms and the s2z
   * parametrization.
   */
  vector partial_center_scale(vector center, vector tau_group, vector g) {
    vector[num_elements(center)] scale;
    for (j in 1 : num_elements(center)) {
      if (center[j] == 0)
        scale[j] = 1;
      else if (center[j] == 1)
        scale[j] = tau_group[j] / g[j];
      else
        scale[j] = (1 - center[j]) + center[j] * tau_group[j] / g[j];
    }
    return scale;
  }
  vector partial_center_loc(vector center, real location, real anchor,
                            vector group_location, vector g) {
    return center .* ((location - anchor - group_location) ./ g);
  }
  vector partial_center_effect(vector coordinate, vector center_location,
                               vector center, vector tau_group, vector g) {
    vector[num_elements(center)] effect;
    for (j in 1 : num_elements(center)) {
      if (center[j] == 0)
        effect[j] = tau_group[j] * (coordinate[j] - center_location[j]);
      else if (center[j] == 1)
        effect[j] = g[j] * (coordinate[j] - center_location[j]);
      else
        effect[j] = tau_group[j] * (coordinate[j] - center_location[j])
                    / ((1 - center[j]) + center[j] * tau_group[j] / g[j]);
    }
    return effect;
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
  
  // random-effect model parametrization: 0 = centered, 1 = non-centered,
  // 2 = partial
  int<lower=0, upper=2> re_param;
  
  /*
   * Per-group partial-centering fraction. This is consulted only when
   * re_param = 2; 0 is the non-centered model parametrization and 1 the
   * centered model parametrization.
   */
  vector<lower=0, upper=1>[n_groups] re_center;
  
  // guesses on the parameter location and scales; tau_raw_guess is one
  // location/scale pair per tau stratum, mirroring beta_raw_guess's existing
  // per-coefficient shape
  array[2] vector[mX] beta_raw_guess;
  array[2] vector[n_tau_strata] tau_raw_guess;

  // per-group physical-scale numeraire for the partial-centering affine map
  // (partial_center_scale() / partial_center_loc() / partial_center_effect()
  // below); consulted only when re_param = 2 and max(re_center) > 0. Any
  // positive value is an exact reparametrization (see the functions' doc
  // comment), so this is populated even when unused.
  vector<lower=0>[n_groups] group_scale_guess;
  // per-group posterior location guess for the same affine map; 0 recovers
  // the previous (location-free) behavior exactly. Scaled by re_center
  // inside partial_center_loc()/the s2z branch below, so it is inert
  // whenever re_center = 0 (any value is fine when unused).
  vector[n_groups] group_location_guess;
  
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
   * recovered exactly (not approximately) in generated quantities.
   *
   * Only applicable when
   *   - the random effects are normal (a vector of iid Student-t values is not
   *     spherically symmetric, so the shift is not data-free),
   *   - there is a single tau stratum (with heterogeneous scales the exact
   *     decomposition is the precision-weighted one, which needs a basis
   *     rebuilt per gradient evaluation and is therefore not
   *     implemented), and
   *   - the model has an overall intercept which can absorb the shift.
   * Otherwise the legacy parametrization is used unchanged.
   *
   * use_s2z is a user escape hatch rather than a correctness condition: the
   * two parametrizations describe the same model, so switching it off changes
   * only the sampling geometry. It exists so that a user who hits pathological
   * behaviour in the new geometry can fall back without downgrading.
   */
  int s2z = use_s2z && (re_dist == 0) && (n_tau_strata == 1) && has_intercept;
  int center_partial = (re_param == 2) && (max(re_center) > 0);
  int center_ncp = (re_param == 1) || ((re_param == 2) && !center_partial);
  // number of sampled group coordinates: J-1 free subspace coordinates under
  // s2z, J group effects otherwise
  int n_re = s2z ? n_groups - 1 : n_groups;
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
   * that the intercept can be folded into the group effects. Both the s2z path
   * and the conventional partial-centering parametrization build theta from X
   * directly and must NOT wire the intercept through X_param -- the
   * parametrization already carries it in the location m_j.
   *
   * The treatment-contrast guard below still executes for s2z fits, but can
   * never trigger for them: the s2z precondition already rules out
   * X[i,1] != 1.
   * That is not a user-visible behaviour change: it fires only when
   * X[i,1] != 1, which the s2z precondition above already rules out.
   */
  X_param = X;
  if (center_ncp || center_partial || s2z) {
    if (center_partial)
      print("parametrization: Partial");
    else if (center_ncp)
      print("parametrization: Non-Centered");
    else
      print("parametrization: Centered");
  } else {
    print("parametrization: Centered");
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
}
transformed parameters {
  vector[H] theta;
  /*
   * The SAMPLED coefficients. Under s2z beta_param[1] is the sampled intercept
   * alpha = beta[1] + mean(eps); otherwise beta_param is the coefficient
   * vector itself. The super-population `beta` is produced in generated
   * quantities and is the only one reported by default.
   */
  vector[mX] beta_param;
  vector[n_tau_strata] tau;
  // s2z partial-centering parametrization: the quadratic form and the
  // log-Jacobian that together replace the density statement on xi_eta.
  real s2z_quad = 0;
  real s2z_log_det = 0;
  
  beta_param = beta_raw_guess[1] + beta_raw_guess[2] .* beta_raw;
  
  // fixed tau distribution ignores raw_tau
  if (tau_prior_dist == -1) 
    tau = tau_prior_stan[1];
  else 
    tau = exp(tau_raw_guess[1] + tau_raw_guess[2] .* tau_raw);
  
  // expand random effect to groups in loop for performance reasons
  if (s2z) {
    /*
     * Group effects inside the zero-sum subspace; marginal sd is
     * tau * sqrt(1 - 1/J), which is exactly the distribution of
     * eps - mean(eps). Do NOT rescale by (1 - 1/J)^-1/2.
     *
     * Partial centering interpolates *linearly in the scale*, matching brms:
     *
     *   sc = (1 - c) + c * tau / g,   w = tau * (Q xi + c .* loc / g) ./ sc,
     *   re = w - mean(w)
     *
     * with g = group_scale_guess the per-group physical-scale numeraire the
     * centered path already uses (shared or heterogeneous, see
     * partial_center_scale()'s doc comment) and loc = group_location_guess a
     * per-group posterior location guess (0 recovers the previous,
     * location-free behavior exactly). c = 0 gives sc = 1 and
     * re = tau * Q xi, the non-centered
     * model parametrization, unconditionally on loc (the c .* loc term
     * vanishes exactly at c = 0);
     * c = 1 gives sc = tau/g and re = g .* (Q xi) + loc, reproducing the
     * existing centered model parametrization pointwise per group plus an
     * exact group-specific shift. Both existing model
     * parametrizations are therefore reproduced exactly.
     *
     * The re-projection `- mean(w)` matters for a *heterogeneous* c, g or
     * loc: dividing (or shifting) a zero-sum vector non-uniformly across
     * groups leaves the subspace, and re-projecting is what makes any
     * per-group loc exact rather than only a shared shift (which the
     * re-projection would cancel entirely, recovering the original
     * marginalized-common-shift model).
     */
    vector[n_groups] qxi = Q * xi_eta;
    vector[n_groups] re;
    
    if (!center_partial) {
      re = center_ncp ? tau[1] * qxi : beta_raw_guess[2, 1] * qxi;
    } else {
      vector[n_groups] tau_group = rep_vector(tau[1], n_groups);
      vector[n_groups] sc = partial_center_scale(re_center, tau_group,
                                                 group_scale_guess);
      vector[n_groups] center_location = -re_center
                                          .* (group_location_guess
                                              ./ group_scale_guess);
      vector[n_groups] w = partial_center_effect(
        qxi, center_location, re_center, tau_group,
        group_scale_guess
      );
      re = w - mean(w);
      /*
       * The induced map on the subspace coordinates is xi -> Q' D Q xi with
       * D = diag(tau ./ sc); the projection drops out because Q'(I - 11'/J) =
       * Q'. Its determinant is available in closed form,
       *
       *   det(Q' diag(d) Q) = (prod d_j) (sum 1/d_j) / J,
       *
       * for *any* positive diagonal d (homogeneous or not; this was
       * verified numerically for heterogeneous d as well), so no
       * factorization is needed and the parametrization stays
       * O(J). Writing
       * the target as the quadratic form in `re` rather than in `xi_eta`
       * absorbs the (J-1) log(tau) of the prior against the same term in the
       * determinant, which is why only mean(sc) survives below -- computed
       * directly from the (per-group) `sc` vector rather than the
       * shared-g-only algebraic shortcut this replaced.
       */
      s2z_quad = dot_self(re) / square(tau[1]);
      s2z_log_det = -sum(log(sc)) + log(mean(sc));
    }

    // note: X, not X_param -- the intercept column must be kept here, and
    // beta_param[1] is alpha, which is exactly what the data see
    for (h in 1 : H) {
      theta[h] = X[h] * beta_param + re[group_index[h]];
    }
  } else if (center_partial) {
    /*
     * Conventional partial-centering model parametrization, see
     * partial_center_scale() / partial_center_loc(). The physical group
     * effect is
     * tau_j * (xi_j - m_j) / sc_j, which reduces to tau_j * xi_j at c = 0 and
     * to g_j * xi_j - (beta_param[1] - anchor) + loc_j at c = 1, where
     * loc_j = group_location_guess[j] is an exact per-group posterior
     * location guess (0 recovers the previous behavior exactly).
     *
     * Note: X, not X_param. The model parametrization carries the intercept
     * through the location m_j, so the intercept column must be kept even when
     * the legacy centered model parametrization zeroed it in X_param. The
     * dedicated centered and non-centered paths below remain division-free.
     */
    vector[n_groups] tau_group = tau[tau_strata_gindex];
    real location = has_intercept ? beta_param[1] : 0;
    real anchor = has_intercept ? beta_raw_guess[1, 1] : 0;
    vector[n_groups] loc = partial_center_loc(re_center, location, anchor,
                                              group_location_guess,
                                              group_scale_guess);
    vector[n_groups] re = partial_center_effect(
      xi_eta, loc, re_center, tau_group, group_scale_guess
    );

    for (h in 1 : H) {
      theta[h] = X[h] * beta_param + re[group_index[h]];
    }
  } else if (center_ncp) {
    if (n_tau_strata == 1) {
      // most common case of just one stratum which simplifies things
      // and in ncp mode
      for (h in 1 : H) {
        theta[h] = X_param[h] * beta_param + xi_eta[group_index[h]] * tau[1];
      }
    } else {
      for (h in 1 : H) {
        theta[h] = X_param[h] * beta_param
                   + xi_eta[group_index[h]]
                     * tau[tau_strata_gindex[group_index[h]]];
      }
    }
  } else {
    for (h in 1 : H) {
      theta[h] = X_param[h] * beta_param + beta_raw_guess[1, 1]
                 + beta_raw_guess[2, 1] * xi_eta[group_index[h]];
    }
  }
}
model {
  if (s2z) {
    // free subspace coordinates; centered and non-centered are the same
    // target, no Jacobian is needed since the density is stated on the
    // sampled variable itself
    if (!center_partial) {
      if (center_ncp) {
        xi_eta ~ std_normal();
      } else {
        xi_eta ~ normal(0, tau[1] / beta_raw_guess[2, 1]);
      }
    } else {
      /*
       * Under partial centering the density is stated on the *physical* group
       * effects and the parametrization enters through the log-Jacobian.
       * c = 0 reduces to `xi_eta ~ std_normal()` and c = 1 to
       * `xi_eta ~ normal(0, tau/g)` exactly, so the two branches above are
       * special cases of this one rather than separate parametrizations.
       */
      target += -0.5 * s2z_quad + s2z_log_det;
    }
    
    // widened intercept prior, alpha | tau ~ N(m1, s1^2 + tau^2/J).
    // The `~` form drops only genuinely constant terms; the -log(sd) here
    // depends on tau and is therefore retained, which is what the tau
    // marginal needs. A log_prob test over a tau grid asserts this.
    beta_param[1] ~ normal(beta_prior_stan[1, 1],
                           hypot(beta_prior_stan[2, 1], tau[1] * inv_sqrt_J));
    
    // remaining coefficients keep their original priors
    if (mX > 1) {
      beta_param[2 : mX] ~ normal(beta_prior_stan[1][2 : mX],
                                  beta_prior_stan[2][2 : mX]);
    }
  } else {
    if (center_partial) {
      /*
       * Conventional partial-centering model parametrization. The density is
       * stated on the sampled coordinate itself with the exact location-scale
       * parameters of the parametrization, so it already contains the
       * Jacobian tau_j / sc_j of the map to the physical effects and no
       * `target +=` correction is needed.
       * c = 0 reduces to the non-centered branch and c = 1 to the centered one
       * below, which is why those remain plain fast paths rather than separate
       * models.
       */
      vector[n_groups] tau_group = tau[tau_strata_gindex];
      real location = has_intercept ? beta_param[1] : 0;
      real anchor = has_intercept ? beta_raw_guess[1, 1] : 0;
      vector[n_groups] sc = partial_center_scale(re_center, tau_group,
                                                 group_scale_guess);
      vector[n_groups] loc = partial_center_loc(re_center, location, anchor,
                                                group_location_guess,
                                                group_scale_guess);
      if (re_dist == 0)
        xi_eta ~ normal(loc, sc);
      if (re_dist == 1)
        xi_eta ~ student_t(re_dist_t_df, loc, sc);
    } else if (center_ncp) {
      // standardized random effect distribution (aka Matt trick)
      if (re_dist == 0) 
        xi_eta ~ normal(0, 1);
      if (re_dist == 1) 
        xi_eta ~ student_t(re_dist_t_df, 0, 1);
    } else {
      // random effect distribution
      real location = has_intercept ? beta_param[1] : 0;
      real anchor = has_intercept ? beta_raw_guess[1, 1] : 0;
      if (re_dist == 0) 
        xi_eta ~ normal((location - anchor) / beta_raw_guess[2, 1],
                        tau[tau_strata_gindex] / beta_raw_guess[2, 1]);
      if (re_dist == 1) 
        xi_eta ~ student_t(re_dist_t_df,
                           (location - anchor) / beta_raw_guess[2, 1],
                           tau[tau_strata_gindex] / beta_raw_guess[2, 1]);
    }

    // assign priors to coefficients
    beta_param ~ normal(beta_prior_stan[1], beta_prior_stan[2]);
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
    target += sum(tau_raw_guess[2] .* tau_raw);
  
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
  // super-population coefficients; under s2z the intercept is recovered below
  vector[mX] beta = beta_param;
  real theta_pred;
  real theta_resp_pred;
  
  if (s2z) {
    /*
     * Recovery of the marginalized common shift,
     *   mean(eps) | alpha, tau ~ N(r^2 (alpha - m1), (s1 r)^2),
     * taken from the RNG stream rather than from an extra sampled parameter.
     * Nothing about the fit -- metric, step size, treedepth, divergences,
     * thinning -- can touch it, so beta[1] is an exact iid draw from its
     * conditional given the retained draw. Written via hypot and the ratio r
     * so that neither s1^2 nor tau^2/J materializes.
     */
    real m1 = beta_prior_stan[1, 1];
    real s1 = beta_prior_stan[2, 1];
    real sd_a = tau[1] * inv_sqrt_J; // sd of the common shift
    real sd_alpha = hypot(s1, sd_a); // widened intercept prior sd
    // r of the design doc, renamed since `r` is the binomial response data
    real r_shift = sd_a / sd_alpha; // in (0, 1)
    real abar = normal_rng(r_shift * r_shift * (beta_param[1] - m1),
                           s1 * r_shift);
    beta[1] = beta_param[1] - abar;
  }

  // make intercept only prediction; drawn after the recovery so that
  // theta_pred stays jointly consistent with beta[1] within a draw
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
