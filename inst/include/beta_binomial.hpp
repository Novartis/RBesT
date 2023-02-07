// [[Rcpp::depends(rstan)]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/arr.hpp>
#include <iostream>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export(rng = false, name="dBetaBinomial")]]
NumericVector dBetaBinomial(IntegerVector r, IntegerVector n, NumericVector alpha, NumericVector beta, bool log = false) {
  std::size_t size = std::max(std::max(r.size(), n.size()), std::max(alpha.size(), beta.size()));
  if(size > 1) {
    if(r.size() == 1) r = NumericVector(size, r[0]);
    if(n.size() == 1) n = NumericVector(size, n[0]);
    if(alpha.size() == 1) alpha = NumericVector(size, alpha[0]);
    if(beta.size() == 1) beta = NumericVector(size, beta[0]);
  }
  if (size != r.size() || size != n.size() || size != alpha.size() || size != beta.size())
    stop("Length of input arguments must match.");
  NumericVector value(size);
  for(std::size_t i=0; i!=size; ++i) {
    value[i] = stan::math::beta_binomial_lpmf(r[i], n[i], alpha[i], beta[i]);
  }
  if(!log)
    return(exp(value));
  return(value);
}

// [[Rcpp::export(rng = false, name=".pBetaBinomial")]]
NumericVector pBetaBinomial(IntegerVector r, IntegerVector n, NumericVector alpha, NumericVector beta, bool lower_tail = true, bool log_p = false) {
  std::size_t size = std::max(std::max(r.size(), n.size()), std::max(alpha.size(), beta.size()));
  if(size > 1) {
    if(r.size() == 1) r = NumericVector(size, r[0]);
    if(n.size() == 1) n = NumericVector(size, n[0]);
    if(alpha.size() == 1) alpha = NumericVector(size, alpha[0]);
    if(beta.size() == 1) beta = NumericVector(size, beta[0]);
  }
  if (size != r.size() || size != n.size() || size != alpha.size() || size != beta.size())
    stop("Length of input arguments must match.");
  NumericVector value(size);
  if(lower_tail) {
    for(std::size_t i=0; i!=size; ++i) {
      try {
        // for r==0 Stan does not return the correct values for the
        // boundaries
        value[i] = r[i] == 0 ?
          stan::math::beta_binomial_lpmf(r[i], n[i], alpha[i], beta[i]) :
          stan::math::beta_binomial_lcdf(r[i], n[i], alpha[i], beta[i]) ;
      } catch(const std::domain_error& e) {
        // this happens if the F32 does not converge. Then we get the
        // requested value brute force.
        std::vector<double> terms(r[i] + 1);
        for(int j=0; j <= r[i]; ++j)
          terms[j] = stan::math::beta_binomial_lpmf(j, n[i], alpha[i], beta[i]);
        value[i] = stan::math::log_sum_exp(terms);
      }
    }
  } else {
    for(std::size_t i=0; i!=size; ++i) {
      try {
        value[i] = r[i] == 0 ?
          stan::math::log_diff_exp(0, stan::math::beta_binomial_lpmf(0, n[i], alpha[i], beta[i])) :
          stan::math::beta_binomial_lccdf(r[i], n[i], alpha[i], beta[i]) ;
      } catch(const std::domain_error& e) {
        // this happens if the F32 does not converge. Then we get the
        // requested value brute force.
        std::vector<double> terms(n[i] - r[i]);
        for(int j=r[i] + 1; j <= n[i]; ++j)
          terms[j-r[i]-1] = stan::math::beta_binomial_lpmf(j, n[i], alpha[i], beta[i]);
        value[i] = stan::math::log_sum_exp(terms);
      }
    }
  }
  if(!log_p)
    return(exp(value));
  return(value);
}
