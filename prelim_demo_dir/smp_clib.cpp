/**
 * This is the code used for the R -> C++ conversion
 * Given the ft function is not implemented,  
 */

#include "RcppArmadillo.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat n_id_matrix;
int global_n;


// [[Rcpp::export]]
RcppExport SEXP createIdentityMatrix(SEXP n_){
  double n = Rcpp::as<double>(n_);
  n_id_matrix = arma::eye<arma::mat>(n, n);
  
  global_n = n;
  return Rcpp::wrap(0);
}


// [[Rcpp::export]]
RcppExport SEXP firstpasFt( SEXP pft_, SEXP s_)
{
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat tmp1 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat tmp2 = arma::solve(n_id_matrix * tmp1, n_id_matrix);
  arma::mat ans = (1.0 / s) * pft * tmp1 * tmp2;
  
  return Rcpp::wrap(ans);
}


// [[Rcpp::export]]
RcppExport SEXP firstpasft( SEXP pft_)
{
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  
  arma::mat tmp1 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat tmp2 = arma::solve(n_id_matrix * tmp1, n_id_matrix);
  arma::mat ans = pft * tmp1 * tmp2;
  
  return Rcpp::wrap(ans);
}

// [[Rcpp::export]]
RcppExport SEXP transprobt( SEXP pft_, SEXP s_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat j = arma::ones(global_n, global_n);
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat ans = 1.0 / s * tmp2 * (n_id_matrix - n_id_matrix * (pft * j));

  //TODO: Look into minor offset on answer...  
  return Rcpp::wrap(ans);
}

// [[Rcpp::export]]
RcppExport SEXP smp_firstpass_cdf( SEXP pft_, SEXP s_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat tmp3 = arma::solve(n_id_matrix * tmp2, n_id_matrix);
  
  arma::mat ans = (1.0 / s) * pft * tmp2 * tmp3;
  return Rcpp::wrap(ans);
}

// [[Rcpp::export]]
RcppExport SEXP smp_firstpass_pdf( SEXP pft_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat tmp3 = arma::solve(n_id_matrix * tmp2, n_id_matrix);
  
  arma::mat ans = pft * tmp2 * tmp3;
  return Rcpp::wrap(ans);
}

// [[Rcpp::export]]
RcppExport SEXP smp_transprob_t( SEXP pft_, SEXP s_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat j = arma::ones(global_n, global_n);
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat ans = 1.0 / s * tmp2 * (n_id_matrix - n_id_matrix * (pft * j));
  
  //TODO: Look into minor offset on answer...  
  return Rcpp::wrap(ans);
}

// [[Rcpp::export]]
RcppExport SEXP smp_expected_visits_to_state( SEXP pft_, SEXP s_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat ans = (1 / s) * (tmp2 - n_id_matrix);
  
  //TODO: Look into minor offset on answer...  
  return Rcpp::wrap(ans);
}

// [[Rcpp::export]]
RcppExport SEXP smp_expected_time_in_state( SEXP pft_, SEXP s_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat j = arma::ones(global_n, global_n);
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat ans = (1.0 / (s * s)) * tmp2 * (n_id_matrix - n_id_matrix * (pft * j));
  
  //TODO: Look into minor offset on answer...  
  return Rcpp::wrap(ans);
}


// [[Rcpp::export]]
RcppExport SEXP smp_trans_to_state_0_times( SEXP pft_, SEXP s_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat j = arma::ones(global_n, global_n);
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat tmp3 = arma::solve(n_id_matrix * tmp2, n_id_matrix);
  arma::mat g = pft * tmp2 * tmp3;
  arma::mat ans = (1.0 / s) * (j - g);
  
  //TODO: Look into minor offset on answer...  
  return Rcpp::wrap(ans);
}


// [[Rcpp::export]]
RcppExport SEXP smp_trans_to_state_1_times( SEXP pft_, SEXP s_){
  arma::mat pft = Rcpp::as<arma::mat>(pft_);
  double s = Rcpp::as<double>(s_);
  
  arma::mat j = arma::ones(global_n, global_n);
  arma::mat tmp2 = arma::solve(n_id_matrix - pft, n_id_matrix);
  arma::mat tmp3 = arma::solve(n_id_matrix * tmp2, n_id_matrix);
  arma::mat g = pft * tmp2 * tmp3;
  arma::mat ans = (1.0 / s) * (j - g * (j * (n_id_matrix * g)));
  
  //TODO: Look into minor offset on answer...  
  return Rcpp::wrap(ans);
}



/*

 SMP.prob.trans.to.state.1.times <- function(s) {
 J <- matrix(1,ncol=n,nrow=n)
 tmp1 <- p*ft(s)  
 tmp2 <- solve(Id-tmp1)
 tmp3 <- solve(Id*tmp2)
 g <- tmp1%*%tmp2%*%tmp3
 return(1/s*(J-g*(J%*%(Id*g)))) 
 }
*/

