#include "RcppArmadillo.h"

{
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat ans = X * Y;
}

RcppExport SEXP rArmSub( SEXP X_, SEXP Y_ )
{
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat ans = X - Y;
}

RcppExport SEXP rArmAdd( SEXP X_, SEXP Y_ )
{
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat ans = X + Y;
}

RcppExport SEXP rArmDiv( SEXP X_, SEXP Y_ )
{
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat ans = X/Y;
}

RcppExport SEXP rArmInv( SEXP X_ )
{
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat ans = arma::inv(arma::sympd(X));
}

RcppExport SEXP rArmSolve( SEXP X_, SEXP Y_)
{
  arma::mat X = Rcpp::as<arma::mat>(X_);
  arma::mat Y = Rcpp::as<arma::mat>(Y_);
  arma::mat ans = arma::solve(X, Y);
}