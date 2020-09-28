#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

// [[Rcpp::export]]
SEXP armaMatMult(arma::mat A, arma::mat B){
  arma::mat C = A * B;

  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}
