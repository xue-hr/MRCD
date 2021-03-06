// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// armaMatMult
SEXP armaMatMult(arma::mat A, arma::mat B);
RcppExport SEXP _MRCD_armaMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(armaMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenMatMult
SEXP eigenMatMult(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _MRCD_eigenMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}
// eigenMapMatMult
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B);
RcppExport SEXP _MRCD_eigenMapMatMult(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(eigenMapMatMult(A, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MRCD_armaMatMult", (DL_FUNC) &_MRCD_armaMatMult, 2},
    {"_MRCD_eigenMatMult", (DL_FUNC) &_MRCD_eigenMatMult, 2},
    {"_MRCD_eigenMapMatMult", (DL_FUNC) &_MRCD_eigenMapMatMult, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRCD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
