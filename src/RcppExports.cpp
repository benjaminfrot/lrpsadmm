// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcppeigen_fit_lrps
const Rcpp::GenericVector rcppeigen_fit_lrps(const Eigen::MatrixXd& Sigma, Eigen::MatrixXd& A, Eigen::MatrixXd& S, Eigen::MatrixXd& L, Eigen::MatrixXd& U, const Eigen::MatrixXd& Zeros, const double l1, const double l2, const double mu, const int max_iter, const double rel_tol, const double abs_tol, const int print_every, const int has_zeros);
RcppExport SEXP _lrpsadmm_rcppeigen_fit_lrps(SEXP SigmaSEXP, SEXP ASEXP, SEXP SSEXP, SEXP LSEXP, SEXP USEXP, SEXP ZerosSEXP, SEXP l1SEXP, SEXP l2SEXP, SEXP muSEXP, SEXP max_iterSEXP, SEXP rel_tolSEXP, SEXP abs_tolSEXP, SEXP print_everySEXP, SEXP has_zerosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type S(SSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type L(LSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type U(USEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Zeros(ZerosSEXP);
    Rcpp::traits::input_parameter< const double >::type l1(l1SEXP);
    Rcpp::traits::input_parameter< const double >::type l2(l2SEXP);
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< const double >::type rel_tol(rel_tolSEXP);
    Rcpp::traits::input_parameter< const double >::type abs_tol(abs_tolSEXP);
    Rcpp::traits::input_parameter< const int >::type print_every(print_everySEXP);
    Rcpp::traits::input_parameter< const int >::type has_zeros(has_zerosSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_fit_lrps(Sigma, A, S, L, U, Zeros, l1, l2, mu, max_iter, rel_tol, abs_tol, print_every, has_zeros));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_lrpsadmm_rcppeigen_fit_lrps", (DL_FUNC) &_lrpsadmm_rcppeigen_fit_lrps, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_lrpsadmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
