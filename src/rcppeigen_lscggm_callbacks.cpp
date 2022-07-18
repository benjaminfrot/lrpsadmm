// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
const Rcpp::List rcppeigen_Schur(const Eigen::MatrixXd A) {
  Eigen::RealSchur<Eigen::MatrixXd> SchurA(A);
  Eigen::MatrixXd R = SchurA.matrixT();
  Eigen::MatrixXd U = SchurA.matrixU();
  
  return Rcpp::List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("U") = U
  );  
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
const Rcpp::List rcppeigen_solve_sylvester(const Eigen::MatrixXd A,
                                           const Eigen::MatrixXd R,
                                           const Eigen::MatrixXd U,
                                           const Eigen::MatrixXd B,
                                           const Eigen::MatrixXd S,
                                           const Eigen::MatrixXd V,
                                           const Eigen::MatrixXd Q) {
  // Eigen::RealSchur<Eigen::MatrixXd> SchurB(B);
  //  Eigen::MatrixXd S = SchurB.matrixT();
  //  Eigen::MatrixXd V = SchurB.matrixU();
  Eigen::MatrixXd F = (U.adjoint() * Q) * V;
  
  Eigen::MatrixXd Y =
    Eigen::internal::matrix_function_solve_triangular_sylvester(R, S, F);
  
  Eigen::MatrixXd X = ((U * Y) * V.adjoint()).real();
  Eigen::MatrixXd Q_calc = A * X + X * B;
  const double diff = (Q - Q_calc).norm();
  return Rcpp::List::create(
    Rcpp::Named("AZX") = X,
    Rcpp::Named("diff") = diff
  ); 
}