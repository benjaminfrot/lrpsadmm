// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
using namespace Eigen;

double l1_norm(const Eigen::MatrixXd& X) {
  return X.lpNorm<1>();
}

double logdet(const Eigen::MatrixXd& X) {
  double logdet = 0.0;
  const int p = X.cols();
  
  SelfAdjointEigenSolver<MatrixXd> es(X);
  for (int i=0; i<p; i++) {
    if(abs(es.eigenvalues()[i]) < 1e-08) {
      logdet = logdet + log(1e-16); 
    } else {
      logdet = logdet + log(es.eigenvalues()[i]);
    }
  }
  
  return logdet;
}

// [[Rcpp::depends(RcppEigen)]]
const Eigen::MatrixXd rcppeigen_update_A(const Eigen::MatrixXd& Sigma,
                                         const Eigen::MatrixXd& S,
                                         const Eigen::MatrixXd& L,
                                         const Eigen::MatrixXd& U,
                                         const double mu) 
{
  const int p = Sigma.cols();
  const Eigen::MatrixXd X1 = mu * (S - L) - Sigma - U;
  const Eigen::MatrixXd X2 = X1.pow(2) + 4 * mu * MatrixXd::Identity(p, p);
  const Eigen::MatrixXd sqrtX2 = X2.sqrt();
  
  return (X1 + sqrtX2) / (2.0 * mu);
}

// [[Rcpp::depends(RcppEigen)]]
const Eigen::MatrixXd rcppeigen_update_S(const Eigen::MatrixXd& A,
                                         const Eigen::MatrixXd& L,
                                         const Eigen::MatrixXd& U,
                                         const double l1,
                                         const double mu,
                                         const Eigen::MatrixXd& Zeros,
                                         const int has_zeros) 
{
  const MatrixXd X1 = A + L + (U / mu);
  // Shrinkage
  MatrixXd X2 = (X1.array().abs() - (l1 / mu)).matrix();
  X2 = ((0.5 * (X2.array().sign() + 1)) * X2.array()).matrix();
  
  // Enforce the prior knowledge on with Zeros:
  MatrixXd S;
  if (has_zeros > 0) {
    S = (Zeros.array() * X2.array() * X1.array().sign()).matrix();
  } else {
    S = (X2.array() * X1.array().sign()).matrix();
  }
  return S;
}

// [[Rcpp::depends(RcppEigen)]]
const Eigen::MatrixXd rcppeigen_update_L(const Eigen::MatrixXd& A,
                                         const Eigen::MatrixXd& S,
                                         const Eigen::MatrixXd& U,
                                         const double l2,
                                         const double mu) 
{
  const MatrixXd X1 = S - A - (U / mu);
  const int p = A.cols();
  MatrixXd shrunk_eigvals = MatrixXd::Identity(p, p);
  SelfAdjointEigenSolver<MatrixXd> es(X1);
  for (int i=0; i<p; i++) {
    shrunk_eigvals(i,i) = fmax(es.eigenvalues()[i] - (l2 / mu), 0.0);
  }
  const MatrixXd L = es.eigenvectors() * shrunk_eigvals * es.eigenvectors().transpose(); 
  
  return L;
}

// [[Rcpp::depends(RcppEigen)]]
const Eigen::MatrixXd rcppeigen_update_U(const Eigen::MatrixXd& A,
                                         const Eigen::MatrixXd& S,
                                         const Eigen::MatrixXd& L,
                                         const Eigen::MatrixXd& U,
                                         const double mu) 
{
  const MatrixXd U2 = U + mu * (A - S + L);
  
  return U2;
}
// [[Rcpp::depends(RcppEigen)]]
double rcppeigen_obj_fun(const Eigen::MatrixXd& Sigma, 
                         const Eigen::MatrixXd& A,
                         const Eigen::MatrixXd& S,
                         const Eigen::MatrixXd& L,
                         double l1,
                         double l2) {
  // We compute
  // -log det A + Tr(Sigma * (A)) + l1 ||S||_1 + l2 Tr(L)
  return (Sigma * A).trace() - 
    logdet(A) + 
    l1 * l1_norm(S) +
    l2 * L.trace();
}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
const Rcpp::GenericVector rcppeigen_fit_lrps(const Eigen::MatrixXd& Sigma,
                                             Eigen::MatrixXd& A,
                                             Eigen::MatrixXd& S,
                                             Eigen::MatrixXd& L,
                                             Eigen::MatrixXd& U,
                                             const Eigen::MatrixXd& Zeros,
                                             const double l1,
                                             const double l2,
                                             const double mu,
                                             const int max_iter,
                                             const double rel_tol,
                                             const double abs_tol,
                                             const int print_every,
                                             const int has_zeros
) {
  Rcpp::NumericVector r_norms;
  Rcpp::NumericVector s_norms;
  Rcpp::NumericVector eps_duals;
  Rcpp::NumericVector eps_pris;
  Rcpp::NumericVector objvals;
  int termcode = -1;
  double obj_fun;
  double s_norm;
  double r_norm;
  double eps_pri;
  double eps_dual;
  const int p = Sigma.cols();
  
  if (print_every > 0) {
    Rcpp::Rcout << "#Iter | Objfun | snorm | eps_dual | rnorm | eps_pri\n";
  }
  
  for (int i=0; i<max_iter; i++) {
    
    A = rcppeigen_update_A(Sigma, S, L, U, mu);
    MatrixXd newS = rcppeigen_update_S(A, L, U, l1, mu, Zeros, has_zeros);
    if (newS.diagonal().any() == 0) {
      termcode = -2;
      break;
    }
    MatrixXd newL = rcppeigen_update_L(A, newS, U, l2, mu);
    
    U = rcppeigen_update_U(A, newS, newL, U, mu);
    
    // Diagnostics 
    obj_fun = rcppeigen_obj_fun(Sigma, A, newS, newL, l1, l2);
    r_norm = (A - (newS - newL)).norm(); //Frob. norm
    s_norm = (mu * ((newS - newL) - (S - L))).norm();
    eps_pri = p * abs_tol + rel_tol * fmax(A.norm(), (newS-newL).norm());
    eps_dual = p * abs_tol + rel_tol * (mu * U).norm();
    
    objvals.push_back(obj_fun);
    s_norms.push_back(s_norm);
    r_norms.push_back(r_norm);
    eps_pris.push_back(eps_pri);
    eps_duals.push_back(eps_dual);
    
    if (print_every > 0) {
      if (i % print_every == 0) {
        Rcpp::Rcout << i;
        Rcpp::Rcout << " | ";
        Rcpp::Rcout << obj_fun;
        Rcpp::Rcout << " | ";
        Rcpp::Rcout << s_norm;
        Rcpp::Rcout << " | ";
        Rcpp::Rcout << eps_dual;
        Rcpp::Rcout << " | ";
        Rcpp::Rcout << r_norm;
        Rcpp::Rcout << " | ";
        Rcpp::Rcout << eps_pri;
        Rcpp::Rcout << "\n";
      }
    }
    
    if ((s_norm < eps_dual) && (r_norm < eps_pri)) {
      termcode = 0;
      break;
    }
    
    // Update S and L
    S = newS;
    L = newL;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("objvals") = objvals,
    Rcpp::Named("snorms") = s_norms,
    Rcpp::Named("rnorms") = r_norms,
    Rcpp::Named("epspris") = eps_pris,
    Rcpp::Named("epsduals") = eps_duals,
    Rcpp::Named("S") = S,
    Rcpp::Named("L") = L,
    Rcpp::Named("U") = U,
    Rcpp::Named("termcode") = termcode
  );
}