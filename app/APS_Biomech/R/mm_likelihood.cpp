#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::ColPivHouseholderQR;
using namespace Rcpp;

// [[Rcpp::export(rng=false)]]
List mm_likelihood(Map<VectorXd> theta, Map<MatrixXd> X, Map<MatrixXd> Z, Map<VectorXd> y) {
  
  // Constructs the loglikelihood of the mixed model representation of
  // a pspline with fixed effects X and random effects Z.
  // given theta with theta[1] being log transformed variance of random effects...
  
  // For use in maximisation of the profile likelihood of a MM representation.
  // From Wood 2017 pp. 79
  
  // untransform parameters...
  double sigma_b = exp(theta[0]);
  double sigma = exp(theta[1]);
  double tau_b = 1/(sigma_b * sigma_b);
  double tau = 1/(sigma*sigma);
  
  // extract dimensions...
  const int n =  y.size();
  const int pr = Z.cols();
  const int pf = X.cols();
  const int K = pr+pf;
  
  const double pi = atan(1)*4;
  
  MatrixXd X1(n, pf+pr);
  X1 << X, Z;
  VectorXd ipsi(K);
  for(int i = 0; i<pf; i++){
    ipsi[i] = 0;
  }
  for(int i = 0; i<pr; i++){
    ipsi[i+pf] = tau_b;
  }

  
  MatrixXd mipsi = ipsi.asDiagonal();
  
  MatrixXd A(K, K);
  A = ((X1.transpose() * X1) * tau) + mipsi;
  
  VectorXd b(K);
  b = (X1.transpose() * y) * tau;

  
  
  VectorXd b_hat(K);
  b_hat = A.colPivHouseholderQr().solve(b);  // Solve Ax = b using LU decomposition.

  
  
  MatrixXd ZpZ(pr, pr);
  ZpZ = Z.transpose()*Z;
  
  
  MatrixXd ipsi_short;
  ipsi_short = MatrixXd::Zero(pr, pr) ;
  
  for(int i = 0; i<pr; i++){
    ipsi_short(i, i) = tau_b;
  }
  MatrixXd B(pr, pr);
  B = (ZpZ * tau) + ipsi_short;
  
  MatrixXd L(B.llt().matrixL());
// ## compute log|Z’Z/sigma^2 + I/sigma.b^2|...
//   ldet <- sum(log(diag(chol(crossprod(Z) / sigma^2 +
//     diag(ipsi[-(1:pf)])))))
  double ldet = 0;
  for(int i = 0; i<pr; i++){
    ldet += log(L(i,i));
  }


  double ll = (-((y - X1*b_hat).array().pow(2)).sum()*tau -
    ((b_hat.array().pow(2) * ipsi.array()).sum()) - 
    (n * log(sigma*sigma)) - 
    pr*log(sigma_b*sigma_b) - 
    2*ldet - 
    n*log(2*pi))/2;
  // VectorXd y_hat(n);
  // y_hat = X1 * b_hat;
  // l <- (-sum((y - X1 %*% b1)^2) / sigma^2 - sum(b1^2 * ipsi) -
  //   n * log(sigma^2) - pr * log(sigma.b^2) - 2 * ldet - n * log(2 * pi)) / 2
  List ret;
  ret["ll"] = -1*ll;
  ret["alp"] = b_hat.head(pf);
  ret["gam"] = b_hat.tail(pr);
  return ret;
  
}


// Testing R Block
// 
// 
// /*** R
// theta = c(log(0.2), log(0.1))
// 
// source("library.R")
// set.seed(1)
// t = seq(0, 1, length.out = 100)
// y = rnorm(100)
// print("Generating BSpline Basis...")
// deg <- 5
// pDeg <- 2
// nK <- 10
// B <- basis(t, nK, deg = deg)
// P <- penalty(dim(B$matrix)[2], pDeg)
// MM <- basis_mm(P)
// X = B$matrix %*% MM$X
// Z = B$matrix %*% MM$Z
// print("Cpp ll")
// print(mm_likelihood(theta, X, Z, y))
// print("R ll")
// print(MM_likelihood(theta, X, Z, y))
// 
// */
