#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

// C++ Functions:
// ScrBeta() = C equivalent to scr.b(): scrore function of beta;
// FsrBeta() = C equivalent to fsr.b(): observed Fisher of beta;
// Lambda() = C equivalent to lambda0(): self-consistent updating function for lambda_k's;
// SWE() = C equivalent to sw(): sandwich variance-covariance estimator;
// PLAC_TI() = wrapper function to calculate the PLAC esimator and its SWE;

// Input data:
// Z = mapped covariate matrix;
// X = mapped matrix containing As, Ys, Ds;
// W = mapped ordered unique observed event times;
// Ind1(i,k) = I(A_i <= w_k <= X_i);
// Ind2(i,k) = I(w_k <= A_i);
// Dn = mapped number of events at each w_k (typically they are all one's);
// (b,h) = parametr values at which the function to be evaluated;

// Parameters:
// n = sample size;
// m = number of distinct observed events;
// p = number of total covariates;

// Intermediates:
// UC = conditional likelihood score function(s);
// UP = pairwise likelihood score functions(s);
// JC = conditional likelihood observed Fisher information;
// JP = pairwise likelihood observed Fisher information;
// Ind2ij = I(w_k <= A_i) - I(w_k <= A_j);
// H_diff = \Lambda(A_i) - \Lambda(A_j);
// rij, Rij = (log) generalized odds ratio;
// Q1ij, Q2ij = \int Q1(t) d\Lambda and \int Q2(t) d\Lambda;

Eigen::MatrixXd Beta(Eigen::Map<Eigen::MatrixXd> Z,
                        Eigen::Map<Eigen::MatrixXd> X,
                        Eigen::Map<Eigen::MatrixXd> Ind1,
                        Eigen::Map<Eigen::MatrixXd> Ind2,
                        Eigen::VectorXd b,
                        Eigen::VectorXd h){
  const int n(X.rows()), p(b.size());
  double H_diff, Rij;
  VectorXd UC_b(p), UP_b(p), ebz(n), Q1ij(p);
  UC_b.fill(0); UP_b.fill(0);
  MatrixXd JC_b(p,p), JP_b(p,p), zebz(p, n);
  JC_b.fill(0); JP_b.fill(0);
  ebz = (Z.adjoint() * b).array().exp();
  zebz = ebz.adjoint().replicate(p,1).array() * Z.array();
  for(int i = 0; i < n; ++i){
    UC_b += (X(i,2) - h.dot(Ind1.col(i)) * ebz(i)) * Z.col(i);
    JC_b += h.dot(Ind1.col(i)) * ebz(i) * (Z.col(i) * Z.col(i).adjoint());
    for(int j = 0; j < n; ++j){
      H_diff = h.dot(Ind2.col(i) - Ind2.col(j));
      Rij = exp((ebz(i) - ebz(j)) * H_diff);
      Q1ij = (zebz.col(i) - zebz.col(j)) * H_diff;
      UP_b += -(Rij * Q1ij) / (1 + Rij);
      JP_b += Rij / pow(1 + Rij, 2) * Q1ij * Q1ij.adjoint() +
              Rij / (1 + Rij) * H_diff *
              (Z.col(i) * Z.col(i).adjoint() * ebz(i) -
               Z.col(j) * Z.col(j).adjoint() * ebz(j));
    }
  }
  return (JC_b / n + JP_b / n / (n - 1)).inverse() * (UC_b / n + UP_b / (n * (n - 1)));
}

VectorXd Lambda(Eigen::Map<Eigen::MatrixXd> Z,
                Eigen::Map<Eigen::MatrixXd> X,
                Eigen::Map<Eigen::MatrixXd> Ind1,
                Eigen::Map<Eigen::MatrixXd> Ind2,
                Eigen::Map<Eigen::ArrayXd> Dn,
                Eigen::VectorXd b,
                Eigen::VectorXd h){
  const int n(X.rows()), m(h.size());
  double Rij;
  VectorXd UC_Lambda(m), UP_Lambda(m), Ind2ij(m), ebz(n);
           UC_Lambda.fill(0); UP_Lambda.fill(0);
  ebz = (Z.adjoint() * b).array().exp();
  for(int i = 0; i < n; ++i){
    UC_Lambda += ebz(i) * Ind1.col(i);
    for(int j = 0; j < n; ++j){
      Ind2ij = Ind2.col(i) - Ind2.col(j);
      Rij = exp((ebz(i) - ebz(j)) * h.dot(Ind2ij));
      UP_Lambda += (ebz(i) - ebz(j)) / (1 + 1 / Rij) * Ind2ij;
    }
  }
  return (Dn / (UC_Lambda / n + UP_Lambda / n / (n - 1)).array()).matrix() / n;
}

Eigen::MatrixXd SWE(Eigen::Map<Eigen::MatrixXd> Z,
                    Eigen::Map<Eigen::MatrixXd> X,
                    Eigen::Map<Eigen::ArrayXd> W,
                    Eigen::Map<Eigen::MatrixXd> Ind1,
                    Eigen::Map<Eigen::MatrixXd> Ind2,
                    Eigen::VectorXd b,
                    Eigen::VectorXd h) {
  const int n(X.rows()), m(h.size()), p(b.size());
  MatrixXd UCs(p + m, n), UPs(p + m, n),
           JC(p + m, p + m), JP(p + m, p + m),
           V(p + m, p + m), J(p + m, p + m), J_(p + m, p + m),
           zebz(p,n), Zi2exbZi(p, p);
           UCs.fill(0); UPs.fill(0); JC.fill(0); JP.fill(0);
  VectorXd Q1ij(p), ebz(n);
  double Rij;
  ebz = (Z.adjoint() * b).array().exp();
  zebz = ebz.adjoint().replicate(p,1).array() * Z.array();
  for(int i = 0; i < n; ++i){
    Zi2exbZi = zebz.col(i)*Z.col(i).adjoint();
    UCs.col(i).head(p) = (X(i,2) - h.dot(Ind1.col(i))*ebz(i)) * Z.col(i);
    UCs.col(i).tail(m) = (W == X(i,1)).select(X(i,2)*h.cwiseInverse(),0).matrix()
                         -ebz(i)*Ind1.col(i);
    JC.topLeftCorner(p,p) += h.dot(Ind1.col(i)) * zebz.col(i) * Z.col(i).adjoint();
    JC.bottomLeftCorner(m,p) +=  Ind1.col(i) * zebz.col(i).adjoint();
    JC.bottomRightCorner(m,m).diagonal() += (W == X(i, 1)).select(X(i, 2) *
                                             h.cwiseInverse().array().pow(2),0).matrix();
    for(int j = 0; j < n; ++j){
      Rij = exp((ebz(i) - ebz(j)) * h.dot(Ind2.col(i) - Ind2.col(j)));
      Q1ij = (zebz.col(i) - zebz.col(j)) * h.dot(Ind2.col(i) - Ind2.col(j));
      UPs.col(i).head(p) += -Q1ij / (1 + 1 / Rij);
      UPs.col(i).tail(m) += -(ebz(i) - ebz(j)) * (Ind2.col(i) - Ind2.col(j)) / (1 + 1 / Rij);
      if(i < j){
        JP.topLeftCorner(p, p) += Rij / pow(1 + Rij, 2) * Q1ij * Q1ij.adjoint() +
                         Rij / (1 + Rij) * h.dot(Ind2.col(i) - Ind2.col(j)) *
                         (Zi2exbZi - zebz.col(j) * Z.col(j).adjoint());
        JP.bottomLeftCorner(m,p) += (Rij * h.dot(Ind2.col(i) - Ind2.col(j)) *
                                    (ebz(i) - ebz(j)) / pow(1 + Rij, 2) + Rij / (1 + Rij)) *
                                    (Ind2.col(i) - Ind2.col(j)) * (zebz.col(i) - zebz.col(j)).adjoint();
        JP.bottomRightCorner(m,m) += (Rij * pow((ebz(i) - ebz(j)) / (1 + Rij), 2)) *
                                     (Ind2.col(i) - Ind2.col(j)) * (Ind2.col(i) - Ind2.col(j)).adjoint();
      }
    }
  }
  J = JC + 2 * JP / (n - 1);
  J.topRightCorner(p, m) = J.bottomLeftCorner(m, p).adjoint();
  V = UCs * UCs.adjoint() + 4 * n / pow(n - 1, 3) * UPs * UPs.adjoint();
  J_ = J.inverse();
  return J_ * V * J_;
}

//' C++ Function for Solving the PLAC Estimator.
//' (with time-invariant convariates only)
//'
//' @param Z matrix for all the covariates history.
//' @param X the response matrix (As, Xs, Ds).
//' @param W the ordered observed event times.
//' @param Ind1 risk-set indicators.
//' @param Ind2 truncation pair indicators.
//' @param Dn number of ties at each observed event time.
//' @param b initial values of the regression coefficients.
//' @param h initial values of the baseline hazard function.
//' @param K maximal iteration number, the default is \code{K = 100}.
//' @return list of model fitting results for both conditional approach and
//' the PLAC estimator.
//' @export
// [[Rcpp::export]]
List PLAC_TI(Eigen::Map<Eigen::MatrixXd> Z,
          Eigen::Map<Eigen::MatrixXd> X,
          Eigen::Map<Eigen::ArrayXd> W,
          Eigen::Map<Eigen::MatrixXd> Ind1,
          Eigen::Map<Eigen::MatrixXd> Ind2,
          Eigen::Map<Eigen::ArrayXd> Dn,
          Eigen::VectorXd b,
          Eigen::VectorXd h,
          int K = 100){
  const int m(h.size()), p(b.size());
  int k(0);
  double diff(99);
  VectorXd b_hat = b, h_hat = h, b_new, h_new, Diff(p + m);
  MatrixXd swe(p + m, p + m);
  while(diff > 0.0005 and k < K){
    h_new = Lambda(Z, X, Ind1, Ind2, Dn, b_hat, h_hat);
    b_new = b_hat + Beta(Z, X, Ind1, Ind2, b_hat, h_new);
    Diff << (b_new - b_hat),
            (h_new - h_hat);
    diff = Diff.cwiseAbs().maxCoeff();
    b_hat = b_new;
    h_hat = h_new;
    k++;
  }
  Rcout << k << " Iterations" << endl;
  swe = SWE(Z, X, W, Ind1, Ind2, b_hat, h_hat);
  VectorXd se_b_hat(p),
           H_hat = VectorXd::Zero(m + 1),
           se_H_hat = VectorXd::Zero(m + 1);
  se_b_hat = swe.topLeftCorner(p, p).diagonal().cwiseSqrt();
  H_hat.tail(m) = h_hat;
  for(int i = 1; i <= m; ++i){
    H_hat(i) += H_hat(i - 1);
    se_H_hat(i) = sqrt(swe.block(p, p, i, i).sum());
  }
  return List::create(Named("b.hat") = b_hat,
                      Named("h0.hat") = h_hat,
                      Named("swe") = swe,
                      Named("se.b.hat") = se_b_hat,
                      Named("H0.hat") = H_hat,
                      Named("se.H0.hat") = se_H_hat,
                      Named("iter") = k);
}
