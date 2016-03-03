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

// SgInd() = C equivalent to sg.ind(): at-risk processes;
// PwInd() = C modification to pw.ind(): truncation processes;
// ScrBeta() = C equivalent to scr.b(): scrore function of beta;
// FsrBeta() = C equivalent to fsr.b(): observed Fisher of beta;
// Lambda() = C equivalent to lambda0(): self-consistent updating function for lambda_k's;
// SWE() = C equivalent to sw(): sandwich variance-covariance estimator;
// PLAC() = wrapper function to calculate the PLAC esimator and its SWE;

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
// Lambda_diff = \Lambda(A_i) - \Lambda(A_j);
// rij, Rij = (log) generalized odds ratio;
// Q1ij, Q2ij = \int Q1(t) d\Lambda and \int Q2(t) d\Lambda;

Eigen::VectorXd ScrBeta(Eigen::Map<Eigen::MatrixXd> Z,
                 Eigen::Map<Eigen::MatrixXd> X,
                 Eigen::Map<Eigen::MatrixXd> Ind1,
                 Eigen::Map<Eigen::MatrixXd> Ind2,
                 Eigen::VectorXd b,
                 Eigen::VectorXd h){
  const int n(Z.rows()), p(Z.cols());
  VectorXd UC_b(p), UP_b(p), Q1ij;
           UC_b.fill(0); UP_b.fill(0);
  double Lambda_diff, Rij, exbZi;
  for(int i = 0; i < n; ++i){
    exbZi = exp(b.dot(Z.row(i)));
    UC_b += (X(i,2) - h.dot(Ind1.row(i)) * exbZi) * Z.row(i);
    for(int j = 0; j < n; ++j){
      Lambda_diff = h.dot(Ind2.row(i) - Ind2.row(j));
      Rij = exp((exbZi - exp(b.dot(Z.row(j)))) * Lambda_diff);
      Q1ij = (Z.row(i)*exbZi - Z.row(j)*exp(b.dot(Z.row(j)))) * Lambda_diff;
      UP_b += -(Rij*Q1ij)/(1+Rij);
    }
  }
  return UC_b/n + UP_b/(n*(n-1));
}

Eigen::MatrixXd FsrBeta(Eigen::Map<Eigen::MatrixXd> Z,
                 Eigen::Map<Eigen::MatrixXd> X,
                 Eigen::Map<Eigen::MatrixXd> Ind1,
                 Eigen::Map<Eigen::MatrixXd> Ind2,
                 Eigen::VectorXd b,
                 Eigen::VectorXd h){
  const int n(Z.rows()), m(Ind1.cols()), p(Z.cols());
  MatrixXd JC_b(p,p),    JP_b(p,p);
           JC_b.fill(0); JP_b.fill(0);
  double Lambda_diff, Rij, exbZi, exbZj;
  VectorXd Q1ij;
  for(int i = 0; i < n; ++i){
    exbZi = exp(b.dot(Z.row(i)));
    JC_b += h.dot(Ind1.row(i))*exbZi/n * (Z.row(i).adjoint()*Z.row(i));
    for(int j = 0; j < n; ++j){
      exbZj = exp(b.dot(Z.row(j)));
      Lambda_diff = h.dot(Ind2.row(i) - Ind2.row(j));
      Rij = exp((exbZi - exbZj) * Lambda_diff);
      Q1ij = (Z.row(i)*exbZi - Z.row(j)*exbZj) * Lambda_diff;
      JP_b += Rij/pow(1+Rij,2)*Q1ij*Q1ij.adjoint() +
              Rij/(1+Rij)*Lambda_diff*
              (Z.row(i).adjoint()*Z.row(i)*exbZi -
               Z.row(j).adjoint()*Z.row(j)*exbZj);
    }
  }
  return JC_b+JP_b/n/(n-1);
}

VectorXd Lambda(Eigen::Map<Eigen::MatrixXd> Z,
                Eigen::Map<Eigen::MatrixXd> X,
                Eigen::Map<Eigen::MatrixXd> Ind1,
                Eigen::Map<Eigen::MatrixXd> Ind2,
                Eigen::Map<Eigen::ArrayXd> Dn,
                Eigen::VectorXd b,
                Eigen::VectorXd h){
  const int n(Z.rows()), m(Ind1.cols()), p(Z.cols());
  VectorXd UC_Lambda(m), UP_Lambda(m), Ind2ij;
           UC_Lambda.fill(0); UP_Lambda.fill(0);
  double Rij, exbZi;
  for(int i = 0; i < n; ++i){
    exbZi = exp(b.dot(Z.row(i))); // to save some computation time in the j-loop
    UC_Lambda += exbZi*Ind1.row(i);
    for(int j = 0; j < n; ++j){
      Ind2ij = Ind2.row(i) - Ind2.row(j);
      Rij = exp((exbZi - exp(b.dot(Z.row(j)))) * h.dot(Ind2ij));
      UP_Lambda += (exbZi - exp(b.dot(Z.row(j)))) / (1+1/Rij)*Ind2ij;
    }
  }
  return (Dn/(UC_Lambda/n+UP_Lambda/n/(n-1)).array()).matrix()/n;
}

Eigen::MatrixXd SWE(Eigen::Map<Eigen::MatrixXd> Z,
                    Eigen::Map<Eigen::MatrixXd> X,
                    Eigen::Map<Eigen::ArrayXd> W,
                    Eigen::Map<Eigen::MatrixXd> Ind1,
                    Eigen::Map<Eigen::MatrixXd> Ind2,
                    Eigen::VectorXd b,
                    Eigen::VectorXd h) {
  const int n(X.rows()), m(h.size()), p(b.size());
  MatrixXd UCs(n,p+m), UPs(n,p+m), JC(p+m,p+m), JP(p+m,p+m),
           V(p+m,p+m), J(p+m,p+m), J_(p+m,p+m), Zi2exbZi(p,p);
           UCs.fill(0); UPs.fill(0); JC.fill(0); JP.fill(0);
  VectorXd Q1ij(p);
  double Rij, exbZi, exbZj;
  // Meat Matrix
  for(int i = 0; i < n; ++i){
    exbZi = exp(b.dot(Z.row(i)));
    Zi2exbZi = Z.row(i).adjoint()*Z.row(i)*exbZi;
    UCs.row(i).head(p) = (X(i,2) - h.dot(Ind1.row(i))*exbZi) * Z.row(i);
    UCs.row(i).tail(m) = (W == X(i,1)).select(X(i,2)*h.cwiseInverse(),0).matrix()
                         -exbZi*Ind1.row(i).adjoint();
    JC.topLeftCorner(p,p) += h.dot(Ind1.row(i))*exbZi*(Z.row(i).adjoint()*Z.row(i));
    JC.bottomLeftCorner(m,p) += exbZi * Ind1.row(i).adjoint()*Z.row(i);
    JC.bottomRightCorner(m,m).diagonal() += (W == X(i,1)).select(X(i,2)*h.cwiseInverse().array().pow(2),0).matrix();
    for(int j = 0; j < n; ++j){
      exbZj = exp(b.dot(Z.row(j)));
      Rij = exp((exbZi - exbZj) * h.dot(Ind2.row(i) - Ind2.row(j)));
      Q1ij = (Z.row(i)*exbZi - Z.row(j)*exbZj) * h.dot(Ind2.row(i) - Ind2.row(j));
      UPs.row(i).head(p) += -Q1ij/(1+1/Rij);
      UPs.row(i).tail(m) += -(exbZi - exbZj)*(Ind2.row(i) - Ind2.row(j))/(1+1/Rij);
      if(i < j){
        JP.topLeftCorner(p,p) += Rij/pow(1+Rij,2)*Q1ij*Q1ij.adjoint() +
                         Rij/(1+Rij)*h.dot(Ind2.row(i) - Ind2.row(j))*(Zi2exbZi -
                         Z.row(j).adjoint()*Z.row(j)*exbZj);
        JP.bottomLeftCorner(m,p) += (Rij*h.dot(Ind2.row(i) - Ind2.row(j))*(exbZi - exbZj)/pow(1+Rij,2) + Rij/(1+Rij)) *
                            (Ind2.row(i) - Ind2.row(j)).adjoint() * (Z.row(i)*exbZi - Z.row(j)*exbZj);
        JP.bottomRightCorner(m,m) += (Rij*pow((exbZi-exbZj)/(1+Rij),2)) *
                                     (Ind2.row(i) - Ind2.row(j)).adjoint() * (Ind2.row(i) - Ind2.row(j));
      }
    }
  }
  J = JC + 2 * JP/(n-1);
  J.topRightCorner(p,m) = J.bottomLeftCorner(m,p).adjoint();
  V = UCs.adjoint()*UCs + 4*n/pow(n-1,3)*UPs.adjoint()*UPs;
  J_ = J.inverse();
  return J_*V*J_;
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
  const int n(Z.rows()), m(h.size()), p(b.size());
  int k(0);
  double diff(99);
  VectorXd b_hat = b, h_hat = h, b_new, h_new, Diff(p+m);
  MatrixXd swe(p+m,p+m);
  while(diff > 0.0005 and k < K){
    h_new = Lambda(Z,X,Ind1,Ind2,Dn,b_hat,h_hat);
    b_new = b_hat + FsrBeta(Z,X,Ind1,Ind2,b_hat,h_new).inverse()*ScrBeta(Z,X,Ind1,Ind2,b_hat,h_new);
    Diff << (b_new-b_hat),
            (h_new-h_hat);
    diff = Diff.cwiseAbs().maxCoeff();
    b_hat = b_new;
    h_hat = h_new;
    k++;
  }
  Rcout << k << " Iterations" << endl;
  swe = SWE(Z,X,W,Ind1,Ind2,b_hat,h_hat);
  VectorXd se_b_hat(p),
           H_hat = VectorXd::Zero(m+1),
           se_H_hat = VectorXd::Zero(m+1);
  se_b_hat = swe.topLeftCorner(p,p).diagonal().cwiseSqrt();
  H_hat.tail(m) = h_hat;
  for(int i = 1; i <= m; ++i){
    H_hat(i) += H_hat(i-1);
    se_H_hat(i) = sqrt(swe.block(p,p,i,i).sum());
  }
  return List::create(Named("b.hat") = b_hat,
                      Named("h0.hat") = h_hat,
                      Named("swe") = swe,
                      Named("se.b.hat") = se_b_hat,
                      Named("H0.hat") = H_hat,
                      Named("se.H0.hat") = se_H_hat,
                      Named("iter") = k);
}
