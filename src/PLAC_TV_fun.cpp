#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>

using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;

// (0): Prelim functions: SgInd(), PwInd() and TvInd():

// SgInd() = C equivalent to sg.ind(): at-risk processes;
// PwInd() = C modification to pw.ind(): truncation processes;
// TvInd() = Function to expand the time-varying covariate at each w_k;

// (1): Point estimation functions: BetaTv() and LambdaTv():

// BetaTv() = compute the scrore/information of beta and update beta with
//            one Newton-Raphson step;
// LambdaTv() = C equivalent to lambda0(): self-consistent updating function
//              for lambda_k's;

// (2): Variance estimation function:

// SWE_Tv() = C equivalent to sw(): sandwich variance-covariance estimator;

// (3): A function wrapper of all the above:

// PLAC_Tv() = wrapper function to calculate the PLAC esimator and its SWE;

// Input data:

// ZFV_ = pxn; time-varying covariates at ones own observed event time;
// Z = (mp)xn matrix of time-varying covariates at all distinct event times;
// X = mapped matrix containing As, Ys, Ds;
// W = mapped ordered unique observed event times;
// Ind1(k,i) = I(A_i <= w_k <= X_i);
// Ind2(k,i) = I(w_k <= A_i);
// Dn = mapped number of events at each w_k (typically they are all 1s);
// (b,h) = parametr values at which the function to be evaluated;

// Size parameters:

// n = sample size;
// m = number of distinct observed events;
// p = number of total covariates;

// Intermediate quantities:

// ebz  = a mxn matrix of exp(b.dot(Z));
// zebz = a (mp)xn matrix of Z*exp(b.dot(Z));
// hzebz = a (mp)xn matrix of h*Z*exp(b.dot(Z));

// rij = log-generalized odds ratio for pair (i, j);
// Rij  = generalized odds ratio for pair (i, j);
// Q0ij = a m-vector of Q0ij(t) at all w_k;
// Q1ij = a pxm-matrix of Q1ij(t) at all w_k;
// Q2ij = a pxp matrix of \int Q2ij(t) d\Lambda;

// UC = conditional likelihood score function(s);
// UP = pairwise likelihood score functions(s);
// JC = conditional likelihood observed Fisher information;
// JP = pairwise likelihood observed Fisher information;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd SgInd1(Eigen::Map<Eigen::MatrixXd> X,
                       Eigen::Map<Eigen::ArrayXd> W){

    const int n(X.rows()), m(W.size());
	MatrixXd out(m,n);

	for(int i = 0; i < n; ++i) {

        out.col(i) = (W > X(i,0) and W <= X(i,1)).cast<double>().matrix();

	}

	return out;

}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd PwInd1(Eigen::Map<Eigen::MatrixXd> X,
                       Eigen::Map<Eigen::ArrayXd> W){

    const int n(X.rows()), m(W.size());
	MatrixXd out(m,n);

	for(int i = 0; i < n; ++i) {

		out.col(i) = (W <= X(i,0)).cast<double>().matrix();

	}

	return out;

}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd TvInd1(Eigen::Map<Eigen::VectorXd> Zv,
                       Eigen::Map<Eigen::ArrayXd> W){

    const int n(Zv.size()), m(W.size());
    MatrixXd out(m,n);

	for(int i = 0; i < n; ++i) {

		out.col(i) = (W > Zv(i)).cast<double>().matrix();

	}

	return out;
}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd BetaTv2(Eigen::Map<Eigen::MatrixXd> ZFV_,
                        Eigen::Map<Eigen::MatrixXd> Z,
                        Eigen::Map<Eigen::MatrixXd> X,
                        Eigen::Map<Eigen::MatrixXd> Ind1,
                        Eigen::Map<Eigen::MatrixXd> Ind2,
                        Eigen::VectorXd b,
                        Eigen::VectorXd h){

    const int n(X.rows()), m(h.size()), p(b.size());

    double rij, Rij;

    VectorXd UC_b(p), UP_b(p), Q1ij(p), Zik(p), Zjk(p);

    UP_b.fill(0);

    MatrixXd ebz(m,n), hzebz(m*p,n), JC_b(p,p), JP_b(p,p), Q2ij(p,p);

    JC_b.fill(0); JP_b.fill(0);

    UC_b = ZFV_ * X.col(2);

    for(int k = 0; k < m; ++k){

        ebz.row(k) = (b.adjoint() * Z.middleRows(k*p,p)).array().exp();
        hzebz.middleRows(k*p,p) = h(k) * Z.middleRows(k*p,p).array() *
                                 ebz.row(k).replicate(p,1).array();

    }

    for(int i = 0; i < n; ++i){

        for(int k = 0; k < m; ++k){

            if(Ind1(k,i) == 1){

                UC_b += -hzebz.block(k*p,i,p,1);
                JC_b +=  hzebz.block(k*p,i,p,1)*Z.block(k*p,i,p,1).adjoint();

            }

        }

        for(int j = 0; j < n; ++j){

            rij = 0;
            Q1ij.fill(0);
            Q2ij.fill(0);

            for(int k = 0; k < m; ++k){

                if(Ind2(k,i) == Ind2(k,j)){

                    continue;

                }else if(Ind2(k,i) < Ind2(k,j)){

                    rij += h(k) * (ebz(k,j) - ebz(k,i));
                    Q1ij += hzebz.block(k*p,j,p,1) - hzebz.block(k*p,i,p,1);
                    Q2ij += hzebz.block(k*p,j,p,1)*Z.block(k*p,j,p,1).adjoint() -
                            hzebz.block(k*p,i,p,1)*Z.block(k*p,i,p,1).adjoint();

                }else{

                    rij += h(k) * (ebz(k,i) - ebz(k,j));
                    Q1ij += hzebz.block(k*p,i,p,1) - hzebz.block(k*p,j,p,1);
                    Q2ij += hzebz.block(k*p,i,p,1)*Z.block(k*p,i,p,1).adjoint() -
                            hzebz.block(k*p,j,p,1)*Z.block(k*p,j,p,1).adjoint();

                }


            }

            Rij = exp(rij);


            UP_b += -(Rij*Q1ij)/(1+Rij);
            JP_b += Rij/pow(1+Rij,2)*Q1ij*Q1ij.adjoint() + Rij/(1+Rij)*Q2ij;

        }

    }

    return (JC_b/n + JP_b/n/(n-1)).inverse() * (UC_b/n + UP_b/(n*(n-1)));

}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::VectorXd LambdaTv2(Eigen::Map<Eigen::MatrixXd> Z,
                          Eigen::Map<Eigen::MatrixXd> X,
                          Eigen::Map<Eigen::MatrixXd> Ind1,
                          Eigen::Map<Eigen::MatrixXd> Ind2,
                          Eigen::Map<Eigen::ArrayXd> Dn,
                          Eigen::VectorXd b,
                          Eigen::VectorXd h){

    const int n(X.rows()), m(h.size()), p(b.size());
    VectorXd UC_Lambda(m), UP_Lambda(m), Q0ij(m);
             UC_Lambda.fill(0); UP_Lambda.fill(0);
    MatrixXd ebz(m,n);

    for(int k = 0; k < m; ++k){

        ebz.row(k) = (b.adjoint() * Z.middleRows(k*p,p)).array().exp();

    }

    UC_Lambda = ebz.cwiseProduct(Ind1).rowwise().sum();

    for(int i = 0; i < n; ++i){

        for(int j = 0; j < n; ++j){

            for(int k = 0; k < m; ++k){

                if(Ind2(k,i) == Ind2(k,j)){

                    Q0ij(k) = 0;

                }else if(Ind2(k,i) < Ind2(k,j)){

                    Q0ij(k) = ebz(k,j) - ebz(k,i);

                }else{

                    Q0ij(k) = ebz(k,i) - ebz(k,j);

                }

            }

            UP_Lambda += Q0ij/(1+1/exp(h.dot(Q0ij)));

        }

    }

    return (Dn/(UC_Lambda/n+UP_Lambda/n/(n-1)).array()).matrix()/n;

}
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Eigen::MatrixXd SWE_Tv2(Eigen::Map<Eigen::MatrixXd> ZFV_,
                        Eigen::Map<Eigen::MatrixXd> Z,
                        Eigen::Map<Eigen::MatrixXd> X,
                        Eigen::Map<Eigen::ArrayXd> W,
                        Eigen::Map<Eigen::MatrixXd> Ind1,
                        Eigen::Map<Eigen::MatrixXd> Ind2,
                        Eigen::VectorXd b,
                        Eigen::VectorXd h) {

    const int n(X.rows()), m(h.size()), p(b.size());
    double Rij, R1ij, R2ij;
    VectorXd Q0ij(m), HQ1ij(p);
    MatrixXd ebz(m,n), zebz(m*p,n),
             UCs(p+m,n),  UPs(p+m,n), JC(p+m,p+m), JP(p+m,p+m),
             V(p+m,p+m),  J(p+m,p+m), J_(p+m,p+m),
             Q1ij(p,m), Q2ij(p,p);

    UCs.fill(0); UPs.fill(0); JC.fill(0);  JP.fill(0);

    UCs.topRows(p) = ZFV_.cwiseProduct(X.col(2).adjoint().replicate(p,1));

    for(int k = 0; k < m; ++k){

        ebz.row(k) = (b.adjoint() * Z.middleRows(k*p,p)).array().exp();
        zebz.middleRows(k*p,p) = Z.middleRows(k*p,p).array() *
                                 ebz.row(k).replicate(p,1).array();

    }

    for(int i = 0; i < n; ++i){

        UCs.col(i).tail(m) = (W == X(i,1)).select(X(i,2)*h.cwiseInverse(),0).matrix();

        for(int k = 0; k < m; ++k){

            if(Ind1(k,i) == 0){

                continue;

            }else{

                // Zik << Z.block(k*p,i,p,1);  // a column vector: Z_i at all obs. event times;

                UCs.col(i).head(p) += - zebz.block(k*p,i,p,1) * h(k);

                UCs(p+k,i) += - ebz(k,i);

                JC.topLeftCorner(p,p) += zebz.block(k*p,i,p,1) *
                                         Z.block(k*p,i,p,1).adjoint() * h(k);

                JC.block(0,p+k,p,1) += zebz.block(k*p,i,p,1);

            }


        }

        JC.diagonal().tail(m) += (W == X(i,1)).select(X(i,2)*h.cwiseInverse().array().pow(2),0).matrix();

        for(int j = 0; j < n; ++j){

            Q0ij.fill(0); Q1ij.fill(0); Q2ij.fill(0);

            for(int k = 0; k < m; ++k){

                if(Ind2(k,i) == Ind2(k,j)){

                   continue;

                }else if(Ind2(k,i) < Ind2(k,j)){

                    Q0ij(k) = (ebz(k,j) - ebz(k,i));
                    Q1ij.col(k) = (zebz.block(k*p,j,p,1) - zebz.block(k*p,i,p,1));
                    Q2ij += (zebz.block(k*p,j,p,1)*Z.block(k*p,j,p,1).adjoint() -
                             zebz.block(k*p,i,p,1)*Z.block(k*p,i,p,1).adjoint()) * h(k);

                }else{

                    Q0ij(k) = (ebz(k,i) - ebz(k,j));
                    Q1ij.col(k) = (zebz.block(k*p,i,p,1) - zebz.block(k*p,j,p,1));
                    Q2ij += (zebz.block(k*p,i,p,1)*Z.block(k*p,i,p,1).adjoint() -
                             zebz.block(k*p,j,p,1)*Z.block(k*p,j,p,1).adjoint()) * h(k);

                }

            }

            Rij = exp(h.dot(Q0ij));
            R1ij = Rij/(1+Rij);
            R2ij = Rij/pow(1+Rij,2);

            HQ1ij = Q1ij * h;

            UPs.col(i).head(p) += - HQ1ij * R1ij;
            UPs.col(i).tail(m) += - Q0ij * R1ij;

            if(i < j){

                JP.topLeftCorner(p,p)     += R2ij*HQ1ij*HQ1ij.adjoint() + R1ij*Q2ij;
                JP.topRightCorner(p,m)    += R2ij*HQ1ij*Q0ij.adjoint()  + R1ij*Q1ij;
                JP.bottomRightCorner(m,m) += R2ij*Q0ij*Q0ij.adjoint();

            }

        }

    }

    J = JC + 2*JP/(n-1);

    J.bottomLeftCorner(m,p) = J.topRightCorner(p,m).adjoint();

    V = UCs*UCs.adjoint() + 4*n/pow(n-1,3)*UPs*UPs.adjoint();

    J_ = J.inverse();

    return J_ * V * J_;

}

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List PLAC_Tv3(Eigen::Map<Eigen::MatrixXd> Z,
              Eigen::Map<Eigen::MatrixXd> ZFV_,
              Eigen::Map<Eigen::MatrixXd> X,
              Eigen::Map<Eigen::ArrayXd> W,
              Eigen::Map<Eigen::MatrixXd> Ind1,
              Eigen::Map<Eigen::MatrixXd> Ind2,
              Eigen::Map<Eigen::ArrayXd> Dn,
              Eigen::VectorXd b,
              Eigen::VectorXd h,
              int K = 100){

    const int n(X.rows()), m(h.size()), p(b.size());

    int k(0);
    double diff(99);
    VectorXd b_hat = b, h_hat = h, b_new, h_new, Diff(p+m);
    MatrixXd swe(p+m,p+m);

    while(diff > 0.0005 and k < K){

        h_new = LambdaTv2(Z, X, Ind1, Ind2, Dn, b_hat, h_hat);

        b_new = b_hat + BetaTv2(ZFV_, Z, X, Ind1, Ind2, b_hat, h_new);

        Diff << (b_new-b_hat),
                (h_new-h_hat);

        diff = Diff.cwiseAbs().maxCoeff();

        b_hat = b_new;
        h_hat = h_new;


        k++;

    }

    cout << k << " Iterations\n" << endl;

    swe = SWE_Tv2(ZFV_, Z, X, W, Ind1, Ind2, b_hat, h_hat);

    VectorXd se_b_hat(p),
             H_hat = VectorXd::Zero(m+1),
             se_H_hat = VectorXd::Zero(m+1);

    se_b_hat = swe.diagonal().head(p).cwiseSqrt();

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
