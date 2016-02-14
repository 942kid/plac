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


// [[Rcpp::export]]
MatrixXd SgInd1(Map<MatrixXd> X, Map<ArrayXd> W){
    
    const int n(X.rows()), m(W.size());
	MatrixXd out(m,n);
    
	for(int i = 0; i < n; ++i) {
        
        out.col(i) = (W > X(i,0) and W <= X(i,1)).cast<double>().matrix();
        
	}
    
	return out;
    
}

// [[Rcpp::export]]
MatrixXd PwInd1(Map<MatrixXd> X, Map<ArrayXd> W){
    
    const int n(X.rows()), m(W.size());
	MatrixXd out(m,n);
    
	for(int i = 0; i < n; ++i) {
        
		out.col(i) = (W <= X(i,0)).cast<double>().matrix();
        
	}
    
	return out;
    
}

// [[Rcpp::export]]
MatrixXd TvInd2(Map<VectorXd> PD,
                Map<VectorXd> Zv, 
                Map<ArrayXd> W){
    
    const int n(Zv.size()), m(W.size());
    MatrixXd out(m,n);
    
	for(int i = 0; i < n; ++i) {
        
		out.col(i) = (W <= Zv(i)).cast<double>().matrix() * PD(i);
        
	}
    
	return out;
}

// [[Rcpp::export]]
MatrixXd TvInd3(Map<VectorXd> PD,
                Map<VectorXd> Zv, 
                Map<ArrayXd> W){
    
    const int n(Zv.size()), m(W.size());
    MatrixXd out(m,n);
    
    for(int i = 0; i < n; ++i) {
        
        out.col(i) = (W <= Zv(i)).cast<double>().matrix() * PD(i)
                    -(W >  Zv(i)).cast<double>().matrix();
        
    }
    
    return out;
}

// [[Rcpp::export]]
MatrixXd BetaTvR1(Map<MatrixXd> ZF,
                 Map<MatrixXd> ZFV_,
                 Map<MatrixXd> Z,
                 Map<MatrixXd> X, 
                 Map<MatrixXd> Ind1,
                 Map<MatrixXd> Ind2,
                 VectorXd b, 
                 VectorXd h){
                     
    const int n(X.rows()), m(h.size()), p(b.size());
    
    double Rij;
    
    VectorXd UC_b(p), UP_b(p), ebzf(n), Ha(n), Q1ij(p);
    
    UP_b.fill(0);
    
    MatrixXd ebz(m,n), hzebz(m*p,n), JC_b(p,p), JP_b(p,p);
    
    JC_b.fill(0); JP_b.fill(0);
             
    UC_b = ZFV_ * X.col(2);
    
    // the exp(b^T * ZF): n x 1 vector;
    ebzf = (ZF.adjoint() * b).array().exp();
    // the cumulative baseline hazard function evaluated at A[i]'s;
    Ha = Ind2.adjoint() * h;
    
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
            
            Q1ij.fill(0);
            // Lambda_diff = h.dot(Ind2.row(i) - Ind2.row(j));
            Rij = exp((ebzf(i) - ebzf(j)) * (Ha(i) - Ha(j)));
            // Note: the last element is still zero (derivative for bv);
            Q1ij = (ZF.col(i)*ebzf(i) - ZF.col(j)*ebzf(j)) * (Ha(i) - Ha(j));
            UP_b += -(Rij*Q1ij)/(1+Rij);
            
            JP_b += Rij/pow(1+Rij,2)*Q1ij*Q1ij.adjoint() + 
                    Rij/(1+Rij)*(Ha(i) - Ha(j))*
                    (ZF.col(i)*ZF.col(i).adjoint()*ebzf(i) - 
                     ZF.col(j)*ZF.col(j).adjoint()*ebzf(j));
            
        }

    }
    
    return (JC_b/n + JP_b/n/(n-1)).inverse() * (UC_b/n + UP_b/(n*(n-1)));
}

// [[Rcpp::export]]
VectorXd LambdaTvR1(Map<MatrixXd> ZF,
                   Map<MatrixXd> Z,
                   Map<MatrixXd> X, 
                   Map<MatrixXd> Ind1,
                   Map<MatrixXd> Ind2,
                   Map<ArrayXd> Dn,
                   VectorXd b, 
                   VectorXd h){
                    
    const int n(X.rows()), m(h.size()), p(b.size());
    VectorXd UC_Lambda(m), UP_Lambda(m), ebzf(n), Ha(n), Ind2ij(m);
    double Rij;
    
    UC_Lambda.fill(0); UP_Lambda.fill(0);
    
    MatrixXd ebz(m,n);
    
    // the exp(bf^T * ZF): n x 1 vector;
    ebzf = (ZF.adjoint() * b).array().exp();
    // the cumulative baseline hazard function evaluated at A[i]'s;
    Ha = Ind2.adjoint() * h;
    
    for(int k = 0; k < m; ++k){
        
        ebz.row(k) = (b.adjoint() * Z.middleRows(k*p,p)).array().exp();
        
    }
    
    UC_Lambda = ebz.cwiseProduct(Ind1).rowwise().sum();
    
    for(int i = 0; i < n; ++i){
                
        for(int j = 0; j < n; ++j){
            
            Ind2ij = Ind2.col(i) - Ind2.col(j);
            
            Rij = exp((ebzf(i) - ebzf(j)) * (Ha(i) - Ha(j)));
            
            UP_Lambda += (ebzf(i) - ebzf(j)) * Ind2ij / (1+1/Rij);
            
            
        }
        
    }
    
    return (Dn/(UC_Lambda/n+UP_Lambda/n/(n-1)).array()).matrix()/n;
    
}

// [[Rcpp::export]]
MatrixXd SWE_TvR1(Map<MatrixXd> ZF,
                 Map<MatrixXd> ZFV_,
                 Map<MatrixXd> Z,
                 Map<MatrixXd> X,
                 Map<ArrayXd> W,
                 Map<MatrixXd> Ind1,
                 Map<MatrixXd> Ind2,
                 VectorXd b, 
                 VectorXd h) {
                 
    const int n(X.rows()), m(h.size()), p(b.size());
    double Rij;
    VectorXd ebzf(n), Ha(n), Q1ij(p);
    MatrixXd ebz(m,n), zebz(m*p,n),
             UCs(p+m,n),  UPs(p+m,n), JC(p+m,p+m), JP(p+m,p+m), 
             V(p+m,p+m),  J(p+m,p+m), J_(p+m,p+m);        
    
    // the exp(bf^T * ZF): n x 1 vector;
    ebzf = (ZF.adjoint() * b).array().exp();
    // the cumulative baseline hazard function evaluated at A[i]'s;
    Ha = Ind2.adjoint() * h;
    
    UCs.fill(0); UPs.fill(0); JC.fill(0);  JP.fill(0);
    
    // Conditional scores correspond to betas
    UCs.topRows(p) = ZFV_.cwiseProduct(X.col(2).adjoint().replicate(p,1));
    
    for(int k = 0; k < m; ++k){
        
        ebz.row(k) = (b.adjoint() * Z.middleRows(k*p,p)).array().exp();
        zebz.middleRows(k*p,p) = Z.middleRows(k*p,p).cwiseProduct(ebz.row(k).replicate(p,1));
        
    }
    
    for(int i = 0; i < n; ++i){
        
        UCs.col(i).tail(m) = (W == X(i,1)).select(X(i,2)*h.cwiseInverse(),0).matrix();
        
        for(int k = 0; k < m; ++k){
            
            if(Ind1(k,i) == 0){
                
                continue;
                
            }else{
                
                UCs.col(i).head(p) += - zebz.block(k*p,i,p,1) * h(k);
            
                UCs(p+k,i) += - ebz(k,i);
                
                JC.topLeftCorner(p,p) += zebz.block(k*p,i,p,1) * 
                                         Z.block(k*p,i,p,1).adjoint() * h(k);
                
                JC.block(0,p+k,p,1) += zebz.block(k*p,i,p,1);
                
            }
            
            
        }
        
        JC.diagonal().tail(m) += (W == X(i,1)).select(X(i,2)*h.cwiseInverse().array().pow(2),0).matrix();    
        
        for(int j = 0; j < n; ++j){
            
            Q1ij.fill(0);
            
            Rij = exp((ebzf(i) - ebzf(j)) * (Ha(i) - Ha(j)));
            Q1ij = (ZF.col(i)*ebzf(i) - ZF.col(j)*ebzf(j)) * (Ha(i) - Ha(j)); // a column;
            
            
            UPs.col(i).head(p) += -Q1ij/(1+1/Rij);
            UPs.col(i).tail(m) += -(ebzf(i) - ebzf(j))*(Ind2.col(i) - Ind2.col(j))/(1+1/Rij);
            
            if(i < j){
                
                JP.topLeftCorner(p,p) += Rij/pow(1+Rij,2)*Q1ij*Q1ij.adjoint() + 
                                          Rij/(1+Rij)*(Ha(i) - Ha(j)) *
                                         (ZF.col(i)*ZF.col(i).adjoint()*ebzf(i) - 
                                          ZF.col(j)*ZF.col(j).adjoint()*ebzf(j));
                
                JP.topRightCorner(p,m) += (Rij*(Ha(i)-Ha(j))*(ebzf(i)-ebzf(j))/pow(1+Rij,2) + Rij/(1+Rij)) *
                                            (ZF.col(i)*ebzf(i) - ZF.col(j)*ebzf(j)) * 
                                            (Ind2.col(i) - Ind2.col(j)).adjoint();
                
                JP.bottomRightCorner(m,m) += (Rij*pow((ebzf(i)-ebzf(j))/(1+Rij),2)) * 
                                             (Ind2.col(i) - Ind2.col(j)) * 
                                             (Ind2.col(i) - Ind2.col(j)).adjoint();
                
            }
            
        }  
        
    }
    
    J = JC + 2*JP/(n-1);
    
    J.bottomLeftCorner(m,p) = J.topRightCorner(p,m).adjoint();
    
    V = UCs*UCs.adjoint() + 4*n/pow(n-1,3)*UPs*UPs.adjoint();

    J_ = J.inverse();
    
    return J_ * V * J_;
    
}


// [[Rcpp::export]]
List PLAC_TvR1(Map<MatrixXd> ZF,
             Map<MatrixXd> ZFV_,
             Map<MatrixXd> Z,
             Map<MatrixXd> X,
             Map<ArrayXd> W,
             Map<MatrixXd> Ind1,
             Map<MatrixXd> Ind2,
             Map<ArrayXd> Dn,
             VectorXd b, 
             VectorXd h,
             int K = 100){
                 
    const int n(X.rows()), m(h.size()), p(b.size());

    int k(0);
    double diff(99);
    VectorXd b_hat = b, h_hat = h, b_new, h_new, Diff(p+m);
    MatrixXd swe(p+m,p+m);
    
    while(diff > 0.0005 and k < K){
       
        h_new = LambdaTvR1(ZF, Z, X, Ind1, Ind2, Dn, b_hat, h_hat);
           
        b_new = b_hat + BetaTvR1(ZF, ZFV_, Z, X, Ind1, Ind2, b_hat, h_new);
        
        Diff << (b_new-b_hat),
                (h_new-h_hat);
        
        diff = Diff.cwiseAbs().maxCoeff();

        b_hat = b_new;
        h_hat = h_new;
        
        
        
        k++;
        
    }
    
    cout << k << " Iterations\n" << endl;
    
    swe = SWE_TvR1(ZF, ZFV_, Z, X, W, Ind1, Ind2, b_hat, h_hat);
    
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
