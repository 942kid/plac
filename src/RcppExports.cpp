// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// SgInd1
Eigen::MatrixXd SgInd1(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::ArrayXd> W);
RcppExport SEXP plac_SgInd1(SEXP XSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    __result = Rcpp::wrap(SgInd1(X, W));
    return __result;
END_RCPP
}
// PwInd1
Eigen::MatrixXd PwInd1(Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::ArrayXd> W);
RcppExport SEXP plac_PwInd1(SEXP XSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    __result = Rcpp::wrap(PwInd1(X, W));
    return __result;
END_RCPP
}
// TvInd1
Eigen::MatrixXd TvInd1(Eigen::Map<Eigen::VectorXd> Zv, Eigen::Map<Eigen::ArrayXd> W);
RcppExport SEXP plac_TvInd1(SEXP ZvSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Zv(ZvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    __result = Rcpp::wrap(TvInd1(Zv, W));
    return __result;
END_RCPP
}
// LambdaTv2
Eigen::VectorXd LambdaTv2(Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Ind1, Eigen::Map<Eigen::MatrixXd> Ind2, Eigen::Map<Eigen::ArrayXd> Dn, Eigen::VectorXd b, Eigen::VectorXd h);
RcppExport SEXP plac_LambdaTv2(SEXP ZSEXP, SEXP XSEXP, SEXP Ind1SEXP, SEXP Ind2SEXP, SEXP DnSEXP, SEXP bSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind1(Ind1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind2(Ind2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type Dn(DnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    __result = Rcpp::wrap(LambdaTv2(Z, X, Ind1, Ind2, Dn, b, h));
    return __result;
END_RCPP
}
// SWE_Tv2
Eigen::MatrixXd SWE_Tv2(Eigen::Map<Eigen::MatrixXd> ZFV_, Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::ArrayXd> W, Eigen::Map<Eigen::MatrixXd> Ind1, Eigen::Map<Eigen::MatrixXd> Ind2, Eigen::VectorXd b, Eigen::VectorXd h);
RcppExport SEXP plac_SWE_Tv2(SEXP ZFV_SEXP, SEXP ZSEXP, SEXP XSEXP, SEXP WSEXP, SEXP Ind1SEXP, SEXP Ind2SEXP, SEXP bSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZFV_(ZFV_SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind1(Ind1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind2(Ind2SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    __result = Rcpp::wrap(SWE_Tv2(ZFV_, Z, X, W, Ind1, Ind2, b, h));
    return __result;
END_RCPP
}
// PLAC_Tv3
List PLAC_Tv3(Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::MatrixXd> ZFV_, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::ArrayXd> W, Eigen::Map<Eigen::MatrixXd> Ind1, Eigen::Map<Eigen::MatrixXd> Ind2, Eigen::Map<Eigen::ArrayXd> Dn, Eigen::VectorXd b, Eigen::VectorXd h, int K);
RcppExport SEXP plac_PLAC_Tv3(SEXP ZSEXP, SEXP ZFV_SEXP, SEXP XSEXP, SEXP WSEXP, SEXP Ind1SEXP, SEXP Ind2SEXP, SEXP DnSEXP, SEXP bSEXP, SEXP hSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZFV_(ZFV_SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind1(Ind1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind2(Ind2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type Dn(DnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    __result = Rcpp::wrap(PLAC_Tv3(Z, ZFV_, X, W, Ind1, Ind2, Dn, b, h, K));
    return __result;
END_RCPP
}
// TvInd2
Eigen::MatrixXd TvInd2(Eigen::Map<Eigen::VectorXd> PD, Eigen::Map<Eigen::VectorXd> Zv, Eigen::Map<Eigen::ArrayXd> W);
RcppExport SEXP plac_TvInd2(SEXP PDSEXP, SEXP ZvSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type PD(PDSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Zv(ZvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    __result = Rcpp::wrap(TvInd2(PD, Zv, W));
    return __result;
END_RCPP
}
// TvInd3
Eigen::MatrixXd TvInd3(Eigen::Map<Eigen::VectorXd> PD, Eigen::Map<Eigen::VectorXd> Zv, Eigen::Map<Eigen::ArrayXd> W);
RcppExport SEXP plac_TvInd3(SEXP PDSEXP, SEXP ZvSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type PD(PDSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type Zv(ZvSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    __result = Rcpp::wrap(TvInd3(PD, Zv, W));
    return __result;
END_RCPP
}
// BetaTvR1
Eigen::MatrixXd BetaTvR1(Eigen::Map<Eigen::MatrixXd> ZF, Eigen::Map<Eigen::MatrixXd> ZFV_, Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Ind1, Eigen::Map<Eigen::MatrixXd> Ind2, Eigen::VectorXd b, Eigen::VectorXd h);
RcppExport SEXP plac_BetaTvR1(SEXP ZFSEXP, SEXP ZFV_SEXP, SEXP ZSEXP, SEXP XSEXP, SEXP Ind1SEXP, SEXP Ind2SEXP, SEXP bSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZF(ZFSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZFV_(ZFV_SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind1(Ind1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind2(Ind2SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    __result = Rcpp::wrap(BetaTvR1(ZF, ZFV_, Z, X, Ind1, Ind2, b, h));
    return __result;
END_RCPP
}
// LambdaTvR1
Eigen::VectorXd LambdaTvR1(Eigen::Map<Eigen::MatrixXd> ZF, Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::MatrixXd> Ind1, Eigen::Map<Eigen::MatrixXd> Ind2, Eigen::Map<Eigen::ArrayXd> Dn, Eigen::VectorXd b, Eigen::VectorXd h);
RcppExport SEXP plac_LambdaTvR1(SEXP ZFSEXP, SEXP ZSEXP, SEXP XSEXP, SEXP Ind1SEXP, SEXP Ind2SEXP, SEXP DnSEXP, SEXP bSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZF(ZFSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind1(Ind1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind2(Ind2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type Dn(DnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    __result = Rcpp::wrap(LambdaTvR1(ZF, Z, X, Ind1, Ind2, Dn, b, h));
    return __result;
END_RCPP
}
// SWE_TvR1
Eigen::MatrixXd SWE_TvR1(Eigen::Map<Eigen::MatrixXd> ZF, Eigen::Map<Eigen::MatrixXd> ZFV_, Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::ArrayXd> W, Eigen::Map<Eigen::MatrixXd> Ind1, Eigen::Map<Eigen::MatrixXd> Ind2, Eigen::VectorXd b, Eigen::VectorXd h);
RcppExport SEXP plac_SWE_TvR1(SEXP ZFSEXP, SEXP ZFV_SEXP, SEXP ZSEXP, SEXP XSEXP, SEXP WSEXP, SEXP Ind1SEXP, SEXP Ind2SEXP, SEXP bSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZF(ZFSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZFV_(ZFV_SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind1(Ind1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind2(Ind2SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    __result = Rcpp::wrap(SWE_TvR1(ZF, ZFV_, Z, X, W, Ind1, Ind2, b, h));
    return __result;
END_RCPP
}
// PLAC_TvR1
List PLAC_TvR1(Eigen::Map<Eigen::MatrixXd> ZF, Eigen::Map<Eigen::MatrixXd> ZFV_, Eigen::Map<Eigen::MatrixXd> Z, Eigen::Map<Eigen::MatrixXd> X, Eigen::Map<Eigen::ArrayXd> W, Eigen::Map<Eigen::MatrixXd> Ind1, Eigen::Map<Eigen::MatrixXd> Ind2, Eigen::Map<Eigen::ArrayXd> Dn, Eigen::VectorXd b, Eigen::VectorXd h, int K);
RcppExport SEXP plac_PLAC_TvR1(SEXP ZFSEXP, SEXP ZFV_SEXP, SEXP ZSEXP, SEXP XSEXP, SEXP WSEXP, SEXP Ind1SEXP, SEXP Ind2SEXP, SEXP DnSEXP, SEXP bSEXP, SEXP hSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZF(ZFSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type ZFV_(ZFV_SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type X(XSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type W(WSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind1(Ind1SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Ind2(Ind2SEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::ArrayXd> >::type Dn(DnSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    __result = Rcpp::wrap(PLAC_TvR1(ZF, ZFV_, Z, X, W, Ind1, Ind2, Dn, b, h, K));
    return __result;
END_RCPP
}