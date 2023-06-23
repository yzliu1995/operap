// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "operap_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// SumL
double SumL(const VectorXd& v);
RcppExport SEXP _operap_SumL(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(SumL(v));
    return rcpp_result_gen;
END_RCPP
}
// getMin
VectorXd getMin(const VectorXd& v);
RcppExport SEXP _operap_getMin(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(getMin(v));
    return rcpp_result_gen;
END_RCPP
}
// quadprogSolveR
List quadprogSolveR(const MatrixXd& Dmat, const VectorXd& dvec, const MatrixXd& Amat);
RcppExport SEXP _operap_quadprogSolveR(SEXP DmatSEXP, SEXP dvecSEXP, SEXP AmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type Dmat(DmatSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type dvec(dvecSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Amat(AmatSEXP);
    rcpp_result_gen = Rcpp::wrap(quadprogSolveR(Dmat, dvec, Amat));
    return rcpp_result_gen;
END_RCPP
}
// sum_risk
double sum_risk(const MatrixXd& indset, const VectorXd& value_vec, const VectorXd& ind_vec, double p);
RcppExport SEXP _operap_sum_risk(SEXP indsetSEXP, SEXP value_vecSEXP, SEXP ind_vecSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type indset(indsetSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type value_vec(value_vecSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type ind_vec(ind_vecSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_risk(indset, value_vec, ind_vec, p));
    return rcpp_result_gen;
END_RCPP
}
// logPL_C
double logPL_C(const MatrixXd& X, const MatrixXd& indset, const VectorXd& cen, const VectorXd& y, const VectorXd& beta, const string type);
RcppExport SEXP _operap_logPL_C(SEXP XSEXP, SEXP indsetSEXP, SEXP cenSEXP, SEXP ySEXP, SEXP betaSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type indset(indsetSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(logPL_C(X, indset, cen, y, beta, type));
    return rcpp_result_gen;
END_RCPP
}
// iter_pen
VectorXd iter_pen(const MatrixXd& Dmat0, const VectorXd& dvec0, const MatrixXd& Apo, const double lambda, const VectorXd& w, const int nb, const VectorXd& inibeta, const int maxiter, const double eps, const bool trace);
RcppExport SEXP _operap_iter_pen(SEXP Dmat0SEXP, SEXP dvec0SEXP, SEXP ApoSEXP, SEXP lambdaSEXP, SEXP wSEXP, SEXP nbSEXP, SEXP inibetaSEXP, SEXP maxiterSEXP, SEXP epsSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type Dmat0(Dmat0SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type dvec0(dvec0SEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Apo(ApoSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type inibeta(inibetaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(iter_pen(Dmat0, dvec0, Apo, lambda, w, nb, inibeta, maxiter, eps, trace));
    return rcpp_result_gen;
END_RCPP
}
// val_pen
double val_pen(const MatrixXd Apo, const VectorXd& w, const int nb, const VectorXd beta);
RcppExport SEXP _operap_val_pen(SEXP ApoSEXP, SEXP wSEXP, SEXP nbSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd >::type Apo(ApoSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(val_pen(Apo, w, nb, beta));
    return rcpp_result_gen;
END_RCPP
}
// roundvec
VectorXd roundvec(const VectorXd& beta, const double eps);
RcppExport SEXP _operap_roundvec(SEXP betaSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(roundvec(beta, eps));
    return rcpp_result_gen;
END_RCPP
}
// lasso_tree_single
VectorXd lasso_tree_single(const MatrixXd& X, const MatrixXd& indset, const MatrixXd& Apo, const int nb, const VectorXd& w, const VectorXd& cen, const VectorXd& y, const double lambda, const VectorXd& inibeta, const int maxiter, const double eps, const bool trace, const bool fullA, const string type);
RcppExport SEXP _operap_lasso_tree_single(SEXP XSEXP, SEXP indsetSEXP, SEXP ApoSEXP, SEXP nbSEXP, SEXP wSEXP, SEXP cenSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP inibetaSEXP, SEXP maxiterSEXP, SEXP epsSEXP, SEXP traceSEXP, SEXP fullASEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type indset(indsetSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Apo(ApoSEXP);
    Rcpp::traits::input_parameter< const int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type inibeta(inibetaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const bool >::type fullA(fullASEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(lasso_tree_single(X, indset, Apo, nb, w, cen, y, lambda, inibeta, maxiter, eps, trace, fullA, type));
    return rcpp_result_gen;
END_RCPP
}
// lasso_tree_multi
MatrixXd lasso_tree_multi(const MatrixXd& X, const MatrixXd& indset, const MatrixXd& Apo, const int nb, const VectorXd& w, const VectorXd& cen, const VectorXd& y, const VectorXd& lambda, const VectorXd& inibeta, const int maxiter, const double eps, const bool trace, const string type, const bool fullA);
RcppExport SEXP _operap_lasso_tree_multi(SEXP XSEXP, SEXP indsetSEXP, SEXP ApoSEXP, SEXP nbSEXP, SEXP wSEXP, SEXP cenSEXP, SEXP ySEXP, SEXP lambdaSEXP, SEXP inibetaSEXP, SEXP maxiterSEXP, SEXP epsSEXP, SEXP traceSEXP, SEXP typeSEXP, SEXP fullASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type indset(indsetSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Apo(ApoSEXP);
    Rcpp::traits::input_parameter< const int >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type inibeta(inibetaSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const bool >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type fullA(fullASEXP);
    rcpp_result_gen = Rcpp::wrap(lasso_tree_multi(X, indset, Apo, nb, w, cen, y, lambda, inibeta, maxiter, eps, trace, type, fullA));
    return rcpp_result_gen;
END_RCPP
}
// quadprog_solveR
List quadprog_solveR(const MatrixXd& Dmat, const VectorXd& dvec, const MatrixXd& Amat, const VectorXd& bvec);
RcppExport SEXP _operap_quadprog_solveR(SEXP DmatSEXP, SEXP dvecSEXP, SEXP AmatSEXP, SEXP bvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type Dmat(DmatSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type dvec(dvecSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Amat(AmatSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type bvec(bvecSEXP);
    rcpp_result_gen = Rcpp::wrap(quadprog_solveR(Dmat, dvec, Amat, bvec));
    return rcpp_result_gen;
END_RCPP
}
// ED
List ED(const MatrixXd& mat);
RcppExport SEXP _operap_ED(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(ED(mat));
    return rcpp_result_gen;
END_RCPP
}
// nNv
VectorXd nNv(const VectorXd& v);
RcppExport SEXP _operap_nNv(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(nNv(v));
    return rcpp_result_gen;
END_RCPP
}
// CnearPD
MatrixXd CnearPD(const MatrixXd& mat);
RcppExport SEXP _operap_CnearPD(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(CnearPD(mat));
    return rcpp_result_gen;
END_RCPP
}
// Sum
double Sum(const VectorXd& v);
RcppExport SEXP _operap_Sum(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(Sum(v));
    return rcpp_result_gen;
END_RCPP
}
// Min
double Min(const VectorXd& v);
RcppExport SEXP _operap_Min(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(Min(v));
    return rcpp_result_gen;
END_RCPP
}
// wInf
bool wInf(const VectorXd& v);
RcppExport SEXP _operap_wInf(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(wInf(v));
    return rcpp_result_gen;
END_RCPP
}
// Round
VectorXd Round(const VectorXd& v, double eps);
RcppExport SEXP _operap_Round(SEXP vSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(Round(v, eps));
    return rcpp_result_gen;
END_RCPP
}
// countDistinct
int countDistinct(const VectorXd& v);
RcppExport SEXP _operap_countDistinct(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(countDistinct(v));
    return rcpp_result_gen;
END_RCPP
}
// cumSum
VectorXd cumSum(const VectorXd& v);
RcppExport SEXP _operap_cumSum(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(cumSum(v));
    return rcpp_result_gen;
END_RCPP
}
// revCumSum
VectorXd revCumSum(const VectorXd& v);
RcppExport SEXP _operap_revCumSum(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(revCumSum(v));
    return rcpp_result_gen;
END_RCPP
}
// gs
VectorXd gs(const double& a, const double& b, const int& c);
RcppExport SEXP _operap_gs(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int& >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(gs(a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// logLK
double logLK(const MatrixXd& Z, const VectorXd& cen, const VectorXd& beta, const VectorXd& y, const MatrixXd& cov, const VectorXd& theta, const std::string type, const bool withCov);
RcppExport SEXP _operap_logLK(SEXP ZSEXP, SEXP cenSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP covSEXP, SEXP thetaSEXP, SEXP typeSEXP, SEXP withCovSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type withCov(withCovSEXP);
    rcpp_result_gen = Rcpp::wrap(logLK(Z, cen, beta, y, cov, theta, type, withCov));
    return rcpp_result_gen;
END_RCPP
}
// derivatives
List derivatives(const MatrixXd& Z, const VectorXd& cen, const VectorXd& beta, const VectorXd& y, const MatrixXd& cov, const VectorXd& theta, const int seed, const bool withCov, const std::string type);
RcppExport SEXP _operap_derivatives(SEXP ZSEXP, SEXP cenSEXP, SEXP betaSEXP, SEXP ySEXP, SEXP covSEXP, SEXP thetaSEXP, SEXP seedSEXP, SEXP withCovSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool >::type withCov(withCovSEXP);
    Rcpp::traits::input_parameter< const std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(derivatives(Z, cen, beta, y, cov, theta, seed, withCov, type));
    return rcpp_result_gen;
END_RCPP
}
// IRLSO
VectorXd IRLSO(const VectorXd& coeff0, const VectorXd& theta0, const VectorXd& cen, const VectorXd& y, const MatrixXd& Z, const MatrixXd& cov, const MatrixXd& Cnstrn, const VectorXd& coeffLambda, const int& seed, const bool& withCov, const string type);
RcppExport SEXP _operap_IRLSO(SEXP coeff0SEXP, SEXP theta0SEXP, SEXP cenSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP covSEXP, SEXP CnstrnSEXP, SEXP coeffLambdaSEXP, SEXP seedSEXP, SEXP withCovSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type coeff0(coeff0SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Cnstrn(CnstrnSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type coeffLambda(coeffLambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type withCov(withCovSEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(IRLSO(coeff0, theta0, cen, y, Z, cov, Cnstrn, coeffLambda, seed, withCov, type));
    return rcpp_result_gen;
END_RCPP
}
// IRLSA
List IRLSA(const VectorXd& coeff, const VectorXd& theta, const VectorXd& cen, const VectorXd& y, const MatrixXd& Z, const MatrixXd& cov, const MatrixXd& Cnstrn, const VectorXd& coeffLambda, const double& eps, const int& maxiter, const int& seed, const bool& withCov, const string type);
RcppExport SEXP _operap_IRLSA(SEXP coeffSEXP, SEXP thetaSEXP, SEXP cenSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP covSEXP, SEXP CnstrnSEXP, SEXP coeffLambdaSEXP, SEXP epsSEXP, SEXP maxiterSEXP, SEXP seedSEXP, SEXP withCovSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type coeff(coeffSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Cnstrn(CnstrnSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type coeffLambda(coeffLambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const bool& >::type withCov(withCovSEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(IRLSA(coeff, theta, cen, y, Z, cov, Cnstrn, coeffLambda, eps, maxiter, seed, withCov, type));
    return rcpp_result_gen;
END_RCPP
}
// IRLSP
VectorXd IRLSP(const VectorXd& coeff0, const VectorXd& theta0, const VectorXd& cen, const VectorXd& y, const MatrixXd& Z, const MatrixXd& cov, const MatrixXd& LCnstrn, const MatrixXd& qCnstrn, const int& mCnstrn, const double& minCoeff, const double& maxCoeff, const double& Lambda, const bool& withCov, const int& seed, const string type);
RcppExport SEXP _operap_IRLSP(SEXP coeff0SEXP, SEXP theta0SEXP, SEXP cenSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP covSEXP, SEXP LCnstrnSEXP, SEXP qCnstrnSEXP, SEXP mCnstrnSEXP, SEXP minCoeffSEXP, SEXP maxCoeffSEXP, SEXP LambdaSEXP, SEXP withCovSEXP, SEXP seedSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type coeff0(coeff0SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type theta0(theta0SEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type LCnstrn(LCnstrnSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type qCnstrn(qCnstrnSEXP);
    Rcpp::traits::input_parameter< const int& >::type mCnstrn(mCnstrnSEXP);
    Rcpp::traits::input_parameter< const double& >::type minCoeff(minCoeffSEXP);
    Rcpp::traits::input_parameter< const double& >::type maxCoeff(maxCoeffSEXP);
    Rcpp::traits::input_parameter< const double& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type withCov(withCovSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(IRLSP(coeff0, theta0, cen, y, Z, cov, LCnstrn, qCnstrn, mCnstrn, minCoeff, maxCoeff, Lambda, withCov, seed, type));
    return rcpp_result_gen;
END_RCPP
}
// IRLSPT
VectorXd IRLSPT(const VectorXd& coeff, const VectorXd& theta, const VectorXd& cen, const VectorXd& y, const MatrixXd& Z, const MatrixXd& cov, const MatrixXd& LCnstrn, const MatrixXd& qCnstrn, const int& mCnstrn, const double& minCoeff, const double& maxCoeff, const double& eps, const int& maxiter, const bool& withCov, const int& seed, const string type);
RcppExport SEXP _operap_IRLSPT(SEXP coeffSEXP, SEXP thetaSEXP, SEXP cenSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP covSEXP, SEXP LCnstrnSEXP, SEXP qCnstrnSEXP, SEXP mCnstrnSEXP, SEXP minCoeffSEXP, SEXP maxCoeffSEXP, SEXP epsSEXP, SEXP maxiterSEXP, SEXP withCovSEXP, SEXP seedSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type coeff(coeffSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type LCnstrn(LCnstrnSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type qCnstrn(qCnstrnSEXP);
    Rcpp::traits::input_parameter< const int& >::type mCnstrn(mCnstrnSEXP);
    Rcpp::traits::input_parameter< const double& >::type minCoeff(minCoeffSEXP);
    Rcpp::traits::input_parameter< const double& >::type maxCoeff(maxCoeffSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type withCov(withCovSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(IRLSPT(coeff, theta, cen, y, Z, cov, LCnstrn, qCnstrn, mCnstrn, minCoeff, maxCoeff, eps, maxiter, withCov, seed, type));
    return rcpp_result_gen;
END_RCPP
}
// IRLSPAT
List IRLSPAT(const VectorXd& coeff, const VectorXd& theta, const VectorXd& cen, const VectorXd& y, const MatrixXd& Z, const MatrixXd& cov, const MatrixXd& LCnstrn, const MatrixXd& qCnstrn, const int& mCnstrn, const double& minCoeff, const double& maxCoeff, const double& eps, const int& maxiter, const bool& withCov, const int& seed, const string type);
RcppExport SEXP _operap_IRLSPAT(SEXP coeffSEXP, SEXP thetaSEXP, SEXP cenSEXP, SEXP ySEXP, SEXP ZSEXP, SEXP covSEXP, SEXP LCnstrnSEXP, SEXP qCnstrnSEXP, SEXP mCnstrnSEXP, SEXP minCoeffSEXP, SEXP maxCoeffSEXP, SEXP epsSEXP, SEXP maxiterSEXP, SEXP withCovSEXP, SEXP seedSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd& >::type coeff(coeffSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type cen(cenSEXP);
    Rcpp::traits::input_parameter< const VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type cov(covSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type LCnstrn(LCnstrnSEXP);
    Rcpp::traits::input_parameter< const MatrixXd& >::type qCnstrn(qCnstrnSEXP);
    Rcpp::traits::input_parameter< const int& >::type mCnstrn(mCnstrnSEXP);
    Rcpp::traits::input_parameter< const double& >::type minCoeff(minCoeffSEXP);
    Rcpp::traits::input_parameter< const double& >::type maxCoeff(maxCoeffSEXP);
    Rcpp::traits::input_parameter< const double& >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const int& >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const bool& >::type withCov(withCovSEXP);
    Rcpp::traits::input_parameter< const int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< const string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(IRLSPAT(coeff, theta, cen, y, Z, cov, LCnstrn, qCnstrn, mCnstrn, minCoeff, maxCoeff, eps, maxiter, withCov, seed, type));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_operap_SumL", (DL_FUNC) &_operap_SumL, 1},
    {"_operap_getMin", (DL_FUNC) &_operap_getMin, 1},
    {"_operap_quadprogSolveR", (DL_FUNC) &_operap_quadprogSolveR, 3},
    {"_operap_sum_risk", (DL_FUNC) &_operap_sum_risk, 4},
    {"_operap_logPL_C", (DL_FUNC) &_operap_logPL_C, 6},
    {"_operap_iter_pen", (DL_FUNC) &_operap_iter_pen, 10},
    {"_operap_val_pen", (DL_FUNC) &_operap_val_pen, 4},
    {"_operap_roundvec", (DL_FUNC) &_operap_roundvec, 2},
    {"_operap_lasso_tree_single", (DL_FUNC) &_operap_lasso_tree_single, 14},
    {"_operap_lasso_tree_multi", (DL_FUNC) &_operap_lasso_tree_multi, 14},
    {"_operap_quadprog_solveR", (DL_FUNC) &_operap_quadprog_solveR, 4},
    {"_operap_ED", (DL_FUNC) &_operap_ED, 1},
    {"_operap_nNv", (DL_FUNC) &_operap_nNv, 1},
    {"_operap_CnearPD", (DL_FUNC) &_operap_CnearPD, 1},
    {"_operap_Sum", (DL_FUNC) &_operap_Sum, 1},
    {"_operap_Min", (DL_FUNC) &_operap_Min, 1},
    {"_operap_wInf", (DL_FUNC) &_operap_wInf, 1},
    {"_operap_Round", (DL_FUNC) &_operap_Round, 2},
    {"_operap_countDistinct", (DL_FUNC) &_operap_countDistinct, 1},
    {"_operap_cumSum", (DL_FUNC) &_operap_cumSum, 1},
    {"_operap_revCumSum", (DL_FUNC) &_operap_revCumSum, 1},
    {"_operap_gs", (DL_FUNC) &_operap_gs, 3},
    {"_operap_logLK", (DL_FUNC) &_operap_logLK, 8},
    {"_operap_derivatives", (DL_FUNC) &_operap_derivatives, 9},
    {"_operap_IRLSO", (DL_FUNC) &_operap_IRLSO, 11},
    {"_operap_IRLSA", (DL_FUNC) &_operap_IRLSA, 13},
    {"_operap_IRLSP", (DL_FUNC) &_operap_IRLSP, 15},
    {"_operap_IRLSPT", (DL_FUNC) &_operap_IRLSPT, 16},
    {"_operap_IRLSPAT", (DL_FUNC) &_operap_IRLSPAT, 16},
    {NULL, NULL, 0}
};

RcppExport void R_init_operap(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
