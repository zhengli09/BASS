// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// Gibbs
Rcpp::List Gibbs(const Rcpp::NumericMatrix X, Rcpp::NumericMatrix xy, const Rcpp::IntegerVector Ns, const int C, const int R, Rcpp::String initMethod, const double kappa, const double alpha0, const double a, const double b, const int k, const int warmUp, const int numSamples, Rcpp::String betaEstApproach, const double betaIn, const double betaMax, const double epsilon, const int M, const int B, const Rcpp::NumericMatrix NHC);
RcppExport SEXP _BASS_Gibbs(SEXP XSEXP, SEXP xySEXP, SEXP NsSEXP, SEXP CSEXP, SEXP RSEXP, SEXP initMethodSEXP, SEXP kappaSEXP, SEXP alpha0SEXP, SEXP aSEXP, SEXP bSEXP, SEXP kSEXP, SEXP warmUpSEXP, SEXP numSamplesSEXP, SEXP betaEstApproachSEXP, SEXP betaInSEXP, SEXP betaMaxSEXP, SEXP epsilonSEXP, SEXP MSEXP, SEXP BSEXP, SEXP NHCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type xy(xySEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< const int >::type C(CSEXP);
    Rcpp::traits::input_parameter< const int >::type R(RSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type initMethod(initMethodSEXP);
    Rcpp::traits::input_parameter< const double >::type kappa(kappaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    Rcpp::traits::input_parameter< const int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const int >::type warmUp(warmUpSEXP);
    Rcpp::traits::input_parameter< const int >::type numSamples(numSamplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::String >::type betaEstApproach(betaEstApproachSEXP);
    Rcpp::traits::input_parameter< const double >::type betaIn(betaInSEXP);
    Rcpp::traits::input_parameter< const double >::type betaMax(betaMaxSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type NHC(NHCSEXP);
    rcpp_result_gen = Rcpp::wrap(Gibbs(X, xy, Ns, C, R, initMethod, kappa, alpha0, a, b, k, warmUp, numSamples, betaEstApproach, betaIn, betaMax, epsilon, M, B, NHC));
    return rcpp_result_gen;
END_RCPP
}
// avgNeighbors
double avgNeighbors(const double r, const Rcpp::NumericMatrix& coord);
RcppExport SEXP _BASS_avgNeighbors(SEXP rSEXP, SEXP coordSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type coord(coordSEXP);
    rcpp_result_gen = Rcpp::wrap(avgNeighbors(r, coord));
    return rcpp_result_gen;
END_RCPP
}
// testFastBetaEst
Rcpp::List testFastBetaEst(const Rcpp::IntegerVector x, const Rcpp::NumericMatrix coord, const double r, double beta, const double epsilon, const int S, const Rcpp::NumericMatrix NHC);
RcppExport SEXP _BASS_testFastBetaEst(SEXP xSEXP, SEXP coordSEXP, SEXP rSEXP, SEXP betaSEXP, SEXP epsilonSEXP, SEXP SSEXP, SEXP NHCSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const int >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type NHC(NHCSEXP);
    rcpp_result_gen = Rcpp::wrap(testFastBetaEst(x, coord, r, beta, epsilon, S, NHC));
    return rcpp_result_gen;
END_RCPP
}
// testPottsSampling
Rcpp::List testPottsSampling(const Rcpp::NumericMatrix coord, const double r, const double beta, const int R, const int M);
RcppExport SEXP _BASS_testPottsSampling(SEXP coordSEXP, SEXP rSEXP, SEXP betaSEXP, SEXP RSEXP, SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< const double >::type r(rSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const int >::type R(RSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(testPottsSampling(coord, r, beta, R, M));
    return rcpp_result_gen;
END_RCPP
}
// testPostBeta
Rcpp::List testPostBeta(const Rcpp::IntegerVector x, const Rcpp::NumericMatrix coord, const int R, double beta, const double epsilon, const int M, const int S);
RcppExport SEXP _BASS_testPostBeta(SEXP xSEXP, SEXP coordSEXP, SEXP RSEXP, SEXP betaSEXP, SEXP epsilonSEXP, SEXP MSEXP, SEXP SSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix >::type coord(coordSEXP);
    Rcpp::traits::input_parameter< const int >::type R(RSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const int >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int >::type S(SSEXP);
    rcpp_result_gen = Rcpp::wrap(testPostBeta(x, coord, R, beta, epsilon, M, S));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BASS_Gibbs", (DL_FUNC) &_BASS_Gibbs, 20},
    {"_BASS_avgNeighbors", (DL_FUNC) &_BASS_avgNeighbors, 2},
    {"_BASS_testFastBetaEst", (DL_FUNC) &_BASS_testFastBetaEst, 7},
    {"_BASS_testPottsSampling", (DL_FUNC) &_BASS_testPottsSampling, 5},
    {"_BASS_testPostBeta", (DL_FUNC) &_BASS_testPostBeta, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_BASS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}