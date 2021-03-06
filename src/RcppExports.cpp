// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// parma
Rcpp::NumericVector parma(const arma::mat& X, arma::vec pik, int nthreads, double EPS);
RcppExport SEXP _tinyFF_parma(SEXP XSEXP, SEXP pikSEXP, SEXP nthreadsSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    rcpp_result_gen = Rcpp::wrap(parma(X, pik, nthreads, EPS));
    return rcpp_result_gen;
END_RCPP
}
// cube
IntegerVector cube(NumericVector prob, NumericMatrix Xbal);
RcppExport SEXP _tinyFF_cube(SEXP probSEXP, SEXP XbalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    rcpp_result_gen = Rcpp::wrap(cube(prob, Xbal));
    return rcpp_result_gen;
END_RCPP
}
// lcube
IntegerVector lcube(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal);
RcppExport SEXP _tinyFF_lcube(SEXP probSEXP, SEXP XspreadSEXP, SEXP XbalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xspread(XspreadSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    rcpp_result_gen = Rcpp::wrap(lcube(prob, Xspread, Xbal));
    return rcpp_result_gen;
END_RCPP
}
// flightphase
NumericVector flightphase(NumericVector prob, NumericMatrix Xbal);
RcppExport SEXP _tinyFF_flightphase(SEXP probSEXP, SEXP XbalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    rcpp_result_gen = Rcpp::wrap(flightphase(prob, Xbal));
    return rcpp_result_gen;
END_RCPP
}
// landingphase
IntegerVector landingphase(NumericVector prob, NumericVector probflight, NumericMatrix Xbal);
RcppExport SEXP _tinyFF_landingphase(SEXP probSEXP, SEXP probflightSEXP, SEXP XbalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probflight(probflightSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    rcpp_result_gen = Rcpp::wrap(landingphase(prob, probflight, Xbal));
    return rcpp_result_gen;
END_RCPP
}
// lcubeflightphase
NumericVector lcubeflightphase(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal);
RcppExport SEXP _tinyFF_lcubeflightphase(SEXP probSEXP, SEXP XspreadSEXP, SEXP XbalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xspread(XspreadSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    rcpp_result_gen = Rcpp::wrap(lcubeflightphase(prob, Xspread, Xbal));
    return rcpp_result_gen;
END_RCPP
}
// lcubelandingphase
IntegerVector lcubelandingphase(NumericVector prob, NumericVector probflight, NumericMatrix Xspread, NumericMatrix Xbal);
RcppExport SEXP _tinyFF_lcubelandingphase(SEXP probSEXP, SEXP probflightSEXP, SEXP XspreadSEXP, SEXP XbalSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type probflight(probflightSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xspread(XspreadSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    rcpp_result_gen = Rcpp::wrap(lcubelandingphase(prob, probflight, Xspread, Xbal));
    return rcpp_result_gen;
END_RCPP
}
// cubestratified
IntegerVector cubestratified(NumericVector prob, NumericMatrix Xbal, IntegerVector integerStrata);
RcppExport SEXP _tinyFF_cubestratified(SEXP probSEXP, SEXP XbalSEXP, SEXP integerStrataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type integerStrata(integerStrataSEXP);
    rcpp_result_gen = Rcpp::wrap(cubestratified(prob, Xbal, integerStrata));
    return rcpp_result_gen;
END_RCPP
}
// lcubestratified
IntegerVector lcubestratified(NumericVector prob, NumericMatrix Xspread, NumericMatrix Xbal, IntegerVector integerStrata);
RcppExport SEXP _tinyFF_lcubestratified(SEXP probSEXP, SEXP XspreadSEXP, SEXP XbalSEXP, SEXP integerStrataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xspread(XspreadSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Xbal(XbalSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type integerStrata(integerStrataSEXP);
    rcpp_result_gen = Rcpp::wrap(lcubestratified(prob, Xspread, Xbal, integerStrata));
    return rcpp_result_gen;
END_RCPP
}
// farma
arma::vec farma(const arma::mat& X, arma::vec pik, unsigned int nthreads, double EPS);
RcppExport SEXP _tinyFF_farma(SEXP XSEXP, SEXP pikSEXP, SEXP nthreadsSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nthreads(nthreadsSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    rcpp_result_gen = Rcpp::wrap(farma(X, pik, nthreads, EPS));
    return rcpp_result_gen;
END_RCPP
}
// flightphase_arma_2
arma::vec flightphase_arma_2(const arma::mat& X, arma::vec pik, double EPS);
RcppExport SEXP _tinyFF_flightphase_arma_2(SEXP XSEXP, SEXP pikSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    rcpp_result_gen = Rcpp::wrap(flightphase_arma_2(X, pik, EPS));
    return rcpp_result_gen;
END_RCPP
}
// parallelffphase_2
arma::vec parallelffphase_2(const arma::mat X, const arma::vec pik);
RcppExport SEXP _tinyFF_parallelffphase_2(SEXP XSEXP, SEXP pikSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type pik(pikSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelffphase_2(X, pik));
    return rcpp_result_gen;
END_RCPP
}
// qfromw
NumericMatrix qfromw(NumericVector& w, int& n);
RcppExport SEXP _tinyFF_qfromw(SEXP wSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type w(wSEXP);
    Rcpp::traits::input_parameter< int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(qfromw(w, n));
    return rcpp_result_gen;
END_RCPP
}
// onestep3
void onestep3(arma::vec& pik, arma::vec u, double EPS);
RcppExport SEXP _tinyFF_onestep3(SEXP pikSEXP, SEXP uSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type u(uSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    onestep3(pik, u, EPS);
    return R_NilValue;
END_RCPP
}
// put01
void put01(arma::uword& done, arma::uword& size, arma::uvec& index, arma::vec& pik);
RcppExport SEXP _tinyFF_put01(SEXP doneSEXP, SEXP sizeSEXP, SEXP indexSEXP, SEXP pikSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uword& >::type done(doneSEXP);
    Rcpp::traits::input_parameter< arma::uword& >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< arma::uvec& >::type index(indexSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type pik(pikSEXP);
    put01(done, size, index, pik);
    return R_NilValue;
END_RCPP
}
// flightphase_arma4
Rcpp::NumericVector flightphase_arma4(Rcpp::NumericMatrix Xr, Rcpp::NumericVector pikr, double EPS);
RcppExport SEXP _tinyFF_flightphase_arma4(SEXP XrSEXP, SEXP pikrSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Xr(XrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pikr(pikrSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    rcpp_result_gen = Rcpp::wrap(flightphase_arma4(Xr, pikr, EPS));
    return rcpp_result_gen;
END_RCPP
}
// greet
void greet();
RcppExport SEXP _tinyFF_greet() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    greet();
    return R_NilValue;
END_RCPP
}
// test2
std::vector<double> test2();
RcppExport SEXP _tinyFF_test2() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(test2());
    return rcpp_result_gen;
END_RCPP
}
// test
int test();
RcppExport SEXP _tinyFF_test() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(test());
    return rcpp_result_gen;
END_RCPP
}
// tinyFF
void tinyFF(const arma::mat& X, arma::vec& pik, double EPS);
RcppExport SEXP _tinyFF_tinyFF(SEXP XSEXP, SEXP pikSEXP, SEXP EPSSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type pik(pikSEXP);
    Rcpp::traits::input_parameter< double >::type EPS(EPSSEXP);
    tinyFF(X, pik, EPS);
    return R_NilValue;
END_RCPP
}
// parallelVectorSum
double parallelVectorSum(NumericVector x);
RcppExport SEXP _tinyFF_parallelVectorSum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelVectorSum(x));
    return rcpp_result_gen;
END_RCPP
}
// mainW
int mainW();
RcppExport SEXP _tinyFF_mainW() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(mainW());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tinyFF_parma", (DL_FUNC) &_tinyFF_parma, 4},
    {"_tinyFF_cube", (DL_FUNC) &_tinyFF_cube, 2},
    {"_tinyFF_lcube", (DL_FUNC) &_tinyFF_lcube, 3},
    {"_tinyFF_flightphase", (DL_FUNC) &_tinyFF_flightphase, 2},
    {"_tinyFF_landingphase", (DL_FUNC) &_tinyFF_landingphase, 3},
    {"_tinyFF_lcubeflightphase", (DL_FUNC) &_tinyFF_lcubeflightphase, 3},
    {"_tinyFF_lcubelandingphase", (DL_FUNC) &_tinyFF_lcubelandingphase, 4},
    {"_tinyFF_cubestratified", (DL_FUNC) &_tinyFF_cubestratified, 3},
    {"_tinyFF_lcubestratified", (DL_FUNC) &_tinyFF_lcubestratified, 4},
    {"_tinyFF_farma", (DL_FUNC) &_tinyFF_farma, 4},
    {"_tinyFF_flightphase_arma_2", (DL_FUNC) &_tinyFF_flightphase_arma_2, 3},
    {"_tinyFF_parallelffphase_2", (DL_FUNC) &_tinyFF_parallelffphase_2, 2},
    {"_tinyFF_qfromw", (DL_FUNC) &_tinyFF_qfromw, 2},
    {"_tinyFF_onestep3", (DL_FUNC) &_tinyFF_onestep3, 3},
    {"_tinyFF_put01", (DL_FUNC) &_tinyFF_put01, 4},
    {"_tinyFF_flightphase_arma4", (DL_FUNC) &_tinyFF_flightphase_arma4, 3},
    {"_tinyFF_greet", (DL_FUNC) &_tinyFF_greet, 0},
    {"_tinyFF_test2", (DL_FUNC) &_tinyFF_test2, 0},
    {"_tinyFF_test", (DL_FUNC) &_tinyFF_test, 0},
    {"_tinyFF_tinyFF", (DL_FUNC) &_tinyFF_tinyFF, 3},
    {"_tinyFF_parallelVectorSum", (DL_FUNC) &_tinyFF_parallelVectorSum, 1},
    {"_tinyFF_mainW", (DL_FUNC) &_tinyFF_mainW, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_tinyFF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
