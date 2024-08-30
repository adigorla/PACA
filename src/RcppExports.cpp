// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_CCA
Rcpp::List cpp_CCA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo, bool normalize, int verbosity);
RcppExport SEXP _PACA_cpp_CCA(SEXP XoSEXP, SEXP YoSEXP, SEXP normalizeSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Xo(XoSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Yo(YoSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_CCA(Xo, Yo, normalize, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// cpp_selectK
Rcpp::List cpp_selectK(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, bool normalize, double threshold, int verbosity);
RcppExport SEXP _PACA_cpp_selectK(SEXP XSEXP, SEXP YSEXP, SEXP normalizeSEXP, SEXP thresholdSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_selectK(X, Y, normalize, threshold, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// cpp_PACA
Rcpp::List cpp_PACA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo, int k, bool normalize, bool retCCA, int verbosity);
RcppExport SEXP _PACA_cpp_PACA(SEXP XoSEXP, SEXP YoSEXP, SEXP kSEXP, SEXP normalizeSEXP, SEXP retCCASEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Xo(XoSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Yo(YoSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    Rcpp::traits::input_parameter< bool >::type retCCA(retCCASEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_PACA(Xo, Yo, k, normalize, retCCA, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// cpp_autoPACA
Rcpp::List cpp_autoPACA(const Eigen::MatrixXd& Xo, const Eigen::MatrixXd& Yo, bool normalize, bool retCCA, double threshold, int verbosity);
RcppExport SEXP _PACA_cpp_autoPACA(SEXP XoSEXP, SEXP YoSEXP, SEXP normalizeSEXP, SEXP retCCASEXP, SEXP thresholdSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Xo(XoSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Yo(YoSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    Rcpp::traits::input_parameter< bool >::type retCCA(retCCASEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_autoPACA(Xo, Yo, normalize, retCCA, threshold, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// cpp_prcomp
Rcpp::List cpp_prcomp(const Eigen::MatrixXd& X, bool center, bool scale, Rcpp::Nullable<int> rank, Rcpp::Nullable<double> tol, int verbosity);
RcppExport SEXP _PACA_cpp_prcomp(SEXP XSEXP, SEXP centerSEXP, SEXP scaleSEXP, SEXP rankSEXP, SEXP tolSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type center(centerSEXP);
    Rcpp::traits::input_parameter< bool >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prcomp(X, center, scale, rank, tol, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// cpp_rPACA
Rcpp::List cpp_rPACA(const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, int k, int niter, int batch, int rank, bool normalize, int verbosity);
RcppExport SEXP _PACA_cpp_rPACA(SEXP XSEXP, SEXP YSEXP, SEXP kSEXP, SEXP niterSEXP, SEXP batchSEXP, SEXP rankSEXP, SEXP normalizeSEXP, SEXP verbositySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type batch(batchSEXP);
    Rcpp::traits::input_parameter< int >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    Rcpp::traits::input_parameter< int >::type verbosity(verbositySEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_rPACA(X, Y, k, niter, batch, rank, normalize, verbosity));
    return rcpp_result_gen;
END_RCPP
}
// normalizeCPP
Eigen::MatrixXd normalizeCPP(Eigen::MatrixXd& x, bool inplace);
RcppExport SEXP _PACA_normalizeCPP(SEXP xSEXP, SEXP inplaceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type inplace(inplaceSEXP);
    rcpp_result_gen = Rcpp::wrap(normalizeCPP(x, inplace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_PACA_cpp_CCA", (DL_FUNC) &_PACA_cpp_CCA, 4},
    {"_PACA_cpp_selectK", (DL_FUNC) &_PACA_cpp_selectK, 5},
    {"_PACA_cpp_PACA", (DL_FUNC) &_PACA_cpp_PACA, 6},
    {"_PACA_cpp_autoPACA", (DL_FUNC) &_PACA_cpp_autoPACA, 6},
    {"_PACA_cpp_prcomp", (DL_FUNC) &_PACA_cpp_prcomp, 6},
    {"_PACA_cpp_rPACA", (DL_FUNC) &_PACA_cpp_rPACA, 8},
    {"_PACA_normalizeCPP", (DL_FUNC) &_PACA_normalizeCPP, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_PACA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
