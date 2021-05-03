// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// llr1d_cpp
Eigen::VectorXd llr1d_cpp(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& X_pred, int kcode, double h, Eigen::VectorXd wt);
RcppExport SEXP _llrRcpp_llr1d_cpp(SEXP XSEXP, SEXP YSEXP, SEXP X_predSEXP, SEXP kcodeSEXP, SEXP hSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type X_pred(X_predSEXP);
    Rcpp::traits::input_parameter< int >::type kcode(kcodeSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(llr1d_cpp(X, Y, X_pred, kcode, h, wt));
    return rcpp_result_gen;
END_RCPP
}
// llr2d_cpp
Eigen::VectorXd llr2d_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X_pred, int kcode, Eigen::VectorXd h, Eigen::VectorXd wt);
RcppExport SEXP _llrRcpp_llr2d_cpp(SEXP XSEXP, SEXP YSEXP, SEXP X_predSEXP, SEXP kcodeSEXP, SEXP hSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X_pred(X_predSEXP);
    Rcpp::traits::input_parameter< int >::type kcode(kcodeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(llr2d_cpp(X, Y, X_pred, kcode, h, wt));
    return rcpp_result_gen;
END_RCPP
}
// llr_cpp
Eigen::VectorXd llr_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X_pred, int kcode, Eigen::VectorXd h, const Eigen::VectorXd& wt);
RcppExport SEXP _llrRcpp_llr_cpp(SEXP XSEXP, SEXP YSEXP, SEXP X_predSEXP, SEXP kcodeSEXP, SEXP hSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X_pred(X_predSEXP);
    Rcpp::traits::input_parameter< int >::type kcode(kcodeSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type h(hSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(llr_cpp(X, Y, X_pred, kcode, h, wt));
    return rcpp_result_gen;
END_RCPP
}
// bin1d_cpp
Rcpp::List bin1d_cpp(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, int bins, const Eigen::VectorXd& wt);
RcppExport SEXP _llrRcpp_bin1d_cpp(SEXP XSEXP, SEXP YSEXP, SEXP binsSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(bin1d_cpp(X, Y, bins, wt));
    return rcpp_result_gen;
END_RCPP
}
// bin2d_cpp
Rcpp::List bin2d_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, int bins, const Eigen::VectorXd& wt);
RcppExport SEXP _llrRcpp_bin2d_cpp(SEXP XSEXP, SEXP YSEXP, SEXP binsSEXP, SEXP wtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type bins(binsSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type wt(wtSEXP);
    rcpp_result_gen = Rcpp::wrap(bin2d_cpp(X, Y, bins, wt));
    return rcpp_result_gen;
END_RCPP
}
// llrt_cpp
Eigen::VectorXd llrt_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::MatrixXd& X_pred, const Eigen::VectorXd& wt, int method, int kcode, double epsilon, const Eigen::VectorXd& h, int N_min);
RcppExport SEXP _llrRcpp_llrt_cpp(SEXP XSEXP, SEXP YSEXP, SEXP X_predSEXP, SEXP wtSEXP, SEXP methodSEXP, SEXP kcodeSEXP, SEXP epsilonSEXP, SEXP hSEXP, SEXP N_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X_pred(X_predSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type kcode(kcodeSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type N_min(N_minSEXP);
    rcpp_result_gen = Rcpp::wrap(llrt_cpp(X, Y, X_pred, wt, method, kcode, epsilon, h, N_min));
    return rcpp_result_gen;
END_RCPP
}
// tgcv_cpp
Rcpp::List tgcv_cpp(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const Eigen::VectorXd& wt, int method, int kcode, double epsilon, const Eigen::MatrixXd& bw, int N_min);
RcppExport SEXP _llrRcpp_tgcv_cpp(SEXP XSEXP, SEXP YSEXP, SEXP wtSEXP, SEXP methodSEXP, SEXP kcodeSEXP, SEXP epsilonSEXP, SEXP bwSEXP, SEXP N_minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< int >::type method(methodSEXP);
    Rcpp::traits::input_parameter< int >::type kcode(kcodeSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type N_min(N_minSEXP);
    rcpp_result_gen = Rcpp::wrap(tgcv_cpp(X, Y, wt, method, kcode, epsilon, bw, N_min));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_llrRcpp_llr1d_cpp", (DL_FUNC) &_llrRcpp_llr1d_cpp, 6},
    {"_llrRcpp_llr2d_cpp", (DL_FUNC) &_llrRcpp_llr2d_cpp, 6},
    {"_llrRcpp_llr_cpp", (DL_FUNC) &_llrRcpp_llr_cpp, 6},
    {"_llrRcpp_bin1d_cpp", (DL_FUNC) &_llrRcpp_bin1d_cpp, 4},
    {"_llrRcpp_bin2d_cpp", (DL_FUNC) &_llrRcpp_bin2d_cpp, 4},
    {"_llrRcpp_llrt_cpp", (DL_FUNC) &_llrRcpp_llrt_cpp, 9},
    {"_llrRcpp_tgcv_cpp", (DL_FUNC) &_llrRcpp_tgcv_cpp, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_llrRcpp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
