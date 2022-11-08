// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _BDS_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpp_backwards
NumericMatrix rcpp_backwards(NumericVector lambda, NumericVector mu, double eta, NumericVector u0, double bl, int nsteps);
RcppExport SEXP _BDS_rcpp_backwards(SEXP lambdaSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP u0SEXP, SEXP blSEXP, SEXP nstepsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u0(u0SEXP);
    Rcpp::traits::input_parameter< double >::type bl(blSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_backwards(lambda, mu, eta, u0, bl, nsteps));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_get_descendants
IntegerVector rcpp_get_descendants(IntegerMatrix edge, int node_idx);
RcppExport SEXP _BDS_rcpp_get_descendants(SEXP edgeSEXP, SEXP node_idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< int >::type node_idx(node_idxSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_get_descendants(edge, node_idx));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_postorder
int rcpp_postorder(NumericVector lambda, NumericVector mu, double eta, IntegerVector po, IntegerMatrix edge, NumericVector branch_lengths, int rootnode);
RcppExport SEXP _BDS_rcpp_postorder(SEXP lambdaSEXP, SEXP muSEXP, SEXP etaSEXP, SEXP poSEXP, SEXP edgeSEXP, SEXP branch_lengthsSEXP, SEXP rootnodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type po(poSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type edge(edgeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type branch_lengths(branch_lengthsSEXP);
    Rcpp::traits::input_parameter< int >::type rootnode(rootnodeSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_postorder(lambda, mu, eta, po, edge, branch_lengths, rootnode));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BDS_rcpp_hello_world", (DL_FUNC) &_BDS_rcpp_hello_world, 0},
    {"_BDS_rcpp_backwards", (DL_FUNC) &_BDS_rcpp_backwards, 6},
    {"_BDS_rcpp_get_descendants", (DL_FUNC) &_BDS_rcpp_get_descendants, 2},
    {"_BDS_rcpp_postorder", (DL_FUNC) &_BDS_rcpp_postorder, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_BDS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}