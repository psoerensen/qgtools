#include <Rcpp.h>
using namespace Rcpp;

#include "RcppHeaders.h"

// -----------------------------------------------------------------------------
// LOW-LEVEL interface using .Call()
// Demonstrates how R interacts with compiled code directly
// -----------------------------------------------------------------------------

extern "C" SEXP mtgrsbed_core_R(
    SEXP fileSEXP,
    SEXP nSEXP,
    SEXP clsSEXP,
    SEXP afSEXP,
    SEXP mSEXP,
    SEXP ntSEXP,
    SEXP scaleSEXP,
    SEXP b_snpSEXP,
    SEXP grsSEXP,
    SEXP nthreadsSEXP,
    SEXP MGSEXP,
    SEXP JBSEXP,
    SEXP TBSEXP
) {
  const char* file = CHAR(STRING_ELT(fileSEXP, 0));
  int n = INTEGER(nSEXP)[0];
  int* cls = INTEGER(clsSEXP);
  double* af = REAL(afSEXP);
  int m = INTEGER(mSEXP)[0];
  int nt = INTEGER(ntSEXP)[0];
  bool scale = LOGICAL(scaleSEXP)[0];
  double* b_snp = REAL(b_snpSEXP);
  double* grs = REAL(grsSEXP);
  int nthreads = INTEGER(nthreadsSEXP)[0];
  int MG = INTEGER(MGSEXP)[0];
  int JB = INTEGER(JBSEXP)[0];
  int TB = INTEGER(TBSEXP)[0];

  try {
    mtgrsbed_core(file, n, cls, af, m, nt, scale,
                  b_snp, grs, nthreads, MG, JB, TB);
  } catch (const std::exception& e) {
    Rf_error(e.what());
  }

  return grsSEXP;
}
