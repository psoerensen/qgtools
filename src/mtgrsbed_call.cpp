#include <R.h>
#include <Rinternals.h>
#include <vector>

#include "mtgrsbed_core.h"
#include "mtgrsbed_core_f.h"

static void build_b_snp(const double* S, int m, int nt, std::vector<double>& b_snp) {
  for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
      b_snp[(size_t)i * nt + t] = S[i + (size_t)t * m];
    }
  }
}

static SEXP make_output(const std::vector<double>& grs_flat, int n, int nt) {
  SEXP out = PROTECT(allocMatrix(REALSXP, n, nt));
  double* out_ptr = REAL(out);

  for (int t = 0; t < nt; ++t) {
    for (int j = 0; j < n; ++j) {
      out_ptr[j + (size_t)t * n] = grs_flat[(size_t)j * nt + t];
    }
  }

  UNPROTECT(1);
  return out;
}

extern "C" SEXP mtgrsbed_c(SEXP fileSEXP,
                          SEXP nSEXP,
                          SEXP clsSEXP,
                          SEXP afSEXP,
                          SEXP scaleSEXP,
                          SEXP SSEXP,
                          SEXP nthreadsSEXP) {
  const char* file = CHAR(STRING_ELT(fileSEXP, 0));
  int n = INTEGER(nSEXP)[0];
  int* cls = INTEGER(clsSEXP);
  double* af = REAL(afSEXP);
  int scale = LOGICAL(scaleSEXP)[0];
  int nthreads = INTEGER(nthreadsSEXP)[0];

  SEXP dim = getAttrib(SSEXP, R_DimSymbol);
  int m  = INTEGER(dim)[0];
  int nt = INTEGER(dim)[1];
  const double* S = REAL(SSEXP);

  std::vector<double> b_snp((size_t)m * nt);
  std::vector<double> grs_flat((size_t)n * nt, 0.0);

  build_b_snp(S, m, nt, b_snp);

  mtgrsbed_core(
    file, n, cls, af, m, nt,
    scale != 0,
    b_snp.data(), grs_flat.data(),
    nthreads, 64, 1024, 32
  );

  return make_output(grs_flat, n, nt);
}

extern "C" SEXP mtgrsbed_f(SEXP fileSEXP,
                          SEXP nSEXP,
                          SEXP clsSEXP,
                          SEXP afSEXP,
                          SEXP scaleSEXP,
                          SEXP SSEXP,
                          SEXP nthreadsSEXP) {
  const char* file = CHAR(STRING_ELT(fileSEXP, 0));
  int n = INTEGER(nSEXP)[0];
  int* cls = INTEGER(clsSEXP);
  double* af = REAL(afSEXP);
  int scale = LOGICAL(scaleSEXP)[0];
  int nthreads = INTEGER(nthreadsSEXP)[0];

  SEXP dim = getAttrib(SSEXP, R_DimSymbol);
  int m  = INTEGER(dim)[0];
  int nt = INTEGER(dim)[1];
  const double* S = REAL(SSEXP);

  std::vector<double> b_snp((size_t)m * nt);
  std::vector<double> grs_flat((size_t)n * nt, 0.0);

  build_b_snp(S, m, nt, b_snp);

  mtgrsbed_core_f(
    file, n, cls, af, m, nt,
    scale != 0,
    b_snp.data(), grs_flat.data(),
    nthreads, 64, 1024, 32
  );

  return make_output(grs_flat, n, nt);
}
