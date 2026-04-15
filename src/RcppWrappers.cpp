// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
using namespace Rcpp;

#include "RcppHeaders.h"

// -----------------------------------------------------------------------------
// HIGH-LEVEL interface (user-facing)
// Handles:
//   - input validation
//   - memory layout conversion
//   - exception handling
// -----------------------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix mtgrsbed_matrix(
    std::string file,
    int n,
    const std::vector<int>& cls,
    const std::vector<double>& af,
    bool scale,
    NumericMatrix S,   // m x nt (SNP x trait)
    int nthreads = 1,
    int MG = 64,
    int JB = 1024,
    int TB = 32
) {
  const int m  = S.nrow();
  const int nt = S.ncol();

  if ((int)cls.size() != m) {
    stop("cls.size() must match number of rows in S");
  }

  if ((int)af.size() != m) {
    stop("af.size() must match number of rows in S");
  }

  // R matrix is column-major: S[i + t*m]
  const double* s = S.begin();

  // Convert to SNP-major layout: b_snp[i * nt + t]
  std::vector<double> b_snp((size_t)m * nt);

  for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
      b_snp[(size_t)i * nt + t] = s[i + (size_t)t * m];
    }
  }

  // Output: individual-major layout
  std::vector<double> grs_flat((size_t)n * nt, 0.0);

  try {
    mtgrsbed_core(
      file.c_str(),
      n,
      cls.data(),
      af.data(),
      m,
      nt,
      scale,
      b_snp.data(),
      grs_flat.data(),
      nthreads,
      MG,
      JB,
      TB
    );
  } catch (const std::exception& e) {
    stop(e.what());
  }

  // Convert to R matrix (n x nt)
  NumericMatrix out(n, nt);

  for (int t = 0; t < nt; ++t) {
    for (int j = 0; j < n; ++j) {
      out(j, t) = grs_flat[(size_t)j * nt + t];
    }
  }

  return out;
}
