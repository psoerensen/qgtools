// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
  void mtgrsbed_core_f(
      const char* file,
      int n,
      const int* cls,
      const double* af,
      int m,
      int nt,
      bool scale,
      const double* b_snp,
      double* grs_flat,
      int nthreads,
      int MG,
      int JB,
      int TB
  );
}


// [[Rcpp::export]]
Rcpp::NumericMatrix mtgrsbed_matrix_f(
    std::string file,
    int n,
    const std::vector<int>& cls,
    const std::vector<double>& af,
    bool scale,
    Rcpp::NumericMatrix S,   // m x nt
    int nthreads = 1,
    int MG = 64,
    int JB = 1024,
    int TB = 32
) {

  const int m  = S.nrow();
  const int nt = S.ncol();

  // --- input checks ---
  if ((int)cls.size() != m) {
    Rcpp::stop("cls.size() must match number of rows in S");
  }

  if ((int)af.size() != m) {
    Rcpp::stop("af.size() must match number of rows in S");
  }

  // --- convert R matrix (column-major) → SNP-major ---
  const double* s = S.begin();

  std::vector<double> b_snp((size_t)m * nt);

  for (int t = 0; t < nt; ++t) {
    for (int i = 0; i < m; ++i) {
      b_snp[(size_t)i * nt + t] = s[i + (size_t)t * m];
    }
  }

  // --- allocate output ---
  std::vector<double> grs_flat((size_t)n * nt, 0.0);

  // --- Fortran expects pointers to scalars ---
  int n_ = n;
  int m_ = m;
  int nt_ = nt;
  int nthreads_ = nthreads;
  int MG_ = MG;
  int JB_ = JB;
  int TB_ = TB;

  // logical(c_bool) → bool is usually OK, but safest:
  bool scale_ = scale;

  // --- call Fortran ---
  mtgrsbed_core_f(
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

  // --- convert back to R matrix ---
  Rcpp::NumericMatrix out(n, nt);

  for (int t = 0; t < nt; ++t) {
    for (int j = 0; j < n; ++j) {
      out(j, t) = grs_flat[(size_t)j * nt + t];
    }
  }

  return out;
}
