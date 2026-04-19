#pragma once

#ifdef __cplusplus
extern "C" {
#endif

  void lib_mtgrsbed_core_c(
      const char* file,
      int n,
      const int* cls,
      const double* af,
      int m,
      int nt,
      int scale,
      const double* b_snp,
      double* grs_flat,
      int nthreads,
      int MG,
      int JB,
      int TB
  );

#ifdef __cplusplus
}
#endif
