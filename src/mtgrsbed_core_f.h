#pragma once

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif
