#include "mtgrsbed_core.h"
#include "mtgrsbed_core_c.h"

extern "C" {

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
  ) {
    if (MG <= 0) MG = 64;
    if (JB <= 0) JB = 1024;
    if (TB <= 0) TB = 32;

    mtgrsbed_core(
      file,
      n,
      cls,
      af,
      m,
      nt,
      static_cast<bool>(scale),
      b_snp,
      grs_flat,
      nthreads,
      MG,
      JB,
      TB
    );
  }

}
