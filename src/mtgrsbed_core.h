#pragma once

void mtgrsbed_core(
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
    int MG = 64,
    int JB = 1024,
    int TB = 32
);
