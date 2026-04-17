#include <cmath>
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

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
) {
  FILE* file_stream = std::fopen(file, "rb");
  if (!file_stream) {
    throw std::runtime_error("Could not open BED file.");
  }

  const size_t nbytes = (n + 3) / 4;

  // Raw BED bytes for a block of MG markers
  unsigned char* block_buffer =
    (unsigned char*) std::malloc((size_t)MG * nbytes);
  if (!block_buffer) {
    std::fclose(file_stream);
    throw std::runtime_error("Failed to allocate BED block buffer.");
  }

  // Shared maps for current marker block
  std::vector<double> map0(MG), map1(MG), map2(MG), map3(MG);

#ifdef _OPENMP
  if (nthreads > 0) {
    omp_set_dynamic(0);
    omp_set_num_threads(nthreads);
  }
#endif

#pragma omp parallel
{
  for (int i0 = 0; i0 < m; i0 += MG) {
    const int imax = std::min(i0 + MG, m);
    const int mlen = imax - i0;

#pragma omp single
{
  for (int ii = 0; ii < mlen; ++ii) {
    const int i = i0 + ii;

    const long long offset =
      (long long)(cls[i] - 1) * (long long)nbytes + 3LL;

    if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
      std::free(block_buffer);
      std::fclose(file_stream);
      throw std::runtime_error("fseek failed.");
    }

    unsigned char* buf_i = block_buffer + (size_t)ii * nbytes;

    const size_t nbytes_read =
      std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);

    if (nbytes_read != nbytes) {
      std::free(block_buffer);
      std::fclose(file_stream);
      throw std::runtime_error("Error reading BED data: nbytes_read != nbytes.");
    }

    const double p = af[i];
    if (scale) {
      const double denom = std::sqrt(2.0 * p * (1.0 - p));
      if (denom <= 0.0 || !std::isfinite(denom)) {
        map0[ii] = 0.0;
        map1[ii] = 0.0;
        map2[ii] = 0.0;
        map3[ii] = 0.0;
      } else {
        map0[ii] = (2.0 - 2.0 * p) / denom;
        map1[ii] = 0.0;                    // missing
        map2[ii] = (1.0 - 2.0 * p) / denom;
        map3[ii] = (-2.0 * p) / denom;
      }
    } else {
      map0[ii] = 2.0;
      map1[ii] = -2.0 * p;                // mean-imputed centered coding
      map2[ii] = 1.0;
      map3[ii] = 0.0;
    }
  }
} // implicit barrier here

#pragma omp for schedule(static)
for (int j0 = 0; j0 < n; j0 += JB) {
  const int jmax = std::min(j0 + JB, n);

  const int byte0 = j0 >> 2;
  const int byte1 = (jmax + 3) >> 2;

  for (int t0 = 0; t0 < nt; t0 += TB) {
    const int tmax = std::min(t0 + TB, nt);
    const int tlen = tmax - t0;

    // Loop over markers in current block
    for (int ii = 0; ii < mlen; ++ii) {
      const int i = i0 + ii;
      const unsigned char* buf_i = block_buffer + (size_t)ii * nbytes;
      const double* bi = &b_snp[(size_t)i * nt + t0];

      const double m0 = map0[ii];
      const double m1 = map1[ii];
      const double m2 = map2[ii];
      const double m3 = map3[ii];

      for (int kb = byte0; kb < byte1; ++kb) {
        unsigned char buf_k = buf_i[kb];
        const int jbase = kb << 2;

        for (int pos = 0; pos < 4; ++pos) {
          const int j = jbase + pos;
          if (j < j0 || j >= jmax || j >= n) {
            buf_k >>= 2;
            continue;
          }

          double xj;
          switch (buf_k & 3u) {
          case 0u: xj = m0; break;   // BED 00 -> 2 copies
          case 1u: xj = m1; break;   // BED 01 -> missing
          case 2u: xj = m2; break;   // BED 10 -> 1 copy
          default: xj = m3; break;   // BED 11 -> 0 copies
          }
          buf_k >>= 2;

          double* gj = &grs_flat[(size_t)j * nt + t0];

#pragma omp simd
          for (int t = 0; t < tlen; ++t) {
            gj[t] += bi[t] * xj;
          }
        }
      }
    }
  }
}
  }
}

std::free(block_buffer);
std::fclose(file_stream);
}

// void mtgrsbed_core(
//     const char* file,
//     int n,
//     const int* cls,
//     const double* af,
//     int m,
//     int nt,
//     bool scale,
//     const double* b_snp,
//     double* grs_flat,
//     int nthreads,
//     int MG = 64,
//     int JB = 1024,
//     int TB = 32
// ) {
//   FILE* file_stream = std::fopen(file, "rb");
//   if (!file_stream) {
//     throw std::runtime_error("Could not open BED file.");
//   }
//
//   const size_t nbytes = (n + 3) / 4;
//
//   unsigned char* block_buffer =
//     (unsigned char*) std::malloc((size_t)MG * nbytes);
//   if (!block_buffer) {
//     std::fclose(file_stream);
//     throw std::runtime_error("Failed to allocate BED block buffer.");
//   }
//
//   std::vector<double> map0(MG), map1(MG), map2(MG), map3(MG);
//
//   // Shared error handling
//   bool error = false;
//   const char* error_msg = nullptr;
//
// #ifdef _OPENMP
//   if (nthreads > 0) {
//     omp_set_num_threads(nthreads);
//   }
// #endif
//
// #pragma omp parallel
// {
//   for (int i0 = 0; i0 < m; i0 += MG) {
//     const int imax = std::min(i0 + MG, m);
//     const int mlen = imax - i0;
//
// #pragma omp single
// {
//   for (int ii = 0; ii < mlen; ++ii) {
//     if (error) continue;  // Early skip if error occurred
//
//     const int i = i0 + ii;
//
//     const long long offset =
//       (long long)(cls[i] - 1) * (long long)nbytes + 3LL;
//
//     // if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
//     //   error = true;
//     //   error_msg = "fseek failed.";
//     //   continue;
//     // }
//
//     if (std::fseek(file_stream, offset, SEEK_SET) != 0) {
// #pragma omp atomic write
//       error = true;
//
// #pragma omp critical
// {
//   if (!error_msg) error_msg = "fseek failed.";
// }
// continue;
//     }
//
//     unsigned char* buf_i = block_buffer + (size_t)ii * nbytes;
//
//     const size_t nbytes_read =
//       std::fread(buf_i, sizeof(unsigned char), nbytes, file_stream);
//
//     // if (nbytes_read != nbytes) {
//     //   error = true;
//     //   error_msg = "Error reading BED data: nbytes_read != nbytes.";
//     //   continue;
//     // }
//
//     if (nbytes_read != nbytes) {
// #pragma omp atomic write
//       error = true;
//
// #pragma omp critical
// {
//   if (!error_msg) error_msg = "Error reading BED data: nbytes_read != nbytes.";
// }
// continue;
//     }
//     const double p = af[i];
//     if (scale) {
//       const double denom = std::sqrt(2.0 * p * (1.0 - p));
//       if (denom <= 0.0 || !std::isfinite(denom)) {
//         map0[ii] = 0.0;
//         map1[ii] = 0.0;
//         map2[ii] = 0.0;
//         map3[ii] = 0.0;
//       } else {
//         map0[ii] = (2.0 - 2.0 * p) / denom;
//         map1[ii] = 0.0;
//         map2[ii] = (1.0 - 2.0 * p) / denom;
//         map3[ii] = (-2.0 * p) / denom;
//       }
//     } else {
//       map0[ii] = 2.0;
//       map1[ii] = -2.0 * p;
//       map2[ii] = 1.0;
//       map3[ii] = 0.0;
//     }
//   }
// } // implicit barrier
//
// #pragma omp for schedule(static)
// for (int j0 = 0; j0 < n; j0 += JB) {
//   if (error) continue;  // Skip work if error
//
//   const int jmax = std::min(j0 + JB, n);
//
//   const int byte0 = j0 >> 2;
//   const int byte1 = (jmax + 3) >> 2;
//
//   for (int t0 = 0; t0 < nt; t0 += TB) {
//     const int tmax = std::min(t0 + TB, nt);
//     const int tlen = tmax - t0;
//
//     for (int ii = 0; ii < mlen; ++ii) {
//       const int i = i0 + ii;
//       const unsigned char* buf_i = block_buffer + (size_t)ii * nbytes;
//       const double* bi = &b_snp[(size_t)i * nt + t0];
//
//       const double m0 = map0[ii];
//       const double m1 = map1[ii];
//       const double m2 = map2[ii];
//       const double m3 = map3[ii];
//
//       for (int kb = byte0; kb < byte1; ++kb) {
//         unsigned char buf_k = buf_i[kb];
//         const int jbase = kb << 2;
//
//         for (int pos = 0; pos < 4; ++pos) {
//           const int j = jbase + pos;
//           if (j < j0 || j >= jmax || j >= n) {
//             buf_k >>= 2;
//             continue;
//           }
//
//           double xj;
//           switch (buf_k & 3u) {
//           case 0u: xj = m0; break;
//           case 1u: xj = m1; break;
//           case 2u: xj = m2; break;
//           default: xj = m3; break;
//           }
//           buf_k >>= 2;
//
//           double* gj = &grs_flat[(size_t)j * nt + t0];
//
// #pragma omp simd
//           for (int t = 0; t < tlen; ++t) {
//             gj[t] += bi[t] * xj;
//           }
//         }
//       }
//     }
//   }
// }
//   }
// }
//
// // Clean up AFTER parallel region
// std::free(block_buffer);
// std::fclose(file_stream);
//
// // Throw error AFTER parallel region
// if (error) {
//   throw std::runtime_error(error_msg);
// }
// }
