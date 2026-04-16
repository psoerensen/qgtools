# .onLoad <- function(libname, pkgname) {
#   lib <- switch(Sys.info()[["sysname"]],
#                 "Linux"   = "linux/libqgtools.so",
#                 "Darwin"  = "mac/libqgtools.dylib",
#                 "Windows" = "windows/libqgtools.dll"
#   )
#   dyn.load(system.file("libs", lib, package = pkgname))
# }

.with_threads <- function(n, mode, expr) {
  old <- Sys.getenv(c(
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_DYNAMIC"
  ), unset = NA)

  on.exit({
    do.call(Sys.setenv, as.list(old))
  }, add = TRUE)

  if (mode == "omp") {
    Sys.setenv(
      OMP_NUM_THREADS = n,
      OMP_DYNAMIC = "FALSE",
      OPENBLAS_NUM_THREADS = 1,
      MKL_NUM_THREADS = 1
    )
  }

  if (mode == "blas") {
    Sys.setenv(
      OMP_NUM_THREADS = 1,
      OPENBLAS_NUM_THREADS = n,
      MKL_NUM_THREADS = n
    )
  }

  force(expr)
}

# example usage
# mtgrsbed_matrix <- function(..., nthreads = 1) {
#
#   .with_threads(nthreads, "omp", {
#     .Call(...)
#   })
#
# }


# my_blas_function <- function(..., nthreads = 1) {
#
#   .with_threads(nthreads, "blas", {
#     # BLAS-heavy code
#   })
#
# }
