#' Set number of threads for qgtools
#' @param n number of threads
#' @export
qgtools_set_threads <- function(n) {
  options(qgtools.threads = as.integer(n))
}


#' Show threading configuration for qgtools
#'
#' Prints information about detected cores and thread settings
#' for OpenMP and BLAS backends.
#'
#' @export
qgtools_threads_info <- function() {

  # Detect cores
  slurm <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA)
  omp   <- Sys.getenv("OMP_NUM_THREADS", unset = NA)
  mkl   <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
  blas  <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  dyn   <- Sys.getenv("OMP_DYNAMIC", unset = NA)

  detected <- parallel::detectCores(logical = FALSE)

  qgthreads <- getOption("qgtools.threads", NA)

  cat("==== qgtools threading info ====\n\n")

  cat("System:\n")
  cat("  Physical cores (detectCores):", detected, "\n")
  cat("  SLURM_CPUS_PER_TASK        :", slurm, "\n\n")

  cat("qgtools:\n")
  cat("  qgtools.threads option     :", qgthreads, "\n\n")

  cat("OpenMP:\n")
  cat("  OMP_NUM_THREADS            :", omp, "\n")
  cat("  OMP_DYNAMIC                :", dyn, "\n\n")

  cat("BLAS:\n")
  cat("  MKL_NUM_THREADS            :", mkl, "\n")
  cat("  OPENBLAS_NUM_THREADS       :", blas, "\n\n")

  # Interpretation -----------------------------------------------------

  cat("Interpretation:\n")

  if (!is.na(slurm) && !is.na(omp) && as.integer(omp) > as.integer(slurm)) {
    cat("  ⚠ OpenMP threads exceed SLURM allocation\n")
  }

  if (!is.na(omp) && !is.na(mkl) && as.integer(omp) > 1 && as.integer(mkl) > 1) {
    cat("  ⚠ OpenMP and MKL both multithreaded (oversubscription risk)\n")
  }

  if (!is.na(omp) && !is.na(blas) && as.integer(omp) > 1 && as.integer(blas) > 1) {
    cat("  ⚠ OpenMP and OpenBLAS both multithreaded (oversubscription risk)\n")
  }

  if (!is.na(omp) && omp == "1") {
    cat("  ℹ OpenMP running single-threaded\n")
  }

  if (!is.na(mkl) && mkl == "1" && !is.na(blas) && blas == "1") {
    cat("  ✔ BLAS threading disabled (safe for OpenMP)\n")
  }

  cat("\n================================\n")

  invisible(NULL)
}

.qgtools_threads <- function() {

  # 1. SLURM (most important)
  n <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))

  # 2. fallback: environment variables (sometimes set manually)
  if (is.na(n)) {
    n <- as.integer(Sys.getenv("OMP_NUM_THREADS", unset = NA))
  }

  # 3. fallback: local machine
  if (is.na(n) || n <= 0) {
    n <- parallel::detectCores(logical = FALSE)
  }

  n
}

.with_threads <- function(mode = c("omp", "blas"), n = NULL, expr) {
  mode <- match.arg(mode)
  if (is.null(n)) n <- .qgtools_threads()

  old <- Sys.getenv(c(
    "OMP_NUM_THREADS",
    "OMP_DYNAMIC",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS"
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


#options(qgtools.threads = 8)
#mtgrsbed_matrix(..., nthreads = 2)

# example usage
# mtgrsbed_matrix <- function(..., nthreads = 1) {
#
#   .with_threads(nthreads, "omp", {
#     .Call(...)
#   })
#
# }

# mtgrsbed_matrix <- function(
#     file, n, cls, af, scale, S,
#     nthreads = NULL, ...
# ) {
#
#   .with_threads("omp", nthreads, {
#     .Call(`_qgtools_mtgrsbed_matrix`, file, n, cls, af, scale, S, ...)
#   })
#
# }

# my_blas_function <- function(x, nthreads = NULL) {
#
#   .with_threads("blas", nthreads, {
#     # BLAS-heavy computation
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
