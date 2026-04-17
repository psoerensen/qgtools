#' Set number of threads for qgtools
#' @param n number of threads
#' @export
qgtools_set_threads <- function(n = NULL,
                                backend = c("openmp", "blas", "auto"),
                                verbose = TRUE) {

  backend <- match.arg(backend)

  # Detect threads
  if (is.null(n)) {
    n <- getOption("qgtools.threads", NA_integer_)
    if (is.na(n)) {
      n <- parallel::detectCores(logical = FALSE)
    }
  }

  # Normalize
  n <- as.integer(n)
  if (is.na(n) || n < 1) n <- 1

  # Save in options
  options(qgtools.threads = n)

  # Detect OS
  sys <- Sys.info()[["sysname"]]

  # Auto backend selection
  if (backend == "auto") {
    # Simple heuristic
    backend <- if (n <= 2) "blas" else "openmp"
  }

  # -----------------------------
  # Apply threading configuration
  # -----------------------------

  if (backend == "openmp") {

    Sys.setenv(
      OMP_NUM_THREADS = n,
      OMP_DYNAMIC = "FALSE",
      OPENBLAS_NUM_THREADS = 1,
      MKL_NUM_THREADS = 1
    )

  } else if (backend == "blas") {

    Sys.setenv(
      OMP_NUM_THREADS = 1,
      OMP_DYNAMIC = "FALSE",
      OPENBLAS_NUM_THREADS = n,
      MKL_NUM_THREADS = n
    )

  }

  # macOS warning
  if (sys == "Darwin" && backend == "openmp") {
    warning("OpenMP may not be available on macOS; performance may be reduced.")
  }

  if (verbose) {
    cat("qgtools threading configured:\n")
    cat("  backend :", backend, "\n")
    cat("  threads :", n, "\n")
  }

  invisible(list(
    backend = backend,
    threads = n
  ))
}

#qgtools_set_threads <- function(n) {
#  options(qgtools.threads = as.integer(n))
#}


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

qgtools_run_openmp <- function(expr) {
  old <- Sys.getenv(c("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"))

  on.exit(do.call(Sys.setenv, as.list(old)), add = TRUE)

  qgtools_set_threads(getOption("qgtools.threads"), backend = "openmp", verbose = FALSE)

  force(expr)
}

qgtools_run_blas <- function(expr) {
  old <- Sys.getenv(c("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"))

  on.exit(do.call(Sys.setenv, as.list(old)), add = TRUE)

  qgtools_set_threads(getOption("qgtools.threads"), backend = "blas", verbose = FALSE)

  force(expr)
}

# options(qgtools.threads = 24)
#
# # OpenMP kernel
# qgtools_run_openmp({
#   mtgrsbed_matrix(...)
# })
#
# # MKL/BLAS computation
# qgtools_run_blas({
#   crossprod(X)
# })

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

# ext <- switch(sys,
#               Linux   = "so",
#               Darwin  = "dylib",
#               Windows = "dll"
# )
#
# dyn.load(paste0("libgrsbed_cpp.", ext))


# run_openmp <- function(...) {
#   old <- Sys.getenv(c("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"))
#   on.exit(do.call(Sys.setenv, as.list(old)), add = TRUE)
#
#   Sys.setenv(
#     OMP_NUM_THREADS = getOption("qgtools.threads", 1),
#     MKL_NUM_THREADS = 1,
#     OPENBLAS_NUM_THREADS = 1
#   )
#
#   .Call("myomp_function", ...)
# }
#
# run_mkl <- function(...) {
#   old <- Sys.getenv(c("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS"))
#   on.exit(do.call(Sys.setenv, as.list(old)), add = TRUE)
#
#   Sys.setenv(
#     OMP_NUM_THREADS = 1,
#     MKL_NUM_THREADS = getOption("qgtools.threads", 1),
#     OPENBLAS_NUM_THREADS = 1
#   )
#
#   .Call("mymkl_function", ...)
# }
#

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
