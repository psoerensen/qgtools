# mtgrsbed <- function(file, n, cls, af, S,
#                      scale = TRUE,
#                      nthreads = 1) {
#
#   lib_path <- getOption("qgtools.lib")
#
#   # try external plugin
#   if (!is.null(lib_path) && file.exists(lib_path)) {
#
#     if (!"libqgtools" %in% names(getLoadedDLLs())) {
#       dyn.load(lib_path, DLLname = "libqgtools")
#     }
#
#     if (is.loaded("lib_mtgrsbed_c", PACKAGE = "libqgtools")) {
#       return(.Call("lib_mtgrsbed_c",
#                    file, n, cls, af, scale, S, nthreads,
#                    PACKAGE = "libqgtools"))
#     }
#   }
#
#   # fallback
#   .Call("mtgrsbed_c",
#         file, n, cls, af, scale, S, nthreads,
#         PACKAGE = "qgtools")
# }
#
mtgrsbed <- function(file, n, cls, af, S,
                     scale = TRUE,
                     nthreads = 1,
                     backend = c("auto", "internal_cpp", "internal_fortran",
                                 "external_cpp", "external_fortran"),
                     lib_path = getOption("qgtools.lib")) {


  cls <- as.integer(cls)
  n   <- as.integer(n)
  af  <- as.numeric(af)
  nthreads <- as.integer(nthreads)

  backend <- match.arg(backend)

  # -------------------------------
  # Load external library if needed
  # -------------------------------
  if (backend %in% c("auto", "external_cpp", "external_fortran")) {
    if (!is.null(lib_path) && file.exists(lib_path)) {
      if (!"libqgtools" %in% names(getLoadedDLLs())) {
        dyn.load(lib_path)
      }
    }
  }

  # -------------------------------
  # Auto backend
  # -------------------------------
  if (backend == "auto") {

    ok_version <- FALSE

    if ("libqgtools" %in% names(getLoadedDLLs()) &&
        isTRUE(is.loaded("libqgtools_version", PACKAGE = "libqgtools"))) {

      ver <- try(.Call("libqgtools_version", PACKAGE = "libqgtools"), silent = TRUE)

      ok_version <- !inherits(ver, "try-error") &&
        identical(as.integer(ver), 1L)
    }

    if (ok_version &&
        isTRUE(is.loaded("lib_mtgrsbed_c", PACKAGE = "libqgtools"))) {

      return(.Call("lib_mtgrsbed_c",
                   file, n, cls, af, scale, S, nthreads,
                   PACKAGE = "libqgtools"))
    }

    return(.Call("mtgrsbed_c",
                 file, n, cls, af, scale, S, nthreads,
                 PACKAGE = "qgtools"))
  }

  # -------------------------------
  # Manual backend selection
  # -------------------------------
  switch(
    backend,

    internal_cpp =
      .Call("mtgrsbed_c",
            file, n, cls, af, scale, S, nthreads,
            PACKAGE = "qgtools"),

    internal_fortran =
      .Call("mtgrsbed_f",
            file, n, cls, af, scale, S, nthreads,
            PACKAGE = "qgtools"),

    external_cpp = {
      if (!isTRUE(is.loaded("lib_mtgrsbed_c", PACKAGE = "libqgtools"))) {
        stop("External C++ backend not available")
      }
      .Call("lib_mtgrsbed_c",
            file, n, cls, af, scale, S, nthreads,
            PACKAGE = "libqgtools")
    },

    external_fortran = {
      if (!isTRUE(is.loaded("lib_mtgrsbed_f", PACKAGE = "libqgtools"))) {
        stop("External Fortran backend not available")
      }
      .Call("lib_mtgrsbed_f",
            file, n, cls, af, scale, S, nthreads,
            PACKAGE = "libqgtools")
    }
  )
}

# mtgrsbed <- function(file, n, cls, af, S,
#                      scale = TRUE,
#                      nthreads = getOption("qgtools.threads", 1),
#                      backend = c("cpp", "fortran", "shared_cpp", "shared_fortran")) {
#
#   backend <- match.arg(backend)
#
#   mtgrsbed_matrix_backend(
#     file, n, cls, af, scale, S,
#     nthreads = nthreads,
#     backend = backend
#   )
# }

# mtgrsbed <- function(file, n, cls, af, S,
#                      scale = TRUE,
#                      nthreads = 1,
#                      backend = c("auto", "cpp", "fortran",
#                                  "cpp_shared", "fortran_shared"),
#                      lib_path = getOption("qgtools.lib")) {
#
#   backend <- match.arg(backend)
#
#   # -------------------------------
#   # Load external lib if needed
#   # -------------------------------
#   if (grepl("_shared", backend) || backend == "auto") {
#
#     if (!is.null(lib_path) && file.exists(lib_path)) {
#
#       if (!"libqgtools" %in% names(getLoadedDLLs())) {
#         dyn.load(lib_path, DLLname = "libqgtools")  # 🔥 important
#       }
#     }
#   }
#
#   # -------------------------------
#   # Auto backend
#   # -------------------------------
#   if (backend == "auto") {
#
#     if ("libqgtools" %in% names(getLoadedDLLs())) {
#
#       # Prefer shared C backend if available
#       if (is.loaded("lib_mtgrsbed_c", PACKAGE = "libqgtools")) {
#
#         return(.Call("lib_mtgrsbed_c",
#                      file, n, cls, af, scale, S, nthreads,
#                      PACKAGE = "libqgtools"))
#       }
#     }
#
#     backend <- "cpp"
#   }
#
#   # -------------------------------
#   # Dispatch
#   # -------------------------------
#   switch(backend,
#
#          cpp = .Call("mtgrsbed_c",
#                      file, n, cls, af, scale, S, nthreads,
#                      PACKAGE = "qgtools"),
#
#          fortran = .Call("mtgrsbed_f",
#                          file, n, cls, af, scale, S, nthreads,
#                          PACKAGE = "qgtools"),
#
#          cpp_shared = .Call("lib_mtgrsbed_c",
#                             file, n, cls, af, scale, S, nthreads,
#                             PACKAGE = "libqgtools"),
#
#          fortran_shared = .Call("lib_mtgrsbed_f",
#                                 file, n, cls, af, scale, S, nthreads,
#                                 PACKAGE = "libqgtools")
#   )
# }
