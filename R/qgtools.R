.onLoad <- function(libname, pkgname) {

  if (is.null(getOption("qgtools.threads"))) {

    n <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))

    if (is.na(n) || n <= 0) {
      n <- parallel::detectCores(logical = FALSE)
    }

    options(qgtools.threads = n)
  }

  # Internal thread manager ---------------------------------------------
  if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Darwin") {
    packageStartupMessage(
      "OpenMP may be unavailable on macOS; performance may be reduced."
    )
  }
  #   lib <- switch(Sys.info()[["sysname"]],
  #                 "Linux"   = "linux/libqgtools.so",
  #                 "Darwin"  = "mac/libqgtools.dylib",
  #                 "Windows" = "windows/libqgtools.dll"
  #   )
  #   dyn.load(system.file("libs", lib, package = pkgname))

}


