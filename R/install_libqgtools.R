install_libqgtools <- function(
    path = NULL,
    persist = TRUE,
    verbose = TRUE
) {

  # --------------------------------------------------
  # Locate demo source
  # --------------------------------------------------
  demo_src <- system.file("demo/shared_lib", package = "qgtools")

  if (demo_src == "") {
    stop("Could not find shared_lib demo directory")
  }

  # --------------------------------------------------
  # Create working directory
  # --------------------------------------------------
  if (is.null(path)) {
    path <- file.path(tempdir(), "qgtools_lib")
  }

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  if (verbose) {
    message("Installing libqgtools into: ", path)
  }

  # --------------------------------------------------
  # Copy source
  # --------------------------------------------------
  file.copy(demo_src, path, recursive = TRUE)

  build_dir <- file.path(path, "shared_lib")

  # --------------------------------------------------
  # Build
  # --------------------------------------------------
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)

  setwd(build_dir)

  if (verbose) message("Building shared library...")

  system("make clean", ignore.stdout = !verbose)
  status <- system("make", ignore.stdout = !verbose)

  if (status != 0) {
    stop("Compilation failed. Ensure make, g++, and gfortran are available.")
  }

  # --------------------------------------------------
  # Detect library file
  # --------------------------------------------------
  sys <- Sys.info()[["sysname"]]

  lib_file <- switch(sys,
                     Windows = "libqgtools.dll",
                     Darwin  = "libqgtools.dylib",
                     Linux   = "libqgtools.so",
                     stop("Unsupported OS")
  )

  lib_path <- file.path(build_dir, lib_file)

  if (!file.exists(lib_path)) {
    stop("Shared library not found after build")
  }

  # --------------------------------------------------
  # Load library
  # --------------------------------------------------
  dlls <- getLoadedDLLs()

  if ("libqgtools" %in% names(dlls)) {
    try(dyn.unload(dlls[["libqgtools"]]$path), silent = TRUE)
  }

  dyn.load(lib_path)

  # --------------------------------------------------
  # Verify installation
  # --------------------------------------------------
  ok <- isTRUE(is.loaded("lib_mtgrsbed_c")) &&
    isTRUE(is.loaded("lib_mtgrsbed_f")) &&
    isTRUE(is.loaded("libqgtools_version"))

  if (!ok) {
    stop("Library loaded but symbols not found")
  }

  # --------------------------------------------------
  # Store path for future sessions
  # --------------------------------------------------
  if (persist) {
    options(qgtools.lib = lib_path)

    if (verbose) {
      message("Stored path in options(qgtools.lib)")
    }
  }

  if (verbose) {
    message("✔ libqgtools installed successfully")
  }

  invisible(lib_path)
}
