# --------------------------------------------------
# Windows: Build and test shared library
# --------------------------------------------------

demo_src <- system.file("demo/shared_lib", package = "qgtools")
demo_root <- file.path(tempdir(), "qgtools_demo")

dir.create(demo_root, showWarnings = FALSE)

# copy demo
file.copy(demo_src, demo_root, recursive = TRUE)

demo_dir <- file.path(demo_root, "shared_lib")
setwd(demo_dir)

cat("\nFiles before build:\n")
print(list.files(recursive = TRUE))

cat("Building shared library...\n")
system("make clean")
status <- system("make")

if (status != 0) {
  stop("Build failed")
}

cat("\nFiles after build:\n")
print(list.files())

# -------------------------------
# Load DLL
# -------------------------------
lib_path <- file.path(demo_dir, "libqgtools.dll")

if (!file.exists(lib_path)) {
  stop("libqgtools.dll not found")
}

# Unload if already loaded (safe)
dlls <- getLoadedDLLs()
if ("libqgtools" %in% names(dlls)) {
  try(dyn.unload(dlls[["libqgtools"]]$path), silent = TRUE)
}

dyn.load(lib_path)

# -------------------------------
# Check symbols
# -------------------------------
cat("\nChecking symbols:\n")

cat("lib_mtgrsbed_c:",
    isTRUE(is.loaded("lib_mtgrsbed_c", PACKAGE = "libqgtools")), "\n")

cat("lib_mtgrsbed_f:",
    isTRUE(is.loaded("lib_mtgrsbed_f", PACKAGE = "libqgtools")), "\n")

cat("libqgtools_version:",
    isTRUE(is.loaded("libqgtools_version", PACKAGE = "libqgtools")), "\n")



