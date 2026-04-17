# --------------------------------------------------
# Windows: Build and test shared library
# --------------------------------------------------

# Path to demo directory (adjust if needed)
demo_src <- system.file("demo/shared_lib", package = "qgtools")
demo_dir <- file.path(tempdir(), "qgtools_demo")

dir.create(demo_dir, showWarnings = FALSE)
file.copy(demo_src, demo_dir, recursive = TRUE)

setwd(demo_dir)


cat("Building shared library (Windows)...\n")
system("make clean")
system("make")

cat("\nFiles after build:\n")
print(list.files())

# Load DLL
lib_path <- file.path(demo_dir, "libqgtools.dll")

# Unload if already loaded
if ("libqgtools" %in% names(getLoadedDLLs())) {
  dyn.unload(lib_path)
}

dyn.load(lib_path)

# Check symbols
cat("\nChecking symbols:\n")
cat("mtgrsbed_core_c:", is.loaded("mtgrsbed_core_c"), "\n")
cat("mtgrsbed_core_f:", is.loaded("mtgrsbed_core_f"), "\n")

cat("\n Windows build OK\n")
