# --------------------------------------------------
# Unix/Linux/macOS: Build and test shared library
# --------------------------------------------------

demo_src <- system.file("demo/shared_lib", package = "qgtools")
demo_root <- file.path(tempdir(), "qgtools_demo")

dir.create(demo_root, showWarnings = FALSE)

# Copy demo
file.copy(demo_src, demo_root, recursive = TRUE)

demo_dir <- file.path(demo_root, "shared_lib")
setwd(demo_dir)

cat("\nFiles before build:\n")
print(list.files(recursive = TRUE))

# --------------------------------------------------
# Build
# --------------------------------------------------
cat("\nBuilding shared library...\n")

system("make clean")
status <- system("make")

if (status != 0) {
  stop("Build failed")
}

cat("\nFiles after build:\n")
print(list.files())

# --------------------------------------------------
# Detect correct library extension
# --------------------------------------------------
sys <- Sys.info()[["sysname"]]

lib_file <- if (sys == "Darwin") {
  "libqgtools.dylib"
} else {
  "libqgtools.so"
}

lib_path <- file.path(demo_dir, lib_file)

if (!file.exists(lib_path)) {
  stop(paste("Shared library not found:", lib_file))
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
# Check symbols
# --------------------------------------------------
cat("\nChecking symbols:\n")

cat("lib_mtgrsbed_c:",
    isTRUE(is.loaded("lib_mtgrsbed_c")), "\n")

cat("lib_mtgrsbed_f:",
    isTRUE(is.loaded("lib_mtgrsbed_f")), "\n")

cat("libqgtools_version:",
    isTRUE(is.loaded("libqgtools_version")), "\n")

# --------------------------------------------------
# Sanity test
# --------------------------------------------------
cat("\nTesting version call:\n")
print(.Call("libqgtools_version"))

cat("\n Unix build OK\n")
