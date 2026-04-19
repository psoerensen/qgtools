# ==================================================
# qgtools demo: internal vs external backend
# ==================================================

# --------------------------------------------------
# Load package
# --------------------------------------------------
library(qgtools)

cat("\nqgtools demo started\n")

# --------------------------------------------------
# Load example data
# --------------------------------------------------
Glist <- readRDS(
  file = "C:/Users/au223366/Dropbox/Mohammed/Simulations/Glist_EUR_10k.RDS"
)

# --------------------------------------------------
# Define inputs
# --------------------------------------------------
file <- Glist$bedfiles   # PLINK .bed file
n <- Glist$n             # number of individuals
m <- Glist$mchr          # number of SNPs
nt <- 3                  # number of traits

# SNP metadata
cls <- sample(1:m, m)
af  <- runif(m, 0.05, 0.5)

# SNP x trait effects
S <- matrix(rnorm(m * nt), nrow = m, ncol = nt)

cat("\nInput prepared\n")

# ==================================================
# INTERNAL BACKENDS
# ==================================================

cat("\nRunning internal C++ backend...\n")

grs_cpp <- qgtools:::mtgrsbed(
  file, n, cls, af, S,
  nthreads = 4,
  backend = "internal_cpp"
)

cat("Done\n")

cat("\nRunning internal Fortran backend...\n")

grs_fortran <- qgtools:::mtgrsbed(
  file, n, cls, af, S,
  nthreads = 4,
  backend = "internal_fortran"
)

cat("Done\n")

# --------------------------------------------------
# Verify internal consistency
# --------------------------------------------------
diff_internal <- max(abs(grs_cpp - grs_fortran))

cat("\nInternal consistency check:\n")
cat("max |cpp - fortran| =", diff_internal, "\n")

if (diff_internal != 0) {
  stop("Internal backends do not match!")
}

# ==================================================
# INSTALL EXTERNAL LIBRARY
# ==================================================

cat("\nInstalling libqgtools...\n")

qgtools:::install_libqgtools()

cat("External library installed\n")

# ==================================================
# EXTERNAL BACKEND
# ==================================================

cat("\nRunning external C++ backend...\n")

grs_external <- qgtools:::mtgrsbed(
  file, n, cls, af, S,
  nthreads = 4,
  backend = "external_cpp"
)

cat("Done\n")

# --------------------------------------------------
# Verify external correctness
# --------------------------------------------------
diff_external <- max(abs(grs_cpp - grs_external))

cat("\nExternal consistency check:\n")
cat("max |cpp - external| =", diff_external, "\n")

if (diff_external != 0) {
  stop("External backend does not match internal!")
}

# ==================================================
# SUCCESS
# ==================================================

cat("\n✔ All checks passed\n")
cat("✔ External backend works correctly\n")
cat("✔ Plugin system is functional\n")

