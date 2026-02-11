cat("Running test_prior_summary.R\n")

library(qgtools)

## ------------------------------------------------------------
## Dummy phenotype data (not actually used by prior)
## ------------------------------------------------------------
pheno <- data.frame(
  id = 1:10,
  BW = rnorm(10)
)

## ------------------------------------------------------------
## Kernel-backed prior
## ------------------------------------------------------------
PED <- makePEDlist(
  fnPED  = "pedigree.txt",
  method = "S-D-NonInbred"
)

p_kernel <- prior(
  variable     = "id",
  traits       = "BW",
  kernel       = PED,
  distribution = iw(df = 4, S = matrix(1)),
  start        = 1
)

cat("\n--- Kernel-backed prior ---\n")
print(summary(p_kernel))


## ------------------------------------------------------------
## Feature-backed prior
## ------------------------------------------------------------
G <- makeGlist(
  fnBED = "chr.bed",
  fnBIM = "chr.bim",
  fnFAM = "chr.fam"
)

features <- makeFeatureSource(
  geno = G
)

featureSets <- list(
  set1 = 1:100,
  set2 = 101:200
)

p_feature <- prior(
  variable      = "id",
  traits        = "BW",
  featureMatrix = features,
  featureSets   = featureSets,
  distribution  = bayesC(
    pi    = c(0.95, 0.05),
    gamma = c(0, 0.001)
  )
)

cat("\n--- Feature-backed prior ---\n")
print(summary(p_feature))


## ------------------------------------------------------------
## Residual prior
## ------------------------------------------------------------
p_resid <- prior(
  variable     = "residual",
  traits       = "BW",
  distribution = iw(df = 4, S = matrix(1))
)

cat("\n--- Residual prior ---\n")
print(summary(p_resid))


cat("\nAll prior summary tests completed successfully.\n")
