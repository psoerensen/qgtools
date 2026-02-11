cat("Running 05_summary_vcs.R\n")

## ------------------------------------------------------------
## Data (in-memory)
## ------------------------------------------------------------
pheno <- data.frame(
  id  = 1:10,
  BW  = rnorm(10),
  sex = factor(rep(c("M", "F"), 5))
)

data <- makeDataSource(
  source = pheno,
  id     = "id",
  roles  = list(
    BW  = "trait",
    sex = "fixed",
    id  = "id"
  )
)

## ------------------------------------------------------------
## Model
## ------------------------------------------------------------
formulas <- list(
  BW = BW ~ (1 | id)
)

## ------------------------------------------------------------
## Features
## ------------------------------------------------------------
G <- makeGlist(
  fnBED = "chr.bed",
  fnBIM = "chr.bim",
  fnFAM = "chr.fam"
)

featureMatrix <- makeFeatureSource(
  geno = G
)

featureSets1 <- list(
  set1 = 1:1000,
  set2 = 1001:2000
)

featureSets2 <- list(
  set1 = 2001:3000,
  set2 = 3001:4000
)

## ------------------------------------------------------------
## Variance components
## ------------------------------------------------------------

vcs <- varcomp(
  animal = vc("id", "BW", kernel = PED),
  marker1 = vc("id", "BW", featureMatrix = featureMatrix, featureSets = featureSets1),
  marker2 = vc("id", "BW", featureMatrix = featureMatrix, featureSets = featureSets2),
  residual = vc("residual", "BW")
)

summary(vcs)
