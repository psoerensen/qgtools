cat("Running 03_feature_reml.R\n")

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
featureMatrix <- makeGlist(
  fnBED = "chr.bed",
  fnBIM = "chr.bim",
  fnFAM = "chr.fam"
)

featureSets <- list(
  set1 = 1:1000,
  set2 = 1001:2000
)

## ------------------------------------------------------------
## Variance components
## ------------------------------------------------------------
vcs <- list(
  marker = vc(
    variable      = "id",
    traits        = "BW",
    featureMatrix = featureMatrix,
    featureSets   = featureSets
  ),
  residual = vc(
    variable = "residual",
    traits   = "BW"
  )
)

## ------------------------------------------------------------
## Export + validate
## ------------------------------------------------------------
bundle_path <- tempfile("bundle_")
dir.create(bundle_path)

export_model_bundle(
  data     = data,
  features = features,
  formulas = formulas,
  vcs      = vcs,
  task     = "reml",
  path     = bundle_path
)

cat("✓ Feature-backed REML workflow passed\n")
