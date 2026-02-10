cat("Running 03_feature_reml.R\n")

## ------------------------------------------------------------
## Data
## ------------------------------------------------------------
data <- makeDataSource(
  source = "pheno.csv",
  id     = "id",
  roles  = list(
    BW = "trait",
    id = "id"
  )
)

## ------------------------------------------------------------
## Model
## ------------------------------------------------------------
formulas <- list(
  BW = BW ~ (1 | marker)
)

## ------------------------------------------------------------
## Features
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
  set1 = 1:1000,
  set2 = 1001:2000
)

## ------------------------------------------------------------
## Variance components
## ------------------------------------------------------------
vcs <- list(
  marker = vc(
    variable      = "marker",
    traits        = "BW",
    featureMatrix = features,
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
