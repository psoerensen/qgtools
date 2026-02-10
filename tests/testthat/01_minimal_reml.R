cat("Running 01_minimal_reml.R\n")

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
  BW = BW ~ sex + (1 | id)
)

## ------------------------------------------------------------
## Variance components
## ------------------------------------------------------------
vcs <- list(
  id = vc(
    variable = "id",
    traits   = "BW",
    kernel   = iid_kernel("id")
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
  features = NULL,
  formulas = formulas,
  vcs      = vcs,
  task     = "reml",
  path     = bundle_path
)

cat("✓ Minimal REML workflow passed\n")
