## ============================================================
## Test script for validate_bundle_json()
## ============================================================

library(qgtools)

cat("Running validate_bundle_json test\n")

## ------------------------------------------------------------
## 1. Minimal DataSource
## ------------------------------------------------------------
data <- makeDataSource(
  source = data.frame(
    id  = 1:10,
    BW  = rnorm(10),
    sex = factor(rep(c("F", "M"), 5))
  ),
  id    = "id",
  roles = list(
    BW  = "trait",
    sex = "fixed",
    id  = "id"
  )
)

data_spec <- qgtools:::as_list.DataSource(data)

## ------------------------------------------------------------
## 2. Minimal model formulas
## ------------------------------------------------------------
formulas <- list(
  BW = BW ~ sex + (1 | id)
)

model_spec <- as_list.model(
  formulas = formulas,
  task     = "reml"
)

## ------------------------------------------------------------
## 3. Kernel specification
## ------------------------------------------------------------
K_id <- iid_kernel("id")
K_id$meta$id <- "IID_id"

vcs <- list(
  id = vc(
    variable = "id",
    traits   = "BW",
    kernel   = K_id
  ),
  residual = vc(
    variable = "residual",
    traits   = "BW"
  )
)

kernel_spec  <- qgtools:::as_list.kernels(vcs)
varcomp_spec <- qgtools:::as_list.vcs(vcs)

## ------------------------------------------------------------
## 4. Run validation (REML case)
## ------------------------------------------------------------
cat("Validating REML bundle...\n")

validate_bundle_json(
  data_spec    = data_spec,
  model_spec   = model_spec,
  kernel_spec  = kernel_spec,
  varcomp_spec = varcomp_spec
)

cat("✔ REML bundle validation passed\n\n")

## ------------------------------------------------------------
## 5. Negative test: missing kernel
## ------------------------------------------------------------
cat("Running negative test: missing kernel\n")

bad_vcs <- list(
  id = list(
    variable   = "id",
    traits     = "BW",
    covariance = list(type = "kernel", kernel = "NON_EXISTENT"),
    structure  = "unstructured"
  )
)

bad_varcomp_spec <- list(
  type = "VarCompSpec",
  components = bad_vcs
)

tryCatch(
  validate_bundle_json(
    data_spec    = data_spec,
    model_spec   = model_spec,
    kernel_spec  = kernel_spec,
    varcomp_spec = bad_varcomp_spec
  ),
  error = function(e) {
    cat("✔ Caught expected error:\n")
    cat(e$message, "\n\n")
  }
)

## ------------------------------------------------------------
## 6. Bayesian test (feature-free, kernel-based)
## ------------------------------------------------------------
cat("Validating Bayesian bundle...\n")

priors <- list(
  id = prior(
    variable     = "id",
    traits       = "BW",
    kernel       = K_id,
    distribution = iw(df = 4, S = matrix(1))
  ),
  residual = prior(
    variable     = "residual",
    traits       = "BW",
    distribution = iw(df = 4, S = matrix(1))
  )
)

prior_spec <- qgtools:::as_list.prior(priors)

validate_bundle_json(
  data_spec   = data_spec,
  model_spec  = model_spec,
  kernel_spec = kernel_spec,
  prior_spec  = prior_spec
)


cat("✔ Bayesian bundle validation passed\n")

cat("All validate_bundle_json tests completed successfully\n")
