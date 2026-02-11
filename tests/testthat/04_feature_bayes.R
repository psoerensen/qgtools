cat("Running 04_feature_bayes.R\n")

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

featureSets <- list(
  set1 = 1:1000,
  set2 = 1001:2000
)

## ------------------------------------------------------------
## Priors
## ------------------------------------------------------------
priors <- list(
  marker = prior(
    variable      = "id",
    traits        = "BW",
    featureMatrix = featureMatrix,
    featureSets   = featureSets,
    distribution  = bayesC(
      pi    = c(0.95, 0.05),
      gamma = c(0, 0.001)
    )
  ),
  residual = prior(
    variable     = "residual",
    traits       = "BW",
    distribution = iw(4, matrix(1))
  )
)

## ------------------------------------------------------------
## JSON export
## ------------------------------------------------------------

validate_bundle_json(
  data_spec   = qgtools:::as_list.DataSource(data),
  model_spec  = qgtools:::as_list.model(formulas = formulas, task = "bayes"),
  kernel_spec = qgtools:::as_list.kernels(vcs),
  prior_spec  = qgtools:::as_list.priors(priors)
)


json <- qgtools:::as_json.data(data)
cat(json)

json <- qgtools:::as_json.features(features)
cat(json)

json <- qgtools:::as_json.priors(priors)
cat(json)

json <- qgtools:::as_json.model(formulas, "bayes")
cat(json)

bundle_path <- tempfile("bundle_")
dir.create(bundle_path)

writeLines(qgtools:::as_json.data(data),      file.path(bundle_path, "data.json"))
writeLines(qgtools:::as_json.features(features), file.path(bundle_path, "features.json"))
writeLines(qgtools:::as_json.priors(priors),  file.path(bundle_path, "priors.json"))
writeLines(qgtools:::as_json.model(formulas, "bayes"), file.path(bundle_path, "model.json"))

cat("✓ Feature-backed Bayesian workflow passed\n")
