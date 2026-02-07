library(qgtools)
# Mouse example
data <- makeDatalist(
  source = system.file(
    "examples", "mouse", "mouse.csv",
    package = "qgtools"
  ),
  format = "CSV",
  id     = "id",
  sep     = ";"
)

PED <- makePEDlist(
  fnPED = system.file(
    "examples", "mouse", "pedigree.csv",
    package = "qgtools"
  ),
  format = "SIMPLE"
)

vcs <- list(
  animal = vc(
    index  = "id",
    traits = c("BW"),
    kernel = PED,
    start  = 1.0
  ),
  residual = vc(
    index  = "Residual",
    traits = c("BW"),
    start  = 1.0
  )
)

formulas <- list(
  BW = BW ~ sex + reps + (1 | id) + (1 | id)
)
cat(as_json.vcs(vcs))
cat(as_json.data(data))
cat(as_json.kernels(vcs))

## Authoring phase (R)

data_spec    <- as_list.Datalist(data)
task <- "reml"
model_spec   <- as_list.model(formulas, task)
kernel_spec  <- as_list.kernels(vcs)
varcomp_spec <- as_list.vcs(vcs)

validate_bundle_reml(
  data_spec    = data_spec,
  model_spec   = model_spec,
  varcomp_spec     = varcomp_spec,
  kernel_spec = kernel_spec
)

# Variance specification in the tradional linear mixed model
vcs <- list(
  animal = vc(
    index     = "id",
    traits    = c("BW", "Gl"),
    kernel    = PED,
    structure = "unstructured"
  ),
  residual = vc(
    index  = "Residual",
    traits = c("BW", "Gl"),
    kernel = iid_kernel()
  )
)

# Variance specification in the Bayesian linear mixed model
priors <- list(
  animal = prior(
    index = "id",
    traits = c("BW", "Gl"),
    kernel    = PED,
    distribution = iw(df = 4, S = diag(2)),
    start = diag(2)
  ),
  residual = prior(
    index = "Residual",
    traits = c("BW", "Gl"),
    distribution = iw(df = 4, S = diag(2)),
    start = diag(2)
  )
)


str(varcomp_spec)
str(kernel_spec)
str(model_spec)

formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | id),
  Gl = Gl ~ sex + reps + (1 | dam) + (1 | id)
)
cat(as_json.model(formulas, task))

formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)
task <- "reml"
as_list.model(formulas, task)

fit <- gfit(
  BW ~ sex + reps + (1 | id),
  data = data,
  vcs  = vcs,
  task = "reml"
)
