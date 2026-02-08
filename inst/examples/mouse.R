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


formulas <- list(
  BW = BW ~ sex + reps + (1 | marker)
)


# Define feature matrix
M <- featureMatrix(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
)

mSets <- list(1:1000, 1001:2000)

# Define prior distribution for marker features used in bayes
marker <- prior(
  variable    = "marker",
  traits   = c("BW", "Gl"),
  features = M,
  featureSets = mSets,
  distribution = bayesC(
    pi     = beta(95, 5),
    sigma2 = invchisq(df = 4, scale = 0.001),
    start  = list(diag(0.005, 2),diag(0.001, 2))
  )
)

# Define prior distribution for marker features used in bayes
marker <- prior(
  variable    = "marker",
  traits   = c("BW", "Gl"),
  features = M,
  featureSets = mSets,
  distribution = ridge(
    sigma2 = invchisq(df = 4, scale = 0.001),
    start  = list(diag(0.005, 2),diag(0.001, 2))
  )
)

marker <- vc(
  variable = "marker",
  traits = c("BW","Gl"),
  features = M,
  featureSets = mSets,
  start  = list(set1 = diag(0.005,2), set2 = diag(0.001,2))
)

# Define variance component for marker features used in reml/solve
marker = vc(
  index     = "id",
  traits    = c("BW", "Gl"),
  features = M,
  sets = sets,  # or group
  kernel = iid_kernel(),
  start = list(diag(0.001, 2),diag(0.002, 2),....diag(0.003, 2))
)


DNA <- featureMatrix(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
)

RNA <- featureMatrix(
  counts   = "rna_counts.h5",
  samples  = "samples.csv"
)

gSets <- list(
  chr1 = 1:1000,
  chr2 = 1001:2000
)

tSets <- list(
  pathway1 = c("geneA","geneB","geneC"),
  pathway2 = c("geneD","geneE")
)


marker <- vc(
  variable    = "marker",
  traits      = c("BW","Gl"),
  features    = DNA,
  featureSets = gSets,
  start       = list(
    chr1 = diag(0.005,2),
    chr2 = diag(0.001,2)
  )
)

transcript <- vc(
  variable    = "transcript",
  traits      = c("BW","Gl"),
  features    = RNA,
  featureSets = tSets,
  start       = diag(0.01, 2)
)

marker <- prior(
  variable    = "marker",
  traits      = c("BW","Gl"),
  features    = DNA,
  featureSets = gSets,
  distribution = bayesC(
    pi     = beta(95,5),
    sigma2 = invchisq(df=4, scale=0.001)
  )
)

transcript <- prior(
  variable    = "transcript",
  traits      = c("BW","Gl"),
  features    = T,
  distribution = ridge(
    sigma2 = invchisq(df=4, scale=0.01)
  )
)

