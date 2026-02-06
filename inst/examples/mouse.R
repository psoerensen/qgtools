# Mouse example
data <- makeDatalist(
  source = system.file(
    "examples", "mouse", "mouse.csv",
    package = "qgtools"
  ),
  format = "CSV",
  id     = "id"
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
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)
cat(as_json.vcs(vcs))
cat(as_json.data(data))
cat(as_json.kernels(vcs))
task <- "reml"
cat(as_json.model(formulas, task))

fit <- gfit(
  BW ~ sex + reps + (1 | id),
  data = data,
  vcs  = vcs,
  task = "reml"
)
