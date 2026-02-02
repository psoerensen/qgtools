formulas <- list(
  V0 = V0 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +
    (1 | Dam) +
    (1 | Animal),

  V1 = V1 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +
    (1 | Dam) +
    (1 | Animal)
)


IID_LITTER <- list(
  ids  = "L_Dam",
  type = "iid"
)
class(IID_LITTER) <- "IIDlist"


IID_RESID <- list(
  ids  = "Animal",
  type = "iid"
)
class(IID_RESID) <- "IIDlist"


PEDlist <- list(
  fnPED = "sheep.ped",
  format = "DMU",
  columns = list(
    id   = "ID",
    sire = "SIRE",
    dam  = "DAM",
    sort = "BIRTH"
  ),
  method = 2,              # DMU: non-inbred
  phantom_parents = TRUE,
  inbreeding = FALSE,
  idsP = "Animal"          # <-- REQUIRED
)
class(PEDlist) <- "PEDlist"


Klist <- qg_Klist(
  PED        = PEDlist,
  IID_LITTER = IID_LITTER,
  IID_RESID  = IID_RESID
)


vcov <- list(

  # Environmental: litter within dam
  litter_env = list(
    index     = "L_Dam",
    traits    = c("V0", "V1"),
    kernel    = "IID_LITTER",
    structure = "unstructured"
  ),

  # Maternal additive genetic effect
  maternal_genetic = list(
    index     = "Dam",
    traits    = c("V0", "V1"),
    kernel    = "PED",
    structure = "unstructured"
  ),

  # Direct additive genetic effect
  direct_genetic = list(
    index     = "Animal",
    traits    = c("V0", "V1"),
    kernel    = "PED",
    structure = "unstructured"
  ),

  # Residual
  residual = list(
    index     = "Residual",
    traits    = c("V0", "V1"),
    kernel    = "IID_RESID",
    structure = "unstructured"
  )
)

priors <- list(
  maternal_genetic = list(
    traits = c("V0", "V1"),
    matrix = matrix(
      c(2.0, 0.3,
        0.3, 1.5),
      nrow = 2, byrow = TRUE
    )
  ),

  direct_genetic = list(
    traits = c("V0", "V1"),
    matrix = matrix(
      c(4.0, 0.8,
        0.8, 3.0),
      nrow = 2, byrow = TRUE
    )
  ),

  residual = list(
    traits = c("V0", "V1"),
    matrix = matrix(
      c(5.0, 0.5,
        0.5, 4.0),
      nrow = 2, byrow = TRUE
    )
  )
)

tmp <- tempfile(fileext = ".yaml")

obj_written <- write_qg_yaml(
  formulas = formulas,
  vcov     = vcov,
  Klist    = Klist,
  file     = tmp
)

obj_read <- read_qg_yaml(tmp)
