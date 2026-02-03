# ------------------------------------------------------------------------------
# Variable names used in formulas refer to columns in the data frame / file
# ------------------------------------------------------------------------------

# Define formulas for each trait
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


qg_parse_formulas(formulas)


# ------------------------------------------------------------------------------
# Define kernels (data sources)
# ------------------------------------------------------------------------------

pedigree <- makePEDlist(fnPED = "pedigree.txt")
grm      <- makeGRMlist(fnGRM = "grm.txt")

# ------------------------------------------------------------------------------
# Define variance components
# ------------------------------------------------------------------------------

vcs <- varcomp(

  litter_env = vc(
    index  = "L_Dam",
    traits = c("V0", "V1")
  ),

  maternal_genetic = vc(
    index  = "Dam",
    traits = c("V0", "V1"),
    kernel = ped_kernel(pedigree),
    prior  = prior_iw(df = 4, S = diag(2))
  ),

  direct_genetic = vc(
    index  = "Animal",
    traits = c("V0", "V1"),
    kernel = ped_kernel(pedigree),
    prior  = prior_iw(
      df = 4,
      S  = matrix(c(4, 0.8,
                    0.8, 3), 2, 2, byrow = TRUE)
    )
  ),

  residual = vc(
    index  = "Residual",
    traits = c("V0", "V1"),
    prior  = prior_diag(c(5, 4))
  )
)


random_effects <- list(

  litter_env = vc(
    index     = "L_Dam",
    traits    = c("V0", "V1"),
    kernel    = IIDkernel("litter"),
    structure = "unstructured"
  ),

  maternal_genetic = vc(
    index     = "Dam",
    traits    = c("V0", "V1"),
    kernel    = ped_kernel,
    structure = "unstructured",
    prior     = prior_iw(
      df = 4,
      S  = diag(2)
    )
  ),

  direct_genetic = vc(
    index     = "Animal",
    traits    = c("V0", "V1"),
    kernel    = ped_kernel,
    structure = "unstructured",
    prior     = prior_iw(
      df = 4,
      S  = matrix(c(4.0, 0.8,
                    0.8, 3.0), 2, 2, byrow = TRUE)
    )
  ),

  residual = vc(
    index     = "Residual",
    traits    = c("V0", "V1"),
    kernel    = IIDkernel("residual"),
    structure = "unstructured",
    prior     = prior_iw(
      df = 4,
      S  = diag(2)
    )
  )
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

# Kernels: covariance structures defined over the indexing units of random effects
#          (e.g. animals, individuals, or SNPs), such as pedigree- or
#          genomic-based relationship matrices and LD matrices.

# (Co)variances: covariance structures defined over traits, describing
#                genetic and residual (co)variation between traits.
#                Typical structures include diagonal, unstructured,
#                and structured (e.g. factor-analytic) forms.

# Kernels: covariance structures over the indexing units of a random effect
#          (e.g. animals, individuals, or SNPs)
#          Implemented as lists carrying kernel-specific information:
#            - Genomic relationship matrix (GRMlist)
#            - Pedigree-based relationship matrix (PEDlist)
#            - Linkage disequilibrium (LD) matrix (LDlist)


# (Co)variances: covariance structures over traits
#                (e.g. genetic and residual (co)variances between traits)
#                Common structures include:
#                  - Diagonal
#                  - Unstructured
#                  - Structured (e.g. banded, factor-analytic)


# Define kernels
IID_LITTER <- list(
  ids = "L_Dam",             # IDS link to ID column in data frame/file
  type = "iid",
  extra_argument_1 = NULL,
  extra_argument_2 = NULL,
  extra_argument_3 = NULL
)
class(IID_LITTER) <- "IIDlist"

IID_RESID <- list(
  ids = "L_Dam",             # IDS link to ID column in data frame/file
  type = "iid",
  extra_argument_1 = NULL,
  extra_argument_2 = NULL,
  extra_argument_3 = NULL
)
class(IID_RESID) <- "IIDlist"


GRMlist <- list(
  ids = "Animal",          # IDS link to ID column in data frame/file
  fnPED = "grm.txt",
  format = "DMU",
  extra_argument_1 = NULL,
  extra_argument_2 = NULL,
  extra_argument_3 = NULL
)
class(GRMlist) <- "GRMlist"




# Create PEDlist
PEDlist <- list(
  ids = "Animal", # IDS link to ID column in data frame/file
  fnPED = "pedigree.txt",               # pedigree in ASCII or in binary format
  format = "ASCII",
  method = "S-D-NonInbred"             # options: "S-D-Inbred","S-D-NonInbred","S-MGS-Inbred","S-MGS-NonInbred","S-D-NonInbred-PP"
)
class(PEDlist) <- "PEDlist"

# make pedigree object
makePED <- function(ids=NULL, fnPED=NULL, format="ASCII", method="S-D-NonInbred") {

  PEDlist <- list(
    ids = ids, # IDS link to ID column in data frame/file
    fnPED = fnPED,               # pedigree in ASCII or in binary format
    format = format,
    method = methods             # options: "S-D-Inbred","S-D-NonInbred","S-MGS-Inbred","S-MGS-NonInbred","S-D-NonInbred-PP"
  )
  class(PEDlist) <- "PEDlist"
}


Glist <- list(
  ids = "Animal", # IDS link to ID column in data frame/file
  fnBED = c("plink_chr1.bed","plink_chr1.bed"),     # One or more bedfiles
  fnBIM = c("plink_chr1.bim","plink_chr1.bim"),     # One or more bimfiles
  fnFAM = c("plink_chr1.fam","plink_chr1.fam"),     # One or more famfiles
  format = "PLINK",                                 # redundant?
  rsids = NULL,                                     # optional
  sets = NULL,                                      # optional
  weights =NULL,                                    # optional
  annotation = NULL,                                # optional
  extra_argument_1 = NULL,
  extra_argument_2 = NULL,
  extra_argument_3 = NULL
)
class(GRMlist) <- "GRMlist"


pedigree = list(
  index     = "individual",
  kernel    = "PED",
  structure = "unstructured"
),
polygenic = list(
  index     = "individual",
  kernel    = "GRM",
  structure = "unstructured"
),
genomic = list(
  index     = "individual",
  kernel    = "G",
  structure = "diagonal",
  sets      = NULL
)


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
