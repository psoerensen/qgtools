# ==============================================================================
# Example: Multivariate mixed / Bayesian genomic model specification
#
# This example illustrates how to define:
#   - per-trait model formulas
#   - kernel "recipes" for genetic and marker effects
#   - variance components with explicit priors
#   - task-specific model fitting (REML / solve / Bayes)
#
# Variable names used in formulas refer to columns in the data frame / file.
# ==============================================================================

mouse <- readRDS(url("https://github.com/psoerensen/bgcourse/raw/main/data/mouse.rds"))
head(mouse)

source(system.file("examples", "simulate_data.R", package = "qgtools"))
head(data)
str(data)

# ------------------------------------------------------------------------------
# Define per-trait model formulas
# ------------------------------------------------------------------------------
# Each element in `formulas` corresponds to one trait.
# Random effects are specified using standard Wilkinson syntax.
# These formulas define *which* random effects exist,
# not how their covariance is modeled (that comes later).

formulas <- list(
  M0 = V0 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +   # litter environmental effect
    (1 | Dam)   +   # maternal genetic effect
    (1 | Animal),  # direct genetic / marker effects

  M1 = V1 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +
    (1 | Dam)   +
    (1 | Animal)
)


# ------------------------------------------------------------------------------
# Define kernel objects (data sources / recipes)
# ------------------------------------------------------------------------------
# Kernel objects are *lightweight descriptors*.
# They do NOT load matrices into memory.
# Instead, they describe how covariance structures are constructed
# or accessed by downstream computational backends.

# Pedigree-based additive genetic relationships
PEDlist <- makePEDlist(fnPED = "pedigree.txt")

# Genomic relationship matrix (optional, not used below)
GRMlist <- makeGRMlist(fnGRM = "grm.txt")

# Marker-level genotype container (e.g. PLINK BED/BIM/FAM)
# This provides information about where genotype data live on disk
# and how animals map to markers.
Glist <- makeGlist(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
)

# Optional: a named kernel registry.
# This is NOT required by the model, but can be used internally
# for caching or sharing kernel objects across tasks.
Klist <- list(
  pedigree = PEDlist,
  grm      = GRMlist,
  geno     = Glist
)


# ------------------------------------------------------------------------------
# Optional marker-level metadata
# ------------------------------------------------------------------------------
# These objects may be used by marker-based priors (e.g. BayesC/BayesR).
# They can encode SNP sets, weights, or functional annotations.
# Setting them to NULL means "use default behavior".

sets       <- NULL
weights    <- NULL
annotation <- NULL


# ------------------------------------------------------------------------------
# Define variance components
# ------------------------------------------------------------------------------
# Each variance component:
#   - corresponds to a random effect in the formulas
#   - specifies how covariance is induced (kernel)
#   - optionally specifies a prior on the (co)variance parameters
#
# The `index` must match the grouping factor used in the formulas.

vcs <- list(

  # Environmental litter effect (IID by default)
  litter_env = vc(
    index  = "L_Dam",
    traits = c("V0", "V1")
  ),

  # Maternal additive genetic effect via pedigree
  maternal_genetic = vc(
    index  = "Dam",
    traits = c("V0", "V1"),
    kernel = PEDlist,
    prior  = prior_iw(
      df = 4,
      S  = diag(2)
    )
  ),

  # Marker-based genetic effect (BayesC prior)
  # The random effect is still indexed by Animal,
  # but covariance is induced implicitly via genotypes and LD.
  marker = vc(
    index  = "Animal",
    traits = c("V0", "V1"),
    kernel = Glist,
    prior  = prior_bayesC(
      pi         = c(0.95, 0.05),
      var        = c(0, 1),
      sets       = sets,
      weights    = weights,
      annotation = annotation
    )
  ),

  # Polygenic additive genetic effect via pedigree
  # This can coexist with the marker effect if desired.
  direct_genetic = vc(
    index  = "Animal",
    traits = c("V0", "V1"),
    kernel = PEDlist,
    prior  = prior_iw(
      df = 4,
      S  = matrix(
        c(4.0, 0.8,
          0.8, 3.0),
        nrow = 2,
        byrow = TRUE
      )
    )
  ),

  # Residual variance (IID by default)
  residual = vc(
    index  = "Residual",
    traits = c("V0", "V1"),
    prior  = prior_diag(c(5, 4))
  )
)


# ------------------------------------------------------------------------------
# Fit models using different computational backends
# ------------------------------------------------------------------------------
# The same model specification can be used for:
#   - REML estimation
#   - solving mixed model equations
#   - Bayesian inference
#
# Each backend validates which kernels and priors it supports internally.

reml(
  fm   = formulas,
  vcs  = vcs,
  data = DATAList
)

solve(
  fm   = formulas,
  vcs  = vcs,
  data = DATAList
)

bayes(
  fm   = formulas,
  vcs  = vcs,
  data = DATAList
)





# ==============================================================================
# DMU Example II (dmut2): Two-trait sheep growth model
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Per-trait model formulas
#    (structure only: fixed + random effect declarations)
# ------------------------------------------------------------------------------

formulas <- list(

  M0 = V0 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +   # environmental: litter within dam
    (1 | Dam)   +   # maternal additive genetic
    (1 | Animal),  # direct additive genetic

  M1 = V1 ~ Month + DamAge + Litter + Sex + HY +
    (1 | L_Dam) +
    (1 | Dam)   +
    (1 | Animal)
)

# ------------------------------------------------------------------------------
# 2. Kernel definitions (lightweight covariance recipes)
# ------------------------------------------------------------------------------

# Pedigree includes all 2729 animals (records: 2187)
PEDlist <- makePEDlist(fnPED = "pedigree.txt")

# ------------------------------------------------------------------------------
# 3. Variance component specifications
#    (this is where covariance structure is defined)
# ------------------------------------------------------------------------------

vcs <- list(

  # Environmental litter-within-dam effect
  litter_env = vc(
    index  = "L_Dam",
    traits = c("V0", "V1")
    # IID kernel is implicit
  ),

  # Maternal additive genetic effect (Dam)
  maternal_genetic = vc(
    index  = "Dam",
    traits = c("V0", "V1"),
    kernel = PEDlist,
    prior  = prior_iw(
      df = 4,
      S  = diag(2)
    )
  ),

  # Direct additive genetic effect (Animal)
  direct_genetic = vc(
    index  = "Animal",
    traits = c("V0", "V1"),
    kernel = PEDlist,
    prior  = prior_iw(
      df = 4,
      S  = diag(2)
    )
  ),

  # Residual error
  residual = vc(
    index  = "Residual",
    traits = c("V0", "V1"),
    prior  = prior_diag(c(1, 1))
  )
)

# ------------------------------------------------------------------------------
# 4. Fit the model
# ------------------------------------------------------------------------------

# REML estimation (DMU-style)
reml(
  fm   = formulas,
  vcs  = vcs,
  data = DATA
)

# (Later, the same specification could be used for:)
# bayes(fm = formulas, vcs = vcs, data = DATA)
# solve(fm = formulas, vcs = vcs, data = DATA)



# ==============================================================================
# Random regression animal model (DMU RR example)
# ==============================================================================

# ------------------------------------------------------------------------------
# Formula (single longitudinal trait)
# ------------------------------------------------------------------------------

formulas <- list(

  GH = GH ~
    yob + breed + p_age + td +          # fixed class effects
    age +                               # fixed regression
    (1 | Animal) +                      # permanent environment intercept
    (L1 + L2 | Animal)                  # random regressions (PE + genetic)
)

# ------------------------------------------------------------------------------
# Kernel definitions
# ------------------------------------------------------------------------------

PEDlist <- makePEDlist(fnPED = "pedigree.txt")

# ------------------------------------------------------------------------------
# Variance components
# ------------------------------------------------------------------------------

vcs <- list(

  # --------------------------------------------------------------------------
  # Permanent environmental random regression
  #   G_PE ⊗ I
  # --------------------------------------------------------------------------

  pe_rr = vc(
    index  = "Animal",
    traits = c("L1", "L2"),
    kernel = iid_kernel("pe"),
    prior  = prior_iw(
      df = 4,
      S  = diag(2)
    )
  ),

  # --------------------------------------------------------------------------
  # Additive genetic random regression
  #   G_A ⊗ A
  # --------------------------------------------------------------------------

  genetic_rr = vc(
    index  = "Animal",
    traits = c("L1", "L2"),
    kernel = PEDlist,
    prior  = prior_iw(
      df = 4,
      S  = diag(2)
    )
  ),

  # --------------------------------------------------------------------------
  # Residual
  # --------------------------------------------------------------------------

  residual = vc(
    index  = "Residual",
    traits = "GH",
    prior  = prior_diag(1)
  )
)

# ------------------------------------------------------------------------------
# Fit model
# ------------------------------------------------------------------------------

reml(
  fm   = formulas,
  vcs  = vcs,
  data = DATA
)



# Genotype container
Glist <- makeGlist(bedfiles="chr.bed", bimfiles="chr.bim", famfiles="chr.fam")

# Annotation: a factor/character vector of length m (one label per SNP, in Glist order)
# e.g. "coding", "regulatory", "other"
annot <- annotation  # length = m, names ideally are rsids

# Optional: weights per SNP (length m) and/or per class
w <- weights         # length = m

vcs <- list(

  marker = vc(
    index  = "Animal",
    traits = "GH",
    kernel = Glist,

    prior  = prior_bayesRC(
      annotation = annot,              # SNP -> class
      pi         = list(               # class-specific mixture weights
        coding     = c(0.95, 0.04, 0.01),
        regulatory = c(0.98, 0.015, 0.005),
        other      = c(0.995, 0.004, 0.001)
      ),
      var        = c(0, 0.01, 0.1),    # mixture variance grid (shared)
      weights    = w                   # optional SNP weights
    )
  ),

  residual = vc(
    index  = "Residual",
    traits = "GH",
    prior  = prior_diag(1)
  )
)

bayes(fm = formulas, vcs = vcs, data = DATA)
