qgtools
================

## Overview

**qgtools** is a framework for specifying and fitting mixed and
hierarchical models, with a focus on quantitative genetics,
high-dimensional predictors, and complex covariance structures.

The key idea behind qgtools is to **separate model specification from
model fitting**, allowing the same model to be estimated using different
software backends and estimation paradigms within a unified modeling
framework.

## Core ideas and orthogonal layers

qgtools separates model specification into independent, orthogonal
layers.  
Each layer answers a distinct modeling question and can be modified
without affecting the others.

1.  **Formulas — *what* enters the model**  
    Per-trait formulas define responses, fixed effects, and named random
    or latent effects  
    (e.g. `(1 | id)`, `(1 | marker)`).  
    Formulas describe *which* effects exist, but not how their
    covariance is modeled.

2.  **Kernels and features — *how* effects act on the data**  
    qgtools supports two complementary mechanisms:

    - **Kernels** define *implicit* covariance structures across
      individuals or traits  
      (e.g. pedigree, GRM, spatial or temporal correlation).
    - **Feature matrices** define *explicit* linear predictors with
      potentially very large numbers of coefficients (e.g. genotypes,
      transcriptomics, other omics layers).

    Both kernels and features are reusable, standalone objects that can
    be shared across models and estimation tasks.

3.  **Variance components or priors — *how much* variation is
    attributed**  
    Each effect is paired with either:

    - variance components (`vc()`) for REML or solver-based estimation,
      or  
    - prior distributions (`prior()`) for Bayesian inference.

    For feature-based effects, optional `featureSets` allow effects to
    be decomposed into multiple components with separate variance
    parameters.

4.  **Task — *how* the model is fitted**  
    Estimation is controlled by `gfit()` via `task = "reml"`, `"solve"`,
    or `"bayes"`,  
    without changing the model structure.

5.  **Data — *where* the information comes from**  
    Data are treated as abstract data sources, enabling transparent use
    of in-memory or disk-backed datasets (e.g. PLINK files, HDF5
    matrices).

All layers are specified independently and validated jointly at fit
time.

## Simple example to illustrate the concept

The following example shows a minimal linear mixed model and illustrates
how qgtools separates **model structure**, **covariance specification**,
and **estimation**.

``` r
## ---------------------------------------------------------------
## Model structure (formulas)
## ---------------------------------------------------------------
## Formulas define which effects enter the model,
## but not how covariance is modeled.

formulas <- list(
  BW = BW ~ sex + (1 | id)
)

## ---------------------------------------------------------------
## Covariance specification (kernels + variance components)
## ---------------------------------------------------------------
## The pedigree kernel defines how covariance is induced
## across levels of the 'id' effect.

PED <- makePEDlist(fnPED = "pedigree.txt")

vcs <- list(
  animal = vc(
    variable = "id",
    traits   = "BW",
    kernel   = PED
  ),
  residual = vc(
    variable = "residual",
    traits   = "BW"
  )
)

## ---------------------------------------------------------------
## Model fitting
## ---------------------------------------------------------------
## The estimation task determines how the model is fit,
## without changing the model structure.

fit <- gfit(formulas, data, vcs, task = "reml")
```

The same model structure can be estimated as a Bayesian model by
replacing *vc()* with *prior()* and setting *task=“bayes”*, without
changing the formulas or kernels.

## Model validation and interoperability

Before fitting, qgtools validates the full model bundle  
(data, formulas, kernels or features, and variance components or priors)
to ensure internal consistency.

This validation step checks, for example, that: - all variables
referenced in formulas are defined, - kernels or feature matrices are
compatible with the corresponding effects, - trait dimensions are
consistent across components, and - required variance components or
priors are supplied.

Model specifications can also be exported as structured JSON, allowing
the same model to be executed by external backends, workflow engines, or
non-R/Python environments.

### Residual effects

The residual variance is represented by a component with
`variable = "residual"`.

- The residual effect is **implicit** in model formulas and must **not**
  appear as `(1 | residual)`
- It must be explicitly specified using `vc()` (REML / solver) or
  `prior()` (Bayesian)
- Residual components never require a kernel or features

## Performance and deployment

qgtools is designed as a lightweight R (or Python) interface to
high-performance computational backends.  
Computationally intensive components—such as likelihood evaluation,
large-scale linear algebra, and sampling—are implemented in compiled
languages such as **C++ or Fortran**, while R (or Python) is used for
model specification and orchestration.

This separation provides: - high computational performance, -
scalability to large datasets, and - a clear boundary between model
definition and numerical implementation.

By representing kernels and feature matrices as *lightweight
descriptors*, qgtools avoids unnecessary memory usage and allows
backends to construct covariance structures and linear predictors only
when needed, using in-memory or disk-backed data as appropriate.

For deployment, qgtools can be packaged into **containerized
environments** (e.g. Docker or Singularity), allowing models to be
executed reproducibly on high-performance computing platforms and cloud
infrastructures such as **AWS** and **Azure**. Containers bundle the
required libraries and runtimes and can be deployed on HPC batch systems
or cloud services without exposing source code.

This design enables qgtools to act as a unifying modeling layer while
supporting multiple computational backends and deployment scenarios,
from local workstations to large-scale cloud and HPC environments.

**qgtools** handles large-scale data by taking advantage of:

- multi-core processing using [openMP](https://www.openmp.org/)
- multithreaded matrix operations implemented in BLAS libraries
  (e.g. [OpenBLAS](https://www.openblas.net/),
  [ATLAS](https://math-atlas.sourceforge.net/), or
  [MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library))
- fast and memory-efficient batch processing of feature data stored on
  disk (e.g. genotype data in
  [PLINK](https://www.cog-genomics.org/plink2) bedfiles)

## R and Python interfaces

qgtools is designed with a **language-agnostic core**, allowing both **R
and Python** interfaces to interact with the same underlying
computational backends.

Model specification and orchestration can be performed in either R or
Python, while computationally intensive tasks (e.g. likelihood
evaluation, large-scale linear algebra, sampling) are handled by shared
C++/Fortran libraries.  
This design ensures consistent results across interfaces and enables
flexible deployment in cloud and HPC environments.

The example below shows the *same mixed model* specified in R and
Python.

``` r
## ---------------------------------------------------------------
## R interface
## ---------------------------------------------------------------
## Model structure and covariance specification are expressed
## using the same abstractions as in previous examples.

formulas <- list(
  BW = BW ~ sex + reps + (1 | id)
)

vcs <- list(
  animal = vc(
    variable = "id",
    traits   = "BW",
    kernel   = PED
  ),
  residual = vc(
    variable = "residual",
    traits   = "BW"
  )
)

fit <- gfit(formulas, data, vcs, task = "reml")

# ---------------------------------------------------------------
# Python interface
# ---------------------------------------------------------------
# The same model specification expressed using Python syntax.
# The underlying model and backend execution are identical.

model = Model(
    formulas={
        "BW": "BW ~ sex + reps + (1 | id)"
    },
    vcs=[
        vc(variable="id", traits=["BW"], kernel=PED),
        vc(variable="residual", traits=["BW"])
    ]
)

fit = gfit(model, data, task="reml")
```

Both interfaces construct the same internal model representation and
invoke the same computational backend. Differences between R and Python
are limited to syntax and user-facing language conventions.

## Annotated example: Multivariate mixed / Bayesian genomic model

``` r
## ------------------------------------------------------------------
## Per-trait model formulas
## ------------------------------------------------------------------
## Formulas define *which* effects enter the model,
## not how covariance is modeled.

formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | id),
  Gl = Gl ~ sex + reps + (1 | dam) + (1 | id)
)

## ------------------------------------------------------------------
## Kernel objects (covariance "recipes")
## ------------------------------------------------------------------
## Kernels describe how covariance across index levels is constructed
## or accessed by the backend.

PED <- makePEDlist(fnPED = "pedigree.txt")

Glist <- makeGlist(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
)

## ------------------------------------------------------------------
## Bayesian prior specification
## ------------------------------------------------------------------
## Priors replace variance components in Bayesian models.
## Each prior corresponds to a random effect or latent component.

priors <- list(
  animal = prior(
    index = "id",
    traits = c("BW", "Gl"),
    kernel = PED,
    distribution = iw(df = 4, S = diag(2))
  ),

  marker = prior(
    index = "marker",
    traits = c("BW", "Gl"),
    kernel = Glist,
    distribution = bayesC(pi = c(0.95, 0.05))
  ),

  residual = prior(
    index = "Residual",
    traits = c("BW", "Gl"),
    distribution = iw(df = 4, S = diag(2))
  )
)

## ------------------------------------------------------------------
## Model fitting
## ------------------------------------------------------------------
## The estimation task determines how the model is fit,
## without changing the model specification.

fit <- gfit(formulas, data, priors, task = "bayes")
```

The same model structure can be estimated either as a classical mixed
model (REML / solver) or as a Bayesian hierarchical model by changing
only the variance specification and the estimation task.

## Kernel and feature objects are lightweight descriptors

Kernel and feature objects in **qgtools** are designed as *lightweight
descriptors* rather than containers of explicit model matrices or
parameters. In particular, they do **not** store large covariance
matrices or design matrices in memory by default.

Instead, these objects describe **how model components should be
constructed or accessed** by downstream computational backends:

- **Kernel objects** specify how covariance should be induced across
  levels of an effect  
  (e.g. from pedigree information, precomputed relationship matrices, or
  other structured sources).

- **Feature objects** specify how high-dimensional predictors enter the
  model  
  (e.g. genotype matrices, transcriptomic measurements, or other omics
  features),  
  potentially stored on disk and accessed in a streaming or block-wise
  fashion.

Depending on the estimation task and backend, covariance and linear
predictors may be:

- derived implicitly from pedigree or relationship information,
- computed on the fly from feature matrices without forming explicit
  covariance matrices, or
- accessed from precomputed objects supplied by the user.

This abstraction decouples **statistical model specification** from
**data representation and numerical implementation**. As a result, the
same model formulation and user-facing code can be applied seamlessly
across a wide range of data scales—from small illustrative examples to
very large genomic or multi-omics datasets—without modification.

By deferring matrix construction and data access to backend-specific
implementations, **qgtools** achieves both flexibility and scalability
while preserving a clear and consistent statistical interface.

- **Variance components vs priors**  
  Classical mixed models associate kernels or features with variance
  components (`vc()`),  
  while Bayesian models replace these with explicit prior distributions
  (`prior()`),  
  without changing formulas, kernels, or feature definitions.

#### Prepare input data

``` r
## Data may be provided as an in-memory data frame or as a disk-backed file
data <- makeDatalist(
  source = "data.txt",
  format = "CSV"
)

## Prepare pedigree kernel for additive genetic effects
PED <- makePEDlist(fnPED = "pedigree.txt", method = "S-D-NonInbred")

## Prepare genomic relationship kernel for additive genomic effects
GRM <- makeGRMlist(fnGRM = "grm_inverse.txt", format = "BINARY", grm_type = "G-inverse")
```

#### Linear mixed models: REML and solver-based estimation

``` r
## Single-trait linear mixed model (REML)

## Model formula:
## (1 | id) represents the additive genetic (animal) effect
formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)

## Variance components define covariance structure
vcs <- list(
  dam_env = vc(
    index  = "dam",
    traits = "BW"
  ),

  animal_genetic = vc(
    index  = "id",
    traits = "BW",
    kernel = PED
  ),

  residual = vc(
    index  = "Residual",
    traits = "BW"
  )
)

## Fit the model using REML
fit <- gfit(formulas, data, vcs, task = "reml")

## Multi-trait linear mixed model (REML)

formulas_mt <- list(
  Gl = Gl ~ sex + reps + (1 | dam) + (1 | id),
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)

## Shared variance components induce correlation between traits
## For REML, supplied values are used as starting values
vcs_mt <- list(
  dam_env = vc(
    index  = "dam",
    traits = c("Gl", "BW"),
    start  = matrix(
      c(1.0, 0.2,
        0.2, 2.0),
      nrow = 2,
      byrow = TRUE
    )
  ),

  animal_genetic = vc(
    index  = "id",
    traits = c("Gl", "BW"),
    kernel = PED,
    start  = matrix(
      c(4.0, 0.8,
        0.8, 3.0),
      nrow = 2,
      byrow = TRUE
    )
  ),

  residual = vc(
    index  = "Residual",
    traits = c("Gl", "BW"),
    start  = diag(c(5, 5))
  )
)

## REML estimation
fit_mt <- gfit(formulas_mt, data, vcs_mt, task = "reml")

## Solve mixed model equations (no variance updates)
fit_mt <- gfit(formulas_mt, data, vcs_mt, task = "solve")
```

#### Bayesian hierarchical linear mixed models

``` r
## Single-trait Bayesian linear mixed model

formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)

# Specify prior distributions for selected factors
priors <- list(
  dam_env = prior(
    index = "dam",
    traits = "BW",
    distribution = iw(df = 4, S = 1),
    start = 1
  ),

  animal_genetic = prior(
    index = "id",
    traits = "BW",
    kernel = PED,
    distribution = iw(df = 4, S = 1),
    start = 1
  ),

  residual = prior(
    index = "Residual",
    traits = "BW",
    distribution = iw(df = 4, S = 1),
    start = 1
  )
)

# Estimate parameters using sampling based methods such as Gibbs
fit <- gfit(formulas, data, priors, task = "bayes")


## Multi-trait Bayesian linear mixed model

formulas_mt <- list(
  Gl = Gl ~ sex + reps + (1 | dam) + (1 | id),
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)

priors_mt <- list(
  dam_env = prior(
    index  = "dam",
    traits = c("Gl", "BW"),
    distribution = iw(df = 4, S = diag(2)),
    start  = matrix(
      c(1.0, 0.2,
        0.2, 2.0),
      nrow = 2,
      byrow = TRUE
    )
  ),

  animal_genetic = prior(
    index  = "id",
    traits = c("Gl", "BW"),
    kernel = PED,
    distribution = iw(df = 4, S = diag(2)),
    start  = matrix(
      c(4.0, 0.8,
        0.8, 3.0),
      nrow = 2,
      byrow = TRUE
    )
  ),

  residual = prior(
    index  = "Residual",
    traits = c("Gl", "BW"),
    distribution = iw(df = 4, S = diag(2)),
    start  = diag(c(5, 5))
  )
)

fit_mt <- gfit(formulas_mt, data, priors_mt, task = "bayes")
```

#### Bayesian multi-component linear models with marker-level priors

Marker-level effects differ fundamentally from classical random effects:

- They typically involve **thousands to millions of coefficients**
- Covariance is induced **implicitly** through the feature matrix
  (e.g. genotypes), linkage disequilibrium, and the chosen prior
- Explicit covariance matrices (e.g. GRMs) are usually avoided for
  scalability

In **qgtools**, marker effects are treated as **feature-based model
components**. They are declared in the model formula as a named effect
(e.g. `(1 | marker)`), and linked to genotype data via a
`featureMatrix()`.

Regularization and covariance structure are then specified through
**marker-level priors**, such as `bayesC()`, optionally using
`featureSets` to define multiple variance components over disjoint (or
annotated) subsets of markers.

This approach mirrors how marker effects are handled in software such as
BGLR and BayesR, while preserving a clear separation between:

- **model structure** (formulas and named effects),
- **data representation** (feature matrices and feature sets), and
- **inference** (Bayesian priors or REML/solver-based estimation).

Unlike classical low-dimensional random effects, marker effects are not
associated with explicit covariance matrices at the observation level.
Instead, their covariance is induced implicitly through the feature
matrix and the prior, enabling scalable multi-component and multi-trait
models.

## Annotated example: Multivariate mixed / Bayesian genomic model

``` r
## ------------------------------------------------------------------
## Per-trait model formulas
## ------------------------------------------------------------------
## Formulas declare *which* effects enter the model,
## but do not specify how covariance is modeled.

formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | animal) + (1 | marker),
  Gl = Gl ~ sex + reps + (1 | dam) + (1 | animal) + (1 | marker)
)

## ------------------------------------------------------------------
## Kernel-based effect: pedigree (random effect)
## ------------------------------------------------------------------
## Pedigree effects are naturally expressed via kernels.

PED <- makePEDlist(fnPED = "pedigree.txt")


## ------------------------------------------------------------------
## Feature-based effect: markers (high-dimensional predictors)
## ------------------------------------------------------------------
## Marker effects are defined through a feature matrix

M <- featureMatrix(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
)

## Optional grouping of markers into multiple variance components
featureSets <- list(
  set1 = 1:1000,
  set2 = 1001:2000
)

## ------------------------------------------------------------------
## Bayesian prior specification
## ------------------------------------------------------------------
## Priors replace variance components in Bayesian models.
## Each prior corresponds to a named model component.

priors <- list(
  dam = prior(
    variable     = "dam",
    traits       = c("BW", "Gl"),
    distribution = iw(df = 4, S = diag(1, 2)),
    start        = diag(1, 2)
  ),

  animal = prior(
    variable     = "animal",
    traits       = c("BW", "Gl"),
    kernel       = PED,
    distribution = iw(df = 4, S = diag(1, 2)),
    start        = diag(1, 2)
  ),

  marker = prior(
    variable     = "marker",
    traits       = c("BW", "Gl"),
    features     = M,
    featureSets  = featureSets,
    distribution = bayesC(
      pi     = beta(95, 5),
      sigma2 = invchisq(df = 4, scale = 0.001),
      start  = list(
        set1 = diag(0.005, 2),
        set2 = diag(0.001, 2)
      )
    )
  ),

  residual = prior(
    variable     = "residual",
    traits       = c("BW", "Gl"),
    distribution = iw(df = 4, S = diag(1, 2)),
    start        = diag(1, 2)
  )
)

## ------------------------------------------------------------------
## Model fitting
## ------------------------------------------------------------------
## The estimation task determines how the model is fit,
## without changing the model structure.

fit_bayes <- gfit(
  formulas = formulas,
  data     = data,
  priors   = priors,
  task     = "bayes"
)
```

The same model structure can be estimated using a marker BLUP
formulation by replacing `prior()` with `vc()` and setting
`task = "solve"`. The formulas, kernels, and feature matrices remain
unchanged.

``` r
## ------------------------------------------------------------------
## Marker BLUP (solver-based estimation)
## ------------------------------------------------------------------

vcs <- list(
  dam = vc(
    variable = "dam",
    traits   = c("BW", "Gl"),
    start    = diag(1, 2)
  ),

  animal = vc(
    variable = "animal",
    traits   = c("BW", "Gl"),
    kernel   = PED,
    start    = diag(1, 2)
  ),

  marker = vc(
    variable    = "marker",
    traits      = c("BW", "Gl"),
    features    = M,
    featureSets = featureSets,
    start       = list(
      set1 = diag(0.005, 2),
      set2 = diag(0.001, 2)
    )
  ),

  residual = vc(
    variable = "residual",
    traits   = c("BW", "Gl"),
    start    = diag(1, 2)
  )
)

fit_solve <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "solve"
)
```
