qgtools
================

## Overview

**qgtools** is a framework for specifying and fitting mixed and
hierarchical models, with a focus on quantitative genetics and complex
covariance structures.

The key idea behind qgtools is to **separate model specification from
model fitting**, allowing the same model to be estimated using different
software backends and estimation paradigms within a unified modeling
framework.

## Core ideas

- **Formulas define model structure**  
  Per-trait formulas specify responses, fixed effects, and random
  effects, but not how covariance is modeled.

- **Variance components define covariance**  
  Variance components describe how random effects induce covariance
  across individuals and traits, using kernels and optional priors.

- **Kernels are reusable objects**  
  Pedigree, genomic, and other kernels are standalone objects that can
  be reused across models and fitting tasks.

- **Tasks control estimation**  
  Models are fitted using `gfit()` with `task = "reml"`, `"bayes"`, or
  `"solve"`, without changing the model specification.

- **Scalable data access**  
  Data are treated as a data source, allowing both in-memory and
  disk-backed datasets to be used transparently.

## R and Python interfaces

qgtools is designed with a language-agnostic core, allowing both **R and
Python** interfaces to interact with the same underlying computational
backends.

Model specification and orchestration can be performed in either R or
Python, while computationally intensive tasks are handled by shared
C++/Fortran libraries. This design ensures consistent results across
interfaces and enables flexible deployment in cloud and HPC
environments.

## Annotated example: mMltivariate mixed / Bayesian genomic model

This example illustrates how qgtools separates: - model structure
(formulas), - covariance definitions (kernels), - variance components or
priors, - and estimation task.

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

## Kernel objects are lightweight descriptors

Kernel objects in **qgtools** are designed as *lightweight descriptors*
rather than containers of explicit covariance matrices. In particular,
kernel objects do **not** store covariance matrices in memory.

Instead, a kernel object specifies **how covariance should be
constructed or accessed** by downstream computational backends.
Depending on the use case, covariance information may be derived from:

- pedigree files,
- precomputed genomic relationship matrices (GRMs), or
- genotype data stored directly on disk.

This abstraction decouples the statistical model specification from the
underlying data representation. As a result, the *same model formulation
and user-facing code* can be applied seamlessly across a wide range of
data scales—from small illustrative examples to very large genomic
datasets—without modification.

By deferring covariance construction to backend-specific
implementations, qgtools achieves both flexibility and scalability while
preserving a clear and consistent statistical interface.

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

# Marker-level genotype container (e.g. PLINK BED/BIM/FAM). This provides information about where genotype data live on disk and how animals map to markers.
Glist <- makeGlist(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
)
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

# Specify prior distrubutions for selected factors
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

Marker-level effects differ from classical random effects:

- They typically involve **thousands to millions of parameters**
- Covariance arises implicitly from genotype structure, LD, and prior
  assumptions
- They are rarely written explicitly in mixed-model formulas

In qgtools, marker effects are therefore specified **only through
priors**, using a special `index = "marker"` together with a marker
kernel (e.g. `Glist`).

This mirrors how marker effects are handled in software such as BGLR and
BayesR, while preserving a clean separation between model structure and
inference.

``` r
## Single-trait Bayesian linear regression with marker-level priors

formulas <- list(
  BW = BW ~ sex + reps + (1 | dam)
)

# Optional marker-level metadata
# Used by marker-based priors (e.g. BayesC / BayesR)
# NULL means default behavior
sets       <- NULL
weights    <- NULL
annotation <- NULL

# Prior specification (Bayesian models only)
priors <- list(
  dam_env = prior(
    index = "dam",
    traits = "BW",
    distribution = iw(df = 4, S = 1),
    start = 1
  ),
  marker_genetic = prior(
    index  = "marker",     # marker-level effects
    traits = "BW",
    kernel = Glist,        # genotype / marker kernel
    distribution = bayesC(
      pi         = c(0.95, 0.05),
      variances  = c(0, 1),
      sets       = sets,
      weights    = weights,
      annotation = annotation
    ),
    start = 1
  ),

  residual = prior(
    index  = "Residual",
    traits = "BW",
    distribution = iw(df = 4, S = 1),
    start = 1
  )
)
```

#### Example R vs Python interface

``` r
# R interface
formulas <- list(
  BW = BW ~ sex + reps + (1 | id)
)

vcs <- list(
  animal = vc(index = "id", traits = "BW", kernel = PED)
)

fit <- gfit(formulas, data, vcs, task = "reml")

# Python interface
model = Model(
    formulas={
        "BW": "BW ~ sex + reps + (1 | id)"
    },
    vcs=[
        vc(index="id", traits=["BW"], kernel=PED)
    ]
)

fit = gfit(model, data, task="reml")
```

## Model validation and interoperability

Before fitting, qgtools validates the full model bundle (data, formulas,
kernels, and variance components or priors) to ensure internal
consistency.

Model specifications can also be exported as structured JSON, allowing
the same model to be executed by external backends, workflow engines, or
non-R/Python environments.

## Orthogonal layers

qgtools separates model specification into four orthogonal layers:

1.  **Formulas**  
    Describe *what* effects enter the model (responses, fixed effects,
    random-effect indices).

2.  **Kernels**  
    Describe *how* covariance is induced across levels of a random
    effect (e.g. pedigree, genomic relationships, marker structure).

3.  **Variance components or priors**  
    Describe *how much* variation is attributed to each component:

    - variance components (`vc()`) for REML / solver-based estimation
    - prior distributions (`prior()`) for Bayesian inference

4.  **Task**  
    Determines *how* parameters are estimated (`"reml"`, `"solve"`,
    `"bayes"`), without changing the model structure.

These layers are specified independently but validated jointly before
fitting.

> **Important**
>
> A model uses **either** variance components (`vc()`) **or** priors
> (`prior()`), but never both.
>
> - `vc()` is used for REML and solver-based estimation.
> - `prior()` is used for Bayesian hierarchical models.
>
> The same model formulas and kernels can be reused across paradigms.

### Residual effects

The residual variance is represented by a component with
`index = "Residual"`.

- It is **implicit** in model formulas and must **not** appear as
  `(1 | Residual)`.
- It must be explicitly specified using `vc()` (REML / solver) or
  `prior()` (Bayesian).
- Residual components never require a kernel.

## Performance and deployment

qgtools is designed as a lightweight R (or Python) interface to
high-performance computing backends.  
Computationally intensive components (e.g. likelihood evaluation, large
linear algebra, and sampling) can be implemented in compiled languages
such as **C++ or Fortran**, while R (or Python) is used for model
specification and orchestration.

This separation provides: - high computational performance, -
scalability to large datasets, - and a clear boundary between model
definition and numerical implementation.

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
  [ATLAS](https://math-atlas.sourceforge.net/) or
  [MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library))  
- fast and memory-efficient batch processing of genotype data stored in
  binary files (e.g. [PLINK](https://www.cog-genomics.org/plink2)
  bedfiles)
