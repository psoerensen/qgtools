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

#### Prepare input data

``` r
## Data may be provided as an in-memory data frame or as a disk-backed file
data <- makeDatalist(
  source = "data.txt",
  format = "CSV"
)

## Prepare pedigree kernel for additive genetic effects
PED <- makePEDlist(fnPED = "pedigree.txt", format = "TEXT", method = "S-D-NonInbred")

## Prepare genomic relationship kernel for additive genomic effects
GRM <- makeGRMlist(fnGRM = "grm_inverse.txt", format = "BINARY", grm_type = "G-inverse")

# Marker-level genotype container (e.g. PLINK BED/BIM/FAM)
# This provides information about where genotype data live on disk and how animals map to markers.
Glist <- makeGlist(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
)
```

#### Single and multiple traits examples

``` r
## Single-trait linear mixed model 

## Model formula:
## (1 | id) represents the additive genetic (animal) effect
formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)

## Variance components define how covariance is modeled
vcs <- list(
  dam_env = vc(index = "dam", traits = "BW"),          # dam environmental effect
  animal_genetic = vc(index = "id", traits = "BW",
                       kernel = PED),                  # additive genetic effect
  residual = vc(index = "Residual", traits = "BW")    # residual variance
)

## Fit the model and estimate variance component using REML
fit <- gfit(formulas, data, vcs, task = "reml")

## Multi-trait linear mixed model 

## Same model structure, now for two correlated traits
formulas_mt <- list(
  Gl = Gl ~ sex + reps + (1 | dam) + (1 | id),
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)

## Shared variance components induce correlation between traits
## Shared variance components induce correlation between traits
## For REML, supplied values are used as starting values
vcs_mt <- list(
  dam_env = vc(
    index  = "dam",
    traits = c("Gl", "BW"),
    prior  = prior_start(
      matrix(
        c(1.0, 0.02,
          0.2, 2.0),
        nrow = 2,
        byrow = TRUE
      )
  ),

  animal_genetic = vc(
    index  = "id",
    traits = c("Gl", "BW"),
    kernel = PED,
    prior  = prior_start(
      matrix(
        c(4.0, 0.8,
          0.8, 3.0),
        nrow = 2,
        byrow = TRUE
      )
    )
  ),

  residual = vc(
    index  = "Residual",
    traits = c("Gl", "BW"),
    prior  = prior_diag(c(5, 5))
  )
)

## Fit the model and estimate variance component using REML
fit_mt <- gfit(formulas_mt, data, vcs_mt, task = "reml")
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
infrastructures such as **AWS** and **Azure**. Containers encapsulate
the required libraries and runtimes and can be deployed to batch
systems, Kubernetes clusters, or managed cloud services without exposing
source code.

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
