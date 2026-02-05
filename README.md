qgtools: A Unified Modeling Framework
================

## Overview

**qgtools** provides a unified framework for specifying and fitting
mixed and hierarchical models, with a focus on quantitative genetics,
genomics, and complex covariance structures.

The core idea behind qgtools is to **separate model specification from
model fitting**, allowing the *same model* to be estimated using
different software backends and estimation paradigms within a single,
coherent modeling universe.

**qgtools** handles large-scale data by taking advantage of:

- multi-core processing using [openMP](https://www.openmp.org/)  
- multithreaded matrix operations implemented in BLAS libraries
  (e.g. [OpenBLAS](https://www.openblas.net/),
  [ATLAS](https://math-atlas.sourceforge.net/) or
  [MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library))  
- fast and memory-efficient batch processing of genotype data stored in
  binary files (e.g. [PLINK](https://www.cog-genomics.org/plink2)
  bedfiles)

------------------------------------------------------------------------

## Motivation

Modern quantitative genetics and statistical modeling rely on a diverse
ecosystem of software:

- likelihood-based tools (e.g. REML solvers),
- Bayesian samplers,
- marker-based and pedigree-based methods,
- large-scale and disk-backed data pipelines.

Each tool often comes with its own model syntax, assumptions, and data
requirements, making it difficult to:

- compare methods fairly,
- reuse model specifications,
- scale analyses from small to very large datasets,
- combine classical and modern approaches in a single workflow.

**qgtools addresses this by providing a common, declarative model
layer** that can be fitted using different estimation tasks and
computational engines without redefining the model.

------------------------------------------------------------------------

## Core concepts

### Model formulas define *what* exists

Per-trait model formulas, written using standard Wilkinson syntax,
define:

- response variables (traits),
- fixed effects,
- random effects (via grouping factors).

Formulas specify **which effects exist**, but not how their covariance
is modeled.

------------------------------------------------------------------------

### Variance components define *how* covariance is modeled

Variance components (`vc()` objects) define:

- which random effects induce covariance,
- how covariance is structured across individuals and traits,
- whether covariance is induced explicitly (e.g. pedigree, GRM) or
  implicitly (e.g. marker-based effects),
- optional prior information.

This separation allows covariance assumptions to be modified
independently of the model structure.

------------------------------------------------------------------------

### Kernels are first-class objects

Covariance is induced through **kernel objects**, such as:

- pedigree-based kernels,
- genomic or marker-based kernels,
- identity or IID kernels.

Kernels are standalone, reusable objects that can be shared across
variance components and model fits, enabling consistent modeling across
different analyses and backends.

------------------------------------------------------------------------

### Tasks control *how* the model is fitted

Models are fitted using `gfit()` with a specified **task**:

- `"reml"` – restricted maximum likelihood estimation,
- `"bayes"` – Bayesian inference using priors,
- `"solve"` – fixed-parameter solutions.

The same model specification can be reused across tasks, enabling direct
comparison of estimation paradigms without changing the model.

------------------------------------------------------------------------

### Data as a scalable data source

qgtools treats data as a **data source**, not just an in-memory table.

- Small datasets can be supplied as ordinary data frames.
- Large datasets can be disk-backed or streamed, analogous to how large
  genotype data are handled.

Model fitting accesses only the variables required by the model and does
not assume that all data must reside in memory.

------------------------------------------------------------------------

## Example

``` r
## Pedigree kernel
PED <- makePEDlist(fnPED = mouse_pedigree)

## Single-trait animal model
formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | id)
)

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

fit <- gfit(
  formulas = formulas,
  data     = mouse,
  vcs      = vcs,
  task     = "reml"
)
```
