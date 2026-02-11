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
    or latent effects (e.g. `(1 | id)`).  
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
    or `"bayes"`, without changing the model structure.

5.  **Data — *where* the information comes from**  
    Data are treated as abstract data sources, enabling transparent use
    of in-memory or disk-backed datasets (e.g. PLINK files, HDF5
    matrices).

All layers are specified independently and validated jointly at fit
time.

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

## Simple example to illustrate the concept

The following example shows a minimal linear mixed model and illustrates
how qgtools separates **data input**, **model structure**, **covariance
specification**, and **estimation**.

``` r
## ---------------------------------------------------------------
## Data input
## ---------------------------------------------------------------
## Observational data provide phenotypes and covariates referenced
## in the model formulas.

data <- makeDataSource(
  source = "data.txt",
  format = "CSV"
)

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

PED <- makePEDlist(fnPED = "pedigree.txt", 
                   method="S-D-NonInbred")

vcs <- list(
  animal = vc(
    variable = "id",
    traits   = "BW",
    kernel   = PED,
    start   = 1.0
  ),
  residual = vc(
    variable = "residual",
    traits   = "BW",
    start   = 1.0
  )
)

## ---------------------------------------------------------------
## Model fitting
## ---------------------------------------------------------------
## The estimation task determines how the model is fit,
## without changing the model structure.

fit <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "reml"
)
```

The same model structure can be estimated as a Bayesian model by
replacing *vc()* with *prior()* and setting *task=“bayes”*, without
changing the formulas, kernels or data.

## Model validation and interoperability

Before fitting, qgtools validates the full model bundle (data, formulas,
kernels or features, and variance components or priors) to ensure
internal consistency.

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
high-performance computational backends. Computationally intensive
components—such as likelihood evaluation, large-scale linear algebra,
and sampling—are implemented in compiled languages such as **C++ or
Fortran**, while R (or Python) is used for model specification and
orchestration.

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
    kernel   = PED,
    start   = 1.0
  ),
  residual = vc(
    variable = "residual",
    traits   = "BW",
    start   = 1.0
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
        vc(variable="id", traits=["BW"], kernel=PED, start=[1.0]),
        vc(variable="residual", traits=["BW"], start=[1.0])
    ]
)

fit = gfit(model, data, task="reml")
```

Both interfaces construct the same internal model representation and
invoke the same computational backend. Differences between R and Python
are limited to syntax and user-facing language conventions.

## Data input and preparation

qgtools distinguishes between three broad classes of input data:
standard observational data, structural data used to define covariance,
and high-dimensional feature data. These inputs are prepared
independently and combined later during model specification.

### Standard data: responses and covariates

Standard data include phenotypes, covariates, and other observed
variables referenced directly in model formulas.

Data may be provided either as an in-memory object or as a disk-backed
source.

``` r
## Observational data (phenotypes and covariates)
data <- makeDataSource(
  source = "data.txt",
  format = "CSV"
)
```

### Structural data: kernels

Kernels describe how covariance is induced across levels of a model
component. They are typically used for low-dimensional random effects,
such as pedigree- based genetic effects.

``` r
## Pedigree-based kernel for additive genetic effects
PED <- makePEDlist(
  fnPED  = "pedigree.txt",
  method = "S-D-NonInbred"
)
```

Kernel objects are lightweight descriptors and do not store explicit
covariance matrices unless required by a backend.

Precomputed relationship matrices may also be supplied when appropriate,
but are not required for most workflows.

``` r
## Optional: precomputed genomic relationship matrix
GRM <- makeGRMlist(
  fnGRM    = "grm_inverse.txt",
  format   = "BINARY",
  grm_type = "G-inverse"
)
```

### High-dimensional predictors: feature matrices

Feature matrices define explicit linear predictors with potentially very
large numbers of coefficients. Typical examples include genotype
matrices, transcriptomic measurements, or other omics data.

Feature data are often stored on disk and accessed lazily by the
backend. Here we illustrate this using PLINK based genotype files:

``` r
## Genotype feature matrix (e.g. PLINK BED/BIM/FAM)
featureMatrix <- makeGlist(
  bedfiles = c("chr1.bed","chr2.bed"),
  bimfiles = c("chr1.bim","chr2.bim"),
  famfiles = c("chr1.fam","chr2.fam")
)
```

Optional feature groupings can be supplied to define multiple variance
components or prior structures over subsets of features.

``` r
## Optional grouping of features
featureSets <- list(
  set1 = 1:1000,
  set2 = 1001:2000
)
```

Feature matrices and feature sets are later linked to model components
via *vc()* (REML / solver) or *prior()* (Bayesian inference).

#### Bayesian multi-component linear models with marker-level priors

Marker-level effects differ fundamentally from classical random effects:

- They typically involve **thousands to millions of coefficients**
- Covariance is induced **implicitly** through the feature matrix
  (e.g. genotypes), linkage disequilibrium, and the chosen prior
- Explicit covariance matrices (e.g. GRMs) are usually avoided for
  scalability

In **qgtools**, marker effects are treated as **feature-based model
components**. They are declared in the model formula as a named effect
(e.g. `(1 | id)`), and linked to genotype data via a `featureMatrix()`.

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
## Formulas declare *which* effects enter the model.
## 'animal' represents the additive genetic effect (pedigree-based),
## 'dam' captures a maternal environmental component, and
## 'marker' represents high-dimensional marker effects
## modeled through a feature matrix.

formulas <- list(
  BW = BW ~ sex + reps + (1 | dam) + (1 | animal),
  Gl = Gl ~ sex + reps + (1 | dam) + (1 | animal)
)

## ------------------------------------------------------------------
## Kernel-based effect: pedigree (random effect)
## ------------------------------------------------------------------
## Pedigree effects are naturally expressed via kernels.

## Pedigree-based kernel for additive genetic effects
PED <- makePEDlist(
  fnPED  = "pedigree.txt",
  method = "S-D-NonInbred"
)


## ------------------------------------------------------------------
## Feature-based effect: markers (high-dimensional predictors)
## ------------------------------------------------------------------
## Marker effects are defined through a feature matrix

## Genotype feature matrix (e.g. PLINK BED/BIM/FAM)
featureMatrix <- makeGlist(
  bedfiles = c("chr1.bed","chr2.bed"),
  bimfiles = c("chr1.bim","chr2.bim"),
  famfiles = c("chr1.fam","chr2.fam")
)

## Optional grouping of markers used for multiple variance components
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
    variable     = "animal",
    traits       = c("BW", "Gl"),
    featureMatrix     = featureMatrix,
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

The same model specification can be estimated using a marker BLUP
formulation by replacing `prior()` with `vc()` and setting
`task = "solve"`. The formulas, kernels, and feature matrices remain
unchanged.

Additional worked examples for REML/solve and Bayesian mixed models
(single- and multi-trait) are provided in the vignettes.

## Classical quantitative genetics case studies (DMU examples)

qgtools is designed to express a wide range of classical quantitative
genetic models while cleanly separating model structure from estimation
strategy.

To demonstrate compatibility with established workflows, several
standard examples from **DMU** are reproduced below using qgtools. In
each case, the *same model specification* can be estimated using REML,
solver-based BLUP, or Bayesian inference by changing only the estimation
task.

### DMU-based examples

- **Genotype × feeding system interaction**  
  Multi-trait model with unrelated sires and G×E expressed through
  cross-trait covariance.  
  [View example
  →](https://psoerensen.github.io/qgtools/examples/dmu_examples/dmu_example_1.html)

- **Direct and maternal genetic effects in sheep growth**  
  Multi-trait pedigree model with direct genetic, maternal genetic, and
  environmental effects.  
  [View example
  →](https://psoerensen.github.io/qgtools/examples/dmu_examples/dmu_example_2.html)

- **Random regression model for growth hormone profiles**  
  Longitudinal random regression model with permanent environmental and
  additive genetic effects using Legendre polynomials.  
  [View example
  →](https://psoerensen.github.io/qgtools/examples/dmu_examples/dmu_example_rr.html)

These examples illustrate how complex legacy DMU models translate
directly into qgtools while benefiting from a clearer and more modular
model specification.
