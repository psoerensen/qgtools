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

data <- makeDatalist(
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
data <- makeDatalist(
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
backend.

``` r
## Genotype feature matrix (e.g. PLINK BED/BIM/FAM)
M <- featureMatrix(
  bedfiles = "chr.bed",
  bimfiles = "chr.bim",
  famfiles = "chr.fam"
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
## Formulas declare *which* effects enter the model.
## 'animal' represents the additive genetic effect (pedigree-based),
## 'dam' captures a maternal environmental component, and
## 'marker' represents high-dimensional marker effects
## modeled through a feature matrix.

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

The same model specification can be estimated using a marker BLUP
formulation by replacing `prior()` with `vc()` and setting
`task = "solve"`. The formulas, kernels, and feature matrices remain
unchanged.

Additional worked examples for REML/solve and Bayesian mixed models
(single- and multi-trait) are provided in the vignettes.

## Case study: Genotype × feeding system interaction (DMU example)

This example reproduces a classical DMU analysis of genotype × feeding
system interactions using **qgtools**, while illustrating how the same
model specification can be estimated using REML, solver-based BLUP, or
Bayesian MCMC.

The data consist of performance records on young bulls tested under two
feeding systems. Measurements under the two systems are treated as
**different traits**, allowing genotype × environment interaction to be
expressed through cross-trait covariance.

### Model description

For each trait (feeding system), the model is:

``` math

y = \text{Month} + \text{BreedYear} + \text{age} + \text{SystemBySire} + e
```

- `Month`, `BreedYear`, and `age` are fixed effects (with `age` as a
  covariate)
- `SystemBySire` is a random sire effect, nested within feeding system
- Sires are assumed **unrelated** (no pedigree or genomic kernel)
- Genotype × feeding system interaction is captured via **multi-trait
  covariance** of sire effects

### Data input

``` r
## Observational data containing:
## - ADG_sys1, ADG_sys2 (traits for feeding system 1 and 2)
## - Month, BreedYear, age (covariates)
## - SystemBySire (sire-by-system identifier)

data <- makeDatalist(
  source = "bulls_data.txt",
  format = "CSV"
)

## Each feeding system is treated as a separate trait

formulas <- list(
  ADG_sys1 = ADG_sys1 ~ Month + BreedYear + age + (1 | SystemBySire),
  ADG_sys2 = ADG_sys2 ~ Month + BreedYear + age + (1 | SystemBySire)
)

vcs <- list(
  sire = vc(
    variable = "SystemBySire",
    traits   = c("ADG_sys1", "ADG_sys2"),
    start    = matrix(
      c(1.0, 0.3,
        0.3, 1.2),
      nrow = 2,
      byrow = TRUE
    )
  ),
  residual = vc(
    variable = "residual",
    traits   = c("ADG_sys1", "ADG_sys2"),
    start    = diag(2)
  )
)

fit_reml <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "reml"
)

## Uses the same model and variance components
## but solves the mixed model equations without updating variances

fit_solve <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "solve"
)

priors <- list(
  sire = prior(
    variable     = "SystemBySire",
    traits       = c("ADG_sys1", "ADG_sys2"),
    distribution = iw(df = 4, S = diag(2)),
    start        = diag(2)
  ),
  residual = prior(
    variable     = "residual",
    traits       = c("ADG_sys1", "ADG_sys2"),
    distribution = iw(df = 4, S = diag(2)),
    start        = diag(2)
  )
)

fit_bayes <- gfit(
  formulas = formulas,
  data     = data,
  priors   = priors,
  task     = "bayes"
)
```

## Case study: Direct and maternal genetic effects in sheep growth (DMU Example II)

This example reproduces *DMU Example II*, which analyzes growth traits
in sheep using a multi-trait mixed model with **direct and maternal
additive genetic effects**, environmental random effects, and
pedigree-based covariance.

The dataset consists of field-recorded growth traits with a pedigree
containing both animals and their dams. Two traits are analyzed jointly.

### Model description

For each trait, the model is:

``` math

y = \text{Month}
  + \text{DamAge}
  + \text{Litter}
  + \text{Sex}
  + \text{HY}
  + \text{L\_Dam}
  + \text{Dam\_Ge}
  + \text{A\_Ge}
  + e
```

where:

- `L_Dam` is a random environmental effect of litter within dam
- `Dam_Ge` is the **maternal additive genetic effect**
- `A_Ge` is the **direct additive genetic effect** of the animal
- Genetic effects are correlated across traits through a pedigree
- Residuals are correlated across traits

This corresponds exactly to the model fitted in DMU.

------------------------------------------------------------------------

### Data input

``` r
## Observational data containing:
## - traits: V1, V2 (e.g. weights at 2 and 4 months)
## - fixed effects: Month, Damage, Litter, Sex, HY
## - random-effect identifiers: A (animal), Dam, P (litter within dam)

data <- makeDatalist(
  source = "sheep_growth.txt",
  format = "CSV"
)

## Two-trait model

formulas <- list(
  V1 = V1 ~ Month + Damage + Litter + Sex + HY +
        (1 | L_Dam) + (1 | Dam_Ge) + (1 | A_Ge),

  V2 = V2 ~ Month + Damage + Litter + Sex + HY +
        (1 | L_Dam) + (1 | Dam_Ge) + (1 | A_Ge)
)

## Pedigree used for both direct and maternal genetic effects

PED <- makePEDlist(
  fnPED  = "pedigree.txt",
  method = "S-D-NonInbred"
)

vcs <- list(
  litter_env = vc(
    variable = "L_Dam",
    traits   = c("V1", "V2")
  ),

  maternal_genetic = vc(
    variable = "Dam_Ge",
    traits   = c("V1", "V2"),
    kernel   = PED
  ),

  direct_genetic = vc(
    variable = "A_Ge",
    traits   = c("V1", "V2"),
    kernel   = PED
  ),

  residual = vc(
    variable = "residual",
    traits   = c("V1", "V2")
  )
)

fit_reml <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "reml"
)

## Solve mixed model equations with fixed variance components

fit_solve <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "solve"
)
```

## Case study: Random regression model for GH challenge profile (DMU RR example)

This example reproduces the DMU random regression (RR) analysis of
plasma growth hormone (GH) concentration measured repeatedly over time
after stimulation.

The model describes the GH profile using **normalized Legendre
polynomials** (as covariates `L1`, `L2`, `L3`) and includes both:

- a **permanent environmental (PE)** animal effect, and
- an **additive genetic (animal)** effect using the pedigree
  relationship matrix.

The same animal identifier enters the model twice, once for each
component, exactly as in DMU.

### Data and covariates

Blood samples are taken at planned times after stimulation (with small
deviations), so time is modeled using covariates derived from normalized
time.

The dataset contains (at minimum):

- response: `GH`
- fixed class effects: `yob`, `breed`, `p_age`, `td`
- fixed covariate: `age`
- Legendre covariates: `L1`, `L2`, `L3`
- animal id: `id`
- pedigree file defining the relationship matrix among animals

### Prepare input data and kernels

``` r
## Observational data (repeated records per animal)
data <- makeDatalist(
  source = "gh_profile.txt",
  format = "CSV"
)

## Pedigree kernel for additive genetic effects
PED <- makePEDlist(
  fnPED  = "pedigree.txt",
  method = "S-D-NonInbred"
)

## Two latent components use the same animal id:
## - pe: permanent environmental RR effect (iid across animals)
## - animal: additive genetic RR effect (pedigree kernel)

## In practice, pe_id and animal_id can both be the 'id' column in the data;
## they are separated conceptually by using different variable names below.

formulas <- list(
  GH = GH ~ yob + breed + p_age + td + age +
    (1  | pe)     + (L1 + L2 + L3 | pe) +
    (1  | animal) + (L1 + L2 + L3 | animal)
)

vcs <- list(
  ## Permanent environmental random regression (iid across animals)
  pe = vc(
    variable = "pe",
    traits   = "GH",
    ## iid (no kernel) → corresponds to G_PE ⊗ I
    start    = diag(4) * 0.5
  ),

  ## Additive genetic random regression (pedigree kernel)
  animal = vc(
    variable = "animal",
    traits   = "GH",
    kernel   = PED,
    ## corresponds to G_A ⊗ A
    start    = diag(4) * 1.0
  ),

  ## Residual variance (single trait)
  residual = vc(
    variable = "residual",
    traits   = "GH",
    start    = 1.0
  )
)

fit_reml <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "reml"
)

## Solve mixed model equations with fixed (co)variance components
fit_solve <- gfit(
  formulas = formulas,
  data     = data,
  vc       = vcs,
  task     = "solve"
)

priors <- list(
  pe = prior(
    variable     = "pe",
    traits       = "GH",
    distribution = iw(df = 6, S = diag(4)),
    start        = diag(4)
  ),

  animal = prior(
    variable     = "animal",
    traits       = "GH",
    kernel       = PED,
    distribution = iw(df = 6, S = diag(4)),
    start        = diag(4)
  ),

  residual = prior(
    variable     = "residual",
    traits       = "GH",
    distribution = invchisq(df = 4, scale = 1.0),
    start        = 1.0
  )
)

fit_bayes <- gfit(
  formulas = formulas,
  data     = data,
  priors   = priors,
  task     = "bayes"
)
```
