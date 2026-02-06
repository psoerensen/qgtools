#' Variance component specification
#'
#' Define a variance component (random effect) in a mixed or hierarchical model.
#' Each variance component:
#' \itemize{
#'   \item corresponds to a random effect in the model formulas,
#'   \item is indexed by a grouping factor (\code{index}),
#'   \item uses a kernel to induce covariance across index levels, and
#'   \item may optionally specify starting values and/or prior distributions
#'         for variance or covariance parameters.
#' }
#'
#' Covariance may be induced explicitly (e.g. via pedigree or GRM-based kernels
#' with variance–covariance parameters) or implicitly through model structure,
#' such as marker-based effects where covariance arises from genotypes, LD, and
#' marker-level priors.
#'
#' The \code{index} must match the grouping factor used for the corresponding
#' random effect in the model formula.
#' The \code{index} is the name of the grouping variable and must correspond
#' exactly to the \code{(1 | index)} term in the model formula.
#'
#' @details
#' The interpretation of supplied variance and covariance information depends on
#' the analysis task:
#' \itemize{
#'   \item \strong{REML}: \code{start} values are used as starting values.
#'   \item \strong{Bayesian}: \code{prior} specifies prior distributions.
#'   \item \strong{Solve}: \code{start} values are treated as fixed parameters.
#' }
#'
#' Not all kernels imply explicit variance–covariance matrices. For example,
#' marker-based kernels may induce covariance implicitly without forming such
#' matrices.
#'
#' @param index Character. Grouping variable indexing the random effect.
#' @param traits Character vector of trait names.
#' @param kernel Kernel object defining covariance across index levels.
#' @param structure Character. Declared trait covariance structure
#'   ("unstructured", "diagonal", or "identity"). This parameter is primarily
#'   informational and may be used for validation or default prior construction
#'   in future versions.
#' @param start Optional numeric scalar, vector, or matrix specifying starting
#'   values (or fixed values) for variance or covariance parameters.
#' @param prior Optional prior specification object.
#'
#' @return An object of class \code{"vc"}.
#' @export
vc <- function(index,
               traits,
               kernel = NULL,
               structure = "unstructured",
               start = NULL,
               prior = NULL) {

  ## ---- index ---------------------------------------------------------------
  if (!is.character(index) || length(index) != 1)
    stop("'index' must be a single character string")

  ## ---- traits --------------------------------------------------------------
  if (!is.character(traits) || length(traits) < 1)
    stop("'traits' must be a non-empty character vector")

  ## ---- structure -----------------------------------------------------------
  allowed_structures <- c("unstructured", "diagonal", "identity")
  if (!structure %in% allowed_structures)
    stop(
      "'structure' must be one of: ",
      paste(allowed_structures, collapse = ", ")
    )

  ## ---- kernel --------------------------------------------------------------
  if (is.null(kernel)) {
    kernel <- iid_kernel("default")
  } else {
    kernel <- as_kernel(kernel)
  }

  ## ---- start ---------------------------------------------------------------
  if (!is.null(start) &&
      !(is.numeric(start) || is.matrix(start))) {
    stop("'start' must be numeric (scalar, vector, or matrix) or NULL")
  }

  ## ---- prior ---------------------------------------------------------------
  if (!is.null(prior) && !inherits(prior, "prior"))
    stop("'prior' must be a prior object")

  ## ---- build object --------------------------------------------------------
  structure(
    list(
      index     = index,
      traits    = traits,
      kernel    = kernel,
      structure = structure,
      start     = start,
      prior     = prior
    ),
    class = "vc"
  )
}




new_kernel <- function(type, data = NULL, meta = list()) {
  structure(
    list(
      type = type,
      data = data,
      meta = meta
    ),
    class = c(paste0(type, "_kernel"), "kernel")
  )
}

#' Pedigree-based kernel
#'
#' @param PEDlist A PEDlist object.
#' @export
ped_kernel <- function(PEDlist) {
  if (!inherits(PEDlist, "PEDlist"))
    stop("Expected a PEDlist object")

  new_kernel(
    type = "PED",
    data = PEDlist
  )
}

#' Genomic relationship kernel
#'
#' @param GRMlist A GRMlist object or similar.
#' @export
grm_kernel <- function(GRMlist) {
  if (!inherits(GRMlist, "GRMlist"))
    stop("Expected a GRMlist object")

  new_kernel(
    type = "GRM",
    data = GRMlist
  )
}

#' LD-based marker kernel
#'
#' @param Glist A Glist object.
#' @param sparse Logical. Whether LD is sparse.
#' @export
ld_kernel <- function(LDlist, sparse = TRUE) {
  if (!inherits(LDlist, "LDlist"))
    stop("Expected an LDlist object")

  new_kernel(
    type = "LD",
    data = LDlist,
    meta = list(sparse = sparse)
  )
}

#' IID kernel
#'
#' @param type Character. Label for IID kernel (e.g. "residual", "litter").
#' @export
iid_kernel <- function(type = "iid") {
  new_kernel(
    type = "IID",
    data = NULL,
    meta = list(label = type)
  )
}

#' @export
print.vc <- function(x, ...) {
  cat("Variance component\n")
  cat("  Index:   ", x$index, "\n")
  cat("  Traits:  ", paste(x$traits, collapse = ", "), "\n")
  cat("  Kernel:  ", x$kernel$type, "\n")
  cat("  Structure:", x$structure, "\n")
  if (!is.null(x$prior)) {
    cat("  Prior:   ", x$prior$type, "\n")
  }
  invisible(x)
}


#' Variance component collection
#'
#' Collect and validate a set of variance component specifications.
#' Each component must be created with \code{vc()}.
#'
#' @param ... Named variance component objects of class \code{"vc"}.
#'
#' @return An object of class \code{"varcomp"}.
#' @export
varcomp <- function(x = NULL, ...) {

  ## ---- collect components --------------------------------------------------
  if (!is.null(x)) {
    if (!is.list(x)) {
      stop("If provided, 'x' must be a list of vc objects")
    }
    comps <- x
  } else {
    comps <- list(...)
  }

  ## ---- empty check ---------------------------------------------------------
  if (length(comps) == 0)
    stop("At least one variance component must be specified")

  ## ---- naming --------------------------------------------------------------
  if (is.null(names(comps)) || any(names(comps) == "")) {
    stop("All variance components must be named")
  }

  ## ---- class checks --------------------------------------------------------
  is_vc <- vapply(comps, inherits, logical(1), what = "vc")
  if (any(!is_vc)) {
    bad <- names(comps)[!is_vc]
    stop(
      "All elements must be 'vc' objects. Invalid components: ",
      paste(bad, collapse = ", ")
    )
  }

  ## ---- index sanity (non-fatal, structural) --------------------------------
  indices <- vapply(comps, function(x) x$index, character(1))
  if (any(duplicated(indices))) {
    warning(
      "Multiple variance components share the same index: ",
      paste(unique(indices[duplicated(indices)]), collapse = ", "),
      "\nThis is allowed but should be intentional ",
      "(e.g. multiple kernels on the same index)."
    )
  }

  ## ---- build object --------------------------------------------------------
  structure(comps, class = "varcomp")
}


#' @export
print.varcomp <- function(x, ...) {
  cat("Variance components (", length(x), ")\n", sep = "")
  for (nm in names(x)) {
    cat("\n- ", nm, "\n", sep = "")
    print(x[[nm]])
  }
  invisible(x)
}


# ------------------------------------------------------------------------------
# Coerce kernel specs or kernel objects to kernel objects
# ------------------------------------------------------------------------------

as_kernel <- function(x) {

  ## already a kernel ----------------------------------------------------------
  if (inherits(x, "kernel")) {
    return(x)
  }

  ## kernel specifications -----------------------------------------------------
  if (inherits(x, "PEDlist")) return(ped_kernel(x))
  if (inherits(x, "GRMlist")) return(grm_kernel(x))
  if (inherits(x, "LDlist"))  return(ld_kernel(x))

  stop(
    "Cannot convert object of class ",
    paste(class(x), collapse = "/"),
    " to a kernel"
  )
}

