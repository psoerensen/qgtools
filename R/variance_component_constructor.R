#' Variance component specification
#'
#' Define a variance component (random effect) in a mixed or hierarchical model.
#' A variance component combines:
#'   - an indexing factor (e.g. Animal, Dam, SNP)
#'   - a kernel defining covariance over index levels
#'   - a covariance structure over traits
#'   - an optional prior on the trait covariance
#'
#' @param index Character. Indexing variable for the random effect.
#' @param traits Character vector of trait names.
#' @param kernel Kernel object (e.g. PED, GRM, LD, IID).
#' @param structure Character. Trait covariance structure.
#' @param prior Optional prior specification object.
#'
#' @return An object of class "vc".
#' @export
vc <- function(index,
               traits,
               kernel = NULL,
               structure = "unstructured",
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
varcomp <- function(...) {

  comps <- list(...)

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
      "\nThis is allowed but should be intentional (e.g. multiple kernels on same index)."
    )
  }

  ## ---- build object --------------------------------------------------------
  structure(
    comps,
    class = "varcomp"
  )
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

