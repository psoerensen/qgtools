#' Variance component specification
#'
#' Define a variance component (random effect) in a mixed or hierarchical model.
#'
#' @param variable Character. Latent model variable used in (1 | variable).
#' @param traits Character vector of trait names.
#' @param kernel Optional kernel object defining covariance.
#' @param featureMatrix Optional FeatureMatrix / FeatureSource object driving
#'   a high-dimensional effect. Mutually exclusive with kernel.
#' @param featureSets Optional named list defining subsets of features.
#'   Requires featureMatrix.
#' @param structure Character. Trait covariance structure.
#' @param start Optional numeric, matrix, or list of starting values.
#'
#' @return An object of class "vc".
#' @export
vc <- function(variable,
               traits,
               kernel = NULL,
               featureMatrix = NULL,
               featureSets = NULL,
               structure = "unstructured",
               start = NULL) {

  ## ---- variable ----------------------------------------------------------
  if (!is.character(variable) || length(variable) != 1)
    stop("'variable' must be a single character string")

  ## ---- traits ------------------------------------------------------------
  if (!is.character(traits) || length(traits) < 1)
    stop("'traits' must be a non-empty character vector")

  ## ---- residual special case --------------------------------------------
  if (variable == "residual") {
    if (!is.null(kernel) || !is.null(featureMatrix))
      stop("Residual variance must not be kernel- or feature-backed")
  }

  ## ---- kernel / featureMatrix exclusivity --------------------------------
  if (variable != "residual") {
    if (!xor(is.null(kernel), is.null(featureMatrix))) {
      stop(
        "Specify exactly one of 'kernel' or 'featureMatrix' for variance component '",
        variable, "'"
      )
    }
  }

  ## ---- featureMatrix / featureSets --------------------------------------
  if (!is.null(featureMatrix)) {

    if (!inherits(featureMatrix, c("FeatureMatrix", "FeatureSource", "FeatureBackend")))
      stop("'featureMatrix' must be a FeatureMatrix / FeatureSource object")

    if (!is.null(featureSets) && !is.list(featureSets))
      stop("'featureSets' must be a named list")

  } else if (!is.null(featureSets)) {
    stop("'featureSets' can only be used when 'featureMatrix' is specified")
  }

  ## ---- structure ---------------------------------------------------------
  allowed_structures <- c("unstructured", "diagonal", "identity")
  if (!structure %in% allowed_structures)
    stop(
      "'structure' must be one of: ",
      paste(allowed_structures, collapse = ", ")
    )

  ## ---- kernel ------------------------------------------------------------
  if (!is.null(kernel)) {
    kernel <- as_kernel(kernel)
  }

  ## ---- start -------------------------------------------------------------
  if (!is.null(start) &&
      !(is.numeric(start) || is.matrix(start) || is.list(start))) {
    stop("'start' must be numeric, matrix, list, or NULL")
  }

  ## ---- build object ------------------------------------------------------
  structure(
    list(
      variable       = variable,
      traits         = traits,
      kernel         = kernel,
      featureMatrix  = featureMatrix,
      featureSets    = featureSets,
      structure      = structure,
      start          = start
    ),
    class = "vc"
  )
}


#' Variance component collection
#'
#' Collect and validate a set of variance component specifications.
#'
#' @param ... Named variance component objects of class "vc".
#'
#' @return An object of class "varcomp".
#' @export
varcomp <- function(x = NULL, ...) {

  ## ---- collect components --------------------------------------------------
  comps <- if (!is.null(x)) {
    if (!is.list(x))
      stop("If provided, 'x' must be a list of vc objects")
    x
  } else {
    list(...)
  }

  ## ---- empty check ---------------------------------------------------------
  if (length(comps) == 0)
    stop("At least one variance component must be specified")

  ## ---- naming --------------------------------------------------------------
  if (is.null(names(comps)) || any(names(comps) == "")) {
    stop("All variance components must be named")
  }

  ## ---- class checks --------------------------------------------------------
  is_vc <- vapply(comps, inherits, logical(1), "vc")
  if (any(!is_vc)) {
    bad <- names(comps)[!is_vc]
    stop(
      "All elements must be 'vc' objects. Invalid components: ",
      paste(bad, collapse = ", ")
    )
  }

  ## ---- variable sanity -----------------------------------------------------
  variables <- vapply(comps, function(x) x$variable, character(1))

  if (any(duplicated(variables))) {
    dup <- unique(variables[duplicated(variables)])
    warning(
      "Multiple variance components share the same variable: ",
      paste(dup, collapse = ", "),
      "\nThis is allowed (e.g. feature subsets) but should be intentional."
    )
  }

  ## ---- kernel vs featureMatrix consistency ---------------------------------
  by_var <- split(comps, variables)

  for (v in names(by_var)) {
    group <- by_var[[v]]

    kernels  <- vapply(group, function(x) !is.null(x$kernel), logical(1))
    features <- vapply(group, function(x) !is.null(x$featureMatrix), logical(1))

    if (any(kernels) && any(features)) {
      stop(
        "Variance components for variable '", v,
        "' mix kernel-backed and feature-backed definitions.\n",
        "Use exactly one covariance mechanism per variable."
      )
    }
  }

  ## ---- featureMatrix identity checks --------------------------------------
  feature_vcs <- Filter(function(x) !is.null(x$featureMatrix), comps)

  if (length(feature_vcs) > 1) {
    fm_ids <- vapply(
      feature_vcs,
      function(x) paste(class(x$featureMatrix), length(x$featureMatrix), sep = ":"),
      character(1)
    )

    if (length(unique(fm_ids)) > 1) {
      warning(
        "Multiple feature-backed variance components use different featureMatrix objects.\n",
        "Ensure this is intentional."
      )
    }
  }

  ## ---- build object --------------------------------------------------------
  structure(comps, class = "varcomp")
}

#' @export
print.vc <- function(x, ...) {

  cat("Variance component\n")
  cat("  Variable:  ", x$variable, "\n")
  cat("  Traits:    ", paste(x$traits, collapse = ", "), "\n")

  ## ---- covariance mechanism --------------------------------------------
  if (!is.null(x$kernel)) {
    cat("  Covariance:", "kernel-backed\n")
    cat("    Kernel:  ", x$kernel$type, "\n")
  } else if (!is.null(x$featureMatrix)) {
    cat("  Covariance:", "feature-backed\n")
    cat("    Feature matrix class: ",
        class(x$featureMatrix)[1], "\n")

    if (!is.null(x$sets)) {
      cat("    Feature sets: ",
          paste(names(x$sets), collapse = ", "), "\n")
    }
  } else if (x$variable == "residual") {
    cat("  Covariance:", "residual\n")
  } else {
    cat("  Covariance:", "implicit (iid)\n")
  }

  ## ---- trait structure --------------------------------------------------
  cat("  Structure: ", x$structure, "\n")

  ## ---- starting values --------------------------------------------------
  if (!is.null(x$start)) {
    cat("  Start:     provided\n")
  }

  ## ---- prior ------------------------------------------------------------
  if (!is.null(x$prior)) {
    cat("  Prior:     ", class(x$prior$distribution)[1], "\n")
  }

  invisible(x)
}


#' @export
print.varcomp <- function(x, ...) {

  cat("Variance components (", length(x), ")\n", sep = "")

  for (nm in names(x)) {
    cat("\n----------------------------------\n")
    cat("Component: ", nm, "\n", sep = "")
    print(x[[nm]])
  }

  invisible(x)
}


#' Summarize a variance component
#'
#' Produce a structured summary of a variance component, describing
#' how covariance is induced and which traits and features are involved.
#'
#' @param object An object of class \code{"vc"}.
#' @param ... Unused.
#'
#' @return A named list summarizing the variance component.
#' @export
summary.vc <- function(object, ...) {

  stopifnot(inherits(object, "vc"))

  ## ---- identify covariance mechanism -----------------------------------
  cov_type <- if (!is.null(object$kernel)) {
    "kernel"
  } else if (!is.null(object$featureMatrix)) {
    "feature"
  } else if (object$variable == "residual") {
    "residual"
  } else {
    "iid"
  }

  ## ---- base summary -----------------------------------------------------
  out <- list(
    variable   = object$variable,
    traits     = object$traits,
    covariance = cov_type,
    structure  = object$structure
  )

  ## ---- kernel-backed ----------------------------------------------------
  if (cov_type == "kernel") {
    out$kernel <- list(
      type = object$kernel$type,
      id   = object$kernel$meta$id %||% NA_character_
    )
  }

  ## ---- feature-backed --------------------------------------------------
  if (cov_type == "feature") {
    out$feature <- list(
      class = class(object$featureMatrix)[1],
      n_rows = length(object$featureMatrix$ids),
      n_cols = length(object$featureMatrix$rsids)
    )

    if (!is.null(object$sets)) {
      out$feature$sets <- list(
        names = names(object$sets),
        sizes = vapply(object$sets, length, integer(1))
      )
    }
  }

  ## ---- starting values -------------------------------------------------
  if (!is.null(object$start)) {
    out$start <- list(
      provided = TRUE,
      type = class(object$start)[1]
    )
  } else {
    out$start <- list(provided = FALSE)
  }

  ## ---- prior -----------------------------------------------------------
  if (!is.null(object$prior)) {
    out$prior <- list(
      distribution = class(object$prior$distribution)[1]
    )
  }

  class(out) <- "summary.vc"
  out
}

#' @export
print.summary.vc <- function(x, ...) {

  cat("Variance component summary\n")
  cat("  Variable:   ", x$variable, "\n")
  cat("  Traits:     ", paste(x$traits, collapse = ", "), "\n")
  cat("  Covariance: ", x$covariance, "\n")
  cat("  Structure:  ", x$structure, "\n")

  if (!is.null(x$kernel)) {
    cat("  Kernel:\n")
    cat("    Type:", x$kernel$type, "\n")
    if (!is.na(x$kernel$id))
      cat("    ID:  ", x$kernel$id, "\n")
  }

  if (!is.null(x$feature)) {
    cat("  Feature matrix:\n")
    cat("    Class:", x$feature$class, "\n")
    cat("    Rows: ", x$feature$n_rows, "\n")
    cat("    Cols: ", x$feature$n_cols, "\n")

    if (!is.null(x$feature$sets)) {
      cat("    Feature sets:\n")
      for (i in seq_along(x$feature$sets$names)) {
        cat("      - ",
            x$feature$sets$names[i],
            " (", x$feature$sets$sizes[i], ")\n",
            sep = "")
      }
    }
  }

  if (x$start$provided) {
    cat("  Start values: provided (", x$start$type, ")\n", sep = "")
  }

  if (!is.null(x$prior)) {
    cat("  Prior distribution:", x$prior$distribution, "\n")
  }

  invisible(x)
}


# ------------------------------------------------------------------------------
# Kernel constructors (explicit covariance mechanisms)
# ------------------------------------------------------------------------------

# Internal helper: construct a kernel object
new_kernel <- function(type, data = NULL, meta = list()) {

  if (!is.character(type) || length(type) != 1)
    stop("'type' must be a single character string")

  ## ------------------------------------------------------------------
  ## Ensure kernel has a stable ID
  ## ------------------------------------------------------------------
  if (is.null(meta$id)) {
    meta$id <- paste0(
      type, "_",
      formatC(abs(stats::runif(1)), digits = 6, format = "f")
    )
  }

  structure(
    list(
      type = type,
      data = data,
      meta = meta
    ),
    class = "kernel"
  )
}


#' IID kernel (independent and identically distributed)
#'
#' Used for random intercepts, residuals, litter effects, etc.
#'
#' @param label Optional label (e.g. "residual", "dam")
#' @export
iid_kernel <- function(label = "iid") {

  new_kernel(
    type = "IID",
    data = NULL,
    meta = list(
      id    = paste0("IID_", label),
      label = label
    )
  )
}


#' Pedigree-based kernel
#'
#' @param PEDlist A PEDlist object
#' @export
ped_kernel <- function(PEDlist) {

  if (!inherits(PEDlist, "PEDlist"))
    stop("Expected a PEDlist object")

  new_kernel(
    type = "PED",
    data = PEDlist
  )
}


#' Genomic relationship matrix (GRM) kernel
#'
#' @param GRMlist A GRMlist object
#' @export
grm_kernel <- function(GRMlist) {

  if (!inherits(GRMlist, "GRMlist"))
    stop("Expected a GRMlist object")

  new_kernel(
    type = "GRM",
    data = GRMlist
  )
}


#' LD-based kernel
#'
#' @param LDlist An LDlist object
#' @param sparse Logical; whether LD is sparse
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


# ------------------------------------------------------------------------------
# Coercion helper
# ------------------------------------------------------------------------------

#' Coerce kernel specifications to kernel objects
#'
#' @param x Kernel object or supported backend object
#' @export
as_kernel <- function(x) {

  ## already a kernel
  if (inherits(x, "kernel"))
    return(x)

  ## backend objects
  if (inherits(x, "PEDlist")) return(ped_kernel(x))
  if (inherits(x, "GRMlist")) return(grm_kernel(x))
  if (inherits(x, "LDlist"))  return(ld_kernel(x))

  stop(
    "Cannot convert object of class ",
    paste(class(x), collapse = "/"),
    " to a kernel"
  )
}
