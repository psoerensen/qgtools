#' Variance component specification
#'
#' Define a variance component associated with a latent model variable
#' (random effect) in a mixed or hierarchical model.
#'
#' A variance component:
#' \itemize{
#'   \item references a latent indexing variable (as used in model formulas),
#'   \item attaches exactly one covariance mechanism (kernel or featureMatrix),
#'   \item defines trait-level covariance structure,
#'   \item optionally provides starting values.
#' }
#'
#' Multiple variance components may reference the same \code{variable}.
#' This allows covariance decomposition of a single random effect into
#' multiple components (e.g. pedigree + marker-based variation).
#'
#' For example, additive genetic variation indexed by \code{id}
#' may include:
#' \itemize{
#'   \item a pedigree-based kernel component, and
#'   \item one or more feature-based components.
#' }
#'
#' The residual component must use \code{variable = "residual"} and
#' must not be kernel- or feature-backed.
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

    # if (any(kernels) && any(features)) {
    #   stop(
    #     "Variance components for variable '", v,
    #     "' mix kernel-backed and feature-backed definitions.\n",
    #     "Use exactly one covariance mechanism per variable."
    #   )
    # }
    ## ---- covariance mechanism diagnostics -----------------------------------
    by_var <- split(comps, variables)

    for (v in names(by_var)) {

      group <- by_var[[v]]

      cov_types <- vapply(group, function(x) {
        if (!is.null(x$kernel)) return("kernel")
        if (!is.null(x$featureMatrix)) return("feature")
        if (x$variable == "residual") return("residual")
        return("iid")
      }, character(1))

      if (length(unique(cov_types)) > 1) {
        warning(
          "Variance components for variable '", v,
          "' use multiple covariance mechanisms: ",
          paste(unique(cov_types), collapse = ", "),
          "\nThis is allowed but should be intentional."
        )
      }
    }


  }

  ## ---- featureMatrix identity checks --------------------------------------
  feature_vcs <- Filter(function(x) !is.null(x$featureMatrix), comps)

  if (length(feature_vcs) > 1) {
    # fm_ids <- vapply(
    #   feature_vcs,
    #   function(x) paste(class(x$featureMatrix), length(x$featureMatrix), sep = ":"),
    #   character(1)
    # )
    fm_ids <- vapply(
      feature_vcs,
      function(x) {
        if (!is.null(x$featureMatrix$id)) {
          x$featureMatrix$id
        } else {
          class(x$featureMatrix)[1]
        }
      },
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
    cat("  Covariance: kernel-backed\n")
    cat("    Kernel type: ", x$kernel$type, "\n")
  } else if (!is.null(x$featureMatrix)) {
    cat("  Covariance: feature-backed\n")
    cat("    Feature class: ", class(x$featureMatrix)[1], "\n")

    if (!is.null(x$featureSets)) {
      cat("    Feature sets: ",
          paste(names(x$featureSets), collapse = ", "), "\n")
    }
  } else if (x$variable == "residual") {
    cat("  Covariance: residual\n")
  } else {
    cat("  Covariance: implicit (iid)\n")
  }

  cat("  Structure: ", x$structure, "\n")

  if (!is.null(x$start))
    cat("  Start:     provided\n")

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

  ## ---- covariance type --------------------------------------------------
  cov_type <- if (!is.null(object$kernel)) {
    "kernel"
  } else if (!is.null(object$featureMatrix)) {
    "feature"
  } else if (object$variable == "residual") {
    "residual"
  } else {
    "iid"
  }

  out <- list(
    variable   = object$variable,
    traits     = object$traits,
    covariance = cov_type,
    structure  = object$structure
  )

  ## ---- kernel-backed ----------------------------------------------------
  if (cov_type == "kernel") {

    kernel_id <- if (!is.null(object$kernel$meta$id)) {
      object$kernel$meta$id
    } else {
      NA_character_
    }

    out$kernel <- list(
      type = object$kernel$type,
      id   = kernel_id
    )
  }

  ## ---- feature-backed ---------------------------------------------------
  if (cov_type == "feature") {

    out$feature <- list(
      class       = class(object$featureMatrix)[1],
      id_type     = object$featureMatrix$ids$type %||% NA_character_,
      column_type = object$featureMatrix$rsids$type %||% NA_character_
    )

    if (!is.null(object$featureSets)) {
      out$feature$sets <- list(
        names = names(object$featureSets),
        sizes = vapply(object$featureSets, length, integer(1))
      )
    }
  }

  out$start <- list(
    provided = !is.null(object$start),
    type     = if (!is.null(object$start)) class(object$start)[1] else NA_character_
  )

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
    cat("    Class:      ", x$feature$class, "\n")
    if (!is.na(x$feature$id_type))
      cat("    Row ID type:", x$feature$id_type, "\n")
    if (!is.na(x$feature$column_type))
      cat("    Col ID type:", x$feature$column_type, "\n")

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

  if (x$start$provided)
    cat("  Start values: provided (", x$start$type, ")\n", sep = "")

  invisible(x)
}

#' Summarize a collection of variance components
#'
#' Produce a structured summary of a varcomp object,
#' including per-component summaries and collection-level diagnostics.
#'
#' @param object An object of class "varcomp".
#' @param ... Unused.
#'
#' @return An object of class "summary.varcomp".
#' @export
summary.varcomp <- function(object, ...) {

  stopifnot(inherits(object, "varcomp"))

  comps <- object
  comp_names <- names(comps)

  ## ---- per-component summaries -----------------------------------------
  summaries <- lapply(comps, summary.vc)
  names(summaries) <- comp_names

  ## ---- collection diagnostics -------------------------------------------
  variables <- vapply(comps, function(x) x$variable, character(1))

  covariance_types <- vapply(
    summaries,
    function(x) x$covariance,
    character(1)
  )

  out <- list(
    n_components      = length(comps),
    component_names   = comp_names,
    n_unique_variables = length(unique(variables)),
    duplicated_variables = unique(variables[duplicated(variables)]),
    covariance_types  = unique(covariance_types),
    components        = summaries
  )

  class(out) <- "summary.varcomp"
  out
}


#' @export
print.summary.varcomp <- function(x, ...) {

  cat("Variance component collection summary\n")
  cat("-------------------------------------\n")
  cat("  Number of components: ", x$n_components, "\n")
  cat("  Unique variables:     ", x$n_unique_variables, "\n")

  if (length(x$duplicated_variables) > 0) {
    cat("  Multiple components per variable:\n")
    for (v in x$duplicated_variables) {
      cat("    - ", v, "\n", sep = "")
    }
  }

  cat("  Covariance mechanisms used: ",
      paste(x$covariance_types, collapse = ", "), "\n")

  cat("\nComponents:\n")

  for (nm in x$component_names) {
    cat("\n----------------------------------\n")
    cat("Component: ", nm, "\n", sep = "")
    print(x$components[[nm]])
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
