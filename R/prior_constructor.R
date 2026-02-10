#' Prior specification for Bayesian hierarchical models
#'
#' Define a prior for a latent model component in a Bayesian model.
#' The structure mirrors \code{vc()}, replacing variance estimation
#' with explicit prior distributions.
#'
#' @param variable Character. Latent model variable (e.g. "animal", "marker").
#' @param traits Character vector of trait names.
#' @param kernel Optional kernel object defining covariance.
#' @param featureMatrix Optional FeatureMatrix / FeatureSource object.
#'   Mutually exclusive with \code{kernel}.
#' @param featureSets Optional named list defining feature subsets.
#'   Requires \code{featureMatrix}.
#' @param structure Character. Trait covariance structure
#'   ("unstructured", "diagonal", "identity").
#' @param distribution Prior distribution object
#'   (e.g. \code{iw()}, \code{bayesC()}, \code{bayesR()}).
#' @param start Optional numeric, matrix, or list (per feature set).
#'
#' @return An object of class \code{"prior"}.
#' @export
prior <- function(variable,
                  traits,
                  kernel = NULL,
                  featureMatrix = NULL,
                  featureSets = NULL,
                  structure = "unstructured",
                  distribution,
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
      stop("Residual prior must not be kernel- or feature-backed")
  }

  ## ---- kernel / featureMatrix exclusivity -------------------------------
  if (variable != "residual") {
    if (!xor(is.null(kernel), is.null(featureMatrix))) {
      stop(
        "Specify exactly one of 'kernel' or 'featureMatrix' for prior '",
        variable, "'"
      )
    }
  }

  ## ---- featureMatrix / featureSets --------------------------------------
  if (!is.null(featureMatrix)) {

    if (!inherits(featureMatrix,
                  c("FeatureMatrix", "FeatureSource", "FeatureBackend")))
      stop("'featureMatrix' must be a FeatureMatrix / FeatureSource object")

    if (!is.null(featureSets) && !is.list(featureSets))
      stop("'featureSets' must be a named list")

  } else if (!is.null(featureSets)) {
    stop("'featureSets' can only be used when 'featureMatrix' is provided")
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

  ## ---- distribution ------------------------------------------------------
  if (missing(distribution))
    stop("A prior 'distribution' must be supplied")

  if (!inherits(distribution, "prior_distribution"))
    stop("'distribution' must be a prior distribution object")

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
      distribution   = distribution,
      start          = start
    ),
    class = "prior"
  )
}

#' Summarize a prior specification
#'
#' Produce a structured summary of a prior object, describing
#' how covariance is induced, which traits and features are involved,
#' and which prior distribution is used.
#'
#' @param object An object of class \code{"prior"}.
#' @param ... Unused.
#'
#' @return A named list summarizing the prior specification.
#' @export
summary.prior <- function(object, ...) {

  stopifnot(inherits(object, "prior"))

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
      class  = class(object$featureMatrix)[1],
      n_rows = length(object$featureMatrix$ids),
      n_cols = length(object$featureMatrix$rsids)
    )

    if (!is.null(object$featureSets)) {
      out$feature$sets <- list(
        names = names(object$featureSets),
        sizes = vapply(object$featureSets, length, integer(1))
      )
    }
  }

  ## ---- prior distribution ----------------------------------------------
  out$distribution <- list(
    type   = object$distribution$type,
    params = object$distribution$params
  )

  ## ---- starting values -------------------------------------------------
  if (!is.null(object$start)) {
    out$start <- list(
      provided = TRUE,
      type     = class(object$start)[1]
    )
  } else {
    out$start <- list(provided = FALSE)
  }

  class(out) <- "summary.prior"
  out
}


#' @export
print.summary.prior <- function(x, ...) {

  cat("Prior summary\n")
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
        cat(
          "      - ",
          x$feature$sets$names[i],
          " (", x$feature$sets$sizes[i], ")\n",
          sep = ""
        )
      }
    }
  }

  cat("  Prior distribution:\n")
  cat("    Type:", x$distribution$type, "\n")

  if (x$start$provided) {
    cat("  Start values: provided (", x$start$type, ")\n", sep = "")
  }

  invisible(x)
}


# ------------------------------------------------------------------------------
# Prior constructors
# ------------------------------------------------------------------------------

#' Base class for prior distributions
#'
#' @export
new_prior_distribution <- function(type, params) {

  if (!is.character(type) || length(type) != 1)
    stop("'type' must be a single character string")

  if (!is.list(params))
    stop("'params' must be a list")

  structure(
    list(
      type   = type,
      params = params
    ),
    class = c(type, "prior_distribution")
  )
}

#' Inverse-Wishart prior
#'
#' @param df Degrees of freedom
#' @param S Scale matrix
#'
#' @return An object of class "iw"
#' @export
iw <- function(df, S) {

  if (!is.numeric(df) || length(df) != 1 || df <= 0)
    stop("'df' must be a positive scalar")

  if (!is.matrix(S))
    stop("'S' must be a matrix")

  if (nrow(S) != ncol(S))
    stop("'S' must be a square matrix")

  new_prior_distribution(
    type = "iw",
    params = list(
      df = df,
      S  = unclass(S)
    )
  )
}

#' BayesC prior for marker effects
#'
#' @param pi Mixing proportions (e.g. c(0.95, 0.05))
#' @param gamma Variance scaling parameters
#' @param sets Optional marker sets
#' @param weights Optional marker weights
#' @param annotation Optional functional annotation
#'
#' @return An object of class "bayesC"
#' @export
bayesC <- function(pi,
                   gamma,
                   sets = NULL,
                   weights = NULL,
                   annotation = NULL) {

  if (!is.numeric(pi) || abs(sum(pi) - 1) > 1e-8)
    stop("'pi' must sum to 1")

  if (!is.numeric(gamma))
    stop("'gamma' must be numeric")

  new_prior_distribution(
    type = "bayesC",
    params = list(
      pi         = pi,
      gamma      = gamma,
      sets       = sets,
      weights    = weights,
      annotation = annotation
    )
  )
}

#' BayesR prior for marker effects
#'
#' @param pi Mixing proportions
#' @param variances Vector of component variances
#' @param sets Optional marker sets
#' @param weights Optional marker weights
#' @param annotation Optional functional annotation
#'
#' @return An object of class "bayesR"
#' @export
bayesR <- function(pi,
                   variances,
                   sets = NULL,
                   weights = NULL,
                   annotation = NULL) {

  if (!is.numeric(pi) || abs(sum(pi) - 1) > 1e-8)
    stop("'pi' must sum to 1")

  if (!is.numeric(variances))
    stop("'variances' must be numeric")

  if (length(pi) != length(variances))
    stop("'pi' and 'variances' must have same length")

  new_prior_distribution(
    type = "bayesR",
    params = list(
      pi         = pi,
      variances  = variances,
      sets       = sets,
      weights    = weights,
      annotation = annotation
    )
  )
}
