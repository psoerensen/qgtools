#' Prior specification for Bayesian hierarchical models
#'
#' Define a prior for a model component in a Bayesian linear or mixed model.
#' The structure mirrors \code{vc()}, but replaces variance estimation with
#' explicit prior distributions.
#'
#' Each prior:
#' \itemize{
#'   \item corresponds to a latent model component,
#'   \item is indexed by a grouping factor (\code{index}),
#'   \item optionally uses a kernel to induce covariance across index levels,
#'   \item specifies a prior distribution for variance, covariance, or effects,
#'   \item may include starting values for MCMC initialization.
#' }
#'
#' The \code{index} must correspond either to:
#' \itemize{
#'   \item a grouping factor in the model formula (e.g. \code{(1 | id)}), or
#'   \item a latent index (e.g. \code{"marker"} or \code{"Residual"}).
#' }
#'
#' @details
#' In Bayesian models:
#' \itemize{
#'   \item \code{kernel} defines covariance across rows/levels,
#'   \item \code{structure} defines covariance across traits,
#'   \item \code{distribution} defines the prior on model parameters.
#' }
#'
#' Marker-based priors (e.g. BayesC, BayesR) typically operate on effect sizes
#' rather than explicit variance–covariance matrices. In these cases,
#' covariance is induced implicitly via the kernel and marker design.
#'
#' @param index Character. Grouping variable or latent index.
#' @param traits Character vector of trait names.
#' @param kernel Optional kernel object defining covariance across index levels.
#' @param structure Character. Trait covariance structure
#'   ("unstructured", "diagonal", or "identity").
#' @param distribution Prior distribution object (e.g. \code{iw()}, \code{bayesC()}).
#' @param start Optional numeric scalar, vector, or matrix used as MCMC starting values.
#'
#' @return An object of class \code{"prior"}.
#' @export
prior <- function(index,
                  traits,
                  kernel = NULL,
                  structure = "unstructured",
                  distribution,
                  start = NULL) {

  ## ---- index -------------------------------------------------------------
  if (!is.character(index) || length(index) != 1)
    stop("'index' must be a single character string")

  ## ---- traits ------------------------------------------------------------
  if (!is.character(traits) || length(traits) < 1)
    stop("'traits' must be a non-empty character vector")

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
      !(is.numeric(start) || is.matrix(start))) {
    stop("'start' must be numeric (scalar, vector, or matrix) or NULL")
  }

  ## ---- build object ------------------------------------------------------
  structure(
    list(
      index        = index,
      traits       = traits,
      kernel       = kernel,
      structure    = structure,
      distribution = distribution,
      start        = start
    ),
    class = "prior"
  )
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
