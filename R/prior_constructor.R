
# ------------------------------------------------------------------------------
# Prior constructors
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Internal prior constructor
# ------------------------------------------------------------------------------

new_prior <- function(type, ...) {

  if (!is.character(type) || length(type) != 1)
    stop("'type' must be a single character string")

  obj <- list(
    type   = type,
    params = list(...)
  )

  class(obj) <- c(paste0(type, "_prior"), "prior")
  obj
}


#' Inverse-Wishart prior for trait covariance
#'
#' @param df Degrees of freedom.
#' @param S Scale matrix.
#' @export
prior_iw <- function(df, S) {

  if (!is.numeric(df) || length(df) != 1 || df <= 0)
    stop("'df' must be a positive scalar")

  if (!is.matrix(S))
    stop("'S' must be a matrix")

  if (nrow(S) != ncol(S))
    stop("'S' must be a square matrix")

  if (!isTRUE(all.equal(S, t(S))))
    stop("'S' must be symmetric")

  if (any(diag(S) <= 0))
    stop("Diagonal elements of 'S' must be positive")

  new_prior(
    type = "inverse_wishart",
    df = df,
    S  = S
  )
}


#' Diagonal covariance prior
#'
#' @param var Numeric vector of variances.
#' @export
prior_diag <- function(var) {

  if (!is.numeric(var) || is.matrix(var))
    stop("'var' must be a numeric vector")

  if (any(var <= 0))
    stop("All variances must be positive")

  new_prior(
    type = "diagonal",
    var = as.numeric(var)
  )
}


#' Fixed covariance prior
#'
#' @param Sigma Fixed covariance matrix.
#' @export
prior_fixed <- function(Sigma) {

  if (!is.matrix(Sigma))
    stop("'Sigma' must be a matrix")

  if (nrow(Sigma) != ncol(Sigma))
    stop("'Sigma' must be a square matrix")

  if (!isTRUE(all.equal(Sigma, t(Sigma))))
    stop("'Sigma' must be symmetric")

  new_prior(
    type = "fixed",
    Sigma = Sigma
  )
}


# ------------------------------------------------------------------------------
# BayesC prior (spike-and-slab, explicit zero component)
# ------------------------------------------------------------------------------

prior_bayesC <- function(pi,
                         var,
                         estimate_pi = FALSE) {

  if (!is.numeric(pi) || any(pi < 0))
    stop("'pi' must be a non-negative numeric vector")

  if (abs(sum(pi) - 1) > 1e-8)
    stop("'pi' must sum to 1")

  if (!is.numeric(var))
    stop("'var' must be numeric")

  if (length(pi) != length(var))
    stop("'pi' and 'var' must have the same length")

  if (length(var) != 2)
    stop("BayesC requires exactly two mixture components: spike (0) and slab")

  if (sum(var == 0) != 1)
    stop("BayesC requires exactly one zero-variance component")

  if (any(var < 0))
    stop("'var' must be non-negative")

  new_prior(
    type = "bayesC",
    pi = pi,
    var = var,
    estimate_pi = estimate_pi
  )
}


# ------------------------------------------------------------------------------
# BayesR prior (mixture of normals, explicit zero component)
# ------------------------------------------------------------------------------

prior_bayesR <- function(pi,
                         var,
                         estimate_pi = FALSE) {

  if (!is.numeric(pi) || any(pi < 0))
    stop("'pi' must be a non-negative numeric vector")

  if (abs(sum(pi) - 1) > 1e-8)
    stop("'pi' must sum to 1")

  if (!is.numeric(var))
    stop("'var' must be numeric")

  if (length(pi) != length(var))
    stop("'pi' and 'var' must have the same length")

  if (sum(var == 0) != 1)
    stop("Exactly one element of 'var' must be 0 (spike at zero)")

  if (any(var < 0))
    stop("'var' must be non-negative")

  new_prior(
    type = "bayesR",
    pi = pi,
    var = var,
    estimate_pi = estimate_pi
  )
}


# ------------------------------------------------------------------------------
# Generic prior print method
# ------------------------------------------------------------------------------

#' @export
print.prior <- function(x, ...) {
  cat("Prior specification\n")
  cat("  Type:", x$type, "\n")
  invisible(x)
}


.print_mixture_prior <- function(x, model_name) {

  pi  <- x$params$pi
  var <- x$params$var

  spike <- which(var == 0)

  cat(model_name, "prior\n", sep = "")
  cat("  Mixture components:", length(var), "\n")

  for (k in seq_along(var)) {
    if (k == spike) {
      cat(sprintf(
        "   [%d] spike at 0      : pi = %.3f\n",
        k, pi[k]
      ))
    } else {
      cat(sprintf(
        "   [%d] N(0, %.4g) : pi = %.3f\n",
        k, var[k], pi[k]
      ))
    }
  }

  if (isTRUE(x$params$estimate_pi)) {
    cat("  Mixing proportions: estimated (pi-hyperprior)\n")
  } else {
    cat("  Mixing proportions: fixed\n")
  }

  invisible(x)
}

# ------------------------------------------------------------------------------
# Print BayesC prior
# ------------------------------------------------------------------------------
#' @export
print.bayesC <- function(x, ...) {
  .print_mixture_prior(x, model_name = "BayesC")
}

# ------------------------------------------------------------------------------
# Print BayesR prior
# ------------------------------------------------------------------------------
#' @export
print.bayesR <- function(x, ...) {
  .print_mixture_prior(x, model_name = "BayesR")
}

