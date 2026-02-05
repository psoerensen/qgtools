#' Fit a qgtools model
#'
#' Fit a mixed or hierarchical model defined by per-trait formulas and variance
#' components using a specified estimation task.
#'
#' Model formulas define which fixed and random effects are included for each
#' trait using standard Wilkinson syntax. Variance components (\code{vc()} objects)
#' define how covariance is induced for those random effects via kernels and
#' optional priors.
#'
#' The data supplied to \code{gfit()} may be an in-memory object (e.g. a data frame)
#' or a disk-backed data source. The fitting procedure accesses only the variables
#' required by the model and does not assume that all data must reside in memory.
#'
#' @details
#' The \code{task} argument controls how the model is fitted:
#' \itemize{
#'   \item \strong{reml}: Restricted maximum likelihood estimation.
#'   \item \strong{bayes}: Bayesian inference using priors specified in variance
#'         components.
#'   \item \strong{solve}: Fixed-parameter solution with known variance components.
#' }
#'
#' Model formulas specify which variables and random effects exist, while variance
#' components specify how covariance is modeled. This separation allows the same
#' model specification to be reused across different estimation tasks and data
#' representations.
#'
#' @param formulas Named list of model formulas, one per trait.
#' @param data Data source containing observations, covariates, and grouping
#'   variables. May be an in-memory data frame or a disk-backed data source
#'   supporting variable-wise and row-wise access.
#' @param vcs Named list of variance component objects created with \code{vc()}.
#' @param task Character. Estimation task to use. One of \code{"reml"},
#'   \code{"bayes"}, or \code{"solve"}.
#' @param control Optional list of task-specific control parameters.
#' @param ... Additional arguments passed to task-specific fitting routines.
#'
#' @return
#' An object of class \code{"gfit"} containing the fitted model, parameter
#' estimates, and task-specific output.
#'
#' @examples
#' ## Single-trait and multi-trait animal models using mouse data
#'
#' ## Pedigree kernel
#' PED <- makePEDlist(fnPED = mouse_pedigree)
#'
#' ## Single-trait animal model (body weight)
#' formulas_single <- list(
#'   BW = BW ~ sex + reps + (1 | dam) + (1 | id)
#' )
#'
#' vcs_single <- list(
#'
#'   ## Dam environmental effect
#'   dam_env = vc(
#'     index  = "dam",
#'     traits = "BW"
#'   ),
#'
#'   ## Additive genetic animal effect
#'   animal_genetic = vc(
#'     index  = "id",
#'     traits = "BW",
#'     kernel = PED
#'   ),
#'
#'   ## Residual variance
#'   residual = vc(
#'     index  = "Residual",
#'     traits = "BW"
#'   )
#' )
#'
#' fit_single <- gfit(
#'   formulas = formulas_single,
#'   data     = mouse,
#'   vcs      = vcs_single,
#'   task     = "reml"
#' )
#'
#'
#' ## Multi-trait animal model (growth and body weight)
#' formulas_multi <- list(
#'   Gl = Gl ~ sex + reps + (1 | dam) + (1 | id),
#'   BW = BW ~ sex + reps + (1 | dam) + (1 | id)
#' )
#'
#' vcs_multi <- list(
#'
#'   ## Dam environmental effect (correlated across traits)
#'   dam_env = vc(
#'     index  = "dam",
#'     traits = c("Gl", "BW")
#'   ),
#'
#'   ## Additive genetic animal effect (animal model)
#'   animal_genetic = vc(
#'     index  = "id",
#'     traits = c("Gl", "BW"),
#'     kernel = PED
#'   ),
#'
#'   ## Residual covariance between traits
#'   residual = vc(
#'     index  = "Residual",
#'     traits = c("Gl", "BW")
#'   )
#' )
#'
#' fit_multi <- gfit(
#'   formulas = formulas_multi,
#'   data     = mouse,
#'   vcs      = vcs_multi,
#'   task     = "reml"
#' )
#'
#' @export
gfit <- function(formulas,
                 data,
                 vcs,
                 task = c("reml", "bayes", "solve"),
                 control = list(),
                 ...) {

  ## ---- task ---------------------------------------------------------------
  task <- match.arg(task)

  ## ---- basic validation ---------------------------------------------------
  if (!is.list(formulas) || length(formulas) < 1)
    stop("'formulas' must be a non-empty list of model formulas")

  if (!is.data.frame(data))
    stop("'data' must be a data frame")

  if (!is.list(vcs) || !all(vapply(vcs, inherits, logical(1), "vc")))
    stop("'vcs' must be a list of 'vc' objects")

  ## ---- build internal model representation --------------------------------
  model <- list(
    formulas = formulas,
    data     = data,
    vcs      = vcs
  )

  ## ---- dispatch by task ---------------------------------------------------
  fit <- switch(
    task,
    reml  = fit_reml(model, control = control, ...),
    bayes = fit_bayes(model, control = control, ...),
    solve = fit_solve(model, control = control, ...)
  )

  ## ---- attach metadata ----------------------------------------------------
  fit$task <- task
  class(fit) <- c("gfit", class(fit))

  fit
}

# ------------------------------------------------------------------
# qgtools: Formula parsing utilities
# ------------------------------------------------------------------

#' Parse a list of per-trait formulas into a normalized qgtools model spec
#'
#' @param formulas Named list of formulas, e.g.
#'   list(bmi = bmi ~ age + (time | id))
#' @return Named list with entries: fixed (chr), random (list)
#' @export
qg_parse_formulas <- function(formulas) {
  if (!is.list(formulas) || is.null(names(formulas)) || any(names(formulas) == "")) {
    stop("`formulas` must be a *named* list of formulas (one name per trait).")
  }

  out <- vector("list", length(formulas))
  names(out) <- names(formulas)

  for (trait in names(formulas)) {
    fm <- formulas[[trait]]
    if (!inherits(fm, "formula")) {
      stop("All entries in `formulas` must be formulas.")
    }

    out[[trait]] <- list(
      fixed  = .qg_extract_fixed_terms(fm),
      random = .qg_extract_random_terms(fm)
    )
  }

  out
}

# ------------------------------------------------------------------
# Random-effect parsing
# ------------------------------------------------------------------

.qg_extract_random_terms <- function(fm) {
  rhs <- fm[[3L]]
  bars <- .qg_find_bar_calls(rhs)
  if (length(bars) == 0L) return(list())
  lapply(bars, .qg_parse_bar_call)
}

.qg_find_bar_calls <- function(expr) {
  out <- list()

  walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1L]], as.name("|"))) {
        out[[length(out) + 1L]] <<- x
      } else {
        for (i in seq_along(x)[-1L]) walk(x[[i]])
      }
    }
  }

  walk(expr)
  out
}

.qg_parse_bar_call <- function(bar_call) {
  lhs <- bar_call[[2L]]
  grp <- bar_call[[3L]]

  eff <- .qg_parse_random_lhs(lhs)

  list(
    index     = .qg_deparse1(grp),
    intercept = eff$intercept,
    slopes    = eff$slopes,
    raw       = paste0("(", .qg_deparse1(lhs), " | ", .qg_deparse1(grp), ")")
  )
}

.qg_parse_random_lhs <- function(lhs) {
  intercept <- TRUE
  slopes <- character(0)

  is_const <- function(x, v) is.numeric(x) && length(x) == 1L && x == v

  # Atomic (symbol or constant)
  if (!is.call(lhs)) {
    if (is_const(lhs, 1)) return(list(intercept = TRUE, slopes = character(0)))
    if (is_const(lhs, 0) || is_const(lhs, -1)) {
      return(list(intercept = FALSE, slopes = character(0)))
    }
    return(list(intercept = TRUE, slopes = .qg_deparse1(lhs)))
  }

  # Handle "+"
  if (identical(lhs[[1L]], as.name("+"))) {
    parts <- .qg_flatten_plus(lhs)

    has0  <- any(vapply(parts, is_const, logical(1), 0))
    has1  <- any(vapply(parts, is_const, logical(1), 1))
    hasm1 <- any(vapply(parts, is_const, logical(1), -1))

    intercept <- !(has0 || hasm1)
    if (has1) intercept <- TRUE

    slopes <- vapply(parts, function(p) {
      if (is.numeric(p)) return(NA_character_)
      .qg_deparse1(p)
    }, character(1))

    slopes <- unique(slopes[!is.na(slopes)])
    return(list(intercept = intercept, slopes = slopes))
  }

  # Any other call (poly(), log(), etc.)
  list(intercept = TRUE, slopes = .qg_deparse1(lhs))
}

.qg_flatten_plus <- function(expr) {
  out <- list()
  walk <- function(x) {
    if (is.call(x) && identical(x[[1L]], as.name("+"))) {
      walk(x[[2L]])
      walk(x[[3L]])
    } else {
      out[[length(out) + 1L]] <<- x
    }
  }
  walk(expr)
  out
}

# ------------------------------------------------------------------
# Fixed effects
# ------------------------------------------------------------------

.qg_extract_fixed_terms <- function(fm) {
  rhs <- fm[[3L]]
  rhs_clean <- .qg_replace_bars_with_zero(rhs)
  fm_clean <- as.formula(call("~", fm[[2L]], rhs_clean), env = environment(fm))
  attr(terms(fm_clean), "term.labels")
}

.qg_replace_bars_with_zero <- function(expr) {
  if (is.call(expr)) {
    if (identical(expr[[1L]], as.name("|"))) return(0)
    for (i in 2:length(expr)) expr[[i]] <- .qg_replace_bars_with_zero(expr[[i]])
  }
  expr
}

.qg_deparse1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

qg_match_random_effects <- function(parsed_formulas, varcomp) {

  # Collect random-effect indices from formulas
  formula_indices <- unique(unlist(lapply(parsed_formulas, function(tr) {
    vapply(tr$random, `[[`, character(1), "index")
  })))

  # Collect vc indices
  vc_indices <- vapply(varcomp, function(vc) vc$index, character(1))

  # Missing vc definitions
  missing_vc <- setdiff(formula_indices, vc_indices)
  if (length(missing_vc) > 0) {
    stop(
      "Random effects present in formulas but missing vc() definitions: ",
      paste(missing_vc, collapse = ", ")
    )
  }

  # Unused vc definitions (warning only)
  unused_vc <- setdiff(vc_indices, formula_indices)
  if (length(unused_vc) > 0) {
    warning(
      "vc() definitions not referenced in formulas: ",
      paste(unused_vc, collapse = ", ")
    )
  }

  invisible(TRUE)
}
