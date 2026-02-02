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
