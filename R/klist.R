
# Klist class ---------------------------------------------------------------

#' Construct a Klist object
#'
#' A Klist is a named collection of kernel/data backends (e.g. GRMlist, Glist, PEDlist),
#' with consistent metadata conventions for disk-backed storage.
#'
#' @param ... Named kernel objects. Names are kernel IDs referenced from vcov$<name>$kernel.
#' @param .list Alternatively, provide a named list.
#' @param strict If TRUE (default), fail on validation issues. If FALSE, store issues as an attribute.
#' @return An object of class "Klist".
#' @export
qg_Klist <- function(..., .list = NULL, strict = TRUE) {
  if (!is.null(.list)) {
    if (!is.list(.list)) stop("`.list` must be a list.")
    if (length(list(...)) > 0) stop("Provide either `...` or `.list`, not both.")
    x <- .list
  } else {
    x <- list(...)
  }

  if (is.null(names(x)) || any(names(x) == "")) {
    stop("Klist entries must be named (kernel IDs).")
  }

  issues <- validate_Klist(x, strict = FALSE)

  class(x) <- c("Klist", "list")

  attr(x, "issues") <- issues
  if (strict && length(issues) > 0) {
    stop(paste0("Invalid Klist:\n- ", paste(issues, collapse = "\n- ")))
  }
  x
}

#' Validate a Klist
#' @param Klist A Klist or named list of kernels
#' @param strict If TRUE, error; else return character vector of issues.
#' @keywords internal
validate_Klist <- function(Klist, strict = TRUE) {
  if (inherits(Klist, "Klist")) x <- unclass(Klist) else x <- Klist

  issues <- character(0)

  if (!is.list(x)) issues <- c(issues, "Klist must be a list.")
  if (is.null(names(x)) || any(names(x) == "")) issues <- c(issues, "All kernels in Klist must be named.")

  if (length(issues) > 0) {
    if (strict) stop(paste0("Invalid Klist:\n- ", paste(issues, collapse = "\n- ")))
    return(issues)
  }

  # Validate each kernel object
  for (nm in names(x)) {
    obj <- x[[nm]]
    k_issues <- validate_kernel(obj, kernel_id = nm)
    issues <- c(issues, k_issues)
  }

  if (strict && length(issues) > 0) stop(paste0("Invalid Klist:\n- ", paste(issues, collapse = "\n- ")))
  issues
}

#' Validate a single kernel object
#' @keywords internal
validate_kernel <- function(obj, kernel_id = "<kernel>") {
  issues <- character(0)

  if (!is.list(obj)) {
    return(paste0(kernel_id, ": kernel must be a list-like object."))
  }

  cls <- class(obj)

  # Allow either explicit class OR duck-typed by required fields
  if ("GRMlist" %in% cls) {
    issues <- c(issues, validate_GRMlist(obj, kernel_id))
  } else if ("Glist" %in% cls) {
    issues <- c(issues, validate_Glist(obj, kernel_id))
  } else if ("PEDlist" %in% cls) {
    issues <- c(issues, validate_PEDlist(obj, kernel_id))
  } else if ("Mlist" %in% cls) {
    issues <- c(issues, validate_Mlist(obj, kernel_id))
  } else if ("LDlist" %in% cls) {
    issues <- c(issues, validate_LDlist(obj, kernel_id))
  } else {
    # Heuristic detection
    if (!is.null(obj$fnG) || !is.null(obj$idsG)) {
      issues <- c(issues, validate_GRMlist(obj, kernel_id))
    } else if (!is.null(obj$bedfiles) || !is.null(obj$fnBED)) {
      issues <- c(issues, validate_Glist(obj, kernel_id))
    } else if (!is.null(obj$pedfiles) || !is.null(obj$fnPED) || !is.null(obj$pedigree)) {
      issues <- c(issues, validate_PEDlist(obj, kernel_id))
    } else if (!is.null(obj$ldfiles) || !is.null(obj$fnLD)) {
      issues <- c(issues, validate_LDlist(obj, kernel_id))
    } else {
      issues <- c(issues, paste0(kernel_id, ": unknown kernel type (add a class like 'GRMlist'/'Glist' or required fields)."))
    }
  }

  issues
}
