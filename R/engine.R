
.qg_engines <- new.env(parent = emptyenv())

#' Register a computational engine
#'
#' @param name Engine name.
#' @param fun A function with signature \code{function(compiled, ...)} returning a fit object.
#' @export
qg_register_engine <- function(name, fun) {
  stopifnot(is.character(name), length(name) == 1)
  stopifnot(is.function(fun))
  assign(name, fun, envir = .qg_engines)
  invisible(TRUE)
}

#' List registered engines
#' @export
qg_list_engines <- function() {
  sort(ls(envir = .qg_engines))
}

# Default "mock" engine for testing / prototyping
qg_register_engine("mock", function(compiled, ...) {
  structure(list(compiled = compiled, engine = "mock"), class = "qg_fit")
})
