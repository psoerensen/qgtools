
.require_jsonlite <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for JSON export")
  }
}



#' Convert a vc object to a list suitable for JSON export
#'
#' @param x An object of class "vc"
#' @param name Optional name of the variance component
#'
#' @return A named list
#' @export
as_list.vc <- function(x) {

  if (!inherits(x, "vc"))
    stop("Object must be of class 'vc'")

  if (is.null(x$kernel)) {
    stop("vc must include a 'kernel' (use iid_kernel('default') or supply a kernel object).")
  }

  out <- list(
    index     = x$index,
    traits    = x$traits,
    kernel    = kernel_name(x$kernel),
    structure = x$structure
  )

  if (!is.null(x$start))
    out$start <- unclass(x$start)

  if (!is.null(x$prior))
    out$prior <- as_list(x$prior)
  else
    out$prior <- NULL

  out
}



#' Export a vc object to JSON
#'
#' @param x An object of class "vc"
#' @param pretty Logical. Pretty-print JSON.
#' @param auto_unbox Logical. Unbox scalars.
#'
#' @return A JSON string
#' @export
as_json.vc <- function(x, pretty = TRUE, auto_unbox = TRUE) {

  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required")

  jsonlite::toJSON(
    as_list.vc(x),
    pretty = pretty,
    auto_unbox = auto_unbox,
    null = "null"
  )
}

#' Export a list of variance components to JSON
#'
#' @param vcs Named list of vc objects
#' @param pretty Logical. Pretty-print JSON.
#'
#' @return A JSON string
#' @export
as_json.vcs <- function(vcs, pretty = TRUE) {

  if (!is.list(vcs) || any(!vapply(vcs, inherits, logical(1), "vc")))
    stop("'vcs' must be a named list of vc objects")

  spec <- list(
    type = "VarCompSpec",
    components = lapply(vcs, as_list.vc)
  )

  spec$components <- setNames(
    spec$components,
    names(vcs)
  )


  jsonlite::toJSON(
    spec,
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}



#' Extract kernel name for JSON export
#'
#' @param kernel A kernel object
#'
#' @return Character scalar identifying the kernel type
#' @export
kernel_name <- function(kernel) {

  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for kernel hashing")
  }


  if (is.null(kernel))
    stop("Kernel is NULL")

  if (!is.list(kernel) || is.null(kernel$type))
    stop("Invalid kernel object: missing 'type' field")

  if (!is.null(kernel$meta$id))
    return(kernel$meta$id)

  paste0(kernel$type, "_",
         digest::digest(
           unclass(kernel$data),
           algo = "xxhash64",
           serialize = TRUE
         )
  )
}


#' Convert a kernel object to a list for JSON export
#'
#' @param x Kernel object
#' @return A named list
#' @export
as_list.kernel <- function(x) {

  if (!inherits(x, "kernel"))
    stop("Object must inherit from class 'kernel'")

  out <- list(
    type = x$type
  )

  if (!is.null(x$data)) {

    if (inherits(x$data, "PEDlist")) {
      out$data <- pedlist_to_list(x$data)
    } else {
      ## strip class just in case
      out$data <- unclass(x$data)
    }
  }

  if (!is.null(x$meta)) {
    out$meta <- unclass(x$meta)
  }

  out <- unclass(out)
  stopifnot(
    is.character(out$type),
    is.null(out$data) || is.list(out$data)
  )

  out
}

#' Export kernel specifications to JSON
#'
#' @param vcs Named list of vc objects
#' @param pretty Logical
#' @return JSON string
#' @export
as_json.kernels <- function(vcs, pretty = TRUE) {

  .require_jsonlite()

  kernels <- lapply(vcs, function(vc) vc$kernel)

  kernels <- kernels[!vapply(kernels, is.null, logical(1))]

  names(kernels) <- vapply(kernels, kernel_name, character(1))

  kernels <- kernels[!duplicated(names(kernels))]

  spec <- list(
    type = "KernelSpec",
    kernels = lapply(kernels, as_list.kernel)
  )

  stopifnot(!any(vapply(
    spec$kernels,
    function(k) inherits(k, "PEDlist"),
    logical(1)
  )))


  jsonlite::toJSON(
    spec,
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}


#' Convert a Datalist object to a list for JSON export
#'
#' @param x Datalist object
#' @return A named list
#' @export
as_list.Datalist <- function(x) {

  if (!inherits(x, "Datalist"))
    stop("Object must be of class 'Datalist'")

  list(
    type     = "DataSpec",
    format   = x$format,
    encoding = x$encoding,
    source   = x$source,
    id       = x$id,
    header   = x$header,
    colnames = x$colnames,
    roles    = x$roles,
    missing  = x$missing,
    options  = x$options
  )
}

#' Export a Datalist object to JSON
#'
#' @param data Datalist object
#' @param pretty Logical
#' @return JSON string
#' @export
as_json.data <- function(data, pretty = TRUE) {

  .require_jsonlite()

  jsonlite::toJSON(
    as_list.Datalist(data),
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}

#' Convert model specification to list for JSON export
#'
#' @param formulas Named list of formulas
#' @param task Character ("reml", "bayes", "solve")
#' @param options Optional list
#' @return Named list
#' @export
as_list.model <- function(formulas, task, options = list()) {

  stopifnot(
    is.list(formulas),
    all(vapply(formulas, inherits, logical(1), "formula"))
  )

  parsed <- qg_parse_formulas(formulas)

  roles <- qg_extract_variable_roles(parsed)
  #roles <- qg_normalize_roles(roles)

  list(
    type = "ModelSpec",
    task = toupper(task),
    traits = names(formulas),

    summary = list(
      fixed  = lapply(parsed, function(x) x$fixed),
      random = lapply(parsed, function(x) x$random)
    ),

    variables = qg_roles_to_spec(roles),
    options   = options
  )
}

#' Export model specification to JSON
#'
#' @param formulas Named list of formulas
#' @param task Character
#' @param options Optional list
#' @param pretty Logical
#' @return JSON string
#' @export
as_json.model <- function(formulas, task, options = list(), pretty = TRUE) {

  .require_jsonlite()

  jsonlite::toJSON(
    as_list.model(formulas, task, options),
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}

#' Export full model specification bundle to JSON files
#'
#' @param data Datalist
#' @param formulas Model formulas
#' @param vcs Variance components
#' @param task Character
#' @param path Output directory
#' @export
export_model_bundle <- function(data, formulas, vcs, task, path = ".") {

  if (!dir.exists(path))
    dir.create(path, recursive = TRUE)

  writeLines(as_json.data(data),      file.path(path, "data.json"))
  writeLines(as_json.kernels(vcs),     file.path(path, "kernels.json"))
  writeLines(as_json.vcs(vcs),         file.path(path, "vcs.json"))
  writeLines(as_json.model(formulas, task), file.path(path, "model.json"))

  invisible(path)
}


pedlist_to_list <- function(x) {

  if (!inherits(x, "PEDlist"))
    stop("Expected a PEDlist object")

  list(
    type        = "PED",
    kernel_type = x$kernel_type,
    file        = x$fnPED,
    format      = x$format,
    encoding    = x$encoding,
    method      = x$method,
    method_code = x$method_code,
    columns     = as.list(x$columns)
  )
}


qg_validate_model_data <- function(data, model) {

  vars_needed <- names(model$variables)
  vars_have   <- data$colnames

  missing <- setdiff(vars_needed, vars_have)

  if (length(missing) > 0) {
    stop(
      "Variables required by the model are missing from the data source: ",
      paste(missing, collapse = ", ")
    )
  }

  invisible(TRUE)
}


#' Convert kernel specifications to list form
#'
#' @param vcs Named list of vc objects
#' @return KernelSpec list
#' @export
as_list.kernels <- function(vcs) {

  kernels <- lapply(vcs, function(vc) vc$kernel)
  kernels <- kernels[!vapply(kernels, is.null, logical(1))]

  kernel_ids <- vapply(kernels, kernel_name, character(1))

  keep <- !duplicated(kernel_ids)
  kernels <- kernels[keep]
  kernel_ids <- kernel_ids[keep]

  kernel_list <- lapply(kernels, as_list.kernel)
  names(kernel_list) <- kernel_ids  # ← CRITICAL LINE

  list(
    type    = "KernelSpec",
    kernels = kernel_list
  )
}

# as_list.kernels <- function(vcs) {
#
#   if (!is.list(vcs) || any(!vapply(vcs, inherits, logical(1), "vc"))) {
#     stop("'vcs' must be a named list of vc objects")
#   }
#
#   kernels <- lapply(vcs, function(vc) vc$kernel)
#   kernels <- kernels[!vapply(kernels, is.null, logical(1))]
#
#   kernel_ids <- vapply(kernels, kernel_name, character(1))
#   names(kernels) <- kernel_ids
#
#   # deduplicate identical kernels
#   kernels <- kernels[!duplicated(names(kernels))]
#
#   list(
#     type = "KernelSpec",
#     kernels = lapply(kernels, as_list.kernel)
#   )
# }

#' Convert variance component specifications to a list
#'
#' @param vcs Named list of vc objects
#' @return VarCompSpec list
#' @export
as_list.vcs <- function(vcs) {

  if (!is.list(vcs) || any(!vapply(vcs, inherits, logical(1), "vc"))) {
    stop("'vcs' must be a named list of vc objects")
  }

  list(
    type = "VarCompSpec",
    components = setNames(
      lapply(vcs, as_list.vc),
      names(vcs)
    )
  )
}

