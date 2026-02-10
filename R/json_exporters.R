# ==============================================================================
# JSON EXPORT UTILITIES
# ==============================================================================

.require_jsonlite <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Package 'jsonlite' is required for JSON export")
}

# ==============================================================================
# KERNEL EXPORT
# ==============================================================================

# Extract a stable kernel identifier (must be explicit)
kernel_name <- function(kernel) {
  stopifnot(inherits(kernel, "kernel"))

  if (!is.null(kernel$meta$id))
    return(kernel$meta$id)

  stop("Kernel has no explicit ID (kernel$meta$id must be set).")
}

# Convert a kernel object to JSON-ready list
as_list.kernel <- function(x) {
  stopifnot(inherits(x, "kernel"))

  out <- list(
    id   = kernel_name(x),
    type = x$type
  )

  # Kernel data must be disk-backed specs only
  if (!is.null(x$data)) {
    stopifnot(inherits(x$data, c("PEDlist", "GRMlist", "LDlist")))
    out$data <- list(ref = out$id)
  }

  if (!is.null(x$meta))
    out$meta <- unclass(x$meta)

  out
}

# Collect unique kernels from variance components
as_list.kernels <- function(vcs) {
  stopifnot(is.list(vcs))

  kernels <- Filter(Negate(is.null), lapply(vcs, `[[`, "kernel"))
  ids     <- vapply(kernels, kernel_name, character(1))

  keep    <- !duplicated(ids)
  kernels <- kernels[keep]
  ids     <- ids[keep]

  list(
    type    = "KernelSpec",
    kernels = setNames(lapply(kernels, as_list.kernel), ids)
  )
}

as_json.kernels <- function(vcs, pretty = TRUE) {
  .require_jsonlite()
  jsonlite::toJSON(
    as_list.kernels(vcs),
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}

# ==============================================================================
# VARIANCE COMPONENT EXPORT
# ==============================================================================

# Convert a single vc() object
as_list.vc <- function(x) {
  stopifnot(inherits(x, "vc"))

  # Identify covariance mechanism (mirrors vc() logic)
  covariance <- if (!is.null(x$kernel)) {
    list(type = "kernel", kernel = kernel_name(x$kernel))
  } else if (!is.null(x$featureMatrix)) {
    list(
      type          = "feature",
      featureMatrix = class(x$featureMatrix)[1],
      featureSets   = x$featureSets
    )
  } else if (x$variable == "residual") {
    list(type = "residual")
  } else {
    list(type = "iid")
  }

  out <- list(
    variable   = x$variable,
    traits     = x$traits,
    covariance = covariance,
    structure  = x$structure
  )

  if (!is.null(x$start))
    out$start <- unclass(x$start)

  out
}

# Convert list(vc(), vc(), ...)
as_list.vcs <- function(vcs) {
  stopifnot(
    is.list(vcs),
    all(vapply(vcs, inherits, logical(1), "vc"))
  )

  list(
    type = "VarCompSpec",
    components = setNames(lapply(vcs, as_list.vc), names(vcs))
  )
}

as_json.vcs <- function(vcs, pretty = TRUE) {
  .require_jsonlite()
  jsonlite::toJSON(
    as_list.vcs(vcs),
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}

# ==============================================================================
# PRIOR EXPORT (MIRRORS VARIANCE COMPONENTS)
# ==============================================================================

# Convert a single prior() object
as_list.prior <- function(x) {
  stopifnot(inherits(x, "prior"))

  covariance <- if (!is.null(x$kernel)) {
    list(type = "kernel", kernel = kernel_name(x$kernel))
  } else if (!is.null(x$featureMatrix)) {
    list(
      type          = "feature",
      featureMatrix = class(x$featureMatrix)[1],
      featureSets   = x$featureSets
    )
  } else if (x$variable == "residual") {
    list(type = "residual")
  } else {
    list(type = "iid")
  }

  out <- list(
    variable     = x$variable,
    traits       = x$traits,
    covariance   = covariance,
    structure    = x$structure,
    distribution = as_list(x$distribution)
  )

  if (!is.null(x$start))
    out$start <- unclass(x$start)

  out
}

# Convert list(prior(), prior(), ...)
as_list.priors <- function(priors) {
  stopifnot(
    is.list(priors),
    all(vapply(priors, inherits, logical(1), "prior"))
  )

  list(
    type = "PriorSpec",
    components = setNames(lapply(priors, as_list.prior), names(priors))
  )
}

as_json.priors <- function(priors, pretty = TRUE) {
  .require_jsonlite()
  jsonlite::toJSON(
    as_list.priors(priors),
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}

# ==============================================================================
# DATA & FEATURE EXPORT
# ==============================================================================

# Data source (phenotypes, covariates, etc.)
as_list.DataSource <- function(x) {
  stopifnot(inherits(x, "DataSource"))

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

as_json.data <- function(data, pretty = TRUE) {
  .require_jsonlite()
  jsonlite::toJSON(
    as_list.DataSource(data),
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}

# Feature sources (genomic, transcriptomic, etc.)
as_list.FeatureSource <- function(x) {
  stopifnot(inherits(x, "FeatureSource"))

  list(
    type     = "FeatureSpec",
    features = lapply(x$features, as_list)
  )
}

as_json.features <- function(features, pretty = TRUE) {
  .require_jsonlite()
  jsonlite::toJSON(
    as_list.FeatureSource(features),
    pretty = pretty,
    auto_unbox = TRUE,
    null = "null"
  )
}

# ==============================================================================
# MODEL BUNDLE EXPORT
# ==============================================================================

# Export a complete model specification bundle
export_model_bundle <- function(data,
                                features,
                                formulas,
                                vcs,
                                task,
                                path = ".") {

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  writeLines(as_json.data(data),           file.path(path, "data.json"))
  writeLines(as_json.features(features),   file.path(path, "features.json"))
  writeLines(as_json.kernels(vcs),          file.path(path, "kernels.json"))
  writeLines(as_json.vcs(vcs),              file.path(path, "vcs.json"))
  writeLines(as_json.model(formulas, task), file.path(path, "model.json"))

  invisible(path)
}


# ==============================================================================
# VALIDATE JSON MODEL SPECIFICATION BUNDLE
# ==============================================================================

#' Validate consistency of a model specification bundle
#'
#' Ensures that DataSpec, ModelSpec, KernelSpec, and either VarCompSpec
#' (REML / solver) or PriorSpec (Bayesian) are mutually consistent.
#'
#' @param data_spec DataSpec list
#' @param model_spec ModelSpec list
#' @param kernel_spec KernelSpec list
#' @param varcomp_spec Optional VarCompSpec list
#' @param prior_spec Optional PriorSpec list
#'
#' @export
validate_bundle_json <- function(data_spec,
                                 model_spec,
                                 kernel_spec,
                                 varcomp_spec = NULL,
                                 prior_spec   = NULL) {

  # --------------------------------------------------------------------------
  # 1. BASIC STRUCTURE CHECKS
  # --------------------------------------------------------------------------

  stopifnot(
    is.list(data_spec),  identical(data_spec$type,  "DataSpec"),
    is.list(model_spec), identical(model_spec$type, "ModelSpec"),
    is.list(kernel_spec),identical(kernel_spec$type,"KernelSpec")
  )

  # Exactly one of VarCompSpec or PriorSpec must be supplied
  if (xor(is.null(varcomp_spec), is.null(prior_spec)) == FALSE) {
    stop(
      "Specify exactly one of:\n",
      "  - VarCompSpec (REML / solver)\n",
      "  - PriorSpec   (Bayesian)"
    )
  }

  components <- if (!is.null(varcomp_spec)) {
    stopifnot(identical(varcomp_spec$type, "VarCompSpec"))
    varcomp_spec$components
  } else {
    stopifnot(identical(prior_spec$type, "PriorSpec"))
    prior_spec$components
  }

  # --------------------------------------------------------------------------
  # 2. KERNEL REFERENCE VALIDATION
  # --------------------------------------------------------------------------

  defined_kernels <- names(kernel_spec$kernels)

  referenced_kernels <- unique(unlist(lapply(
    components,
    function(x) {
      cov <- x$covariance
      if (is.null(cov) || cov$type != "kernel") return(NULL)
      cov$kernel
    }
  )))

  missing_kernels <- setdiff(referenced_kernels, defined_kernels)

  if (length(missing_kernels) > 0) {
    stop(
      "Components reference undefined kernels: ",
      paste(missing_kernels, collapse = ", ")
    )
  }

  # --------------------------------------------------------------------------
  # 3. RANDOM EFFECTS VS MODEL FORMULAS
  # --------------------------------------------------------------------------

  # Random-effect variables appearing in formulas
  model_vars <- unique(unlist(lapply(
    model_spec$summary$random,
    function(tr) {
      if (length(tr) == 0) return(character(0))
      vapply(tr, `[[`, character(1), "variable")
    }
  )))

  # Variables defined in VarComp / Prior spec
  spec_vars <- vapply(
    components,
    function(x) x$variable,
    character(1)
  )

  # Variables allowed without explicit formula terms
  allowed_extra <- c("residual", "marker")

  # Missing specs
  missing_specs <- setdiff(model_vars, spec_vars)
  if (length(missing_specs) > 0) {
    stop(
      "Random effects in formulas but missing specification: ",
      paste(missing_specs, collapse = ", ")
    )
  }

  # Unused specs
  unused_specs <- setdiff(spec_vars, c(model_vars, allowed_extra))
  if (length(unused_specs) > 0) {
    warning(
      "Specifications not referenced in formulas: ",
      paste(unused_specs, collapse = ", ")
    )
  }

  # --------------------------------------------------------------------------
  # 4. DATA VARIABLE AVAILABILITY
  # --------------------------------------------------------------------------

  vars_needed <- names(model_spec$variables)
  vars_have   <- if (is.list(data_spec$colnames)) {
    unique(unlist(data_spec$colnames))
  } else {
    data_spec$colnames
  }

  missing_vars <- setdiff(vars_needed, vars_have)

  if (length(missing_vars) > 0) {
    stop(
      "Variables required by the model are missing from the data source: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # --------------------------------------------------------------------------
  # 5. BAYESIAN-ONLY CONSISTENCY CHECKS
  # --------------------------------------------------------------------------

  if (!is.null(prior_spec)) {

    for (nm in names(prior_spec$components)) {
      p <- prior_spec$components[[nm]]

      cov <- p$covariance

      if (p$variable != "residual" &&
          cov$type == "iid") {
        stop(
          "Prior '", nm,
          "' has implicit IID covariance.\n",
          "Bayesian priors must specify kernel- or feature-backed covariance."
        )
      }
    }
  }

  invisible(TRUE)
}

#' .require_jsonlite <- function() {
#'   if (!requireNamespace("jsonlite", quietly = TRUE)) {
#'     stop("Package 'jsonlite' is required for JSON export")
#'   }
#' }
#'
#'
#'
#' #' Convert a vc object to a list suitable for JSON export
#' #'
#' #' @param x An object of class "vc"
#' #' @param name Optional name of the variance component
#' #'
#' #' @return A named list
#' #' @export
#' as_list.vc <- function(x) {
#'
#'   stopifnot(inherits(x, "vc"))
#'
#'   ## ---- identify covariance mechanism -----------------------------------
#'   cov <- if (!is.null(x$kernel)) {
#'     list(
#'       type   = "kernel",
#'       kernel = kernel_name(x$kernel)
#'     )
#'   } else if (!is.null(x$featureMatrix)) {
#'     list(
#'       type          = "feature",
#'       featureMatrix = class(x$featureMatrix)[1],
#'       featureSets   = x$featureSets
#'     )
#'   } else if (x$variable == "residual") {
#'     list(type = "residual")
#'   } else {
#'     list(type = "iid")
#'   }
#'
#'   out <- list(
#'     variable   = x$variable,
#'     traits     = x$traits,
#'     covariance = cov,
#'     structure  = x$structure
#'   )
#'
#'   if (!is.null(x$start))
#'     out$start <- unclass(x$start)
#'
#'   out
#' }
#'
#'
#'
#'
#' #' Export a vc object to JSON
#' #'
#' #' @param x An object of class "vc"
#' #' @param pretty Logical. Pretty-print JSON.
#' #' @param auto_unbox Logical. Unbox scalars.
#' #'
#' #' @return A JSON string
#' #' @export
#' as_json.vc <- function(x, pretty = TRUE, auto_unbox = TRUE) {
#'
#'   if (!requireNamespace("jsonlite", quietly = TRUE))
#'     stop("Package 'jsonlite' is required")
#'
#'   jsonlite::toJSON(
#'     as_list.vc(x),
#'     pretty = pretty,
#'     auto_unbox = auto_unbox,
#'     null = "null"
#'   )
#' }
#'
#' #' Export a list of variance components to JSON
#' #'
#' #' @param vcs Named list of vc objects
#' #' @param pretty Logical. Pretty-print JSON.
#' #'
#' #' @return A JSON string
#' #' @export
#' as_json.vcs <- function(vcs, pretty = TRUE) {
#'
#'   if (!is.list(vcs) || any(!vapply(vcs, inherits, logical(1), "vc")))
#'     stop("'vcs' must be a named list of vc objects")
#'
#'   spec <- list(
#'     type = "VarCompSpec",
#'     components = lapply(vcs, as_list.vc)
#'   )
#'
#'   spec$components <- setNames(
#'     spec$components,
#'     names(vcs)
#'   )
#'
#'
#'   jsonlite::toJSON(
#'     spec,
#'     pretty = pretty,
#'     auto_unbox = TRUE,
#'     null = "null"
#'   )
#' }
#'
#'
#' #' Extract kernel name for JSON export
#' #'
#' #' @param kernel A kernel object
#' #'
#' #' @return Character scalar identifying the kernel type
#' #' @export
#' kernel_name <- function(kernel) {
#'
#'   if (!is.null(kernel$meta$id)) {
#'     return(kernel$meta$id)
#'   }
#'
#'   stop(
#'     "Kernel has no explicit ID. ",
#'     "Provide kernel$meta$id or construct kernels via helper functions."
#'   )
#' }
#'
#'
#' #' Convert a kernel object to a list for JSON export
#' #'
#' #' @param x Kernel object
#' #' @return A named list
#' #' @export
#' as_list.kernel <- function(x) {
#'
#'   if (!inherits(x, "kernel"))
#'     stop("Object must inherit from class 'kernel'")
#'
#'   id <- kernel_name(x)
#'
#'   out <- list(
#'     id   = id,
#'     type = x$type
#'   )
#'
#'   if (!is.null(x$data)) {
#'
#'     if (inherits(x$data, c("PEDlist", "GRMlist", "LDlist"))) {
#'       out$data <- list(ref = id)
#'     } else {
#'       stop(
#'         "Kernel data of class ",
#'         paste(class(x$data), collapse = "/"),
#'         " cannot be serialized. Use disk-backed kernel specifications only."
#'       )
#'     }
#'   }
#'
#'   if (!is.null(x$meta)) {
#'     out$meta <- unclass(x$meta)
#'   }
#'
#'   unclass(out)
#' }
#'
#' #' Export kernel specifications to JSON
#' #'
#' #' @param vcs Named list of vc objects
#' #' @param pretty Logical
#' #' @return JSON string
#' #' @export
#' as_json.kernels <- function(vcs, pretty = TRUE) {
#'
#'   .require_jsonlite()
#'
#'   kernels <- lapply(vcs, function(vc) vc$kernel)
#'
#'   kernels <- kernels[!vapply(kernels, is.null, logical(1))]
#'
#'   names(kernels) <- vapply(kernels, kernel_name, character(1))
#'
#'   kernels <- kernels[!duplicated(names(kernels))]
#'
#'   spec <- list(
#'     type = "KernelSpec",
#'     kernels = lapply(kernels, as_list.kernel)
#'   )
#'
#'   stopifnot(!any(vapply(
#'     spec$kernels,
#'     function(k) inherits(k, "PEDlist"),
#'     logical(1)
#'   )))
#'
#'
#'   jsonlite::toJSON(
#'     spec,
#'     pretty = pretty,
#'     auto_unbox = TRUE,
#'     null = "null"
#'   )
#' }
#'
#'
#' #' Convert a Datalist object to a list for JSON export
#' #'
#' #' @param x Datalist object
#' #' @return A named list
#' #' @export
#' as_list.DataSource <- function(x) {
#'
#'   stopifnot(inherits(x, "DataSource"))
#'
#'   list(
#'     type     = "DataSpec",
#'     format   = x$format,
#'     encoding = x$encoding,
#'     source   = x$source,
#'     id       = x$id,
#'     header   = x$header,
#'     colnames = x$colnames,
#'     roles    = x$roles,
#'     missing  = x$missing,
#'     options  = x$options
#'   )
#' }
#'
#' #' Export a Datalist object to JSON
#' #'
#' #' @param data Datalist object
#' #' @param pretty Logical
#' #' @return JSON string
#' #' @export
#' as_json.data <- function(data, pretty = TRUE) {
#'
#'   .require_jsonlite()
#'
#'   jsonlite::toJSON(
#'     as_list.Datalist(data),
#'     pretty = pretty,
#'     auto_unbox = TRUE,
#'     null = "null"
#'   )
#' }
#'
#' #' Convert model specification to list for JSON export
#' #'
#' #' @param formulas Named list of formulas
#' #' @param task Character ("reml", "bayes", "solve")
#' #' @param options Optional list
#' #' @return Named list
#' #' @export
#' as_list.model <- function(formulas, task, options = list()) {
#'
#'   stopifnot(
#'     is.list(formulas),
#'     all(vapply(formulas, inherits, logical(1), "formula"))
#'   )
#'
#'   parsed <- qg_parse_formulas(formulas)
#'
#'   roles <- qg_extract_variable_roles(parsed)
#'   #roles <- qg_normalize_roles(roles)
#'
#'   list(
#'     type = "ModelSpec",
#'     task = toupper(task),
#'     traits = names(formulas),
#'
#'     summary = list(
#'       fixed  = lapply(parsed, function(x) x$fixed),
#'       random = lapply(parsed, function(x) x$random)
#'     ),
#'
#'     variables = qg_roles_to_spec(roles),
#'     options   = options
#'   )
#' }
#'
#' #' Export model specification to JSON
#' #'
#' #' @param formulas Named list of formulas
#' #' @param task Character
#' #' @param options Optional list
#' #' @param pretty Logical
#' #' @return JSON string
#' #' @export
#' as_json.model <- function(formulas, task, options = list(), pretty = TRUE) {
#'
#'   .require_jsonlite()
#'
#'   jsonlite::toJSON(
#'     as_list.model(formulas, task, options),
#'     pretty = pretty,
#'     auto_unbox = TRUE,
#'     null = "null"
#'   )
#' }
#'
#'
#'
#' #' Convert FeatureSpec to list for JSON export
#' #'
#' #' @param x FeatureSpec object
#' #' @return Named list
#' #' @export
#' as_list.FeatureSource <- function(x) {
#'
#'   stopifnot(inherits(x, "FeatureSource"))
#'
#'   list(
#'     type = "FeatureSpec",
#'     features = lapply(x$features, as_list)
#'   )
#' }
#'
#' #' Export FeatureSpec to JSON
#' #'
#' #' @param features FeatureSpec object
#' #' @param pretty Logical
#' #' @return JSON string
#' #' @export
#' as_json.features <- function(features, pretty = TRUE) {
#'
#'   .require_jsonlite()
#'
#'   jsonlite::toJSON(
#'     as_list.FeatureSpec(features),
#'     pretty = pretty,
#'     auto_unbox = TRUE,
#'     null = "null"
#'   )
#' }
#'
#'
#' #' Export full model specification bundle to JSON files
#' #'
#' #' @param data Datalist
#' #' @param formulas Model formulas
#' #' @param vcs Variance components
#' #' @param task Character
#' #' @param path Output directory
#' #' @export
#' export_model_bundle <- function(data,
#'                                 features,
#'                                 formulas,
#'                                 vcs,
#'                                 task,
#'                                 path = ".") {
#'
#'   if (!dir.exists(path))
#'     dir.create(path, recursive = TRUE)
#'
#'   writeLines(as_json.data(data),          file.path(path, "data.json"))
#'   writeLines(as_json.features(features), file.path(path, "features.json"))
#'   writeLines(as_json.kernels(vcs),        file.path(path, "kernels.json"))
#'   writeLines(as_json.vcs(vcs),            file.path(path, "vcs.json"))
#'   writeLines(as_json.model(formulas, task), file.path(path, "model.json"))
#'
#'   invisible(path)
#' }
#'
#'
#' pedlist_to_list <- function(x) {
#'
#'   if (!inherits(x, "PEDlist"))
#'     stop("Expected a PEDlist object")
#'
#'   list(
#'     type        = "PED",
#'     kernel_type = x$kernel_type,
#'     file        = x$fnPED,
#'     format      = x$format,
#'     encoding    = x$encoding,
#'     method      = x$method,
#'     method_code = x$method_code,
#'     columns     = as.list(x$columns)
#'   )
#' }
#'
#'
#' qg_validate_model_data <- function(data, model) {
#'
#'   vars_needed <- names(model$variables)
#'   vars_have <- if (is.list(data$colnames)) unique(unlist(data$colnames)) else data$colnames
#'
#'
#'   missing <- setdiff(vars_needed, vars_have)
#'
#'   if (length(missing) > 0) {
#'     stop(
#'       "Variables required by the model are missing from the data source: ",
#'       paste(missing, collapse = ", ")
#'     )
#'   }
#'
#'   invisible(TRUE)
#' }
#'
#'
#' #' Convert kernel specifications to list form
#' #'
#' #' @param vcs Named list of vc objects
#' #' @return KernelSpec list
#' #' @export
#' as_list.kernels <- function(vcs) {
#'
#'   kernels <- lapply(vcs, function(vc) vc$kernel)
#'   kernels <- kernels[!vapply(kernels, is.null, logical(1))]
#'
#'   kernel_ids <- vapply(kernels, kernel_name, character(1))
#'
#'   keep <- !duplicated(kernel_ids)
#'   kernels <- kernels[keep]
#'   kernel_ids <- kernel_ids[keep]
#'
#'   kernel_list <- lapply(kernels, as_list.kernel)
#'   names(kernel_list) <- kernel_ids  # ← CRITICAL LINE
#'
#'   list(
#'     type    = "KernelSpec",
#'     kernels = kernel_list
#'   )
#' }
#'
#' # as_list.kernels <- function(vcs) {
#' #
#' #   if (!is.list(vcs) || any(!vapply(vcs, inherits, logical(1), "vc"))) {
#' #     stop("'vcs' must be a named list of vc objects")
#' #   }
#' #
#' #   kernels <- lapply(vcs, function(vc) vc$kernel)
#' #   kernels <- kernels[!vapply(kernels, is.null, logical(1))]
#' #
#' #   kernel_ids <- vapply(kernels, kernel_name, character(1))
#' #   names(kernels) <- kernel_ids
#' #
#' #   # deduplicate identical kernels
#' #   kernels <- kernels[!duplicated(names(kernels))]
#' #
#' #   list(
#' #     type = "KernelSpec",
#' #     kernels = lapply(kernels, as_list.kernel)
#' #   )
#' # }
#'
#' #' Convert variance component specifications to a list
#' #'
#' #' @param vcs Named list of vc objects
#' #' @return VarCompSpec list
#' #' @export
#' as_list.vcs <- function(vcs) {
#'
#'   stopifnot(
#'     is.list(vcs),
#'     all(vapply(vcs, inherits, logical(1), "vc"))
#'   )
#'
#'   list(
#'     type = "VarCompSpec",
#'     components = setNames(
#'       lapply(vcs, as_list.vc),
#'       names(vcs)
#'     )
#'   )
#' }
#'
#'
#'
#' #' Convert a prior object to a list suitable for JSON export
#' #'
#' #' @param x An object of class "prior"
#' #' @return Named list
#' #' @export
#' as_list.prior <- function(x) {
#'
#'   if (!inherits(x, "prior"))
#'     stop("Object must be of class 'prior'")
#'
#'   if (is.null(x$distribution))
#'     stop("prior() must include a 'distribution'")
#'
#'   out <- list(
#'     variable     = x$variable,
#'     traits    = x$traits,
#'     structure = x$structure
#'   )
#'
#'   ## ---- kernel (optional) -----------------------------------------------
#'   if (!is.null(x$kernel)) {
#'     out$kernel <- kernel_name(as_kernel(x$kernel))
#'   } else {
#'     out$kernel <- NULL
#'   }
#'
#'   ## ---- starting values -------------------------------------------------
#'   if (!is.null(x$start))
#'     out$start <- unclass(x$start)
#'
#'   ## ---- distribution ----------------------------------------------------
#'   out$distribution <- as_list(x$distribution)
#'
#'   out
#' }
#'
#'
#' #' Convert prior specifications to list form
#' #'
#' #' @param priors Named list of prior objects
#' #' @return PriorSpec list
#' #' @export
#' as_list.prior <- function(x) {
#'
#'   stopifnot(inherits(x, "prior"))
#'
#'   cov <- if (!is.null(x$kernel)) {
#'     list(
#'       type   = "kernel",
#'       kernel = kernel_name(x$kernel)
#'     )
#'   } else if (!is.null(x$featureMatrix)) {
#'     list(
#'       type          = "feature",
#'       featureMatrix = class(x$featureMatrix)[1],
#'       featureSets   = x$featureSets
#'     )
#'   } else if (x$variable == "residual") {
#'     list(type = "residual")
#'   } else {
#'     list(type = "iid")
#'   }
#'
#'   out <- list(
#'     variable     = x$variable,
#'     traits       = x$traits,
#'     covariance   = cov,
#'     structure    = x$structure,
#'     distribution = as_list(x$distribution)
#'   )
#'
#'   if (!is.null(x$start))
#'     out$start <- unclass(x$start)
#'
#'   out
#' }
#'
#'
#' #' Export prior specifications to JSON
#' #'
#' #' @param priors Named list of prior objects
#' #' @param pretty Logical
#' #' @return JSON string
#' #' @export
#' as_json.priors <- function(priors, pretty = TRUE) {
#'
#'   .require_jsonlite()
#'
#'   jsonlite::toJSON(
#'     as_list.priors(priors),
#'     pretty = pretty,
#'     auto_unbox = TRUE,
#'     null = "null"
#'   )
#' }
#'
#' #' Validate consistency of a model specification bundle
#' #'
#' #' @param data_spec DataSpec list
#' #' @param model_spec ModelSpec list
#' #' @param kernel_spec KernelSpec list
#' #' @param varcomp_spec VarCompSpec list (REML / solver)
#' #' @param prior_spec PriorSpec list (Bayesian)
#' #'
#' #' @export
#' validate_bundle_json <- function(data_spec,
#'                             model_spec,
#'                             kernel_spec,
#'                             varcomp_spec = NULL,
#'                             prior_spec   = NULL) {
#'
#'   ## ------------------------------------------------------------------
#'   ## Basic structure checks
#'   ## ------------------------------------------------------------------
#'   stopifnot(
#'     is.list(data_spec),  identical(data_spec$type,  "DataSpec"),
#'     is.list(model_spec), identical(model_spec$type, "ModelSpec"),
#'     is.list(kernel_spec),identical(kernel_spec$type,"KernelSpec")
#'   )
#'
#'   if (!xor(is.null(varcomp_spec), is.null(prior_spec))) {
#'     stop("Specify exactly one of VarCompSpec (REML/solver) or PriorSpec (Bayesian)")
#'   }
#'
#'   ## ------------------------------------------------------------------
#'   ## Kernel references
#'   ## ------------------------------------------------------------------
#'   defined_kernels <- names(kernel_spec$kernels)
#'
#'   components <- if (!is.null(varcomp_spec)) {
#'     varcomp_spec$components
#'   } else {
#'     prior_spec$components
#'   }
#'
#'   referenced_kernels <- unique(unlist(lapply(
#'     components,
#'     function(x) x$kernel
#'   )))
#'
#'   referenced_kernels <- referenced_kernels[!is.na(referenced_kernels)]
#'
#'   missing_kernels <- setdiff(referenced_kernels, defined_kernels)
#'
#'   if (length(missing_kernels) > 0) {
#'     stop(
#'       "Components reference undefined kernels: ",
#'       paste(missing_kernels, collapse = ", ")
#'     )
#'   }
#'
#'   ## ------------------------------------------------------------------
#'   ## Random-effect indices vs formulas
#'   ## ------------------------------------------------------------------
#'   model_indices <- unique(unlist(lapply(
#'     model_spec$summary$random,
#'     function(tr) {
#'       if (length(tr) == 0) return(character(0))
#'       vapply(tr, `[[`, character(1), "variable")
#'     }
#'   )))
#'
#'   spec_indices <- vapply(
#'     components,
#'     function(x) x$variable,
#'     character(1)
#'   )
#'
#'   ## Marker and Residual are allowed even if not in formula
#'   allowed_extra <- c("residual", "marker")
#'
#'   missing_specs <- setdiff(
#'     model_indices,
#'     spec_indices
#'   )
#'
#'   if (length(missing_specs) > 0) {
#'     stop(
#'       "Random effects in formulas but missing specification: ",
#'       paste(missing_specs, collapse = ", ")
#'     )
#'   }
#'
#'   unused_specs <- setdiff(
#'     spec_indices,
#'     c(model_indices, allowed_extra)
#'   )
#'
#'   if (length(unused_specs) > 0) {
#'     warning(
#'       "Specifications not referenced in formulas: ",
#'       paste(unused_specs, collapse = ", ")
#'     )
#'   }
#'
#'   ## ------------------------------------------------------------------
#'   ## Variable availability
#'   ## ------------------------------------------------------------------
#'   vars_needed <- names(model_spec$variables)
#'   vars_have   <- data_spec$colnames
#'
#'   missing_vars <- setdiff(vars_needed, vars_have)
#'
#'   if (length(missing_vars) > 0) {
#'     stop(
#'       "Variables required by the model are missing from the data source: ",
#'       paste(missing_vars, collapse = ", ")
#'     )
#'   }
#'
#'   ## ------------------------------------------------------------------
#'   ## Bayesian-only checks
#'   ## ------------------------------------------------------------------
#'   if (!is.null(prior_spec)) {
#'
#'     for (nm in names(prior_spec$components)) {
#'       p <- prior_spec$components[[nm]]
#'
#'       if (p$variable != "residual" && is.null(p$kernel)) {
#'         stop(
#'           "Prior '", nm,
#'           "' must specify a kernel (except for residual)."
#'         )
#'       }
#'     }
#'   }
#'
#'   invisible(TRUE)
#' }
