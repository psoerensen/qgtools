
.qg_stop <- function(...) stop(..., call. = FALSE)

.qg_validate_vcov_entry <- function(name, x, Klist_keys = NULL) {
  req <- c("index","traits","kernel","structure")
  missing <- setdiff(req, names(x))
  if (length(missing)) .qg_stop("vcov.", name, " missing fields: ", paste(missing, collapse = ", "))
  if (!is.character(x$index) || length(x$index) != 1) .qg_stop("vcov.", name, ".index must be a single string")
  if (!is.character(x$kernel) || length(x$kernel) != 1) .qg_stop("vcov.", name, ".kernel must be a single string key into Klist")
  if (!is.character(x$structure) || length(x$structure) != 1) .qg_stop("vcov.", name, ".structure must be a single string")
  if (!is.character(x$traits)) .qg_stop("vcov.", name, ".traits must be a character vector")

  if (!is.null(Klist_keys) && !(x$kernel %in% Klist_keys)) {
    .qg_stop("vcov.", name, ".kernel = '", x$kernel, "' not found in Klist (available: ", paste(Klist_keys, collapse = ", "), ")")
  }
  invisible(TRUE)
}

.qg_validate_random_term <- function(rt, path = "random") {
  req <- c("intercept","slopes","index","raw")
  if (!is.list(rt)) .qg_stop(path, " must be a list.")
  missing <- setdiff(req, names(rt))
  if (length(missing)) .qg_stop(path, " missing fields: ", paste(missing, collapse = ", "))
  if (!is.logical(rt$intercept) || length(rt$intercept) != 1L) .qg_stop(path, ".intercept must be TRUE/FALSE")
  if (!is.character(rt$slopes)) .qg_stop(path, ".slopes must be a character vector")
  if (!is.character(rt$index) || length(rt$index) != 1L) .qg_stop(path, ".index must be a single string")
  if (!is.character(rt$raw) || length(rt$raw) != 1L) .qg_stop(path, ".raw must be a single string")
  invisible(TRUE)
}

.qg_validate_spec <- function(spec) {
  if (!is.list(spec)) .qg_stop("Specification must be a list.")
  if (is.null(spec$models)) .qg_stop("Specification must include `models`.")
  if (is.null(spec$vcov)) spec$vcov <- list()
  if (!is.list(spec$models)) .qg_stop("`models` must be a named list of traits.")

  # validate Klist first (so we can validate vcov$kernel)
  if (is.null(spec$Klist)) spec$Klist <- list()
  if (!(is.list(spec$Klist) || inherits(spec$Klist, "Klist"))) .qg_stop("`Klist` must be a list or Klist object.")
  if (!inherits(spec$Klist, "Klist")) {
    # Validate but don't force class; keep backward compatibility.
    validate_Klist(spec$Klist, strict = TRUE)
  } else {
    validate_Klist(spec$Klist, strict = TRUE)
  }

  kkeys <- names(spec$Klist)

  # validate models
  for (tr in names(spec$models)) {
    m <- spec$models[[tr]]
    if (!is.list(m)) .qg_stop("models.", tr, " must be a list.")
    if (is.null(m$fixed)) m$fixed <- character()
    if (is.null(m$random)) m$random <- list()
    if (!is.character(m$fixed)) .qg_stop("models.", tr, ".fixed must be character vector")

    if (is.character(m$random)) {
      # YAML-friendly: allow random terms as strings and coerce here
      m$random <- lapply(m$random, .qg_parse_random_string)
      spec$models[[tr]]$random <- m$random
    }
    if (!is.list(m$random)) .qg_stop("models.", tr, ".random must be a list of random terms.")
    if (length(m$random)) {
      for (i in seq_along(m$random)) {
        .qg_validate_random_term(m$random[[i]], path = paste0("models.", tr, ".random[[", i, "]]"))
      }
    }
  }

  # validate vcov
  if (!is.list(spec$vcov)) .qg_stop("`vcov` must be a list.")
  if (length(spec$vcov)) {
    if (is.null(names(spec$vcov)) || any(names(spec$vcov) == "")) {
      .qg_stop("`vcov` must be a *named* list (component names like 'polygenic', 'pedigree', ...).")
    }
    for (nm in names(spec$vcov)) .qg_validate_vcov_entry(nm, spec$vcov[[nm]], Klist_keys = kkeys)
  }

  spec
}

#' Compile a qgtools specification into an engine-agnostic model object
#'
#' @param models Either a list of formulas (Wilkinson) OR a YAML spec list.
#' @param vcov Optional vcov list (overrides / complements YAML).
#' @param Klist Optional kernel list (overrides / complements YAML).
#' @param data Optional data.frame (can be NULL when engines read from disk).
#' @param datalist Optional list describing data on disk (engine-specific).
#' @return An object of class \code{qg_model}.
#' @export
qg_compile <- function(models,
                       vcov = NULL,
                       Klist = NULL,
                       data = NULL,
                       datalist = NULL) {

  # Case 1: YAML list
  if (is.list(models) && !inherits(models, "formula") && !is.null(models$models)) {
    spec <- models
  } else {
    # Case 2: formulas -> spec
    parsed <- .qg_parse_formulas(models)
    spec <- list(
      models = lapply(parsed, function(p) list(fixed = p$fixed, random = p$random)),
      vcov = list(),
      Klist = list()
    )
  }

  # merge overrides
  if (!is.null(vcov)) spec$vcov <- vcov
  if (!is.null(Klist)) spec$Klist <- Klist
  if (!is.null(data)) spec$data <- list(source = "env", frame = "<provided>", n = NROW(data))
  if (!is.null(datalist)) spec$datalist <- datalist

  spec <- .qg_validate_spec(spec)

  structure(
    list(
      spec = spec,
      models = spec$models,
      vcov = spec$vcov,
      Klist = spec$Klist,
      data = data,
      datalist = datalist
    ),
    class = "qg_model"
  )
}

#' Fit a qgtools model with a chosen engine
#'
#' @param models Either a list of formulas (Wilkinson) or a YAML spec (list) or path to YAML file.
#' @param data Optional data.frame in memory.
#' @param datalist Optional disk-backed data descriptor (engine-specific).
#' @param vcov Variance component / covariance specification list.
#' @param Klist List of kernel objects / references (e.g., Glist/GRMlist/PEDlist).
#' @param engine Name of the computational engine.
#' @param ... Passed to the engine.
#' @export
qg_fit <- function(models,
                   data = NULL,
                   datalist = NULL,
                   vcov = NULL,
                   Klist = NULL,
                   engine = "mock",
                   ...) {

  # allow YAML path
  if (is.character(models) && length(models) == 1 && grepl("\\.ya?ml$", models, ignore.case = TRUE)) {
    models <- read_qg_yaml(models)
  }

  compiled <- qg_compile(models = models, vcov = vcov, Klist = Klist, data = data, datalist = datalist)

  if (!exists(engine, envir = .qg_engines, inherits = FALSE)) {
    .qg_stop("Unknown engine '", engine, "'. Available engines: ", paste(qg_list_engines(), collapse = ", "))
  }
  fun <- get(engine, envir = .qg_engines, inherits = FALSE)
  fit <- fun(compiled, ...)
  if (!inherits(fit, "qg_fit")) {
    fit$engine <- engine
    class(fit) <- unique(c("qg_fit", class(fit)))
  }
  fit
}

# printing
#' @export
print.qg_model <- function(x, ...) {
  cat("<qg_model>\n")
  cat("Traits:", paste(names(x$models), collapse = ", "), "\n")
  cat("VC components:", if (length(x$vcov)) paste(names(x$vcov), collapse = ", ") else "(none)", "\n")
  cat("Klist keys:", if (length(x$Klist)) paste(names(x$Klist), collapse = ", ") else "(none)", "\n")
  invisible(x)
}

#' @export
print.qg_fit <- function(x, ...) {
  cat("<qg_fit> engine:", x$engine, "\n")
  if (!is.null(x$compiled)) print(x$compiled)
  invisible(x)
}
