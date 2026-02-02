# ============================================================
# qgtools: YAML I/O
# ============================================================

#' Read a qgtools YAML specification
#'
#' @param path Path to a YAML file
#' @return A list representing the qgtools model specification
#' @export
read_qg_yaml <- function(path) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read YAML files")
  }
  if (!file.exists(path)) {
    stop("YAML file does not exist: ", path)
  }

  yaml::read_yaml(path)
}

#' Write qgtools model specification to YAML
#'
#' @param formulas Named list of formulas, one per trait
#' @param vcov List specifying trait-level covariance structures
#' @param Klist Optional Klist object (used only for kernel names)
#' @param file Output YAML filename
#'
#' @return Invisibly returns the YAML object
#' @export
write_qg_yaml <- function(formulas, vcov = NULL, Klist = NULL, file) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to write YAML files")
  }

  # --- Parse formulas using canonical parser ---
  models <- qg_parse_formulas(formulas)

  yaml_obj <- list(
    models = models,
    vcov   = vcov,
    kernels = if (!is.null(Klist)) names(Klist) else NULL
  )

  yaml::write_yaml(yaml_obj, file)
  invisible(yaml_obj)
}
