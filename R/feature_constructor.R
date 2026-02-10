#' Create a FeatureSource object
#'
#' Describe one or more high-dimensional feature sources (e.g. genotypes,
#' omics, image-derived features) used in a qgtools model.
#'
#' Each feature source must define:
#'   - row identifiers (e.g. individual IDs)
#'   - column identifiers (e.g. marker IDs, features)
#'
#' These identifiers are used to align features with the DataSource
#' (via row IDs) and with model components such as variance components
#' and feature sets (via column IDs).
#'
#' @param ... Named feature objects (e.g. Glist, FeatureBackend)
#'
#' @return An object of class "FeatureSource"
#' @export
makeFeatureSource <- function(...) {

  features <- list(...)

  if (length(features) == 0)
    stop("At least one feature must be supplied")

  if (is.null(names(features)) || any(names(features) == ""))
    stop("All features must be named")

  is_valid_ids <- function(x) {
    is.character(x) || (is.list(x) && !is.null(x$type))
  }

  for (nm in names(features)) {

    f <- features[[nm]]

    if (!inherits(f, c("Glist", "FeatureBackend"))) {
      stop(
        "Feature '", nm,
        "' must inherit from class 'Glist' or 'FeatureBackend'"
      )
    }

    ## Row identifiers
    if (is.null(f$ids) || !is_valid_ids(f$ids)) {
      stop(
        "Feature '", nm,
        "' must provide 'ids' as a character vector or symbolic descriptor"
      )
    }

    ## Column identifiers
    if (is.null(f$rsids) || !is_valid_ids(f$rsids)) {
      stop(
        "Feature '", nm,
        "' must provide 'rsids' as a character vector or symbolic descriptor"
      )
    }

    ## Uniqueness checks (only if materialized)
    if (is.character(f$ids) && anyDuplicated(f$ids)) {
      stop("Duplicate row IDs detected in feature '", nm, "'")
    }

    if (is.character(f$rsids) && anyDuplicated(f$rsids)) {
      stop("Duplicate column IDs detected in feature '", nm, "'")
    }
  }

  structure(
    list(
      type     = "FeatureSpec",
      features = features
    ),
    class = "FeatureSource"
  )
}

#' @export
as_list.Glist <- function(x) {

  list(
    backend = "PLINK",
    files = list(
      bed = x$fnBED,
      bim = x$fnBIM,
      fam = x$fnFAM
    ),
    format = x$format,
    extra  = x$extra
  )
}


#' Create a Glist object
#'
#' Create a genotype data specification object (\code{Glist}) based on PLINK
#' \code{.bed}, \code{.bim}, and \code{.fam} files. The object describes the
#' location and format of genotype data but contains no model- or
#' phenotype-specific information.
#'
#' Linking genotype records to analysis data is handled at a higher (model)
#' level.
#'
#' @param fnBED Character vector. One or more PLINK \code{.bed} files.
#'
#' @param fnBIM Character vector. One or more PLINK \code{.bim} files.
#'
#' @param fnFAM Character vector. One or more PLINK \code{.fam} files.
#'
#' @param format Character string specifying the genotype file format.
#'   Currently only \code{"PLINK"} is supported.
#'
#' @param ... Additional optional arguments reserved for future extensions.
#'
#' @return
#' An object of class \code{"Glist"} describing a genotype data source.
#'
#' @export
makeGlist <- function(fnBED = NULL,
                      fnBIM = NULL,
                      fnFAM = NULL,
                      format = "PLINK",
                      ...) {

  ## ---- input checks ----------------------------------------------------------
  if (is.null(fnBED) || is.null(fnBIM) || is.null(fnFAM)) {
    stop("fnBED, fnBIM, and fnFAM must be provided.")
  }

  if (!is.character(fnBED) || !is.character(fnBIM) || !is.character(fnFAM)) {
    stop("fnBED, fnBIM, and fnFAM must be character vectors.")
  }

  nfiles <- unique(c(length(fnBED), length(fnBIM), length(fnFAM)))
  if (length(nfiles) > 1) {
    stop("'fnBED', 'fnBIM', and 'fnFAM' must have the same length.")
  }

  format <- toupper(format)
  if (format != "PLINK") {
    stop("Unsupported 'format': ", format)
  }

  ## ---- build object ----------------------------------------------------------
  Glist <- list(
    fnBED  = fnBED,
    fnBIM  = fnBIM,
    fnFAM  = fnFAM,
    format = format,

    ## NEW: symbolic identifiers
    ids = list(
      type = "plink_fam",
      fields = c("FID", "IID")
    ),

    rsids = list(
      type = "plink_bim",
      field = "rsid"
    ),

    extra = list(...)
  )

  class(Glist) <- "Glist"
  Glist
}


#' Validate a Glist object
#'
#' Perform consistency and integrity checks on a \code{Glist} object.
#' This function validates only the genotype file specification.
#' Alignment with phenotype or pedigree data is handled at the model level.
#'
#' @param x An object of class \code{"Glist"}.
#' @param check.file Logical. If \code{TRUE} (default), checks that genotype files exist.
#'
#' @return
#' Invisibly returns \code{TRUE} if validation succeeds.
#'
#' @export
validateGlist <- function(x, check.file = TRUE) {

  ## ---- class check -----------------------------------------------------------
  if (!inherits(x, "Glist")) {
    stop("Object must be of class 'Glist'.")
  }

  ## ---- required fields -------------------------------------------------------
  required_fields <- c("fnBED", "fnBIM", "fnFAM", "format")
  missing_fields  <- setdiff(required_fields, names(x))

  if (length(missing_fields) > 0) {
    stop(
      "Glist is missing required field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  ## ---- format check ----------------------------------------------------------
  if (x$format != "PLINK") {
    stop("Invalid 'format' in Glist: ", x$format)
  }

  ## ---- file length consistency ----------------------------------------------
  if (!(length(x$fnBED) == length(x$fnBIM) &&
        length(x$fnBED) == length(x$fnFAM))) {
    stop("fnBED, fnBIM, and fnFAM must have the same length.")
  }

  ## ---- file existence --------------------------------------------------------
  if (check.file) {
    files <- c(x$fnBED, x$fnBIM, x$fnFAM)
    missing <- files[!file.exists(files)]

    if (length(missing) > 0) {
      stop(
        "The following genotype file(s) do not exist:\n",
        paste(missing, collapse = "\n")
      )
    }
  }

  invisible(TRUE)
}

