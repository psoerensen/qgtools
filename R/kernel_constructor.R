# Kernels: covariance structures over the indexing units of a random effect
#          (e.g. animals, individuals, or SNPs)
#          Implemented as lists carrying kernel-specific information:
#            - Genomic relationship matrix (GRMlist)
#            - Pedigree-based relationship matrix (PEDlist)
#            - Linkage disequilibrium (LD) matrix (LDlist)


#' Create a PEDlist object
#'
#' Pedigree kernels are standalone objects and can be reused across variance
#' components and model fits.
#'
#' Create a pedigree data specification object (\code{PEDlist}) used to define
#' the structure of additive genetic relationships based on pedigree information.
#'
#' By default, \code{makePEDlist()} assumes a simple pedigree format with columns
#' identifying individual, sire, dam, and an optional ordering variable.
#'
#' The \code{PEDlist} object contains no phenotype or model-specific information.
#' Linking pedigree records to analysis data is handled at a higher (model) level.
#'
#' The pedigree file (or data frame) must contain information identifying
#' individuals and their parents. The exact interpretation of columns depends
#' on the selected pedigree \code{format}.
#'
#' For formats that support inbreeding, a sorting variable (e.g. birth date) may
#' be required; otherwise it may be set to 0.
#'
#' When phantom parent grouping is used, unknown parents must be replaced by
#' negative group codes to distinguish them from regular individual IDs.
#'
#' @param fnPED Character. Path to the pedigree file.
#' @param format Character. Pedigree format specification. One of
#'   \code{"DMU"}, \code{"SIMPLE"}, \code{"PLINK"}, or \code{"CUSTOM"}.
#' @param encoding Character. One of \code{"ASCII"} or \code{"BINARY"}.
#'   ASCII input is expected to conform to INTEGER*4 and REAL*4 ranges.
#'   Binary input must follow unformatted Fortran conventions.
#' @param method Character. Pedigree relationship construction method.
#' @param columns Optional named character vector or list defining column
#'   mappings for \code{format = "CUSTOM"}.
#'
#' @examples
#' ## DMU-compatible pedigree (default use case)
#' \dontrun{
#' PED_dmu <- makePEDlist(
#'   fnPED    = "pedigree.txt",
#'   format   = "DMU",
#'   encoding = "ASCII",
#'   method   = "S-D-NonInbred"
#' )
#' }
#'
#' ## Simple text pedigree (ID, sire, dam, optional order)
#' \dontrun{
#' PED_simple <- makePEDlist(
#'   fnPED  = "ped.txt",
#'   format = "SIMPLE"
#' )
#' }
#'
#' ## PLINK .fam pedigree
#' \dontrun{
#' PED_plink <- makePEDlist(
#'   fnPED  = "data.fam",
#'   format = "PLINK"
#' )
#' }
#'
#' ## Custom pedigree schema with user-defined column mapping
#' \dontrun{
#' PED_custom <- makePEDlist(
#'   fnPED  = "ped.txt",
#'   format = "CUSTOM",
#'   columns = list(
#'     id    = "animal",
#'     sire  = "father",
#'     dam   = "mother",
#'     order = "birth_year"
#'   )
#' )
#' }
#'
#' @export
makePEDlist <- function(fnPED = NULL,
                        format   = c("SIMPLE", "DMU", "PLINK", "CUSTOM"),
                        encoding = c("ASCII", "BINARY"),
                        method   = "S-D-NonInbred",
                        columns  = NULL) {

  format   <- match.arg(format)
  encoding <- match.arg(encoding)

  ## ---- basic checks ------------------------------------------------------
  if (!is.null(fnPED) && (!is.character(fnPED) || length(fnPED) != 1)) {
    stop("'fnPED' must be a single character string.")
  }

  if (format == "CUSTOM" && is.null(columns)) {
    stop("'columns' must be supplied when format = 'CUSTOM'.")
  }

  ## ---- method validation ------------------------------------------------
  allowed_methods <- c(
    "S-D-Inbred",
    "S-D-NonInbred",
    "S-MGS-Inbred",
    "S-MGS-NonInbred",
    "S-D-NonInbred-PP"
  )

  if (!method %in% allowed_methods) {
    stop(
      "Invalid 'method'. Must be one of: ",
      paste(allowed_methods, collapse = ", ")
    )
  }

  method_map <- c(
    "S-D-Inbred"       = 1,
    "S-D-NonInbred"    = 2,
    "S-MGS-Inbred"     = 3,
    "S-MGS-NonInbred"  = 4,
    "S-D-NonInbred-PP" = 6
  )

  ## ---- default column semantics -----------------------------------------
  default_columns <- switch(
    format,
    DMU    = c(id = "id", sire = "sire", dam = "dam_or_mgs", order = "order"),
    SIMPLE = c(id = "id", sire = "sire", dam = "dam", order = "order"),
    PLINK  = c(id = "IID", sire = "PID", dam = "MID"),
    CUSTOM = columns
  )

  ## ---- build object ------------------------------------------------------
  structure(
    list(
      kernel_type = "additive",
      fnPED       = fnPED,
      format      = format,
      encoding    = encoding,
      method      = method,
      method_code = unname(method_map[method]),
      columns     = default_columns
    ),
    class = "PEDlist"
  )
}

#' Validate a PEDlist object
#'
#' Perform consistency and integrity checks on a \code{PEDlist} object.
#'
#' @param x An object of class \code{"PEDlist"}.
#' @param check.file Logical. Check that the pedigree file exists.
#' @param check.format Logical. Lightweight structure check (ASCII only).
#' @param check.logic Logical. Optional deep pedigree logic checks (default FALSE).
#'
#' @export
validatePEDlist <- function(x,
                            check.file   = TRUE,
                            check.format = TRUE,
                            check.logic  = FALSE) {

  ## ---- class & structure -------------------------------------------------
  if (!inherits(x, "PEDlist")) {
    stop("Object must be of class 'PEDlist'.")
  }

  required_fields <- c(
    "kernel_type", "fnPED", "format", "encoding",
    "method", "method_code", "columns"
  )

  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop(
      "PEDlist is missing required field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (x$kernel_type != "additive") {
    stop("Invalid 'kernel_type' for PEDlist: ", x$kernel_type)
  }

  ## ---- format & encoding -------------------------------------------------
  valid_formats <- c("DMU", "SIMPLE", "PLINK", "CUSTOM")
  if (!x$format %in% valid_formats) {
    stop(
      "Invalid 'format' in PEDlist: ", x$format,
      ". Must be one of: ", paste(valid_formats, collapse = ", ")
    )
  }

  if (!x$encoding %in% c("ASCII", "BINARY")) {
    stop("Invalid 'encoding' in PEDlist: ", x$encoding)
  }

  ## ---- file existence ----------------------------------------------------
  if (check.file) {
    if (is.null(x$fnPED) || !nzchar(x$fnPED)) {
      stop("'fnPED' is NULL or empty.")
    }
    if (!file.exists(x$fnPED)) {
      stop("Pedigree file does not exist: ", x$fnPED)
    }
  }

  ## ---- lightweight structure check --------------------------------------
  ## Only meaningful for ASCII input
  if (check.format && check.file && x$encoding == "ASCII") {

    ped_head <- tryCatch(
      utils::read.table(
        x$fnPED,
        header = FALSE,
        nrows  = 5,
        stringsAsFactors = FALSE,
        comment.char = ""
      ),
      error = function(e) {
        stop("Failed to read pedigree header: ", conditionMessage(e))
      }
    )

    if (ncol(ped_head) < length(x$columns)) {
      stop(
        "Pedigree file must contain at least ",
        length(x$columns), " columns.\n",
        "Found ", ncol(ped_head), " column(s)."
      )
    }
  }

  ## ---- optional deep logic checks ----------------------------------------
  if (check.logic && check.file && x$encoding == "ASCII") {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' is required for logical pedigree checks.")
    }

    ## Read only ID and parent columns
    ped <- data.table::fread(
      x$fnPED,
      select = 1:3,
      showProgress = FALSE
    )

    if (any(ped[[1]] == ped[[2]] | ped[[1]] == ped[[3]], na.rm = TRUE)) {
      stop("Logical pedigree error: individual appears as its own parent.")
    }
  }

  invisible(TRUE)
}


#' @export
print.PEDlist <- function(x, ...) {
  cat("PEDlist object\n")
  cat("  Kernel type : ", x$kernel_type, "\n", sep = "")
  cat("  File        : ", x$fnPED, "\n", sep = "")
  cat("  Format      : ", x$format, "\n", sep = "")
  cat("  Encoding    : ", x$encoding, "\n", sep = "")
  cat("  Method      : ", x$method,
      " (code ", x$method_code, ")\n", sep = "")
  invisible(x)
}


#' Create a GRMlist object
#'
#' Create a genomic relationship matrix (GRM) specification object used to
#' define additive genetic relationships based on marker data.
#'
#' The \code{GRMlist} object contains no phenotype or model-specific information.
#' Linking the GRM to analysis data is handled at a higher (model) level.
#'
#' @param fnGRM Character string. File name of the GRM (or its inverse).
#' @param format Character string. Format of the GRM file
#'   (\code{"TEXT"} or \code{"BINARY"}).
#' @param grm_type Character string specifying the GRM representation.
#'   Typical values include \code{"G"}, \code{"Ginv"}, or \code{"scaled_G"}.
#'
#' @return An object of class \code{"GRMlist"}.
#'
#' @export
makeGRMlist <- function(fnGRM = NULL,
                        format = "BINARY",
                        grm_type = "G") {

  if (!is.null(fnGRM) && (!is.character(fnGRM) || length(fnGRM) != 1)) {
    stop("'fnGRM' must be a single character string.")
  }

  format <- toupper(format)
  if (!format %in% c("TEXT", "BINARY")) {
    stop("'format' must be either 'TEXT' or 'BINARY'.")
  }

  allowed_types <- c("G", "Ginv", "scaled_G")

  if (!grm_type %in% allowed_types) {
    stop(
      "Invalid 'grm_type'. Must be one of: ",
      paste(allowed_types, collapse = ", ")
    )
  }

  structure(
    list(
      kernel_type = "additive",
      fnGRM       = fnGRM,
      format      = format,
      grm_type    = grm_type
    ),
    class = "GRMlist"
  )
}

#' Create a GRMlist object
#'
#' Create a genomic relationship matrix (GRM) specification object used to
#' define additive genetic relationships based on marker data.
#'
#' The \code{GRMlist} object contains no phenotype or model-specific information.
#' Linking the GRM to analysis data is handled at a higher (model) level.
#'
#' @param fnGRM Character string. File name of the GRM (or its inverse).
#' @param format Character string. Format of the GRM file
#'   (\code{"TEXT"} or \code{"BINARY"}).
#' @param grm_type Character string specifying the GRM representation.
#'   Typical values include \code{"G"}, \code{"Ginv"}, or \code{"scaled_G"}.
#'
#' @return An object of class \code{"GRMlist"}.
#'
#' @export
makeGRMlist <- function(fnGRM = NULL,
                        format = "BINARY",
                        grm_type = "G") {

  if (!is.null(fnGRM) && (!is.character(fnGRM) || length(fnGRM) != 1)) {
    stop("'fnGRM' must be a single character string.")
  }

  format <- toupper(format)
  if (!format %in% c("TEXT", "BINARY")) {
    stop("'format' must be either 'TEXT' or 'BINARY'.")
  }

  allowed_types <- c("G", "Ginv", "scaled_G")

  if (!grm_type %in% allowed_types) {
    stop(
      "Invalid 'grm_type'. Must be one of: ",
      paste(allowed_types, collapse = ", ")
    )
  }

  structure(
    list(
      kernel_type = "additive",
      fnGRM       = fnGRM,
      format      = format,
      grm_type    = grm_type
    ),
    class = "GRMlist"
  )
}


#' Validate a GRMlist object
#'
#' Perform consistency and integrity checks on a \code{GRMlist} object.
#'
#' @param x An object of class \code{"GRMlist"}.
#' @param check.file Logical. Check that the GRM file exists.
#' @param check.format Logical. Lightweight structure check (TEXT only).
#'
#' @return Invisibly returns TRUE if validation succeeds.
#'
#' @export
validateGRMlist <- function(x,
                            check.file   = TRUE,
                            check.format = TRUE) {

  ## ---- class & structure -----------------------------------------------------
  if (!inherits(x, "GRMlist")) {
    stop("Object must be of class 'GRMlist'.")
  }

  required_fields <- c("kernel_type", "fnGRM", "format", "grm_type")
  missing_fields  <- setdiff(required_fields, names(x))

  if (length(missing_fields) > 0) {
    stop(
      "GRMlist is missing required field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (x$kernel_type != "additive") {
    stop("Invalid 'kernel_type' for GRMlist: ", x$kernel_type)
  }

  if (!x$format %in% c("TEXT", "BINARY")) {
    stop("Invalid 'format' in GRMlist: ", x$format)
  }

  ## ---- file existence --------------------------------------------------------
  if (check.file) {
    if (is.null(x$fnGRM) || !nzchar(x$fnGRM)) {
      stop("'fnGRM' is NULL or empty.")
    }
    if (!file.exists(x$fnGRM)) {
      stop("GRM file does not exist: ", x$fnGRM)
    }
  }

  ## ---- lightweight structure check (TEXT only) -----------------------------
  if (check.format && check.file && x$format == "TEXT") {

    grm_head <- tryCatch(
      utils::read.table(
        x$fnGRM,
        header = FALSE,
        nrows  = 5,
        stringsAsFactors = FALSE,
        comment.char = ""
      ),
      error = function(e) {
        stop("Failed to read GRM header: ", conditionMessage(e))
      }
    )

    if (ncol(grm_head) < 3) {
      stop(
        "TEXT GRM file must contain at least 3 columns ",
        "(row, column, value)."
      )
    }
  }

  invisible(TRUE)
}

#' @export
print.GRMlist <- function(x, ...) {
  cat("GRMlist object\n")
  cat("  Kernel type : ", x$kernel_type, "\n", sep = "")
  cat("  File        : ", x$fnGRM, "\n", sep = "")
  cat("  Format      : ", x$format, "\n", sep = "")
  cat("  GRM type    : ", x$grm_type, "\n", sep = "")
  invisible(x)
}



#' Create an LDlist object
#'
#' Create a linkage disequilibrium (LD) matrix specification object used to
#' define covariance between SNP effects.
#'
#' The \code{LDlist} object contains no phenotype or model-specific information.
#' Linking SNPs to genotype or summary statistic data is handled at a higher
#' (model) level.
#'
#' @param fnLD Character string. File name of the LD matrix (or its inverse).
#' @param format Character string. Format of the LD file
#'   (\code{"TEXT"} or \code{"BINARY"}).
#' @param ld_type Character string specifying the LD representation.
#'   Typical values include \code{"R"}, \code{"Rinv"}, \code{"sparse_R"},
#'   or \code{"banded_R"}.
#' @param snp_index Optional character string specifying the SNP indexing
#'   convention (e.g. reference panel or chromosome).
#'
#' @return An object of class \code{"LDlist"}.
#'
#' @export
makeLDlist <- function(fnLD = NULL,
                       format = "BINARY",
                       ld_type = "R",
                       snp_index = NULL) {

  if (!is.null(fnLD) && (!is.character(fnLD) || length(fnLD) != 1)) {
    stop("'fnLD' must be a single character string.")
  }

  format <- toupper(format)
  if (!format %in% c("TEXT", "BINARY")) {
    stop("'format' must be either 'TEXT' or 'BINARY'.")
  }

  allowed_types <- c(
    "R",
    "Rinv",
    "sparse_R",
    "banded_R",
    "block_R"
  )

  if (!ld_type %in% allowed_types) {
    stop(
      "Invalid 'ld_type'. Must be one of: ",
      paste(allowed_types, collapse = ", ")
    )
  }

  structure(
    list(
      kernel_type = "snp",
      fnLD        = fnLD,
      format      = format,
      ld_type     = ld_type,
      snp_index   = snp_index
    ),
    class = "LDlist"
  )
}



#' Validate an LDlist object
#'
#' Perform consistency and integrity checks on an \code{LDlist} object.
#'
#' @param x An object of class \code{"LDlist"}.
#' @param check.file Logical. Check that the LD file exists.
#' @param check.format Logical. Lightweight structure check (TEXT only).
#'
#' @return Invisibly returns TRUE if validation succeeds.
#'
#' @export
validateLDlist <- function(x,
                           check.file   = TRUE,
                           check.format = TRUE) {

  ## ---- class & structure -----------------------------------------------------
  if (!inherits(x, "LDlist")) {
    stop("Object must be of class 'LDlist'.")
  }

  required_fields <- c(
    "kernel_type", "fnLD", "format", "ld_type"
  )

  missing_fields <- setdiff(required_fields, names(x))
  if (length(missing_fields) > 0) {
    stop(
      "LDlist is missing required field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (x$kernel_type != "snp") {
    stop("Invalid 'kernel_type' for LDlist: ", x$kernel_type)
  }

  if (!x$format %in% c("TEXT", "BINARY")) {
    stop("Invalid 'format' in LDlist: ", x$format)
  }

  ## ---- file existence --------------------------------------------------------
  if (check.file) {
    if (is.null(x$fnLD) || !nzchar(x$fnLD)) {
      stop("'fnLD' is NULL or empty.")
    }
    if (!file.exists(x$fnLD)) {
      stop("LD file does not exist: ", x$fnLD)
    }
  }

  ## ---- lightweight structure check (TEXT only) -----------------------------
  if (check.format && check.file && x$format == "TEXT") {

    ld_head <- tryCatch(
      utils::read.table(
        x$fnLD,
        header = FALSE,
        nrows  = 5,
        stringsAsFactors = FALSE,
        comment.char = ""
      ),
      error = function(e) {
        stop("Failed to read LD header: ", conditionMessage(e))
      }
    )

    if (ncol(ld_head) < 3) {
      stop(
        "TEXT LD file must contain at least 3 columns ",
        "(row, column, value)."
      )
    }
  }

  invisible(TRUE)
}


#' @export
print.LDlist <- function(x, ...) {
  cat("LDlist object\n")
  cat("  Kernel type : ", x$kernel_type, "\n", sep = "")
  cat("  File        : ", x$fnLD, "\n", sep = "")
  cat("  Format      : ", x$format, "\n", sep = "")
  cat("  LD type     : ", x$ld_type, "\n", sep = "")
  if (!is.null(x$snp_index)) {
    cat("  SNP index   : ", x$snp_index, "\n", sep = "")
  }
  invisible(x)
}


