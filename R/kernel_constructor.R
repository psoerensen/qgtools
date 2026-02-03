# Kernels: covariance structures over the indexing units of a random effect
#          (e.g. animals, individuals, or SNPs)
#          Implemented as lists carrying kernel-specific information:
#            - Genomic relationship matrix (GRMlist)
#            - Pedigree-based relationship matrix (PEDlist)
#            - Linkage disequilibrium (LD) matrix (LDlist)


#' Create a PEDlist object
#'
#' Create a pedigree data specification object (\code{PEDlist}) used to define
#' the structure of additive genetic relationships based on pedigree information.
#'
#' The \code{PEDlist} object contains no phenotype or model-specific information.
#' Linking pedigree records to analysis data is handled at a higher (model) level.
#'
#' @export
makePEDlist <- function(fnPED = NULL,
                        format = "ASCII",
                        method = "S-D-NonInbred") {

  if (!is.null(fnPED) && (!is.character(fnPED) || length(fnPED) != 1)) {
    stop("'fnPED' must be a single character string.")
  }

  format <- toupper(format)
  if (!format %in% c("ASCII", "BINARY")) {
    stop("'format' must be either 'ASCII' or 'BINARY'.")
  }

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

  structure(
    list(
      kernel_type = "additive",
      fnPED       = fnPED,
      format      = format,
      method      = method,
      method_code = unname(method_map[method]),
      columns     = c("id", "sire", "dam_or_mgs", "order")
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

  ## ---- class & structure -----------------------------------------------------
  if (!inherits(x, "PEDlist")) {
    stop("Object must be of class 'PEDlist'.")
  }

  required_fields <- c(
    "kernel_type", "fnPED", "format",
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

  if (!x$format %in% c("ASCII", "BINARY")) {
    stop("Invalid 'format' in PEDlist: ", x$format)
  }

  ## ---- file existence --------------------------------------------------------
  if (check.file) {
    if (is.null(x$fnPED) || !nzchar(x$fnPED)) {
      stop("'fnPED' is NULL or empty.")
    }
    if (!file.exists(x$fnPED)) {
      stop("Pedigree file does not exist: ", x$fnPED)
    }
  }

  ## ---- lightweight structure check (safe for large files) -------------------
  if (check.format && check.file && x$format == "ASCII") {

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

  ## ---- optional deep logic checks (opt-in only) ------------------------------
  if (check.logic && check.file && x$format == "ASCII") {

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' is required for logical pedigree checks.")
    }

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
#'   (\code{"ASCII"} or \code{"BINARY"}).
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
  if (!format %in% c("ASCII", "BINARY")) {
    stop("'format' must be either 'ASCII' or 'BINARY'.")
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
#'   (\code{"ASCII"} or \code{"BINARY"}).
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
  if (!format %in% c("ASCII", "BINARY")) {
    stop("'format' must be either 'ASCII' or 'BINARY'.")
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
#' @param check.format Logical. Lightweight structure check (ASCII only).
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

  if (!x$format %in% c("ASCII", "BINARY")) {
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

  ## ---- lightweight structure check (ASCII only) -----------------------------
  if (check.format && check.file && x$format == "ASCII") {

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
        "ASCII GRM file must contain at least 3 columns ",
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
#'   (\code{"ASCII"} or \code{"BINARY"}).
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
  if (!format %in% c("ASCII", "BINARY")) {
    stop("'format' must be either 'ASCII' or 'BINARY'.")
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
#' @param check.format Logical. Lightweight structure check (ASCII only).
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

  if (!x$format %in% c("ASCII", "BINARY")) {
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

  ## ---- lightweight structure check (ASCII only) -----------------------------
  if (check.format && check.file && x$format == "ASCII") {

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
        "ASCII LD file must contain at least 3 columns ",
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
    extra  = list(...)
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
