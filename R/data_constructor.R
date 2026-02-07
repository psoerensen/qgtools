#' Create a data source specification object
#'
#' Create a \code{Datalist} object describing how phenotype and covariate data
#' are provided to a model. The data may be supplied as an in-memory data frame
#' or as disk-backed files in a variety of common formats.
#'
#' qgtools requires that all data variables have explicit names. Column names may
#' be obtained either from headers in the data source or supplied directly by the
#' user.
#'
#' For formats that do not include headers (e.g. DMU-style input), column names
#' must be provided explicitly using the \code{colnames} argument. When headers
#' are present, they are used by default unless overridden.
#'
#' Internally, all variables are accessed by name rather than by position. This
#' ensures robust model specification, avoids ambiguity, and allows the same
#' model to be applied consistently across different data formats.
#'
#' Variable roles (e.g. traits, covariates, identifiers) may be specified
#' explicitly using the \code{roles} argument or inferred later from model
#' formulas. Conflicting or ambiguous specifications result in an error.
#'
#' For DMU-style input, data may be provided in either ASCII (text) or binary
#' format. ASCII input is expected to conform to INTEGER*4 and REAL*4 ranges,
#' while binary input must follow the conventions of unformatted Fortran files
#' on the target system. Binary input places responsibility for value ranges and
#' portability on the user.
#'
#' @param source Data source. Either a data frame (for \code{format = "DATAFRAME"}),
#'   a file path, or a named list of file paths depending on the chosen format.
#'
#' @param format Character. Data format. One of \code{"DATAFRAME"},
#'   \code{"CSV"}, \code{"DMU"}, \code{"PLINK"}, or \code{"PARQUET"}.
#'
#' @param encoding Character. One of \code{"ASCII"} or \code{"BINARY"}.
#'   Only meaningful for \code{format = "DMU"}; ignored for other formats.
#'
#' @param id Optional character. Name of the individual identifier variable.
#'
#' @param roles Optional named list specifying variable roles (e.g. traits,
#'   covariates, weights). If not supplied, roles may be inferred from model
#'   formulas.
#'
#' @param missing Optional missing-value specification. The interpretation
#'   depends on the data format.
#'
#' @param header Logical. Whether the data source contains column headers.
#'   Defaults depend on the selected format.
#'
#' @param colnames Optional character vector or named list of column names.
#'   Required when \code{header = FALSE}.
#'
#' @param ... Additional format-specific options, stored but not interpreted
#'   at construction time.
#'
#' @return An object of class \code{"Datalist"}.
#'
#' @examples
#' ## In-memory data frame
#' data_df <- makeDatalist(
#'   source = mouse,
#'   format = "DATAFRAME"
#' )
#'
#' ## CSV file with header
#' data_csv <- makeDatalist(
#'   source = "mouse.csv",
#'   format = "CSV",
#'   id     = "id"
#' )
#'
#' ## CSV file without header
#' data_csv_noheader <- makeDatalist(
#'   source   = "data.txt",
#'   format   = "CSV",
#'   header   = FALSE,
#'   colnames = c("id", "sex", "reps", "BW")
#' )
#'
#' ## DMU-style input with separate integer and real files
#' data_dmu <- makeDatalist(
#'   source = list(
#'     integer = "data_int.txt",
#'     real    = "data_real.txt"
#'   ),
#'   format   = "DMU",
#'   encoding = "ASCII",
#'   colnames = list(
#'     integer = c("id", "dam", "sex"),
#'     real    = c("BW")
#'   ),
#'   missing = list(
#'     integer = 0,
#'     real    = -999
#'   )
#' )
#'
#' @export
makeDatalist <- function(source,
                         format   = c("DATAFRAME", "CSV", "DMU", "PLINK", "PARQUET"),
                         encoding = c("ASCII", "BINARY"),
                         id = NULL,
                         roles = NULL,
                         missing = NULL,
                         header = NULL,
                         colnames = NULL,
                         ...) {
  ## ---- match arguments --------------------------------------------------
  format   <- match.arg(format)
  encoding <- match.arg(encoding)

  ## ---- basic checks -----------------------------------------------------
  if (is.null(source))
    stop("'source' must be provided")

  if (!is.null(id) && (!is.character(id) || length(id) != 1))
    stop("'id' must be a single character string or NULL")

  if (!is.null(roles) && !is.list(roles))
    stop("'roles' must be a named list or NULL")

  ## ---- defaults for header ----------------------------------------------
  if (is.null(header)) {
    header <- format %in% c("DATAFRAME", "CSV", "PLINK", "PARQUET")
  }

  ## ---- enforce column names ---------------------------------------------
  if (!header && is.null(colnames)) {
    stop(
      "Column names must be supplied via 'colnames' when the data source ",
      "does not contain headers."
    )
  }

  ## ---- format-specific validation ---------------------------------------
  spec <- switch(
    format,

    DATAFRAME = {
      if (!is.data.frame(source))
        stop("For format = 'DATAFRAME', 'source' must be a data frame")
      list(data = source)
    },

    CSV = {
      if (!is.character(source) || length(source) != 1)
        stop("For format = 'CSV', 'source' must be a file path")
      list(file = source)
    },

    DMU = {
      if (!(is.character(source) || is.list(source)))
        stop("For format = 'DMU', 'source' must be a file path or a named list")
      list(files = source)
    },

    PLINK = {
      if (!is.list(source) || !"fam" %in% names(source))
        stop("PLINK source must be a named list including at least 'fam'")
      list(files = source)
    },

    PARQUET = {
      if (!is.character(source) || length(source) != 1)
        stop("For format = 'PARQUET', 'source' must be a file path")
      list(file = source)
    }
  )

  ## ---- resolve column names -----------------------------------------------
  resolved_colnames <- colnames

  if (is.null(resolved_colnames)) {

    resolved_colnames <- switch(
      format,

      DATAFRAME = {
        names(source)
      },

      CSV = {
        if (header) {

          opts <- list(...)

          utils::read.table(
            source,
            header = TRUE,
            nrows  = 0,
            sep    = opts$sep %||% ",",
            quote  = opts$quote %||% "\"",
            comment.char = opts$comment.char %||% "",
            stringsAsFactors = FALSE,
            check.names = FALSE
          ) |> names()

        } else {
          NULL
        }
      },

      PARQUET = {
        # do NOT read full file
        # backend or arrow can resolve later
        NULL
      },

      DMU = {
        # must already be supplied
        colnames
      },

      PLINK = {
        # minimal standard .fam schema
        c("FID", "IID", "PID", "MID", "SEX", "PHENOTYPE")
      }
    )
  }

  ## ---- build object -----------------------------------------------------
  structure(
    list(
      format   = format,
      encoding = encoding,   # only meaningful for DMU
      source   = spec,
      id       = id,
      roles    = roles,
      missing  = missing,
      header   = header,
      colnames = resolved_colnames,
      options  = list(...)
    ),
    class = "Datalist"
  )
}
