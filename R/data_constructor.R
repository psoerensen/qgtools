#' Create a data source specification object
#'
#' qgtools requires that all data variables have explicit names. Column names may
#' be obtained either from headers in the data source or supplied directly by the
#' user.
#'
#' For data formats that do not include headers (e.g. DMU-style input), column
#' names must be provided explicitly using the \code{colnames} argument. When
#' headers are present, they are used by default unless overridden.
#'
#' Internally, all variables are accessed by name rather than by position. This
#' ensures robust model specification, avoids ambiguity, and allows the same
#' model specification to be applied consistently across different data formats.
#'
#' Variable roles (e.g. traits, covariates, identifiers) may be specified
#' explicitly using the \code{roles} argument or inferred from model formulas.
#' Conflicting or ambiguous specifications result in an error.
#'
#'
#' @param source Data source (data.frame, file path, or named list)
#' @param format One of "DATAFRAME", "CSV", "DMU", "PLINK", "PARQUET"
#' @param id Optional name of individual identifier variable
#' @param roles Optional named list of variable roles
#' @param missing Optional missing-value specification
#' @param header Logical. Whether the data source contains column headers.
#'   Defaults depend on format.
#' @param colnames Optional character vector or named list of column names.
#'   Required if headers are not present.
#' @param ... Format-specific options
#'
#' @return An object of class "Datalist"
#' @examples
#' ## In-memory data frame (column names already present)
#' data_df <- makeDatalist(
#'   source = mouse,
#'   format = "DATAFRAME"
#' )
#'
#' ## CSV file with header (default behavior)
#' data_csv <- makeDatalist(
#'   source = "data.csv",
#'   format = "CSV"
#' )
#'
#' ## CSV file without header (column names must be supplied)
#' data_csv_noheader <- makeDatalist(
#'   source   = "data.txt",
#'   format   = "CSV",
#'   header   = FALSE,
#'   colnames = c("id", "sex", "reps", "BW", "Gl")
#' )
#'
#' ## DMU-style input with separate integer and real files
#' data_dmu <- makeDatalist(
#'   source = list(
#'     integer = "data_int.txt",
#'     real    = "data_real.txt"
#'   ),
#'   format   = "DMU",
#'   colnames = list(
#'     integer = c("id", "dam", "sex", "reps"),
#'     real    = c("BW", "Gl")
#'   ),
#'   missing = list(
#'     integer = 0,
#'     real    = -999
#'   )
#' )
#'
#' ## PLINK phenotype and covariate input
#' data_plink <- makeDatalist(
#'   source = list(
#'     fam        = "data.fam",
#'     covariates = "covariates.txt"
#'   ),
#'   format = "PLINK",
#'   id     = "IID"
#' )
#'
#' ## Parquet file (disk-backed, columnar format)
#' data_parquet <- makeDatalist(
#'   source = "data.parquet",
#'   format = "PARQUET"
#' )
#' @export
makeDatalist <- function(source,
                         format = c("DATAFRAME", "CSV", "DMU", "PLINK", "PARQUET"),
                         id = NULL,
                         roles = NULL,
                         missing = NULL,
                         header = NULL,
                         colnames = NULL,
                         ...) {

  format <- match.arg(format)

  ## ---- basic checks ------------------------------------------------------
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

  ## ---- build object ------------------------------------------------------
  structure(
    list(
      format   = format,
      source   = spec,
      id       = id,
      roles    = roles,
      missing  = missing,
      header   = header,
      colnames = colnames,
      options  = list(...)
    ),
    class = "Datalist"
  )
}
