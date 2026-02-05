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
      colnames = colnames,
      options  = list(...)
    ),
    class = "Datalist"
  )
}
