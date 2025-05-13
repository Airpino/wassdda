#' @title DDtable class
#' @description
#' A flexible class to store heterogeneous tabular data where each column may be of a different type,
#' including custom types such as `MatH` (histogram-based distributions) and `MatD` (discrete distributions).
#' This class mimics the behavior of a data.frame, allowing column-wise heterogeneity,
#' with consistency enforced across rows.
#'
#' @slot data A named list of column vectors or matrix-like distribution objects (`MatH`, `MatD`, etc.).
#' @slot col.types A character vector indicating the class of each column (e.g., "distributionH", "distributionD", "numeric", "factor").
#' @slot row.names A character vector with row identifiers.
#' @slot col.names A character vector with column identifiers.
#'
#' @section Validity:
#' - All columns must have the same number of rows.
#' - `col.names` must have the same length as the number of columns.
#' - `row.names` must have the same length as the number of rows.
#'
#' @seealso [MatH-class], [MatD-class]
#' @export
setClass(
  Class = "DDtable",
  slots = c(
    data = "list",
    col.types = "character",
    row.names = "character",
    col.names = "character"
  ),
  validity = function(object) {
    cols <- object@data
    if (length(cols) == 0) return(TRUE)

    # Lunghezze
    lengths <- sapply(cols, function(c) {
      if (is(c, "MatH") || is(c, "MatD")) {
        nrow(c@M)
      } else {
        length(c)
      }
    })

    if (length(unique(lengths)) != 1) {
      return("All columns must have the same number of rows")
    }

    # Verifica nomi
    if (length(object@col.names) != length(cols)) {
      return("Length of 'col.names' must match number of columns")
    }
    if (length(object@row.names) != lengths[1]) {
      return("Length of 'row.names' must match number of rows")
    }

    TRUE
  }
)

#' @title Constructor for DDtable objects
#' @description
#' Creates a new `DDtable` object by combining multiple columns of potentially heterogeneous types,
#' such as distributions (`MatH`, `MatD`), numeric vectors, factors, or characters. All columns must
#' have the same number of rows.
#'
#' @param ... Named arguments representing columns. Supported types include `MatH`, `MatD`, `numeric`, `factor`, and `character`.
#' @param row.names Optional character vector of row names. If not provided, default names like "R1", "R2", ... will be used.
#'
#' @return An object of class \code{"DDtable"}.
#' @examples
#' # Create basic MatH and MatD objects (with default empty distributions)
#' mh <- new("MatH", nrows = 3, ncols = 1)
#' md <- new("MatD", nrows = 3, ncols = 1)
#'
#' # Create some standard columns
#' age <- c(25, 30, 35)
#' sex <- factor(c("M", "F", "M"))
#'
#' # Construct a DDtable with mixed column types
#' dd <- DDtable(hist = mh, disc = md, age = age, sex = sex)
#'
#' @seealso \linkS4class{DDtable}
#' @export
DDtable <- function(..., row.names = NULL) {
  cols <- list(...)
  n <- length(cols)
  if (n == 0) stop("No columns provided")

  col.types <- sapply(cols, function(c) {
    if (is(c, "MatH")) return("distributionH")
    if (is(c, "MatD")) return("distributionD")
    if (is.factor(c)) return("factor")
    if (is.numeric(c)) return("numeric")
    if (is.character(c)) return("character")
    class(c)[1]
  })

  lengths <- sapply(cols, function(c) {
    if (is(c, "MatH") || is(c, "MatD")) {
      nrow(c@M)
    } else {
      length(c)
    }
  })

  if (length(unique(lengths)) != 1) {
    stop("All columns must have the same number of rows")
  }

  nrows <- unique(lengths)
  if (is.null(row.names)) {
    row.names <- paste0("R", seq_len(nrows))
  }

  col.names <- names(cols)
  if (is.null(col.names)) {
    col.names <- paste0("V", seq_len(n))
  }

  new("DDtable",
      data = cols,
      col.types = col.types,
      row.names = row.names,
      col.names = col.names)
}
