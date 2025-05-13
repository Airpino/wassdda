# Classes for Discrete distributions

# A class for histogram-valued data ---------------------------------------
#' Class distributionD.
#' @name distributionD-class
#' @rdname distributionD-class
#' @exportClass distributionD
setClass(
  "distributionD",
  slots = c(
    x = "numeric",
    p = "numeric",
    m = "numeric",
    s = "numeric"
  ),
  validity = function(object) {
    if (length(object@x) != length(object@p)) {
      return("Slots 'x' e 'p' devono avere la stessa lunghezza")
    }
    if (length(object@x) < 1) {
      return("Almeno un valore è richiesto in 'x'")
    }
    if (is.unsorted(object@x, strictly = FALSE)) {
      return("Slot 'x' deve essere ordinato in modo non decrescente")
    }
    if (is.unsorted(object@p, strictly = FALSE)) {
      return("Slot 'p' deve essere non decrescente")
    }
    if (any(object@p < 0 | object@p > 1)) {
      return("Tutti i valori di 'p' devono essere compresi tra 0 e 1")
    }
    if (abs(tail(object@p, 1) - 1) > 1e-8) {
      return("L'ultimo valore di 'p' deve essere uguale a 1")
    }

    return(TRUE)
  }
)
#' Constructor method of distributionH class
#'
#' Class \code{distributionD} defines a histogram object
#'
#' @name distributionD
#' @rdname distributionD-class
#' @aliases initialize,distributionD-method
#' @description Class \code{"distribution"} defines an discrete distributional object
#' The class describes a discrete DDA by means of its cumulative distribution
#' function. The methods are develoved accordingly to the L2 Wasserstein
#' distance between distributions.
#'
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("distributionD", x, p, m, s)}.
#' @param .Object the type ("distributionD")
#' @param x a numeric vector. it is the domain of the distribution .
#' @param p a numeric vector (of the same lenght of x). It is the cumulative distribution function CDF.
#' @param m (optional) a numeric value. Is the mean of the histogram.
#' @param s (optional) a numeric positive value. It is the standard deviation of a histogram.
#' @author Antonio Irpino
#' @references Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI
#' 10.1007/s11634-014-0176-4
#' @keywords classes
#' @examples #---- initialize a distributionH object mydist
#' # from a simple histogram
#' # ----------------------------
#' # | x    |  Prob  | cdf   |
#' # ----------------------------
#' # | 1   |  0.4   | 0.4   |
#' # | 2   |  0.6   | 1.0   |
#' # ----------------------------
#' # | Tot.    |  1.0   | -     |
#' # ----------------------------
#' mydist <- new("distributionD", c(1, 2, 3), c(0.1, 0.5, 1))
#' str(mydist)
#' # OUTPUT
#' # Formal class 'distributionD' [package "wassdda"] with 4 slots
#' #   ..@@ x: num [1:2] 1 2 the values
#' #   ..@@ p: num [1:2] 0.4 1 the cdf
#' #   ..@@ m: num  the mean
#' #   ..@@ s: num  the standard deviation
#' @seealso \code{\link{meanD}} computes the mean. \code{\link{stdD}} computes the standard deviation.
# showClass("distributionD")
setMethod("initialize", "distributionD",
          definition = function(.Object, x = numeric(0), p = numeric(0), m = NA_real_, s = NA_real_) {
            .Object@x <- x
            .Object@p <- p

            if (length(x) > 0 && length(p) > 0) {
              validObject(.Object)

              if (is.na(m)) .Object@m <- meanD(.Object) else .Object@m <- m
              if (is.na(s)) .Object@s <- stdD(.Object) else .Object@s <- s
            } else {
              .Object@m <- m
              .Object@s <- s
            }

            return(.Object)
          }
)

# A class of a matrix of discrete DD ----------------------------
#' Class MatD.
#' @name MatD-class
#' @rdname MatD-class
#' @exportClass MatD
#' @docType class
#' @description   Class \code{MatD} defines a matrix of \code{distributionD} objects
#' @author Antonio Irpino
#' @references Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI
#' 10.1007/s11634-014-0176-4
#' @keywords classes
#' @examples
#'
#' ## ---- create a list of six distributionD objects
#' ListOfDist <- vector("list", 6)
#' ListOfDist[[1]] <- distributionD(c(1, 2, 3), c(0.2, 0.4, 1))
#' ListOfDist[[2]] <- distributionD(c(7, 8, 10, 15), c(0.1, 0.2, 0.7, 1))
#' ListOfDist[[3]] <- distributionD(c(9, 11, 20), c(0.3, 0.5, 1))
#' ListOfDist[[4]] <- distributionD(c(2, 5, 8), c(0.1, 0.3, 1))
#' ListOfDist[[5]] <- distributionD(c(8, 10, 15), c(0.4, 0.75, 1))
#' ListOfDist[[6]] <- distributionD(c(20, 22, 24), c(0.05, 0.12, 1))
#'
#' ## create a MatH object filling it by columns
#' MyMAT <- new("MatD",
#'   nrows = 3, ncols = 2, ListOfDist = ListOfDist,
#'   names.rows = c("I1", "I2", "I3"), names.cols = c("Var1", "Var2"), by.row = FALSE
#' )
#'
#' showClass("MatD")
#' @param .Object the object type "MatD"
#' @param ListOfDist a vector or a list of \code{distributionD} objects
#' @param names.rows a vector or list of strings with thenames of the rows
#' @param names.cols a vector or list of strings with thenames of the columns (variables)

setClass(
  Class = "MatD",
  slots = c(M = "matrix")
)

setMethod("initialize", "MatD",
          function(.Object,
                   nrows = 1, ncols = 1, ListOfDist = NULL,
                   names.rows = NULL, names.cols = NULL,
                   by.row = FALSE) {

            # Creazione matrice di oggetti distributionD
            mat <- matrix(vector("list", nrows * ncols), nrows, ncols)
            for (i in seq_len(nrows)) {
              for (j in seq_len(ncols)) {
                mat[[i, j]] <- new("distributionD")
              }
            }

            .Object@M <- mat

            # Inserimento delle distribuzioni dalla lista
            if (!is.null(ListOfDist)) {
              if (!all(sapply(ListOfDist, function(x) is(x, "distributionD")))) {
                stop("All elements in ListOfDist must be 'distributionD' objects")
              }

              nOBJ <- length(ListOfDist)
              count <- 1
              idx <- if (by.row) expand.grid(i = seq_len(nrows), j = seq_len(ncols), KEEP.OUT.ATTRS = FALSE)
              else expand.grid(j = seq_len(ncols), i = seq_len(nrows), KEEP.OUT.ATTRS = FALSE)[, c("i", "j")]

              for (k in seq_len(nrows * ncols)) {
                i <- idx[k, 1]
                j <- idx[k, 2]
                .Object@M[[i, j]] <- ListOfDist[[count]]
                count <- ifelse(count == nOBJ, 1, count + 1)
              }
            }

            # Nomi righe e colonne
            rownames(.Object@M) <- if (!is.null(names.rows)) names.rows else paste0("I", seq_len(nrows))
            colnames(.Object@M) <- if (!is.null(names.cols)) names.cols else paste0("X", seq_len(ncols))

            validObject(.Object)
            return(.Object)
          }
)


