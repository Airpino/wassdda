# Classes for Discrete distributions

# A class for histogram-valued data ---------------------------------------
#' Class distributionD.
#' @name distributionD-class
#' @rdname distributionD-class
#' @exportClass distributionD
setClass(
  Class = "distributionD",
  representation = representation(
    x = "numeric",
    p = "numeric",
    m = "numeric",
    s = "numeric"
  ),
  # distributionH=
  validity = function(object) {
    if (length(object@x) <= 0) {
      x <- NA
    }
    else {
      if (length(object@p)) {
        nv <- length(object@x)
        p <- c(1:nv) / nv
      }
      if (length(object@x) != length(object@p)) {
        stop("the x and p vectors must be of the same length")
      }
      if (length(object@x) == 1) {
        object@x <- c(object@x, object@x)
        object@p <- c(1)
      }
      nv <- length(object@x)
      if ((length(object@x) > 1) &&
          (min(diff(object@x)) < 0)) {
        print(object@x)
        print(object@p)
        return("the x must be in not descending order")
      }
      if ((length(object@x) > 1) && (min(diff(object@p)) < 0) &&
          (object@p[1] > 1e-14) && (object@p[nv] < 1 - 1e-14)) {
        print(object@x)
        print(object@p)
        return("the p must be in not descending order from 0 to 1")
      }

      return(TRUE)
    }
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
          definition = function(.Object, x = numeric(0), p = numeric(0), m = numeric(0), s = numeric(0)) {
            .Object@x <- x
            .Object@p <- p
            if (length(x) > 0) {
              validObject(.Object)
              if (length(m) == 0) .Object@m <- meanD(.Object) else .Object@m <- m
              if (length(s) == 0) .Object@s <- stdD(.Object) else .Object@s <- s
            }

            return(.Object)
          }
)

# A class of a matrix of discrete DD ----------------------------
#' Class MatH.
#' @name MatH-class
#' @rdname MatH-class
#' @exportClass MatH
#' @docType class
#' @description   Class \code{MatH} defines a matrix of \code{distributionH} objects
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
  representation = representation(M = "matrix"),
)
#' Constructor method for MatD class
#' @name MatD
#' @rdname MatD-class
#' @aliases initialize,MatD-method
setMethod("initialize", "MatD",
          definition = function(.Object,
                                nrows = 1, ncols = 1, ListOfDist = NULL,
                                names.rows = NULL, names.cols = NULL,
                                by.row = FALSE) {
            tt <- list(new("distributionD"))
            .Object@M <- matrix(tt, nrows, ncols)

            if (length(ListOfDist) > 0) {
              nOBJ <- length(ListOfDist)
              if (by.row) {
                count <- 0
                for (i in 1:nrows) {
                  for (j in 1:ncols) {
                    count <- count + 1
                    if (count == nOBJ) count <- 1
                    .Object@M[i, j][[1]] <- ListOfDist[[count]]
                  }
                }
              }
              else {
                count <- 0
                for (j in 1:ncols) {
                  for (i in 1:nrows) {
                    count <- count + 1
                    if (count > nOBJ) count <- 1
                    .Object@M[i, j][[1]] <- ListOfDist[[count]]
                  }
                }
              }
            }

            if (length(names.rows) > 0) {
              count <- 0
              rnames <- vector("list", nrows)

              for (i in 1:nrows) {
                count <- count + 1

                if (count > length(names.rows)) {
                  rnames[[count]] <- paste("I", count, sep = "")
                }
                else {
                  rnames[[count]] <- names.rows[[count]]
                }
              }
              rownames(.Object@M) <- rnames
            }
            else {
              rownames(.Object@M) <- paste("I", 1:nrows, sep = "")
            }

            if (length(names.cols) > 0) {
              count <- 0
              cnames <- vector("list", ncols)
              for (i in 1:ncols) {
                count <- count + 1
                if (count > length(names.cols)) {
                  cnames[[count]] <- paste("X", count, sep = "")
                }
                else {
                  cnames[[count]] <- names.cols[[count]]
                }
                colnames(.Object@M) <- cnames
              }
            } else {
              colnames(.Object@M) <- paste("X", 1:ncols, sep = "")
            }

            return(.Object)
          }
)
