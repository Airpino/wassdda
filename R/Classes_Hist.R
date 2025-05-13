# A class for histogram-valued data ---------------------------------------
#' Class distributionH.
#' @name distributionH-class
#' @rdname distributionH-class
#' @exportClass distributionH
#' @slot x Numeric vector: domain (quantiles or bin edges).
#' @slot p Numeric vector: cumulative distribution function (same length as x).
#' @slot m Numeric: mean of the distribution.
#' @slot s Numeric: standard deviation of the distribution.
setClass(
  "distributionH",
  slots = c(
    x = "numeric",
    p = "numeric",
    m = "numeric",
    s = "numeric"
  ),
  validity = function(object) {
    # Length check
    if (length(object@x) != length(object@p)) {
      return("Slots 'x' and 'p' must have the same length")
    }
    if (length(object@x) <= 1) {
      return("At least two points are required in 'x'")
    }

    # Monotonicity and range checks
    if (is.unsorted(object@x, strictly = FALSE)) {
      return("Slot 'x' must be sorted in non-decreasing order")
    }
    if (is.unsorted(object@p, strictly = FALSE)) {
      return("Slot 'p' must be non-decreasing")
    }
    if (abs(object@p[1]) > 1e-8 || abs(tail(object@p, 1) - 1) > 1e-8) {
      return("Slot 'p' must start at 0 and end at 1")
    }

    TRUE
  }
)


#' Constructor method of distributionH class
#'
#' Class \code{distributionH} defines a histogram object
#'
#' @name distributionH
#' @rdname distributionH-class
#' @aliases initialize,distributionH-method
#' @description Class \code{"distributionH"} desfines an histogram object
#' The class describes a histogram by means of its cumulative distribution
#' function. The methods are develoved accordingly to the L2 Wasserstein
#' distance between distributions.
#'
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("distributionH", x, p, m, s)}.
#' @param .Object the type ("distributionH")
#' @param x a numeric vector. it is the domain of the distribution (i.e. the
#' extremes of bins).
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
#' # | Bins    |  Prob  | cdf   |
#' # ----------------------------
#' # | [1,2)   |  0.4   | 0.4   |
#' # | [2,3]   |  0.6   | 1.0   |
#' # ----------------------------
#' # | Tot.    |  1.0   | -     |
#' # ----------------------------
#' mydist <- new("distributionH", c(1, 2, 3), c(0, 0.4, 1))
#' str(mydist)
#' # OUTPUT
#' # Formal class 'distributionH' [package "HistDAWass"] with 4 slots
#' #   ..@@ x: num [1:3] 1 2 3 the quantiles
#' #   ..@@ p: num [1:3] 0 0.4 1 the cdf
#' #   ..@@ m: num 2.1 the mean
#' #   ..@@ s: num 0.569 the standard deviation
#' @seealso \code{\link{meanH}} computes the mean. \code{\link{stdH}} computes the standard deviation.
# showClass("distributionH")
setMethod("initialize", "distributionH",
          function(.Object, x = numeric(0), p = numeric(0), m = numeric(0), s = numeric(0)) {
            .Object@x <- x
            .Object@p <- p

            if (length(x) > 1 && length(p) == length(x)) {
              # Calcola m e s solo se non forniti
              if (length(m) == 0) .Object@m <- meanH(.Object) else .Object@m <- m
              if (length(s) == 0) .Object@s <- stdH(.Object) else .Object@s <- s
            } else {
              .Object@m <- m
              .Object@s <- s
            }

            validObject(.Object)
            return(.Object)
          }
)

# A class of a matrix of histogram-valued data ----------------------------
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
#' ## ---- create a list of six distributionH objects
#' ListOfDist <- vector("list", 6)
#' ListOfDist[[1]] <- distributionH(c(1, 2, 3), c(0, 0.4, 1))
#' ListOfDist[[2]] <- distributionH(c(7, 8, 10, 15), c(0, 0.2, 0.7, 1))
#' ListOfDist[[3]] <- distributionH(c(9, 11, 20), c(0, 0.5, 1))
#' ListOfDist[[4]] <- distributionH(c(2, 5, 8), c(0, 0.3, 1))
#' ListOfDist[[5]] <- distributionH(c(8, 10, 15), c(0, 0.75, 1))
#' ListOfDist[[6]] <- distributionH(c(20, 22, 24), c(0, 0.12, 1))
#'
#' ## create a MatH object filling it by columns
#' MyMAT <- new("MatH",
#'   nrows = 3, ncols = 2, ListOfDist = ListOfDist,
#'   names.rows = c("I1", "I2", "I3"), names.cols = c("Var1", "Var2"), by.row = FALSE
#' )
#'
#' showClass("MatH")
#' @param .Object the object type "MatH"
#' @param ListOfDist a vector or a list of \code{distributionH} objects
#' @param names.rows a vector or list of strings with thenames of the rows
#' @param names.cols a vector or list of strings with thenames of the columns (variables)

setMethod("initialize", "MatH",
          function(.Object, nrows = 1, ncols = 1, ListOfDist = NULL, names.rows = NULL, names.cols = NULL, by.row = FALSE) {

            .Object@M <- matrix(
              replicate(nrows * ncols, new("distributionH"), simplify = FALSE),
              nrow = nrows, ncol = ncols
            )

            if (length(ListOfDist) > 0) {
              count <- 0
              for (i in 1:nrows) {
                for (j in 1:ncols) {
                  count <- count + 1
                  if (count > length(ListOfDist)) count <- 1
                  if (by.row) {
                    .Object@M[i, j] <- ListOfDist[[count]]
                  } else {
                    .Object@M[j, i] <- ListOfDist[[count]]
                  }
                }
              }
            }

            rownames(.Object@M) <- if (!is.null(names.rows)) {
              rep_len(names.rows, nrows)
            } else {
              paste0("I", seq_len(nrows))
            }

            colnames(.Object@M) <- if (!is.null(names.cols)) {
              rep_len(names.cols, ncols)
            } else {
              paste0("X", seq_len(ncols))
            }

            return(.Object)
          }
)

## Classes for Histogram Time Series ------
## A single distribution with a time stamp ------
#' Class TdistributionH
#'
#' Class \code{TdistributionH} defines a histogram with a time (point or period)
#'
#' @name TdistributionH-class
#' @rdname TdistributionH-class
#' @exportClass TdistributionH
setClass(
  Class = "TdistributionH",
  contains = "distributionH",
  slots = c(
    tstamp = "numeric",
    period = "numeric"
  ),
  validity = function(object) {
    if (length(object@tstamp) != 1) {
      return("The 'tstamp' slot must contain exactly one numeric value")
    }
    if (length(object@period) > 2) {
      return("The 'period' slot may contain at most two values")
    }
    if (length(object@period) == 2 && object@period[1] > object@period[2]) {
      return("The 'period' start must not be greater than its end")
    }
    return(TRUE)
  }
)
## Initialize ------
#' Constructor method of TdistributionH Class
#'
#' @name TdistributionH
#' @rdname TdistributionH-class
#' @aliases initialize,TdistributionH-method
#' @param .Object the type of object ("TdistributionH") a \code{"distributionH"} object with a time reference
#' @param tstamp a numeric value related to  a timestamp
#' @param period a list of two values, the starting time and the ending time (alternative to tstamp if the
#' distribution is observed along a period and not on a timestamp)
#' @param x a vector of increasing values, the domain of the distribution (the same of \code{distributionH} object)
#' @param p a vector of increasing values from 0 to 1,
#' the CDF of the distribution (the same of \code{distributionH} object)
#' @param m a number, the mean of the distribution (the same of \code{distributionH} object)
#' @param s a positive number, the standard deviation of the distribution (the same of \code{distributionH} object)
setMethod("initialize", "TdistributionH",
          function(.Object,
                   tstamp = numeric(1),
                   period = numeric(2),  # preferibile rispetto alla list
                   x = numeric(0),
                   p = numeric(0),
                   m = numeric(0),
                   s = numeric(0)) {

            .Object@x <- x
            .Object@p <- p
            .Object@tstamp <- tstamp
            .Object@period <- period

            # Calcola m e s solo se x è valorizzato
            if (length(x) > 0) {
              .Object@m <- if (length(m) == 0 || is.na(m)) meanH(.Object) else m
              .Object@s <- if (length(s) == 0 || is.na(s)) stdH(.Object) else s
            } else {
              .Object@m <- m
              .Object@s <- s
            }

            validObject(.Object)  # Ultima cosa da fare
            return(.Object)
          }
)
# Coerce a TdistributionH into a distributionH
setAs(
  from = "TdistributionH", to = "distributionH",
  function(from, to) {
    to <- new("distributionH",
      x = from@x,
      p = from@p,
      m = from@m,
      s = from@s
    )
    return(to)
  }
)

## A mulvariate MatH with a time stamp  ----
#' Class TMatH
#'
#' Class \code{TMatH} defines a matrix of histograms, a \code{TMatH} object, with a time (a timepoint or a time window).
#'
#' @name TMatH-class
#' @rdname TMatH-class
#' @exportClass TMatH
setClass(
  Class = "TMatH",
  representation = representation(
    tstamp = "numeric",
    period = "numeric"  # Consigliato: vettore numerico di lunghezza 2
  ),
  contains = "MatH",
  validity = function(object) {
    if (length(object@tstamp) > 1) {
      return("Slot 'tstamp' must contain a single numeric value")
    }
    if (length(object@period) > 2) {
      return("Slot 'period' must contain at most two numeric values (start, end)")
    }
    if ((length(object@period) == 2) && (object@period[1] > object@period[2])) {
      return("Slot 'period': start time is greater than end time")
    }
    TRUE
  }
)
## Initialize TMatH ------
#' Constructor method of TdistributionH Class
#'
#' @name TMatH
#' @rdname TMatH-class
#' @aliases initialize,TMatH-method
#' @param .Object the type of object ("TMatH")
#' @param tstamp a vector of time stamps, numeric.
#' @param period a list of pairs with a vectorof starting time and a vector ofending time.
#' This parameter is used alternatively to \code{tstamp} if the distributions are related to time periods
#' instead of timestamps
#' @param mat a \code{MatH} object
setMethod("initialize", "TMatH",
          definition = function(.Object,
                                tstamp = numeric(0),
                                period = numeric(0),  # usa numeric invece di list
                                mat = new("MatH")) {
            .Object@M <- mat@M  # assegna correttamente la matrice contenuta
            .Object@tstamp <- tstamp
            .Object@period <- period
            validObject(.Object)
            return(.Object)
          }
)

## A Histogram Time Series HTS-----
#' Class HTS
#'
#' Class \code{HTS} defines a histogram time series, i.e. a set of histograms observed along time
#'
#' @name HTS-class
#' @rdname HTS-class
#' @exportClass HTS
setClass(
  Class = "HTS",
  representation = representation(data = "list"),
  validity = function(object) {
    if (length(object@data) > 0) {
      # Verifica tipo del primo elemento
      type <- is(object@data[[1]])[1]
      if (!(type %in% c("TdistributionH", "TMatH"))) {
        return("Slot 'data' must contain only TdistributionH or TMatH objects")
      }

      for (i in 2:length(object@data)) {
        if (is(object@data[[i]])[1] != type) {
          return("All elements in 'data' must be of the same class")
        }

        if (type == "TMatH") {
          current <- object@data[[i]]
          previous <- object@data[[i - 1]]

          if (nrow(current@M) != nrow(previous@M) ||
              ncol(current@M) != ncol(previous@M)) {
            return("All TMatH objects must have the same dimensions")
          }
        }
      }
    }
    TRUE
  }
)
## Initialize HTS ------
#' Constructor method of HTS Class (Histogram Time Series)
#'
#' @name HTS
#' @rdname HTS-class
#' @aliases initialize,HTS-method
#' @param .Object the object type ("HTS") a histogram time series
#' @param epocs the number of histograms (one for each timepoint or period)
#' @param ListOfTimedElements a vector of \code{TdistributionH} objects
setMethod("initialize", "HTS",
          definition = function(.Object, epocs = 1, ListOfTimedElements = list(new("TdistributionH"))) {
            if (!is(ListOfTimedElements, "list")) {
              stop("'ListOfTimedElements' must be a list")
            }

            if (length(ListOfTimedElements) == 0) {
              stop("ListOfTimedElements must contain at least one element")
            }

            replicated_list <- vector("list", epocs)
            for (i in seq_len(epocs)) {
              index <- ((i - 1) %% length(ListOfTimedElements)) + 1
              replicated_list[[i]] <- ListOfTimedElements[[index]]
            }

            .Object@data <- replicated_list
            validObject(.Object)
            return(.Object)
          }
)
