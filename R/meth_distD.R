# Constructor
#' Wrapper function distributionD
#'
#' A discrete DD object can be created also with the function \code{distributionD(...)}, the costructor function for creating an object containing the description of
#' a discrete DD.
#'
#' @name distributionD
#' @rdname distributionD-class
#' @export
#' @return A \code{distributionD} object
#' @examples
#' # or using
#' mydist <- distributionD(x = c(1, 2, 3), p = c(0.1, 0.4, 1))
#' @import methods
distributionD <- function(x = numeric(0), p = numeric(0)) {
  object <- new("distributionD", x = x, p = p)
  return(object)
}

# Get methods -----

#' Method \code{get.m}: the mean of a distribution
#' @name get.m
#' @rdname get.m-methods
#' @exportMethod get.m
setGeneric("get.m", function(object) standardGeneric("get.m"))
#' @rdname get.m-methods
#' @aliases get.m,distributionD-method
#' @description This functon return the mean of a \code{distributionD} object.
#' @param object a \code{distributionD} object
#' @return A numeric value
#' @examples
#' D <- distributionD(x = c(1, 2, 3, 4), p = c(0.1, 0.2, 0.6, 1))
#' get.m(D) # returns the mean of D
setMethod(
  "get.m", "distributionD",
  function(object) {
    if (!is.null(object)) {
      return(object@m)
    } else {
      return(NA)
    }
  }
)

#' Method \code{get.s}: the standard deviation of a distribution
#' @name get.s
#' @rdname get.s-methods
#' @exportMethod get.s
setGeneric("get.s", function(object) standardGeneric("get.s"))

#' @rdname get.s-methods
#' @aliases get.s,distributionD-method
#' @description This functon return the standard deviation of a \code{distributionD} object.
#' @param object a \code{distributionD} object.
#' @return A numeric positive value, the standard deviation.
#' @examples
#' D <- distributionD(x = c(1, 2, 3, 4), p = c(0.1, 0.2, 0.6, 1))
#' get.s(D) # returns the standard deviation of D
setMethod(
  "get.s", "distributionD",
  function(object) {
    if (!is.null(object)) {
      return(object@s)
    } else {
      return(NA)
    }
  }
)
#' Method \code{get.distr}: show the distribution
#' @name get.distr
#' @rdname get.distr-methods
#' @exportMethod get.distr
setGeneric("get.distr", function(object) standardGeneric("get.distr"))

#' @rdname get.distr-methods
#' @aliases get.distr,distributionD-method
#' @description This functon return the cumulative distribution function of a \code{distributionD} object.
#' @param object a \code{distributionD} object.
#' @return A data frame: the first column contains the domain the second the CDF values.
#' @examples
#' D <- distributionD(x = c(1, 2, 3, 4), p = c(0, 0.2, 0.6, 1))
#' get.distr(D) # a data.frame describing the CDF of D
setMethod(
  "get.distr", "distributionD",
  function(object) {
    if (!is.null(object)) {
      MAT <- cbind(object@x, object@p)
      colnames(MAT) <- c("x", "p")
      return(MAT = as.data.frame(MAT))
    } else {
      return(NA)
    }
  }
)
#' Method \code{get.bars}: show the distribution with bins
#' @name get.bars
#' @rdname get.bars-methods
#' @exportMethod get.bars
setGeneric("get.bars", function(object) standardGeneric("get.bars"))

# Get the distribution of a discrete DD
#' @rdname get.bars-methods
#' @aliases get.bars,distributionD-method
#' @description This functon return a data.frame describing the mass frequency of a \code{distributionD} object.
#' @param object a \code{distributionD} object.
#' @return A matrix: the two columns contains the values of DD and the the probablity (or the relative frequency) of the bin.
#' @examples
#' D <- distributionD(x = c(1, 2, 3, 4), p = c(0.1, 0.2, 0.6, 1))
#' get.bars(D) # returns the discrete DD representation of D by a data.frame
setMethod(
  "get.bars", "distributionD",
  function(object) {
    if (!is.null(object)) {
      M <- get.distr(object)
      MAT <- cbind(
        M$x,
        c(M$p[1],diff(M$p))
      )
      colnames(MAT) <- c("x", "p")
      return(MAT = as.data.frame(MAT))
    } else {
      return(NA)
    }
  }
)

# Basic statistics of distributions --------
#' Method \code{meanD}: computes the mean of a distribution
#' @name meanD
#' @rdname meanD-methods
#' @exportMethod meanD
setGeneric("meanD", function(object) standardGeneric("meanD"))
#' Method \code{stdD}: computes the standard deviation of a distribution
#' @name stdD
#' @rdname stdD-methods
#' @exportMethod stdD
setGeneric("stdD", function(object) standardGeneric("stdD"))
#' Method \code{skewD}: computes the skewness of a distribution
#' @name skewD
#' @rdname skewD-methods
#' @exportMethod skewD
setGeneric("skewD", function(object) standardGeneric("skewD"))
#' Method \code{kurtD}: computes the kurthosis of a distribution
#' @name kurtD
#' @rdname kurtD-methods
#' @exportMethod kurtD
setGeneric("kurtD", function(object) standardGeneric("kurtD"))

#' @rdname meanD-methods
#' @aliases meanD,distributionD-method
#' @description Mean of a discrete  DD (First moment of the distribution)
#' @param object a \code{distributionD} object
#' @return the mean of the distribution
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#'
#' ## ---- A mydist distribution ----
#' mydist <- distributionD(x = c(1, 2, 3, 10), p = c(0.05, 0.1, 0.5, 1))
#' ## ---- Compute the mean of mydist ----
#' meanD(mydist)
setMethod(
  "meanD", "distributionD", # Mean of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      m <- sum(object@x*c(object@p[1],diff(object@p)))
      return(m)
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
#' @rdname stdD-methods
#' @aliases stdD,distributionD-method
#' @description Standard deviation of a discrete DD (i.e., the square root of the centered
#' second moment)
#'
#' @param object a \code{distributionD} object
#' @return A value for the standard deviation
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#'
#' ## ---- A mydist distribution ----
#' mydist <- distributionD(x = c(1, 2, 3, 10), p = c(0.05, 0.1, 0.5, 1))
#' ## ---- Compute the standard deviation of mydist ----
#' stdD(mydist)
setMethod(
  "stdD", "distributionD", # Std of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      w<-c(object@p[1],diff(object@p))
      s <- sqrt(sum(object@x^2*w-(object@x*w)^2))
      #   #resu=crwtransform(object)
      #   c=0.5*(object@x[2:length(object@x)]+object@x[1:(length(object@x)-1)])
      #   r=0.5*(object@x[2:length(object@x)]-object@x[1:(length(object@x)-1)]) # resu[[2]]
      #   w=(object@p[2:length(object@p)]-object@p[1:(length(object@p)-1)])#resu[[3]]
      #   std=sqrt(abs(sum(w*c^2+1/3*w*r^2)-(sum(w*c))^2))
      # #  std=as.numeric(std)
      return(s)
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
#' @rdname skewD-methods
#' @aliases skewD,distributionD-method
#' @description Skewness of a discrete DD (using the third standardized moment)
#'
#' @param object a \code{distributionD} object
#' @return A value for the skewness index
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#'
#' ## ---- A mydist distribution ----
#' mydist <- distributionD(x = c(1, 2, 3, 10), p = c(0.05, 0.1, 0.5, 1))
#' ## ---- Compute the skewness of mydist ----
#' skewD(mydist)
setMethod(
  "skewD", "distributionD", # skewness of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      w<-c(object@p[1],diff(object@p))
      m<-sum(object@x*w)
      m3<-sum((object@x-m)^3*w)
      s<-sqrt(sum(object@x^2*w-(object@x*w)^2))
      sk <- m3/(s^3)
      return(sk)
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)
#' @rdname kurtD-methods
#' @aliases kurtD,distributionD-method
#' @description Kurtosis of a discrete DD (using the fourth standardized moment)
#' @param object a \code{distributionD} object
#' @return A value for the kurtosis index, 3 is the kurtosis of a Gaussian
#' distribution
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#'
#' ## ---- A mydist distribution ----
#' mydist <- distributionD(x = c(1, 2, 3, 10), p = c(0.05, 0.1, 0.5, 1))
#' ## ---- Compute the kurtosis of mydist ----
#' kurtD(mydist)
setMethod(
  "kurtD", "distributionD", # Kurtosis of a distr --------
  function(object) {
    if (!is.null(object@x) || !is.null(object@p)) {
      w<-c(object@p[1],diff(object@p))
      m<-sum(object@x*w)
      m4<-sum((object@x-m)^4*w)
      s<-sqrt(sum(object@x^2*w-(object@x*w)^2))
      ku <- m4/(s^4)
      return(ku)
    } else {
      stop("Something wrong, null domain or cdf")
    }
  }
)

##  cancellare

# Overloading of the sum of two distribution according to the L2 w --------
#' Method +
#' @name +
#' @aliases +,distributionD,distributionD-method
#' @description the sum of two distribution according to the L2 Wassserstein
#' @param e1 a \code{distributionD} object or a number
#' @param e2 a \code{distributionD} object or a number
#' @return a \code{distributionD} object
#' @export
#' @docType methods
#' @rdname plus-methods
#'
setMethod(
  "+",
  signature(e1 = "distributionD", e2 = "distributionD"),
  function(e1, e2) {
    if (!identical(e1@p, e2@p)) {
      tmp <- register(e1, e2)
      x <- callGeneric(tmp[[1]]@x, tmp[[2]]@x)
      e1@p <- tmp[[1]]@p
    } else {
      x <- callGeneric(e1@x, e2@x)
    }

    e1@x <- x

    e1@m <- e1@m + e2@m
    e1@s <- stdD(e1)
    return(e1)
    #            OBJ_NEW=new("distributionD",x,tmp[[1]]@p,(tmp[[1]]@m+tmp[[2]]@m))
  }
)


#' Method +
#' @name +
#' @aliases +,numeric,distributionD-method
#' @description the sum of a number and a distribution according to the L2 Wasserstein
#' @export
#' @docType methods
#' @rdname plus-methods

setMethod(
  "+",
  signature(e1 = "numeric", e2 = "distributionD"),
  function(e1, e2) {
    x <- callGeneric(rep(e1, length(e2@x)), e2@x)
    OBJ_NEW <- new("distributionD", x, e2@p, (e1 + e2@m), e2@s)
  }
)
#' Method +
#' @name +
#' @aliases +,distributionD,numeric-method
#' @description the sum of adistribution and a number according to the L2 Wassserstein
#' @export
#' @docType methods
#' @rdname plus-methods
setMethod(
  "+",
  signature(e1 = "distributionD", e2 = "numeric"),
  function(e1, e2) {
    x <- callGeneric(rep(e2, length(e1@x)), e1@x)
    OBJ_NEW <- new("distributionD", x, e1@p, (e2 + e1@m), e1@s)
  }
)
# Overloading of the difference of two distributions according to the L2 wasserstein --------
#' Method -
#' @name minus
#' @aliases -,distributionD,distributionD-method
#' @description the difference of two distribution according to the L2 Wasserstein
#' @export
setMethod(
  "-",
  signature(e1 = "distributionD", e2 = "distributionD"),
  function(e1, e2) {
    if (!identical(e1@p, e2@p)) {
      tmp <- register(e1, e2)
      x <- callGeneric(tmp[[1]]@x, tmp[[2]]@x)
      e1@p <- tmp[[1]]@p
    } else {
      x <- callGeneric(e1@x, e2@x)
    }
    x <- callGeneric(tmp[[1]]@x, tmp[[2]]@x)
    OBJ_NEW <- new("distributionD", x, e1@p)
  }
)
#' Method -
#' @name minus
#' @aliases -,numeric,distributionD-method
#' @param e1 a \code{distributionD} object or a number
#' @param e2 a \code{distributionD} object or a number
#' @description the difference of a number and a distribution according to the L2 Wasserstein
setMethod(
  "-",
  signature(e1 = "numeric", e2 = "distributionD"),
  function(e1, e2) {
    x <- callGeneric(rep(e1, length(e2@x)), e2@x)
    OBJ_NEW <- new("distributionD", x, e2@p)
  }
)
#' Method -
#' @name minus
#' @aliases -,distributionD,numeric-method
#' @description the difference of a distribution and a number according to the L2 Wasserstein
#' @note it may not works properly if the difference is not a distribution
setMethod(
  "-",
  signature(e1 = "distributionD", e2 = "numeric"),
  function(e1, e2) {
    x <- callGeneric(e1@x, rep(e2, length(e1@x)))
    OBJ_NEW <- new("distributionD", x, e1@p, (e1@m - e2), e1@s)
  }
)

# Overloading of the product of a number by a distribution according to the L2 w --------
#' Method *
#' @name *-methods
#' @aliases *,distributionD,distributionD-method
#' @description the product of a number and a distribution according to the L2 Wasserstein
#' @param e1 a \code{distributionD} object or a number
#' @param e2 a \code{distributionD} object or a number
#' @export
setMethod(
  "*",
  signature(e1 = "distributionD", e2 = "distributionD"),
  function(e1, e2) {
    stop("please use dotpW function product between distributions")
  }
)
#' Method *
#' @name *-methods
#' @aliases *,numeric,distributionD-method
#' @description the product of a number and a distribution according to the L2 Wasserstein
setMethod(
  "*",
  signature(e1 = "numeric", e2 = "distributionD"),
  function(e1, e2) {
    x <- callGeneric(rep(e1, length(e2@x)), e2@x)

    e2@x <- x
    e2@p <- e2@p
    e2@m <- e1 * e2@m
    e2@s <- abs(e1) * e2@s
    return(e2)
  }
)
#' Method *
#' @name *-methods
#' @aliases *,distributionD,numeric-method
#' @description the product of a number and a distribution according to the L2 Wasserstein
setMethod(
  "*",
  signature(e1 = "distributionD", e2 = "numeric"),
  function(e1, e2) {
    x <- callGeneric(rep(e2, length(e1@x)), e1@x)
    e1@x <- x
    e1@p <- e1@p
    e1@m <- e2 * e1@m
    e1@s <- abs(e2) * e1@s
    return(e1)
  }
)

# Utilities for single or couples of distributionD --------------------------------
#' Method \code{checkEmptyBins}
#' @name checkEmptyBins
#' @rdname checkEmptyBins-methods
#' @exportMethod checkEmptyBins
setGeneric("checkEmptyBins", function(object) standardGeneric("checkEmptyBins"))
#' Method \code{compQ}
#' @name compQ
#' @rdname compQ-methods
#' @exportMethod compQ
setGeneric("compQ", function(object, p) standardGeneric("compQ"))
#' Method \code{compP}
#' @name compP
#' @rdname compP-methods
#' @exportMethod compP
setGeneric("compP", function(object, q) standardGeneric("compP"))
#' Method \code{register}
#' @name register
#' @rdname register-methods
#' @exportMethod register
setGeneric("register", function(object1, object2) standardGeneric("register"))

#' @rdname register-methods
#' @aliases register,distributionD-method
#' @description Given two \code{distributionD} objects, it returns two equivalent distributions such that
#' they share the same cdf values. This function is useful for computing basic statistics.
#'
#' @param object1 A \code{distributionD} object
#' @param object2 A \code{distributionD} object
#' @return The two \code{distributionD} objects in input sharing the same cdf (the \code{p}
#' slot)
#' @author Antonio Irpino
#' @references Irpino, A., Lechevallier, Y. and Verde, R. (2006): \emph{Dynamic
#' clustering of histograms using Wasserstein metric} In: Rizzi, A., Vichi, M.
#' (eds.) COMPSTAT 2006. Physica-Verlag, Berlin, 869-876.\cr Irpino, A.,Verde,
#' R. (2006): \emph{A new Wasserstein based distance for the hierarchical
#' clustering of histogram symbolic data} In: Batanjeli, V., Bock, H.H.,
#' Ferligoj, A., Ziberna, A. (eds.) Data Science and Classification, IFCS 2006.
#' Springer, Berlin, 185-192.
#' @keywords distribution
#' @examples
#'
#' ## ---- initialize two distributionD objects mydist1 and mydist2
#' mydist1 <- distributionD(c(1, 2, 3), c(0.1, 0.4, 1))
#' mydist2 <- distributionD(c(7, 8, 10, 15), c(0.1, 0.2, 0.7, 1))
#' ## register the two distributions
#' regDist <- register(mydist1, mydist2)
#'
#' ## OUTPUT:
#' ## regDist$[[1]]
#' ## An object of class "distributionD"
#' ## Slot "x": [1] 1.0 1.5 2.0 2.5 3.0
#' ## Slot "p": [1] 0.1 0.2 0.4 0.7 1.0
#' ## ...
#' ## regDist$[[2]]
#' ## An object of class "distributionD"
#' ## Slot "x": [1] 7.0 8.0 8.8 10.0 15.0
#' ## Slot "p": [1] 0.1 0.2 0.4  0.7  1.0
#' ## ...
#' # The REGISTER function ----
setMethod(
  f = "register", signature = c(object1 = "distributionD", object2 = "distributionD"),
  function(object1, object2) {
    res <- REGISTER2(object1, object2)
    return(res)

  }
)
#' @rdname checkEmptyBins-methods
#' @aliases checkEmptyBins,distributionD-method
#' @description The method checking for empty bins in a distribution, i.e. if two cdf consecutive
#' values are equal. In that case a probability value of \code{1e-7} is
#' assigned to the empty bin and the cdf is recomputed. This methods is useful
#' for numerical reasons.
#'
#'
#' @param object a \code{distributionD} object
#' @return A \code{distributionD} object without empty bins
#' @author Antonio Irpino
#' @keywords distribution
#' @examples
#'
#' ## ---- A mydist distribution with an empty bin i.e. two consecutive values of p are equal----
#' mydist <- distributionD(x = c(1, 2, 3, 10), p = c(0, 0.5, 0.5, 1))
#' ## ---- Checks for empty byns and returns the newdist object without empty bins ----
#' newdist <- checkEmptyBins(mydist)
setMethod(
  f = "checkEmptyBins", signature = "distributionD",
  function(object) {
    w <- diff(object@p)
    TOL <- 1e-14
    if (length(which(w <= TOL))) {
      w[which(w < TOL)] <- 10 * TOL
      object@p <- cumsum(w) / sum(w)
    }
    return(object)
  }
)
#' @rdname compQ-methods
#' @aliases compQ,distributionD-method
#' @description Compute the quantile value of a discrete DD for a given probability.
#'
#'
#' @param object an object of \env{distributionD} class
#' @param p a number between 0 and 1
#' @return \deqn{y= F^{-1}(p)=Q(p)} A number that is the quantile of the passed
#' discrete DD \env{object} at level \env{p}.
#' @author Antonio Irpino
#' @examples
#'
#' ## ---- A mydist distribution ----
#' mydist <- distributionD(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' ## ---- Compute the quantile of mydist for different values of p ----
#' y <- compQ(mydist, 0.5) # the median
#' y <- compQ(mydist, 0) # the minimum
#' y <- compQ(mydist, 1) # the maximum
#' y <- compQ(mydist, 0.25) # the first quartile
#' y <- compQ(mydist, 0.9) # the ninth decile
setMethod(
  f = "compQ", signature = c(object = "distributionD", p = "numeric"),
  function(object, p) {
    # %Computes the p-th quantile p=[0,1] of the distribution o1
    # %INPUT  - p  a value in [0,1]
    # %           - o1 a distribution
    # %OUTPUT-  res the computed quantile
    # %example:q=compQ(o1,0.5) returns the median of the
    # % distribution o1

    # Check for errors
    if (p < 0 || p > 1) stop("p must be a value between 0 and 1")

    if (p <= 0) {
      return(q = object@x[1])
    }
    if (p >= 1) {
      return(q = object@x[length(object@x)])
    }

    qua<-e1071::qdiscrete(p, c(object@p[1],diff(object@p)), values = object@x)
    return(qua)
  }
)
#' @rdname compP-methods
#' @aliases compP,distributionD-method
#' @description Compute the cdf probability at a given value for a discrete DD
#'
#' @param object is an object of \env{distributionD} class
#' @param q is a numeric value
#' @return Returns a value between 0 and 1.
#' @keywords distribution
#' @examples
#'
#' ## ---- A mydist distribution ----
#' mydist <- distributionD(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' ## ---- Compute the cfd value for q=5 (not observed) ----
#' p <- compP(mydist, 5)
setMethod(
  f = "compP", signature = c(object = "distributionD", q = "numeric"),
  function(object, q) {
    p<-e1071::pdiscrete(q, c(object@p[1],diff(object@p)), values = object@x)
    return(p)
  }
)
# L2 Wasserstein distance between two distributions and related results ----
#' Method \code{WassSqDistD}
#' @name WassSqDistD
#' @rdname WassSqDistD-methods
#' @exportMethod WassSqDistD
setGeneric("WassSqDistD", function(object1, object2, ...) standardGeneric("WassSqDistD")) # Wasserstein distance between two distributions
#' Method \code{rQQ}
#' @name rQQ
#' @rdname rQQ-methods
#' @exportMethod rQQ
setGeneric("rQQ", function(e1, e2) standardGeneric("rQQ")) # Quantile-Quantile correlation between two distributions

#' @rdname WassSqDistD-methods
#' @aliases WassSqDistD,distributionD-method
#' @description Computes the squared L2 Wasserstein distance between two \code{distributionD} objects.
#' @param object1 is an object of \env{distributionD} class
#' @param object2 is an object of \env{distributionD} class
#' @param ... optional parameters
#' @param details (optional, default=FALSE) is a logical value, if TRUE returns the decomposition of the distance
#' @return
#' If \code{details=FALSE}, the function returns the squared L2 Wasserstein distance.\cr
#' If \code{details=TRUE}, the function returns list containing the squared distance, its
#' decomposition in three parts (position, size and shape) and the correlation coefficient between the quantile functions.
#' @references
#' Irpino, A. and Romano, E. (2007): \emph{Optimal histogram representation of large data sets:
#' Fisher vs piecewise linear approximations}. RNTI E-9, 99-110.\cr
#' Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI 10.1007/s11634-014-0176-4
#' @keywords distribution
#' @examples
#' ## ---- create two distributionD objects ----
#' mydist1 <- distributionD(x = c(1, 2, 3), p = c(0.1, 0.4, 1))
#' mydist2 <- distributionD(x = c(7, 8, 10, 15), p = c(0.1, 0.2, 0.7, 1))
#' # -- compute the squared L2 Waaserstein distance
#' WassSqDistD(mydist1, mydist2)
#' # -- compute the squared L2 Waaserstein distance with details
#' WassSqDistD(mydist1, mydist2, details = TRUE)
setMethod(
  f = "WassSqDistD", signature = c(object1 = "distributionD", object2 = "distributionD"),
  # Computes the L2 Wasserstein squared distance between two distributions
  # INPUT: object1 and object2 - two distributionD objects
  # OUTPUT: A list containing the distance and its decomposition in three parts (position, size and shape)
  function(object1 = object1, object2 = object2, details = FALSE) {
    a<-object1@x
    b<-object2@x
    wa<-c(object1@p[1],diff(object1@p))
    wb<-c(object2@p[1],diff(object2@p))

    m <- length(a)
    n <- length(b)
    if (is.null(wa)) {
      wa <- rep(1, m)
    }
    else {
      wha <- wa > 0
      wa <- wa[wha]
      a <- a[wha]
      m <- length(a)
    }
    if (is.null(wb)) {
      wb <- rep(1, n)
    }
    else {
      whb <- wb > 0
      wb <- wb[whb]
      b <- b[whb]
      n <- length(b)
    }
    ua <- (wa/sum(wa))[-m]
    ub <- (wb/sum(wb))[-n]
    cua <- c(cumsum(ua))
    cub <- c(cumsum(ub))
    arep<- countC(cub,cua)
    brep<- countC(cua,cub)
    aa <- rep(a, times = arep)
    bb <- rep(b, times = brep)
    uu <- sort(c(cua, cub),method="quick")
    uu0 <- c(0, uu)
    uu1 <- c(uu, 1)
    D <- (sum((uu1 - uu0) * (bb - aa)^2))

    if (details) {
      DC <- (object1@m - object2@m)^2
      DS <- (object1@s - object2@s)^2

      DR <- ifelse(abs(D - DC - DS) < 1e-30, yes = 0, no = abs(D - DC - DS))
      rho <- 1-DR/(2*object1@s*object2@s)
      if (rho < 0) rho <- 0
      resu <- c(D, DC, DS, DR, rho)
      names(resu) <- c("SQ_W_dist", "POSITION", "SIZE", "SHAPE", "rQQ")
      return(resu)
    }
    else {
      return(as.numeric(D))
    }
  }
)

#' Method \code{dotpW}
#' @name dotpW
#' @rdname dotpW-methods
#' @exportMethod dotpW
setGeneric("dotpW", function(e1, e2) standardGeneric("dotpW")) # dotproduct from L2 Wasserstein
#' @rdname dotpW-methods
#' @aliases dotpW,distributionD-method
#' @description The dot product of two distributions inducing the L2 Wasserstein metric
#'
#' @param e1 a \code{distributionD} object or a number
#' @param e2 a \code{distributionD} object or a number
#' @return A numeric value
#' @author Antonio Irpino
#' @references  Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI 10.1007/s11634-014-0176-4
#' @keywords distribution
#' @examples
#'
#' ## let's define two distributionD objects
#' mydist1 <- distributionD(x = c(1, 2, 3, 10), p = c(0, 0.1, 0.5, 1))
#' mydist2 <- distributionD(x = c(5, 7, 15), p = c(0, 0.7, 1))
#'
#' ## the dot product between the distributions
#' dotpW(mydist1, mydist2) #---> 39.51429
#'
#' ## the dot product between a distribution and a numeric
#' dotpW(mydist1, 3) #---> 13.2
#' dotpW(3, mydist1) #---> 13.2
#'
#'
#' # DOTPW method -----
setMethod("dotpW",
          signature(e1 = "distributionD", e2 = "distributionD"),
          definition = function(e1, e2) {

            return(c_dotpW(e1, e2))
          }
)
#' @rdname dotpW-methods
#' @aliases dotpW,distributionD-method
#' @description The dot product of a number (considered as an impulse distribution function) and a distribution
setMethod(
  "dotpW",
  signature(e1 = "numeric", e2 = "distributionD"),
  function(e1, e2) {
    dprod <- e1 * e2@m
    return(dprod)
  }
)
#' @rdname dotpW-methods
#' @aliases dotpW,distributionD-method
#' @description The dot product of a distribution and a number (considered as an impulse distribution function).
setMethod(
  "dotpW",
  signature(e1 = "distributionD", e2 = "numeric"),
  function(e1, e2) {
    dprod <- e2 * e1@m
    return(dprod)
  }
)
#' @rdname rQQ-methods
#' @aliases rQQ,distributionD-method
#' @description Quantile-Quantile correlation between two distributions
#' @param e1  A \code{distributionD} object
#' @param e2  A \code{distributionD} object
#' @return Pearson correlation index between quantiles
#' @author Antonio Irpino
#' @references  Irpino, A., Verde, R. (2015) \emph{Basic
#' statistics for distributional symbolic variables: a new metric-based
#' approach} Advances in Data Analysis and Classification, DOI 10.1007/s11634-014-0176-4
#' @examples
#'
#' ## ---- initialize two distributionD object mydist1 and mydist2
#' mydist1 <- distributionD(x = c(1, 2, 3), p = c(0.1, 0.4, 1))
#' mydist2 <- distributionD(x = c(7, 8, 10, 15), p = c(0.05, 0.2, 0.7, 1))
#' ## computes the rQQ
#' rQQ(mydist1, mydist2)

#' @export rQQ
setMethod("rQQ",
          signature(e1 = "distributionD", e2 = "distributionD"),
          definition = function(e1, e2) {
            rQQ <- (dotpW(e1, e2) - e1@m * e2@m) / (e1@s * e2@s)
            return(rQQ)
          }
)
## ---- Show overridding for distributionD and MatH ----
#' Method show for distributionD
#' @name show
#' @rdname show-distributionD-methods
#' @docType methods
#' @aliases show,distributionD-method
#' @description An overriding show function for a \code{distributionD} object. The function returns a representation
#' of the discrete DD, if the number of bins is high the central part of the discrete DD is truncated.
#' @param object a \code{distributionD} object
#' @examples
#' ## ---- initialize a distributionD
#' mydist <- distributionD(x = c(7, 8, 10, 15), p = c(0, 0.2, 0.7, 1))
#' # show the discrete DD
#' mydist
setMethod("show",
          signature(object = "distributionD"),
          definition = function(object) {
            if (length(object@p) > 1) {
              mymat <- matrix(0, length(object@p) - 1, 2)
              if (length(object@p) > 11) {
                cat("Output shows the first five and the last five rows due to eccesive length \n")
                mymat <- matrix(0, 12, 2)
                mymat[1, 1] <- "X"
                mymat[1, 2] <- "p"
                pp<-c(object@p[1],diff(object@p))
                count <- 0
                for (i in 1:5) {
                  count <- count + 1
                  mymat[count + 1, 1] <- format(object@x[i], digits = 5)
                  mymat[count + 1, 2] <- format(pp[i], digits = 4)
                }
                count <- count + 1
                mymat[count + 1, 1] <- paste("...")
                mymat[count + 1, 2] <- paste("...")
                for (i in (length(object@p) - 4):length(object@p)) {
                  count <- count + 1
                  mymat[count + 1, 1] <- format(object@x[i], digits = 5)
                  mymat[count + 1, 2] <- format(pp[i], digits = 4)
                }
                rownames(mymat) <- c(
                  "", paste("v", 1:5, sep = "_"), "...",
                  paste("v", (length(object@p) - 4):(length(object@p) ), sep = "_")
                )
                write.table(format(mymat, justify = "right", digits = 5),
                            row.names = T, col.names = F, quote = F
                )
              }
              else {
                mymat <- matrix(0, length(object@p)+1, 2)
                mymat[1, 1] <- "X"
                mymat[1, 2] <- "p"
                pp<-c(object@p[1],diff(object@p))
                if (length(object@p) > 2) {
                  for (i in 1:length(object@p) ) {
                    mymat[i+1, 1] <- format(object@x[i], digits = 5)
                    mymat[i+1, 2] <- format(pp[i], digits = 4)
                  }
                }
                rownames(mymat) <- c(" ", paste("v", 1:(length(object@p)), sep = "_"))
                write.table(format(mymat, justify = "right"),
                            row.names = T, col.names = F, quote = F
                )
              }
              cat(paste("\n mean = ", format(signif(object@m, 6), digits = 6), "  std  = ", format(signif(object@s, 6), digits = 6), "\n "))
            } else {
              (cat("Empty distributionD\n"))
            }
          }
)




## --- Plot  for distributionD  ----
#' plot for a distributionD object
#' @name plot-distributionD
#' @docType methods
#' @aliases plot,distributionD-method
#' @description A plot function for a \code{distributionD} object. The function returns a representation
#' of the discrete DD.
#' @param x  a \code{distributionD} object
#' @param type (optional) a string describing the type of plot, default="BAR".\cr Other allowed types are
#' \cr"STAOR" stached horizontal \cr"STAVER" stacked vertical
#' \cr"CDF"=Cumulative distribution function, \cr"QF"= quantile function,
#' \cr"HBOXPLOT"=horizontal boxplot, \cr"VBOXPLOT"= vertical boxplot,
#' @param col (optional) a string the color of the plot, default="green".
#' @param border (optional) a string the color of the border of the plot, default="black".
#' @examples
#' ## ---- initialize a distributionD
#' mydist <- distributionD(x = c(7, 8, 10, 15), p = c(0.1, 0.2, 0.7, 1))
#' # show the discrete DD
#' plot(mydist) # plots mydist
#' plot(mydist, type = "BAR", col = "red", border = "blue") # plots mydist
#' plot(mydist, type = "STAOR", col = "red", border = "blue") # plots mydist stacked orizz
#' plot(mydist, type = "STAVER", col = "red", border = "blue") # plots mydist stacked vert
#' plot(mydist, type = "HBOXPLOT") # plots a horizontal boxplot for mydist
#' plot(mydist, type = "VBOXPLOT") # plots a vertical boxplot for mydist
#' plot(mydist, type = "CDF") # plots the cumulative distribution function of mydist
#' plot(mydist, type = "QF") # plots the quantile function of mydist
#' @importFrom utils write.table
#' @export
setMethod(
  "plot",
  signature(x = "distributionD"),
  function(x, type = "BAR", col = "green", border = "black",title=NULL) {
    require(dplyr)
    require(ggplot2)
    require(colorspace)
    require(forcats)
    require(tidyr)
    if(type=="BAR"){
      DF<-data.frame(x=x@x,p=c(x@p[1],diff(x@p)))
      p<-ggplot2::ggplot(DF)+geom_bar(aes(x=x,y=p),stat="identity",fill=col,color=border)+
        ggtitle(title)+xlab((paste("\n mean = ", format(signif(x@m, 6), digits = 6), "  std  = ", format(signif(x@s, 6), digits = 6), "\n ")))
    }
    if(type=="STAVER"){
      DF<-data.frame(x=factor(x@x),p=c(x@p[1],diff(x@p)))%>% mutate(x=fct_rev(x))
      p<-ggplot2::ggplot(DF)+geom_bar(aes(x=1,y=p,fill=x),stat="identity",color=border)+
        ggtitle(title)+
        scale_fill_discrete_diverging(palette = "Blue-Red2", rev = TRUE)+
        guides(fill = guide_legend(reverse = TRUE))+
        xlab((paste("\n mean = ", format(signif(x@m, 6), digits = 6), "  std  = ", format(signif(x@s, 6), digits = 6), "\n ")))
    }
    if(type=="STAOR"){
      DF<-data.frame(x=factor(x@x),p=c(x@p[1],diff(x@p))) %>% mutate(x=fct_rev(x))
      p<-ggplot2::ggplot(DF)+geom_bar(aes(x=1,y=p,fill=x),stat="identity",color=border)+
        ggtitle(title)+coord_flip()+
        scale_fill_discrete_diverging(palette = "Blue-Red2", rev = TRUE)+
        guides(fill = guide_legend(reverse = TRUE))+
        labs(caption=(paste("\n mean = ", format(signif(x@m, 6), digits = 6), "  std  = ", format(signif(x@s, 6), digits = 6), "\n ")))

    }
    if(type=="CDF"){
      DF<-data.frame(x=c(x@x[1],x@x),p=c(0,x@p))
      DF2<-DF %>%
        mutate(type = "cdf") %>%
        bind_rows(DF %>%
                    mutate(type = "prior",
                           p = lag(p))) %>%
        drop_na() %>%
        arrange(x, desc(type))
      p<- ggplot(DF2) +
        geom_segment(aes(x = lag(x), y = lag(p),
                         xend = x, yend = p,
                         lty = type,color = type,alpha=type),linewidth = 1,show.legend = F) +
        geom_point(aes(x, p, fill = type),
                   shape = 21,size=2,show.legend = F) +
        scale_fill_manual(values = c("black", "white"))+
        scale_linetype_manual(values = c("dashed", "solid"))+
        scale_color_manual(values = c("black","black"))+
        scale_alpha_manual(values = c( 0.3,1))+
        theme_bw()+xlab("x")+ylab("p")+labs(subtitle="CDF")
    }
    if(type=="QF"){
      DF<-data.frame(x=c(x@x[1],x@x),p=c(0,x@p))
      DF2<-DF %>%
        mutate(type = "cdf") %>%
        bind_rows(DF %>%
                    mutate(type = "prior",
                           p = lag(p))) %>%
        drop_na() %>%
        arrange(x, desc(type))
      p<- ggplot(DF2) +
        geom_segment(aes(x = lag(x), y = lag(p),
                         xend = x, yend = p,
                         lty = type,alpha = type),color="black",linewidth = 1,show.legend = F) +
        geom_point(aes(x, p, fill = type),
                   shape = 21,size=2,show.legend = F) +
        scale_fill_manual(values = c("white","black"))+
        scale_linetype_manual(values = c("solid","dashed"))+
        scale_alpha_manual(values = c( 1,0.3))+
        theme_bw()+xlab("x")+ylab("p")+labs(subtitle="QF")+coord_flip()
    }
   p
  }
)
