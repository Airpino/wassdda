#' @title Class 'distributionK'
#' @description Represents a univariate kernel density estimate using a cumulative distribution (CDF)
#' @slot x numeric vector of domain values (sorted)
#' @slot p numeric vector of CDF values, starting at 0 and ending at 1
#' @slot m numeric, mean of the distribution
#' @slot s numeric, standard deviation
#' @slot kernel character, type of kernel used (e.g., "gaussian")
#' @slot bw numeric, bandwidth of the kernel estimate
#' @export
setClass(
  "distributionK",
  slots = c(
    x = "numeric",
    p = "numeric",
    m = "numeric",
    s = "numeric",
    kernel = "character",
    bw = "numeric"
  ),
  validity = function(object) {
    if (length(object@x) != length(object@p)) return("x and p must be same length")
    if (length(object@x) <= 1) return("At least two x values required")
    if (is.unsorted(object@x)) return("x must be sorted")
    if (is.unsorted(object@p)) return("p must be non-decreasing")
    if (abs(object@p[1]) > 1e-8 || abs(tail(object@p, 1) - 1) > 1e-8) return("p must start at 0 and end at 1")
    if (!object@kernel %in% c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")) {
      return("Invalid kernel type")
    }
    TRUE
  }
)

#' @title Compute Mean of Kernel CDF
#' @param x numeric vector of domain values
#' @param p numeric vector of CDF values
#' @return numeric mean
#' @export
meanK <- function(x, p) {
  dx <- diff(x)
  dp <- diff(p)
  mid_x <- (x[-length(x)] + x[-1]) / 2
  sum(mid_x * dp)
}

#' @title Compute Standard Deviation of Kernel CDF
#' @param x numeric vector of domain values
#' @param p numeric vector of CDF values
#' @param mean (optional) precomputed mean
#' @return numeric standard deviation
#' @export
stdK <- function(x, p, mean = NULL) {
  if (is.null(mean)) mean <- meanK(x, p)
  dx <- diff(x)
  dp <- diff(p)
  mid_x <- (x[-length(x)] + x[-1]) / 2
  sqrt(sum((mid_x - mean)^2 * dp))
}

#' @title Initializer for distributionK
#' @keywords internal
setMethod("initialize", "distributionK",
          function(.Object, x = numeric(0), p = numeric(0), m = numeric(0), s = numeric(0),
                   kernel = "gaussian", bw = 1) {
            .Object@x <- x
            .Object@p <- p
            .Object@kernel <- kernel
            .Object@bw <- bw

            if (length(x) > 1 && length(p) == length(x)) {
              if (length(m) == 0) .Object@m <- meanK(x, p) else .Object@m <- m
              if (length(s) == 0) .Object@s <- stdK(x, p, .Object@m) else .Object@s <- s
            } else {
              .Object@m <- m
              .Object@s <- s
            }

            validObject(.Object)
            return(.Object)
          }
)

#' @title Class 'MatK': Matrix of distributionK
#' @description A rectangular array of kernel-based distributions
#' @slot M matrix of distributionK objects
#' @slot nrows number of rows
#' @slot ncols number of columns
#' @slot names.rows character vector with row names
#' @slot names.cols character vector with column names
#' @slot by.row logical, whether the matrix was filled row-wise
#' @export
setClass(
  "MatK",
  slots = c(
    M = "matrix",              # matrix of distributionK objects
    nrows = "numeric",
    ncols = "numeric",
    names.rows = "character",
    names.cols = "character",
    by.row = "logical"
  ),
  validity = function(object) {
    if (!all(sapply(object@M, is, "distributionK"))) return("All entries must be distributionK")
    if (nrow(object@M) != object@nrows || ncol(object@M) != object@ncols) return("Dimension mismatch")
    TRUE
  }
)

#' @title Constructor for distributionK objects
#' @description
#' Creates a new object of class \code{distributionK}, representing a kernel density estimate
#' in cumulative form (CDF), along with its mean and standard deviation.
#'
#' @param x A numeric vector of sorted domain values.
#' @param p A numeric vector of CDF values corresponding to \code{x}. Must start at 0 and end at 1.
#' @param kernel A character string indicating the kernel type. Must be one of:
#'   \code{"gaussian"}, \code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"},
#'   \code{"biweight"}, \code{"cosine"}, \code{"optcosine"}. Default is \code{"gaussian"}.
#' @param bw A positive numeric value specifying the bandwidth used in kernel estimation. Default is 1.
#'
#' @return An object of class \code{distributionK}.
#' @export
#'
#' @examples
#' x <- seq(-3, 3, length.out = 100)
#' p <- pnorm(x)
#' dk <- distributionK(x, p, kernel = "gaussian", bw = 0.5)
distributionK <- function(x, p, kernel = "gaussian", bw = 1) {
  allowed_kernels <- c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")
  if (!kernel %in% allowed_kernels) stop("Invalid kernel type.")
  if (!is.numeric(bw) || length(bw) != 1 || bw <= 0) stop("Bandwidth must be a positive number.")
  if (length(x) != length(p)) stop("x and p must have the same length.")

  m <- meanK(x, p)
  s <- stdK(x, p, m)
  new("distributionK", x = x, p = p, m = m, s = s, kernel = kernel, bw = bw)
}


#' @title Constructor for MatK objects
#' @description
#' Builds a rectangular matrix of kernel density distributions (class \code{distributionK}).
#' Can be used to initialize an empty matrix or from a pre-defined list of \code{distributionK} objects.
#'
#' @param x A list of objects of class \code{distributionK}, in row- or column-wise order.
#'   If \code{NULL}, an empty matrix of default \code{distributionK} objects is created.
#' @param nrows Integer, number of rows in the matrix.
#' @param ncols Integer, number of columns in the matrix.
#' @param rownames Optional character vector of row names.
#' @param varnames Optional character vector of column (variable) names.
#' @param by.row Logical. If \code{TRUE}, the matrix is filled row-wise; otherwise, column-wise.
#'
#' @return An object of class \code{MatK}.
#' @export
#'
#' @examples
#' x <- seq(-3, 3, length.out = 100)
#' p <- pnorm(x)
#' dk1 <- distributionK(x, p)
#' dk2 <- distributionK(x, p, kernel = "epanechnikov", bw = 0.3)
#' mat_k <- MatK(x = list(dk1, dk2, dk1, dk2), nrows = 2, ncols = 2, by.row = TRUE)
MatK <- function(x = NULL, nrows = 1, ncols = 1,
                 rownames = NULL, varnames = NULL, by.row = FALSE) {
  # Se x è NULL, crea matrice vuota
  if (is.null(x)) {
    empty_dist <- new("distributionK", x = numeric(2), p = c(0, 1), m = 0, s = 0, kernel = "gaussian", bw = 1)
    mat <- matrix(rep(empty_dist, nrows * ncols), nrow = nrows, ncol = ncols)
  } else {
    if (!all(sapply(x, function(e) inherits(e, "distributionK")))) {
      stop("All elements in x must be of class 'distributionK'")
    }
    if (length(x) != nrows * ncols) {
      stop("Length of x must match nrows * ncols")
    }
    mat <- if (by.row) {
      matrix(x, nrow = nrows, ncol = ncols, byrow = TRUE)
    } else {
      matrix(x, nrow = nrows, ncol = ncols, byrow = FALSE)
    }
  }

  if (is.null(rownames)) rownames <- paste0("R", seq_len(nrows))
  if (is.null(varnames)) varnames <- paste0("V", seq_len(ncols))
  dimnames(mat) <- list(rownames, varnames)

  new("MatK",
      M = mat,
      nrows = nrows,
      ncols = ncols,
      names.rows = rownames,
      names.cols = varnames,
      by.row = by.row)
}


#### METHODS ----------------

#' @title Show method for distributionK
#' @description Displays basic information about a distributionK object
#' @param object An object of class distributionK
#' @export
setMethod("show", "distributionK", function(object) {
  cat("An object of class 'distributionK'\n")
  cat("Kernel: ", object@kernel, "\n")
  cat("Bandwidth: ", object@bw, "\n")
  cat("Mean: ", round(object@m, 4), "\n")
  cat("Std Dev: ", round(object@s, 4), "\n")
  cat("Support: [", min(object@x), ", ", max(object@x), "]\n", sep = "")
})

#' @title Summary method for distributionK
#' @description Summarizes the distributionK object
#' @param object An object of class distributionK
#' @param ... Unused
#' @return A list summary
#' @export
setMethod("summary", "distributionK", function(object, ...) {
  list(
    kernel = object@kernel,
    bandwidth = object@bw,
    mean = object@m,
    sd = object@s,
    x.range = range(object@x),
    length = length(object@x)
  )
})

#' @title Plot method for distributionK using ggplot2
#' @description Plot the CDF of a kernel density distribution using ggplot2
#' @param x An object of class distributionK
#' @param y Not used
#' @param main Title of the plot
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param ... Further graphical parameters (currently unused)
#' @return A ggplot object
#' @import ggplot2
#' @export
setMethod("plot", signature(x = "distributionK", y = "missing"),
          function(x, y, main = "Kernel CDF", xlab = "x", ylab = "F(x)", ...) {
            df <- data.frame(x = x@x, p = x@p)
            ggplot2::ggplot(df, ggplot2::aes(x = x, y = p)) +
              ggplot2::geom_line(color = "steelblue", size = 1) +
              ggplot2::labs(title = main, x = xlab, y = ylab) +
              ggplot2::theme_minimal()
          })
