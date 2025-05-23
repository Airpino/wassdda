% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Met_MatD.R
\name{WH.var.covar}
\alias{WH.var.covar}
\alias{WH.var.covar,MatD-method}
\title{Method WH.var.covar}
\usage{
WH.var.covar(object, ...)

\S4method{WH.var.covar}{MatD}(object, w = numeric(0))
}
\arguments{
\item{object}{a \code{MatD} object}

\item{...}{some optional parameters}

\item{w}{it is possible to add a vector of weights (positive numbers)
having the same size of the rows of the \code{MatD object},
default = equal weight for each row}
}
\value{
a squared \code{matrix} with the (weighted) variance-covariance values
}
\description{
Compute the variance-covariance matrix of a \code{MatD} object, i.e.
a matrix of values consistent with
a set of distributions equipped with a L2 wasserstein metric.
}
\examples{
WH.var.covar(BLOOD)
# generate a set of random weights
RN <- runif(get.MatD.nrows(BLOOD))
WH.var.covar(BLOOD, w = RN)
}
\references{
Irpino, A., Verde, R. (2015) \emph{Basic
statistics for distributional symbolic variables: a new metric-based
approach} Advances in Data Analysis and Classification, DOI
10.1007/s11634-014-0176-4
}
