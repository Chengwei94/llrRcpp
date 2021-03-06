% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/llr.R
\name{llr}
\alias{llr}
\alias{llr-class}
\alias{llr.default}
\alias{llr.bin}
\title{Local linear regression (Regular and binned version)}
\usage{
llr(x, ...)

\method{llr}{default}(
  x,
  y,
  xpred,
  kernel = "epanechnikov",
  bandwidth,
  weight,
  kdtree = FALSE,
  approx = FALSE,
  epsilon = 0.05,
  N_min = 1,
  ...
)

\method{llr}{bin}(x, xpred, kernel = "epanechnikov", bandwidth, ...)
}
\arguments{
\item{x}{bin object}

\item{...}{further arguments to be passed}

\item{y}{a numeric vector of y data corresponding to \code{x}.}

\item{xpred}{a numeric vector or matrix of same dimension as \code{x}. x values to be predicted.}

\item{kernel}{kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic".}

\item{bandwidth}{a numeric vector or single number of same dimension as \code{x}.}

\item{weight}{a numeric vector of \code{length(x)} for weight of each data point.}

\item{kdtree}{boolean flag: If \code{TRUE}, kdtree is used for computation of local linear regression.}

\item{approx}{boolean flag: If \code{TRUE}, kdtree approximation is used . Only used when \code{kdtree = TRUE}.}

\item{epsilon}{margin of error allowed for llr approximation using kdtree. Only used when both \code{kdtree = TRUE} and \code{approx = TRUE}.}

\item{N_min}{minimum number of points stored in the kd-tree. Only used when both \code{kdtree = TRUE} and \code{approx = TRUE}. Currently not in use.}
}
\value{
returns a S3 object of class "llr" containing
\describe{
   \item{x}{sorted numeric vector or matrix of \code{xpred}}
   \item{y}{estimated values corresponding to 'llr' class \code{x}}
}
}
\description{
Local linear regression (Regular and binned version)
}
\examples{
\dontrun{
n <- 1000
x <- seq(0, 10, length.out = n)
x1 <- rnorm(n, 0, 0.2)
y <- sin(x) + x1
w <- rep(1 / n, n)
binned <- bin(x, y, bins = 400, w)
## local linear regression for exact without kdtree
llr_exact <- llr(x, y, x, bandwidth = 0.2, weight = w)
## local linear regression for kdtree exact
llr_kdexact <- llr(x, y, x, bandwidth = 0.2, weight = w, kdtree = TRUE)
## local linear regression for kdtree approximation
llr_kdapprox <- llr(x, y, x, bandwidth = 0.2, weight = w, kdtree = TRUE, approx = TRUE)
## local linear regression for data after binning.
llr_bin <- llr(binned, x, bandwidth = 0.2)
}
}
