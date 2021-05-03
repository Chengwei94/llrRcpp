#' Local linear regression (Regular and binned version)
#' @param x a numeric vector/matrix of x data or a data object(bin) used to select the method
#' @param ... further arguments to be passed
#' @aliases llr-class llr
#' @return returns a S3 object of class "llr" containing
#' \describe{
#'    \item{x}{sorted numeric vector or matrix of \code{xpred}}
#'    \item{y}{estimated values corresponding to 'llr' class \code{x}}
#' }
#' @export
llr <- function(x, ...) UseMethod("llr")

#' @rdname llr
#' @method llr default
#' @param y a numeric vector of y data corresponding to \code{x}.
#' @param xpred a numeric vector or matrix of same dimension as \code{x}. x values to be predicted.
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic".
#' @param bandwidth a numeric vector or single number of same dimension as \code{x}.
#' @param weight a numeric vector of \code{length(x)} for weight of each data point.
#' @param kdtree boolean flag: If \code{TRUE}, kdtree is used for computation of local linear regression.
#' @param approx boolean flag: If \code{TRUE}, kdtree approximation is used . Only used when \code{kdtree = TRUE}.
#' @param epsilon margin of error allowed for llr approximation using kdtree. Only used when both \code{kdtree = TRUE} and \code{approx = TRUE}.
#' @param N_min minimum number of points stored in the kd-tree. Only used when both \code{kdtree = TRUE} and \code{approx = TRUE}. Currently not in use.
#' @examples
#' \dontrun{
#' n <- 1000
#' x <- seq(0, 10, length.out = n)
#' x1 <- rnorm(n, 0, 0.2)
#' y <- sin(x) + x1
#' w <- rep(1 / n, n)
#' binned <- bin(x, y, bins = 400, w)
#' ## local linear regression for exact without kdtree
#' llr_exact <- llr(x, y, x, bandwidth = 0.2, weight = w)
#' ## local linear regression for kdtree exact
#' llr_kdexact <- llr(x, y, x, bandwidth = 0.2, weight = w, kdtree = TRUE)
#' ## local linear regression for kdtree approximation
#' llr_kdapprox <- llr(x, y, x, bandwidth = 0.2, weight = w, kdtree = TRUE, approx = TRUE)
#' ## local linear regression for data after binning.
#' llr_bin <- llr(binned, x, bandwidth = 0.2)
#' }
#' @export
llr.default <- function(x, y, xpred, kernel = "epanechnikov", bandwidth, weight, kdtree = FALSE, approx = FALSE,
                        epsilon = 0.05, N_min = 1, ...) {
  
  kernel <- match.arg(kernel)
  
  switch(kernel,
    epanechnikov = {
      kcode <- 1
    },
    rectangular = {
      kcode <- 2
    },
    triangular = {
      kcode <- 3
    },
    quartic = {
      kcode <- 4
    },
    triweight = {
      kcode <- 5
    },
    tricube = {
      kcode <- 6
    },
    cosine = {
      kcode <- 7
    },
    gauss = {
      kcode <- 21
    },
    logistic = {
      kcode <- 22
    },
    sigmoid = {
      kcode <- 23
    },
    silverman = {
      kcode <- 24
    }
  )
  # helper functions
  normalize <- function(x) {
    return(max(x) - min(x))
  }

  x <- as.matrix(x)
  y <- as.numeric(y)
  xpred <- as.matrix(xpred)
  wt <- as.numeric(weight)
  bandwidth <- as.numeric(bandwidth)

  xpred <- as.matrix(xpred[order(xpred[, 1]), ])
  if (ncol(x) > 1 && ncol(xpred) == 1) {
     xpred <- t(xpred)
  }
  
  if(ncol(x) != ncol(xpred)){
    stop("Dimension of x and dimension of xpred is not the same ")
  }
  
  if (nrow(x) != length(y) || nrow(x) != length(wt)){
    stop('x, y and weight must have the same length')
  }
  
  if(ncol(x) != length(bandwidth)) {
    stop('x and bandwidth should have the same dimension')
  }
  #Ordering of the x, y, wt by the first column
  df <- cbind(x, y, wt)
  df <- df[order(df[, 1]), ]
  x <- as.matrix(df[, 1:ncol(x)])
  y <- as.numeric(df[, ncol(df) - 1])
  wt <- as.numeric(df[, ncol(df)])

  scale <- apply(x, 2, FUN = normalize)
  scale <- round(scale, 3)
  bandwidth <- bandwidth * scale

  if (kdtree) {
    if (!approx) {
      ypred <- llrt_cpp(x, y, xpred, wt, 1, kcode, epsilon, bandwidth, N_min)
    }
    else {
      ypred <- llrt_cpp(x, y, xpred, wt, 2, kcode, epsilon, bandwidth, N_min)
    }
  }
  else {
    if (ncol(x) == 1) {
      ypred <- llr1d_cpp(x, y, xpred, kcode, bandwidth, wt)
    }
    if (ncol(x) == 2) {
      ypred <- llr2d_cpp(x, y, xpred, kcode, bandwidth, wt)
    }
    if (ncol(x) >= 3) {
      ypred <- llr_cpp(x, y, xpred, kcode, bandwidth, wt)
    }
  }
  results <- list("x" = xpred, "fitted" = ypred)
  class(results) <- "llr"
  
  results
}

#' @rdname llr
#' @method llr bin
#' @param x bin object
#' @param xpred a numeric vector or matrix of same dimension as \code{x}. x values to be predicted.
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic".
#' @export
llr.bin <- function(x, xpred, kernel = "epanechnikov", bandwidth, ...) {

  kernel <- match.arg(kernel)
  
  if (!inherits(x, "bin")) {
    stop("Function only works for objects of class bin")
  }
  
  switch(kernel,
    epanechnikov = {
      kcode <- 1
    },
    rectangular = {
      kcode <- 2
    },
    triangular = {
      kcode <- 3
    },
    quartic = {
      kcode <- 4
    },
    triweight = {
      kcode <- 5
    },
    tricube = {
      kcode <- 6
    },
    cosine = {
      kcode <- 7
    },
    gauss = {
      kcode <- 21
    },
    logistic = {
      kcode <- 22
    },
    sigmoid = {
      kcode <- 23
    },
    silverman = {
      kcode <- 24
    }
  )
  # helper functions
  normalize <- function(x) {
    return(max(x) - min(x))
  }
  
  y <- as.numeric(x$y)
  wt <- as.numeric(x$weight)
  x <- as.matrix(x$x)
  
  xpred <- as.matrix(xpred)
  xpred <- xpred[order(xpred[, 1]), ]
  xpred <- as.matrix(xpred)
  bandwidth <- as.numeric(bandwidth)
  # as.matrix transform a 1 x c matrix to c x 1 matrix by default
  if (ncol(x) > 1 && ncol(xpred == 1)) {
   xpred <- t(xpred)
  }
  
  if(ncol(x) != ncol(xpred)){
    stop("Dimension of x and xpred should be the same ")
  }
  
  if (nrow(x) != length(y) || nrow(x) != length(wt)){
    stop('x, y and weight must have the same length')
  }
  
  if(ncol(x) != length(bandwidth)) {
    stop('x and bandwidth should have the same dimension')
  }
  
  if (ncol(x) >= 3){
    stop("Only supports object of 2d")
  }
  
  df <- data.frame(x, y, wt)
  df <- df[order(df[, 1]), ]
  
  x <- as.matrix(df[, 1:ncol(x)])
  y <- as.numeric(df[, ncol(df) - 1])
  wt <- as.numeric(df[, ncol(df)])

  scale <- apply(x, 2, FUN = normalize)
  scale <- round(scale, 3)
  bandwidth <- bandwidth * scale

  if (ncol(x) == 1) {
    ypred <- llr1d_cpp(x, y, xpred, kcode, bandwidth, wt)
  }
  if (ncol(x) == 2) {
    ypred <- llr2d_cpp(x, y, xpred, kcode, bandwidth, wt)
  }
  results <- list("x" = xpred, "fitted" = ypred)
  class(results) <- "llr"
  
  results
}
