#' Bandwidth selection using leave one out cross validation
#' @param x a numeric vector or matrix of x data
#' @param y a numeric vector of y data corresponding to \code{x} data
#' @param weight a numeric vector of \code{length(x)} for weight of each data point.
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic". "sigmoid" and "silverman".
#' @param approx boolean flag: if \code{true} kdtree approximation is used.
#' @param epsilon margin of error allowed for llr approximation using kdtree. Only used when \code{approx = TRUE}.
#' @param N_min minimum number of points stored in the kd-tree. Only used when \code{approx} = TRUE.
#' @param bandwidth a numeric vector or matrix of bandwidth considered for selection.
#' @return returns a single or numeric vector of \code{bandwidth} that gives the smallest mean square error.
#' @examples
#' \dontrun{
#' n <- 1000
#' x <- seq(0, 10, length.out = n)
#' x1 <- rnorm(n, 0, 0.2)
#' y <- sin(x) + x1
#' w <- rep(1 / n, n)
#' bandwidth <- seq(0.02, 0.4, by = 0.01)
#' binned <- bin(x, y, bins = 400, w)
#' ## Bandwidth selection of binned data
#' h_bin <- gcv.llr(binned$x, binned$y, binned$weight, bandwidth = bandwidth)
#' ## Bandwidth selection of exact local linear regression
#' h_exact <- gcv.llr(x, y, w, bandwidth = bandwidth)
#' ## Bandwidth selection of approx local linear regression with kdtree
#' h_kdapprox <- gcv.llr(x, y, w, approx = TRUE, bandwidth = bandwidth)
#' }
#' @export
loocv.llr <- function(x, y, weight, kernel = "epanechnikov", approx = FALSE, epsilon = 0.05,
                    N_min = 1, bandwidth) {
  
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
  wt <- as.numeric(weight)
  bandwidth <- data.frame(bandwidth)
  
  if (ncol(x) > 1 && ncol(bandwidth) == 1) {
    bandwidth <- t(bandwidth)
  }
  
  if (nrow(x) != length(y) || nrow(x) != length(wt)){
    stop('x, y and weight must have the same length')
  }
  
  if(ncol(x) != ncol(bandwidth)) {
    stop('x and bandwidth should have the same dimension')
  }
  
  scale <- apply(x, 2, FUN = normalize)
  scale <- as.numeric(round(scale , 3))
  bandwidth1 <- mapply('*', bandwidth, scale)
  bandwidth1 <- as.matrix(bandwidth1)
  
  if (ncol(x) > 1 && ncol(bandwidth1) == 1) {
    bandwidth1 <- t(bandwidth1)
  }

  if (approx == FALSE) {
    results <- tgcv_cpp(x, y, wt, 1, kcode, epsilon, bandwidth1, N_min)
  }
  else {
    results <- approx_gcv_cpp(x, y, wt, 2, kcode, epsilon, bandwidth1, N_min)
  }
  
   results$bw_opt <- results$bw_opt / scale
   results$bandwidth <- bandwidth 
   
   if (any(results$bw_opt == 0)) {
     stop("Choose a larger range of bandwidth")
   }
  results
}
