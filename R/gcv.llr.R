#' Bandwidth selection using generalized cross validation 
#' @param x a numeric vector or matrix of x data
#' @param y a numeric vector of y data corresponding to \code{x} data
#' @param weight a numeric vector of \code{length(x)} for weight of each data point.
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic". "sigmoid" and "silverman".
#' @param approx boolean flag: if \code{true} kdtree approximation is used.
#' @param epsilon margin of error allowed for llr approximation using kdtree. Only used when \code{approx = TRUE}.
#' @param N_min minimum number of points stored in the kd-tree. Only used when \code{approx} = TRUE.
#' @param bw a numeric vector or matrix of bandwidth considered for selection.
#' @return returns a single or numeric vector of \code{bw} that gives the smallest mean square error.
#' @examples    
#' n <- 1000
#' x <- seq(0,10,length.out = n)
#' x1 <- rnorm(n,0,0.2)
#' y <- sin(x) + x1
#' w <- rep(1/n, n)
#' binned <- bin(x, y, bins=400, w)
#' ## Bandwidth selection of binned data
#' h_bin <- gcv.llr(binned$x, binned$y, binned$weight)
#' ## Bandwidth selection of exact local linear regression 
#' h_exact <- gcv.llr(x, y, w)
#' ## Bandwidth selection of approx local linear regression with kdtree
#' h_kdapprox <- gcv.llr(x, y, w, approx = TRUE) 
#' @export
gcv.llr <- function(x, y, weight, kernel = "epanechnikov", approx = FALSE, epsilon = 0.05,
                    N_min = 1, bw = seq(0.05, 0.4, by = 0.01)){
  
  x <- as.matrix(x)
  y <- as.numeric(y)
  bw <- as.matrix(bw)
  
  switch(kernel, 
         epanechnikov = {kcode <- 1}, 
         rectangular = {kcode <-2},
         triangular = {kcode <- 3}, 
         quartic = {kcode <- 4}, 
         triweight = {kcode <- 5}, 
         tricube = {kcode <- 6}, 
         cosine = {kcode <- 7}, 
         gauss = {kcode <- 21}, 
         logistic = {kcode<- 22},
         sigmoid = {kcode <- 23}, 
         silverman = {kcode <- 24}
  )
  
  #helper functions
  normalize <- function(x)
  {
    return(max(x)- min(x))
  }
  
  scale <- apply(x, 2, FUN = normalize)
  bw <- bw * scale
  
  if(approx == FALSE){ 
    hopt <- tgcv_cpp(x, y, weight, 1, kcode, epsilon, bw, N_min)  
  }
  if(approx == TRUE) { 
    hopt <- tgcv_cpp(x ,y ,weight, 2, kcode, epsilon, bw, N_min)  
  }
  hopt/scale  
}