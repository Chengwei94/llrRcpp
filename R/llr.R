#' Local linear Regression 
#' @param x a numeric vector/matrix of x data or a data object(bin) used to select the method
#' @param ... further arguments to be passed
#' @aliases llr-class llr
#' @return returns a S3 object of class "llr" containing
#' \itemize{
#'    \item {\code{x} sorted numeric vector or matrix of \code{xpred}.}
#'    \item {\code{fitted} estimated values corresponding to \code{x}.}
#' } 
#' @export
llr <- function(x, ...)  UseMethod("llr")

#' @rdname llr 
#' @method llr default
#' @param y a numeric vector of y data corresponding to \code{x}.
#' @param xpred a numeric vector or matrix of same dimension as \code{x}.
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic".
#' @param bw a numeric vector or single number of same dimension as \code{x}.
#' @param weight a numeric vector of \code{length(x)} for weight of each data point.
#' @param kdtree boolean flag: If \code{TRUE}, kdtree is used for computation of local linear regression.
#' @param approx boolean flag: If \code{TRUE}, kdtree approximation is used . Only used when \code{kdtree = TRUE}. 
#' @param epsilon margin of error allowed for llr approximation using kdtree. Only used when both \code{kdtree = TRUE} and \code{approx = TRUE}.
#' @param N_min minimum number of points stored in the kd-tree. Only used when both \code{kdtree = TRUE} and \code{approx = TRUE}.
#' @examples 
#' n <- 1000
#' x <- seq(0,10,length.out = n)
#' x1 <- rnorm(n,0,0.2)
#' y <- sin(x) + x1
#' w <- rep(1/n, n)
#' binned <- bin(x, y, bins=400, w)
#' ## local linear regression for exact without kdtree
#' llr_exact <- llr(x, y, x, bw =0.2, weight = w)
#' ## local linear regression for kdtree exact
#' llr_kdexact <- llr(x, y, x, bw = 0.2, weight = w, kdtree = TRUE)
#' ## local linear regression for kdtree approximation
#' llr_kdapprox <- llr(x, y, x, bw = 0.2, weight = w, kdtree = TRUE, approx = TRUE)
#' ## local linear regression for data after binning.
#' llr_bin <- llr(binned, x , bw = 0.2)
#' @export 
llr.default<- function(x, y, xpred, kernel = "epanechnikov", bw, weight, kdtree = FALSE, approx = FALSE,
                        epsilon = 0.05, N_min =1, ...) { 
  
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
    
  x <- as.matrix(x) 
  y <- as.numeric(y)
  xpred <- as.matrix(xpred)
  
  scale <- apply(x, 2, FUN = normalize)
  bw <- bw * scale
  
  xy <- cbind(x, y, weight)
  xy <- xy[order(xy[,1]),]
  xpred <- as.matrix(xpred[order(xpred[,1]),])
  x <- as.matrix(xy[,1:ncol(x)])
  y <- as.numeric(xy[,ncol(xy)-1])
  wt <- as.numeric(xy[,ncol(xy)])

  
  if (kdtree == TRUE){ 
    if(approx == FALSE){
      ypred <- llrt_cpp(x, y, xpred, wt, 1, kcode, epsilon, bw, N_min)
    }
    else if (approx == TRUE){ 
      ypred <- llrt_cpp(x, y, xpred, wt, 2, kcode, epsilon, bw, N_min)
    }
  }
  
  if (kdtree == FALSE){
    if (ncol(x) == 1) {
      ypred <- llr1d_cpp(x, y, xpred, kcode, bw, wt)
    }
    if (ncol(x) == 2){
      ypred <- llr2d_cpp(x, y, xpred, kcode, bw, wt)
    }
    if (ncol(x) >= 3) { 
      ypred <- llr_cpp(x, y, xpred, kcode, bw, wt)
    }
  }
  
  results <- list ("x" = xpred, "fitted" = ypred)
  class(results) <- "llr"
  return(results)
}

#' @rdname llr 
#' @method llr bin
#' @export 
llr.bin <- function(x, xpred, kernel = "epanechnikov", bw, ... ){
  if (!inherits(x, "bin"))
    stop("function only works for objects of class bin")
  
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

  x_ <- as.matrix(x$x) 
  y_ <- as.numeric(x$y)
  wt_ <- as.numeric(x$xweight)
  xy <- cbind(x_,y_, wt_)
  xy <- xy[order(xy[,1]),]

  xpred <- as.matrix(xpred)
  xpred <- xpred[order(xpred[,1]),]
  x_ <- as.matrix(xy[,1:ncol(x_)])
  y_ <- as.numeric(xy[,ncol(xy)-1])
  wt_ <- as.numeric(xy[,ncol(xy)])
  
  scale <- apply(x_, 2, FUN = normalize)
  bw <- bw * scale

  if (ncol(x_) == 1){
    ypred <- llr1d_cpp(x_, y_, xpred, kcode, bw, wt_)
  }
  if (ncol(x_) == 2) {
    ypred <- llr2d_cpp(x_, y_, xpred, kcode, bw, wt_)
  }
  results <- list ("x" = xpred, "fitted" = ypred)
  class(results) <- "llr"
  return (results)
}




