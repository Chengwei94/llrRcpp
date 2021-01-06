#' Splitting 1D/2D data into bins 
#' @param x numeric vector or matrix of 2 dimensons of x data 
#' @param y numeric vector of y data corresponding to \code{x}
#' @param bins number of bins to split the data into
#' @param wt weight of each observation 
#' @return returns an S3 object of the class 'bin' containing
#' \itemize{
#'    \item {\code{x} a numeric vector of x data after binned}
#'    \item {\code{y} a numeric vector of y data corresponding to binned \code{x}}
#'    \item {\code{weight} weight corresponding to binned \code{x}}
#' }
#' @examples
#' n <- 1000
#' x <- seq(0,10,length.out = n)
#' x1 <- rnorm(n,0,0.2)
#' y <- sin(x) + x1
#' w <- rep(1/n, n)
#' z <- bin(x,y,bins=400, w)
#' @export
bin <-  function(x, y, bins = 400, weight){  
  x <- as.matrix(x)
  y <- as.numeric(y)
  if ((ncol(x) > 2))  stop('x must be of dimension 1 or 2')
  if (!identical(nrow(x),length(y))) stop('x and y must have the same length')
  
  if (ncol(x) == 1){
    r <- bin1d_cpp(x, y, bins, weight)
  }
  
  if (ncol(x) == 2){
    r <- bin2d_cpp(x,y, bins, weight)
  }
  
  class(r) <- "bin" 
  
  index <- which(r$weight != 0)
  r$x <- as.matrix(r$x)
  r$x <- r$x[index,]
  r$y <- r$y[index]
  r$weight <- r$weight[index]
  
  return (r)
}

