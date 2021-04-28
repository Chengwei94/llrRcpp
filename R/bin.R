#' Splitting the data into bins  
#' @param x numeric vector or matrix of 2 dimensons of x data 
#' @param y numeric vector of y data with length corresponding to \code{x}
#' @param bins number of bins to split the data
#' @param weight numeric vector corresponding to weight of each observation 
#' @return returns an S3 object of the class 'bin' containing
#' \describe{
#'    \item{x}{a numeric vector of x data of size \code{bins}}
#'    \item{y}{a numeric vector of y data corresponding to 'bin' class \code{x}}
#'    \item{weight}{weight corresponding to 'bin' class \code{x}}
#'}
#' @examples
#' n <- 1000
#' x <- seq(0,10,length.out = n)
#' x1 <- rnorm(n,0,0.2)
#' y <- sin(x) + x1
#' w <- rep(1/n, n)
#' binned <- bin(x,y,bins=400, w)
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

