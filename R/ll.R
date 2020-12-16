loclin <- function(x, ...)  UseMethod("llr")

#' @title Local linear regression
#' 
#' @description Multivariate local linear regression 
#' 
#' @param X  a matrix of the estimators with the observed values in the last column
#' @param Y  a matrix of the estimators with the observed values in the last column
#' @param method estimation method 
#' @param kernel type of kernel used 
#' @param epsilon level of error allowed
#' @param bw bandwidth for estimation of local linear regression
#' @param N_min number of points stored in the leaf of the tree (only works for approx)
#' 
#' @return h 
#' @export 

llr.default<- function(x, y, xpred, kernel = "epanechnikov", bw, weight){ 
  
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
  scale <- apply(x, 2, FUN = normalize)
  bw <- bw * scale
  
  if (ncol(x) == 1) {
    ypred <- llr1d_cpp(x, y, xpred, kcode, bw, weight)
  }
  else if (ncol(x) == 2){
    ypred <- llr2d_cpp(x, y, xpred, kcode, bw, weight)
  }
  else{ 
    ypred <- llr_cpp(x, y, xpred, kcode, bw, weight)
  }
  class(ypred) <- "llr"
  
  return(ypred)
}

llr.binned <- function(x, xpred, kernel = "epanechnikov", bw){
  if (!inherits(x, "binned"))
    stop("function only works for objects of class binned")
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
  
  x$x <- as.matrix(x$x)
  x$y <- as.numeric(x$y)
  
  scale <- apply(x$x, 2, FUN = normalize)
  bw <- bw * scale
  
  if (ncol(x$x) == 1){
    ypred <- llr1d_cpp (x$x, x$y, xpred, kcode, bw, x$xweight)
  }
  if (ncol(x$x) == 2) {
    # to be field up with bin2d
  }
  return (ypred)
}
# ll <- function(X, Y, method = c("approx", "exact"), kernel = "epanechnikov", 
#                epsilon = 0.05, bw, N_min = 1){
#   
#   method <- match.arg(method) 
#   
#   switch(kernel, 
#          epanechnikov = {kcode <- 1}, 
#          rectangular = {kcode <-2},
#          triangular = {kcode <- 3}, 
#          quartic = {kcode <- 4}, 
#          triweight = {kcode <- 5}, 
#          tricube = {kcode <- 6}, 
#          cosine = {kcode <- 7}, 
#          gauss = {kcode <- 21}, 
#          logistic = {kcode<- 22},
#          sigmoid = {kcode <- 23}, 
#          silverman = {kcode <- 24}
#   )
#   switch(method, 
#          approx = {metd <- 2}, 
#          exact = {metd <- 1}
#   )
#   
#   #helper functions
#   normalize <- function(x)
#   {
#     return(max(x)- min(x))
#   } 
#   
#   scale <- apply(as.matrix(X) ,2, FUN = normalize)
#   bw <- bw * scale
#   X <- as.matrix(X)
#   Xpred <- loclin(X, Y, metd, kcode, epsilon, bw, N_min)
#   Xpred
# }