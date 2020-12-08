#' @title Local linear regression prediction
#'
#' @description Predict values through local linear smoothing
#'
#' @param XY  a matrix of the estimators with the observed values in the last column
#' @param X_pred a vector of predictions for observed values
#' @param method estimation method
#' @param kernel type of kernel used
#' @param epsilon level of error allowed
#' @param bw bandwidth for estimation of local linear regression
#' @param N_min number of points stored in the leaf of the tree (only works for approx)
#'
#' @export
ll_2d <- function(X, Y, X_pred, kernel = 'epanechnikov',
                  bw){
  
  X_pred <- as.matrix(X_pred)
  X <- X[order(X[,1]), ]
  X_pred <- X_pred[order(X_pred[,1]), ]
  X_pred <- as.matrix(X_pred)
                   
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
  
  X <- as.matrix(X)
  scale <- apply(X, 2, FUN = normalize)
  bw <- bw * scale
  y_pred <- predict2dd(X, Y, X_pred, kcode, bw)
  y_pred
}
