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
ll_bin1d <- function(X, Y, X_pred, kernel = 'epanechnikov', bw, bins = 400){

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

  scale <- max(X) - min(X)
  bw <- bw * scale
  y_pred <- bin1d(X, Y, X_pred, kcode, bw, bins)
  y_pred
}



