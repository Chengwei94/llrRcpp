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
ll_predict <- function(XY, X_pred, method = c('approx', 'exact'), kernel = 'epanechnikov', epsilon = 0.05,
                       bw, N_min = 1){
  
  method <- match.arg(method)
  if (ncol(XY) != ncol(X_pred) + 1) stop('Dimensions of predictors in XY must be same dimensions of predictors in X_pred')
  
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
  switch(method,
         approx = {metd <- 2},
         exact = {metd <- 1}
  )
  #helper functions
  normalize <- function(x)
  {
    return(max(x)- min(x))
  }
  
  scale <- apply(as.matrix(XY[,1:ncol(XY)-1]),2, FUN = normalize)
  bw <- bw * scale
  predict_values <- predict(XY, X_pred, metd, kcode, epsilon, bw, N_min)
  predict_values
}