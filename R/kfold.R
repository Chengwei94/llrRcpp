#' @title kfold cross valdation
#'
#' @description choosing the best bandwidth through kfold cross validation
#'
#' @param XY  a matrix of the estimators with the observed values in the last colum
#' @param method estimation method
#' @param kernel type of kernel used
#' @param epsilon level of error allowed
#' @param bw bandwidth for estimation of local linear regression
#' @param N_min number of points stored in the leaf of the tree (only works for approx)
#' @param k number of folds for cv
#'
#' @importFrom dplyr mutate
#'
#' @export
ll.cv <- function(X, Y, method = c("approx", "exact"), kernel = "epanechnikov", epsilon = 0.05, bw, N_min = 1, k = 5){
  
  method <- match.arg(method)

  #helper function  
  ll_predict <- function(X, Y, X_pred, method = c('approx', 'exact'), kernel = 'epanechnikov', epsilon = 0.05,
                         bw, N_min = 1){
    
    method <- match.arg(method)
    
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
    
    X <- as.matrix(X)
    X_pred <- as.matrix(X_pred)
    scale <- apply(X, 2, FUN = normalize)
    bw <- bw * scale
    predict_values <- predict(X, Y, X_pred, metd, kcode, epsilon, bw, N_min)
    predict_values
  }
  
  XY <- cbind.data.frame(X, Y)
  XY <- XY[sample(nrow(XY)),]
  folds <- cut(seq(1, nrow(XY)), breaks= k, labels=FALSE)
  
  MSE_opt <- -1
  h_opt <- -1
  bw <- as.matrix(bw)
  
  for (j in 1:nrow(bw)){
    h <- bw[j]
    SSE <- 0
    for (i in 1:k){
      testIndexes <- which(folds==i,arr.ind=TRUE)
      test <- XY[testIndexes, ]
      train <- XY[-testIndexes, ]
      y_pred <- ll_predict(train$X, train$Y, test$X, method, kernel, epsilon, h, N_min)
      SE <- (test$Y - y_pred)^2
      SSE <- c(SSE, SE)
    }
    MSE <- mean(SSE)
    if(MSE_opt == -1 && MSE> 0){
      MSE_opt <- MSE
      h_opt <- h
    }
    else if (MSE <= MSE_opt) {
      MSE_opt <- MSE
      h_opt <- h
    }
  }
  h_opt
}
