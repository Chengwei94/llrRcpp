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
ll_1d.cv <- function(X, Y, kernel = "epanechnikov", bw, k = 5 ){

  #helper function
  ll_1d <- function(X, Y, X_pred, kernel = 'epanechnikov',
                    bw){
    
    
    X <- sort(X)
    X_pred <- sort(X_pred)
    
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
    
    scale <- max(X) - min(X)
    bw <- bw * scale
    y_pred <- predict1dd(X, Y, X_pred, kcode, bw)
    y_pred
  }
  
  XY <- cbind.data.frame(X, Y)
  XY <- XY[sample(nrow(XY)),]
  folds <- cut(seq(1,nrow(XY)),breaks= k,labels=FALSE)
  
  MSE_opt <- -1
  h_opt <- -1
  for (j in 1:length(bw)){
    h <- bw[j]
    SSE <- 0
    for (i in 1:k){
      testIndexes <- which(folds==i,arr.ind=TRUE)
      test <- XY[testIndexes, ]
      train <- XY[-testIndexes, ]
      y_pred <- ll_1d(train$X, train$Y, test$X, kernel, h)
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