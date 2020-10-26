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
bw_cvkfold <- function(XY, method = c("approx", "exact"), kernel = "epanechnikov", epsilon =0.05, bw, N_min = 1, k = 5 ){
  
  method <- match.arg(method)
  
  #helper function
  partition_data<- function(XY,k){
    set.seed(123)
    XY <- data.frame(XY)
    if (k != nrow(XY)){
      XY <- dplyr::mutate(XY, myfolds = sample(1:k,
                                                size = nrow(XY),
                                                replace = TRUE))
    }
    else {
      XY <- dplyr::mutate(XY, myfolds = sample(1:k,
                                                size = nrow(XY),
                                                replace = FALSE))
    }
    XY
  }
  
  
  ll_predict <- function(XY, X_pred, method = c('approx', 'exact'), kernel = 'epanechnikov', epsilon = 0.05,
                         bw, N_min = 1){
    
    if (ncol(XY) != ncol(X_pred)+1) stop('Dimensions of predictors in XY must be same dimensions of predictors in X_pred')
    
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
  
  XY <- partition_data(XY,k)
  MSE_opt <- -1
  for (j in 1:nrow(bw)){
    h <- bw[j, ]
    MSE <- 0
    for (i in 1:k){
      train <- subset(XY, myfolds != i)
      train <- train[,1:(ncol(train)-1)]
      train <- as.matrix(train)
      test <- subset(XY, myfolds == i)
      Y_val <- test[,ncol(test)-1]
      test <- test[,1:(ncol(test)-2)] # 2 (1 for myfolds, 1 for Y)
      test <- as.matrix(test)
      Y_val_predict <- ll_predict(train, test, method, kernel, epsilon, h, N_min)
      SE <- (Y_val - Y_val_predict)^2
      MSE <- MSE + mean(SE)
    }
    if(MSE_opt == -1 && MSE> 0){
      MSE_opt <- MSE
      h_opt = h
    }
    else if (MSE <= MSE_opt) {
      MSE_opt <- MSE
      h_opt <- h
    }
  }
  h_opt
}