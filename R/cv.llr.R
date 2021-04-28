#' Bandwidth selection using k-fold cross validation 
#' @param x a numeric vector or matrix of x data
#' @param y a numeric vector of y data corresponding to \code{x} data
#' @param weight a numeric vector of \code{length(x)} for weight of each data point.
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic". "sigmoid" and "silverman".
#' @param kdtree boolean flag: If \code{TRUE}, a kdtree is used for computation of local linear regression.
#' @param approx boolean flag: if \code{TRUE} kdtree approximation is used.
#' @param epsilon margin of error allowed for llr approximation using kdtree. Only used when \code{kdtree = TRUE} and \code{approx = TRUE}.
#' @param N_min minimum number of points stored in the kd-tree. Only used when \code{kdtree = TRUE} and \code{approx} = TRUE.
#' @param bw a numeric vector or matrix of bandwidth considered for selection.  
#' @param k number of folds used for cross validation
#' @return returns a single numeric value or vector of \code{bw} that gives the smallest mean square error from cross validation.
#' @examples 
#' n <- 1000
#' x <- seq(0,10,length.out = n)
#' x1 <- rnorm(n,0,0.2)
#' y <- sin(x) + x1
#' w <- rep(1/n, n)
#' binned <- bin(x,y,bins=400, w)
#' bw <- seq(0.02, 0.3, by = 0.02)
#' ## Bandwidth selection of binned data
#' h_bin <- cv.llr(binned$x, binned$y, binned$weight, bw = bw)
#' ## Bandwidth selection of exact local linear regression
#' h_exact <- cv.llr(x, y, w , bw = bw)
#' ## Bandwidth selection of exact local linear regression with kdtree
#' h_kdexact <- cv.llr(x, y, w, kdtree = TRUE, approx = FALSE, bw = bw)
#' ## Bandwidth selection of approx local linear regression with kdtree
#' h_kdapprox <- cv.llr(x, y, w , kdtree = TRUE , approx = TRUE, bw = bw) 
#' @export
cv.llr <- function(x, y, weight, kernel = "epanechnikov", bw, 
                     kdtree = FALSE, approx = FALSE, epsilon = 0.05, N_min = 1, k = 5){
  
  x <- as.matrix(x)
  y <- as.numeric(y) 
  weight <- as.numeric(weight)
  bw <- as.matrix(bw)
  
  xy <- cbind.data.frame(x, y, weight)
  xy <- xy[sample(nrow(xy)),]
  folds <- cut(seq(1,nrow(xy)),breaks= k,labels=FALSE)

  MSE_opt <- -1
  h_opt <- -1
  TSE <- 0 
  for (j in 1:nrow(bw)){
    h <- bw[j,]
    SSE <- 0
    for (i in 1:k){
      testIndexes <- which(folds==i,arr.ind=TRUE)
      test <- xy[testIndexes, ]
      train <- xy[-testIndexes, ]
      test <- test[order(test[,1]), ]
      if (kdtree == FALSE){
        ypred <- llr(train[,1:ncol(x)], train$y, test[,1:ncol(x)], kernel, h, train$weight)
      }
      else { 
        ypred <- llr(train[,1:ncol(x)], train$y, test[,1:ncol(x)], kernel, h, train$weight, kdtree, approx, epsilon, N_min)
      }
      
    SE <- (test$y - ypred$fitted)^2
    SSE <- c(SSE, SE)
    }
    
    MSE <- mean(SSE)
    TSE <- cbind(TSE, MSE)
    
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