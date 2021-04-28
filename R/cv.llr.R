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
  
  calc_sse <- function(h, data, kdtree, approx) {
    mse_tot <- 0
    for (i in 0:(k - 1)) {
      train <- data[data$idx != i, ]
      test <- data[data$idx == i, ]
      test <- test[order(test[,1]),]
      
      if (!kdtree){
          test$y_est <- llr(train[, 1:(ncol(data)-3)], train$y, test[, 1:(ncol(data)-3)], 
                            bw = h, weight = train$weight)$fitted
      }
      else {
        if (approx){
          test$y_est <- llr(train[, 1:(ncol(data)-3)], train$y, test[, 1:(ncol(data)-3)], 
                            kdtree = TRUE, bw = h, weight = train$weight)$fitted
        } 
        else {
          test$y_est <- llr(train[, 1:(ncol(data)-3)], train$y, test[, 1:(ncol(data)-3)], 
                            kdtree = TRUE, approx = TRUE, bw = h, weight = train$weight)$fitted
        }
      }
      
      sse <- sum((test$y_est - test$y)^2, na.rm = TRUE)
      n <- sum(is.finite(test$y_est))
      mse <- sse / n
      mse_tot <- mse_tot + mse
    }
    mse_tot
  }
  x <- as.matrix(x)
  y <- as.numeric(y) 
  weight <- as.numeric(weight)
  bw <- as.matrix(bw)
  
  n <- length(y)
  idx <- sample(n, n) %% k
  df <- data.frame(x, y, weight, idx)
  table <- data.frame(bw)
  table$mse <- apply(table, 1, calc_sse, data = df, kdtree = kdtree , approx = approx)

  h_opt <- table[which.min(table$mse),]
  h_opt <- h_opt[,1:ncol(bw)]
}