#' Bandwidth selection using k-fold cross validation
#' @param x a numeric vector or matrix of x data
#' @param y a numeric vector of y data corresponding to \code{x} data
#' @param weight a numeric vector of \code{length(x)} for weight of each data point.
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic". "sigmoid" and "silverman".
#' @param kdtree boolean flag: If \code{TRUE}, a kdtree is used for computation of local linear regression.
#' @param approx boolean flag: if \code{TRUE} kdtree approximation is used.
#' @param epsilon margin of error allowed for llr approximation using kdtree. Only used when \code{kdtree = TRUE} and \code{approx = TRUE}.
#' @param N_min minimum number of points stored in the kd-tree. Only used when \code{kdtree = TRUE} and \code{approx} = TRUE.
#' @param bandwidth a numeric vector or matrix of bandwidth considered for selection.
#' @param k number of folds used for cross validation
#' @return returns a single numeric value or vector of \code{bandwidth} that gives the smallest mean square error from cross validation.
#' @examples
#' n <- 1000
#' x <- seq(0, 10, length.out = n)
#' x1 <- rnorm(n, 0, 0.2)
#' y <- sin(x) + x1
#' w <- rep(1 / n, n)
#' binned <- bin(x, y, bins = 400, w)
#' bandwidth <- seq(0.02, 0.3, by = 0.02)
#' ## Bandwidth selection of binned data
#' h_bin <- cv.llr(binned$x, binned$y, binned$weight, bandwidth = bandwidth)
#' ## Bandwidth selection of exact local linear regression
#' h_exact <- cv.llr(x, y, w, bandwidth = bandwidth)
#' ## Bandwidth selection of exact local linear regression with kdtree
#' h_kdexact <- cv.llr(x, y, w, kdtree = TRUE, approx = FALSE, bandwidth = bandwidth)
#' ## Bandwidth selection of approx local linear regression with kdtree
#' h_kdapprox <- cv.llr(x, y, w, kdtree = TRUE, approx = TRUE, bandwidth = bandwidth)
#' @export
cv.llr <- function(x, y, weight, kernel = "epanechnikov", bandwidth,
                   kdtree = FALSE, approx = FALSE, epsilon = 0.05, N_min = 1, k = 5) {
  
  kernel <- match.arg(kernel)
  #helper functions
  calc_sse <- function(h, data, kdtree, approx, epsilon) {
    
    sse_tot <- 0
    
    
    for (i in 0:(k - 1)) {
      train <- data[data$idx != i, ]
      test <- data[data$idx == i, ]
      test <- test[order(test[, 1]), ]

      if (!kdtree) {
        test$y_est <- llr(train[, 1:(ncol(data) - 3)], train$y, test[, 1:(ncol(data) - 3)],
          bandwidth = h, weight = train$wt, epsilon = epsilon
        )$fitted
      }
      else {
        if (approx) {
          test$y_est <- llr(train[, 1:(ncol(data) - 3)], train$y, test[, 1:(ncol(data) - 3)],
            kdtree = TRUE, approx = TRUE, bandwidth = h, weight = train$wt, epsilon = epsilon
          )$fitted
        }
        else {
          test$y_est <- llr(train[, 1:(ncol(data) - 3)], train$y, test[, 1:(ncol(data) - 3)],
            kdtree = TRUE, bandwidth = h, weight = train$wt, epsilon = epsilon
          )$fitted
        }
      }
      sse <- sum((test$y_est - test$y)^2)
      sse_tot <- sse_tot + sse
    }
    sse_tot
  }
  
  x <- as.matrix(x)
  y <- as.numeric(y)
  wt <- as.numeric(weight)
  bandwidth <- as.matrix(bandwidth)
  
  if (nrow(x) != length(y) || nrow(x) != length(wt)){
    stop('x, y and weight must have the same length')
  }
  
  if(ncol(x) != ncol(bandwidth)) {
    stop('x and bandwidth should have the same dimension')
  }

  n <- length(y)
  idx <- sample(n, n) %% k
  df <- data.frame(x, y, wt, idx)
  table <- data.frame(bandwidth)
  table$sse <- apply(table, 1, calc_sse, data = df, kdtree = kdtree, approx = approx, epsilon = epsilon)
  table1 <- table[is.finite(table$sse), ]
  h_opt <- table1[which.min(table1$sse), ]
  h_opt <- h_opt[, 1:ncol(bandwidth)]
  h_opt <- as.numeric(h_opt)
  
  if (length(h_opt) == 0){
    stop("Choose a larger range of bandwidth")
  }
  
  results <- list('bw_opt' = h_opt, 'bws' = bandwidth, 'SSE' = table$sse)
}
