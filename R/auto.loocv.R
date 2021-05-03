#' Auto selection of bandwidth using metaheuristic algorithms 
#'
#' @param x a numeric vector or matrix of x data
#' @param y a numeric vector of y data corresponding to \code{x} data
#' @param weight a numeric vector of \code{length(x)} for weight of each data point
#' @param kernel kernel type used; supported are 'epanechnikov', "rectangular", "triangular", "quartic", "triweight", "tricube", "cosine", "gauss", "logistic". "sigmoid" and "silverman".
#' @param bw_range range of bandwidth for the metaheuristic algorithm to search 
#' @param approx boolean flag: if \code{true} kdtree approximation is used.
#' @param epsilon  margin of error allowed for llr approximation using kdtree. Only used when \code{approx = TRUE}.
#' @param control control parameter from R package  metaheuristicOpt
#' @param algorithm algorithm parameter from R pacakage  metaheuristicOpt
#' @examples  
#' \dontrun {
#' n <- 1000
#' x <- runif(n, 0 , 1)
#' y <- sin(x) + rnorm(n, 0, 0.2)
#' w <- rep(1/n , n)
#' ## Optimization using Grasshopper Optimisation Algorithm
#' algorithm <- 'GOA'
#' control <- list(numPopulation = 10, maxIter = 100)
#' h <- autoloocv.llr(x, y, w, control = control , algorithm = algorithm)
#' }
#' @export
autoloocv.llr <- function(x, y, weight, kernel = "epanechnikov", bw_range = c(0.01, 0.5), approx = FALSE, epsilon = 0.05,
                       control = list(numPopulation=15, maxIter=75, Vmax=2, ci=1.49445, cg=1.49445, w=0.729), 
                       algorithm = "PSO", seed = 1) {

  x <- as.matrix(x)
  y <- as.numeric(y)
  weight <- as.numeric(weight)
    
  loss_function <- function(h) {
    loss <- loocv.llr_nowarn(x, y, weight = weight, bandwidth = h, approx = approx, epsilon = epsilon)
    return (loss$SSE)
  }
  
  control <- control
  numvar <- ncol(x)
  rangevar <- matrix(bw_range, nrow =2)
  
  best_SSE <- metaheuristicOpt::metaOpt(loss_function, optimType="MIN", algorithm=algorithm, numvar,
                           rangevar, control, seed)
  
  best_SSE
  return (list('bw_opt' = best_SSE$result, 'SSE' = best_SSE$optimumValue))
}

#Internal function without warning for autocv.llr. Exactly the same function as loocv.llr but without 
#the warning. 
loocv.llr_nowarn <- function(x, y, weight, kernel = "epanechnikov", approx = FALSE, epsilon = 0.05,
                             N_min = 1, bandwidth) {
  
  kernel <- match.arg(kernel)
  
  switch(kernel,
         epanechnikov = {
           kcode <- 1
         },
         rectangular = {
           kcode <- 2
         },
         triangular = {
           kcode <- 3
         },
         quartic = {
           kcode <- 4
         },
         triweight = {
           kcode <- 5
         },
         tricube = {
           kcode <- 6
         },
         cosine = {
           kcode <- 7
         },
         gauss = {
           kcode <- 21
         },
         logistic = {
           kcode <- 22
         },
         sigmoid = {
           kcode <- 23
         },
         silverman = {
           kcode <- 24
         }
  )
  
  # helper functions
  normalize <- function(x) {
    return(max(x) - min(x))
  }
  
  x <- as.matrix(x)
  y <- as.numeric(y)
  wt <- as.numeric(weight)
  bandwidth <- data.frame(bandwidth)
  
  if (ncol(x) > 1 && ncol(bandwidth) == 1) {
    bandwidth <- t(bandwidth)
  }
  
  if (nrow(x) != length(y) || nrow(x) != length(wt)){
    stop('x, y and weight must have the same length')
  }
  
  if(ncol(x) != ncol(bandwidth)) {
    stop('x and bandwidth should have the same dimension')
  }
  
  scale <- apply(x, 2, FUN = normalize)
  scale <- as.numeric(round(scale , 3))
  bandwidth1 <- mapply('*', bandwidth, scale)
  bandwidth1 <- as.matrix(bandwidth1)
  
  if (ncol(x) > 1 && ncol(bandwidth1) == 1) {
    bandwidth1 <- t(bandwidth1)
  }
  
  if (approx == FALSE) {
    results <- tgcv_cpp(x, y, wt, 1, kcode, epsilon, bandwidth1, N_min)
  }
  else {
    results <- tgcv_cpp(x, y, wt, 2, kcode, epsilon, bandwidth1, N_min)
  }
  
  results$bw_opt <- results$bw_opt / scale
  results$bandwidth <- bandwidth 
  
  results
}



