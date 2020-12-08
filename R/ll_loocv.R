#' @title LOOCV bandwidth selection
#' 
#' @description choosing the best bandwidth through LOOCV 
#' 
#' @param XY  a matrix of the estimators with the observed values in the last colum 
#' @param method estimation method 
#' @param kernel type of kernel used 
#' @param epsilon level of error allowed
#' @param bw bandwidth for estimation of local linear regression
#' @param N_min number of points stored in the leaf of the tree (only works for approx)
#' 
#' @export 
bw.gcv <- function(X, Y, method = c('approx', 'exact'), kernel = 'epanechnikov',
                      epsilon = 0.05, bw, N_min = 1){
  
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
  scale <- apply(X, 2, FUN = normalize)
  bw <- bw * scale
  h <- bw_loocv(X, Y, metd, kcode, epsilon, bw, N_min)
  h <- h/scale
}





