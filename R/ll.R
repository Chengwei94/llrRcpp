
#' @title Local linear regression
#' 
#' @description Multivariate local linear regression 
#' 
#' @param XY  a matrix of the estimators with the observed values in the last colum 
#' @param method estimation method 
#' @param kernel type of kernel used 
#' @param epsilon level of error allowed
#' @param bw bandwidth for estimation of local linear regression
#' @param N_min number of points stored in the leaf of the tree (only works for approx)
#' 
#' @return h 
#' @export 
ll <- function(XY, method = c("approx", "exact"), kernel = "epanechnikov", 
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
  
  scale <- apply(as.matrix(XY[,1:ncol(XY)-1]),2, FUN = normalize)
  bw <- bw * scale
  Xpred <- loclin(XY, metd, kcode, epsilon, bw, N_min)
  Xpred
}