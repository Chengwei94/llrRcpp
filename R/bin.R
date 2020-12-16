bin <-  function(x, y, bins = 400){  
  x <- as.matrix(x)
  y <- as.numeric(y)
  if ((ncol(x) > 2))  stop('x must be of dimension 1 or 2')
  if (!identical(nrow(x),length(y))) stop('x and y must have the same length')
  if (ncol(x) == 1){
  r <- bin1d_cpp(x, y, bins)
  }
  if (ncol(x) == 2){
  # to be filled 
  }
  class(r) <- "binned" 
  return (r)
}

