eval_kernel <- function(z){ 
  if (abs(z) < 1){
    y = 3/4*(1-z^2)
  }
  else { 
    y = 0
  }
  y
}

kern_weight <- function(vec_x, h){
  vec_x = vec_x / h 
  y = sapply(vec_x,FUN = eval_kernel) / h 
  y = prod(y)
  y
}

ll_R <- function(XY_mat,h){ 
  X_ones = rep(1,nrow(XY_mat))
  X_mat = XY_mat[,1:ncol(XY_mat)-1]
  Y_mat = XY_mat[,ncol(XY_mat)]
  W = diag(nrow(XY_mat))
  y_pred = rep(0,nrow(XY_mat))
  X_mat = as.matrix(X_mat)
  Y_mat = as.matrix(Y_mat)
  for(i in 1:nrow(XY_mat)){
    X_mat_i = sweep(X_mat,2, X_mat[i,])
    X_mat_i = cbind(X_ones, X_mat_i)
    for(j in 1:nrow(XY_mat)){
      W[j,j] = kern_weight(X_mat_i[j,2:ncol(X_mat_i)],h)
    }
    y_pred[i] = solve(t(X_mat_i) %*% W %*% X_mat_i, t(X_mat_i)%*% W %*%Y_mat)[1]    
  }
  y_pred
}

test_that("Check output of exact",{ 
  n = 300
  x1 = runif(n,0,0.5)
  x2 = runif(n,0,0.3)
  y = sin(x1)+x2
  matrix1 = cbind(x1,x2,y)
  bw = 0.2
  h = c(0.2,0.2)
  expect_equal(ll_R(matrix1, bw), ll(matrix1, method = 'exact', bw = h), tolerance = 0.01)
})

test_that("Check output of approx",{ 
  n = 300
  x1 = runif(n,0,0.5)
  x2 = runif(n,0,0.3)
  y = sin(x1)+x2
  matrix1 = cbind(x1,x2,y)
  bw = 0.2
  h = c(0.2,0.2)
  expect_equal(ll_R(matrix1, bw), ll(matrix1, method = 'approx', epsilon = 0, bw=h), tolerance = 0.01)
})

