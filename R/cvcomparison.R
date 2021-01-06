library(microbenchmark)
library(llrRcpp)

n <- 5000
variance <- 0.3
x <- runif(n, 0, 10)
x1 <- rnorm(n, 0, variance)
y <- sin(x) + x1
w <- rep(1/n, n)
d <- data.frame(x = runif(n))

microbenchmark(
  cv_exact = cv.llr(x, y, w, kdtree = TRUE),
  cv_approx = cv.llr(x,y, w, kdtree = TRUE, approx = FALSE),
  cv_1D = cv.llr(x, y, w),
  gcv_exact = gcv.llr(x, y, w),
  gcv_approx = gcv.llr(x,y, w, approx = TRUE),
  times = 10
)
    
