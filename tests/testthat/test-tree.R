test_that("Check kdtree_exact estimation is same as kdtree_approx(when epsilon = 0) and exact method", {
  n <- 1000
  x <- runif(n,0,10)
  xnew <- seq(0.1, 9.9, by = 0.01)
  x1 <- rnorm(n,0,0.2)
  y <- sin(x) + x1
  w <- rep(1/n, n)
  llr_exact <- llr(x, y, xnew, bandwidth = 0.1, weight = w)
  llr_kdexact <- llr(x, y, xnew, bandwidth = 0.1, weight = w, kdtree = TRUE)
  llr_kdapprox <- llr(x, y, xnew, bandwidth = 0.1, weight = w, kdtree = TRUE, approx = TRUE, epsilon = 0)
  expect_equal(llr_exact$fitted, llr_kdexact$fitted)
  expect_equal(llr_exact$fitted, llr_kdapprox$fitted)
  expect_equal(llr_kdexact$fitted, llr_kdapprox$fitted)
})

test_that("Check kdtree_exact estimation is same as kdtree_approx(when epsilon = 0) and exact method
          under different weights", {
  n <- 1000
  x <- runif(n,0,10)
  xnew <- seq(0.1, 9.9, by = 0.01)
  x1 <- rnorm(n,0,0.2)
  y <- sin(x) + x1
  rand <- runif(n, 0, 1)
  w <- rand/sum(rand)
  llr_exact <- llr(x, y, xnew, bandwidth = 0.1, weight = w)
  llr_kdexact <- llr(x, y, xnew, bandwidth = 0.1, weight = w, kdtree = TRUE)
  llr_kdapprox <- llr(x, y, xnew, bandwidth = 0.1, weight = w, kdtree = TRUE, approx = TRUE, epsilon = 0)
  expect_equal(llr_exact$fitted, llr_kdexact$fitted)
  expect_equal(llr_exact$fitted, llr_kdapprox$fitted)
  expect_equal(llr_kdexact$fitted, llr_kdapprox$fitted)
})

test_that("Check kdtree_exact estimation is same as kdtree_approx(when epsilon = 0)
          and exact method for 2D", {
  n <- 1000
  x1 <- runif(n, 0, 10)
  x2 <- runif(n, 0, 10)
  x3 <- rnorm(n, 0, 0.2)
  y <- sin(x1) + cos(x2) + x3
  x <- cbind(x1, x2)
  xnew1 <- seq(0.1, 9.9, by = 0.2)
  xnew2 <- seq(0.1, 9.9, by = 0.2)
  xnew <- expand.grid(xnew1, xnew2)
  w <- rep(1/n, n)
  h <- c(0.2,0.2)
  llr_exact <- llr(x, y, xnew, bandwidth = h, weight = w)
  llr_kdexact <- llr(x, y, xnew, bandwidth = h, weight = w, kdtree = TRUE)
  llr_kdapprox <- llr(x, y, xnew, bandwidth = h, weight = w, kdtree = TRUE, approx = TRUE, epsilon = 0)
  expect_equal(llr_exact$fitted, llr_kdexact$fitted)
  expect_equal(llr_exact$fitted, llr_kdapprox$fitted)
  expect_equal(llr_kdexact$fitted, llr_kdapprox$fitted)
})

test_that("Check kdtree_exact estimation is same as kdtree_approx(when epsilon = 0)
          and exact method for 3D", {
  n <- 1000
  x1 <- runif(n, 0, 10)
  x2 <- runif(n, 0, 10)
  x3 <- runif(n, 0 ,10)
  x4 <- rnorm(n, 0, 0.2)
  y <- sin(x1) + cos(x2) + sin(x3/2) + x4
  x <- cbind(x1, x2, x3)
  xnew1 <- seq(0.1, 9.9, by = 1)
  xnew2 <- seq(0.1, 9.9, by = 1)
  xnew3 <- seq(0.1, 9.9, by = 1)
  xnew <- expand.grid(xnew1, xnew2, xnew3)
  w <- rep(1/n, n)
  h <- c(0.2, 0.2, 0.2)
  llr_exact <- llr(x, y, xnew, bandwidth = h, weight = w)
  llr_kdexact <- llr(x, y, xnew, bandwidth = h, weight = w, kdtree = TRUE)
  llr_kdapprox <- llr(x, y, xnew, bandwidth = h, weight = w, kdtree = TRUE, approx = TRUE, epsilon = 0)
  expect_equal(llr_exact$fitted, llr_kdexact$fitted)
  expect_equal(llr_exact$fitted, llr_kdapprox$fitted)
  expect_equal(llr_kdexact$fitted, llr_kdapprox$fitted)
})



