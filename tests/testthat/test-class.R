test_that("Check bin split is correct for 1D", {
  n <- 300
  x <- runif(n,0,10)
  x1 <- rnorm(n,0,0.2)
  y <- sin(x) + x1
  w <- rep(1/n, n)
  binned <- bin(x, y, 400, w)
  expect_is(binned, 'bin')
  llr <- llr(binned, xpred = x, bandwidth = 0.1)
  expect_is(llr, 'llr')
})

test_that("Check bin split is correct for 1D", {
  n <- 300
  x <- runif(n,0,10)
  x1 <- rnorm(n,0,0.2)
  y <- sin(x) + x1
  w <- rep(1/n, n)
  llr <- llr(x, y, x, weight = w, bandwidth = 0.1)
  expect_is(llr, 'llr')
})
