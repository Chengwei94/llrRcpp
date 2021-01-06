test_that("Check bin split is correct for 1D", {
  n <- 300
  x <- runif(n,0,10)
  x1 <- rnorm(n,0,0.2)
  y <- sin(x) + x1
  w <- rep(1/n, n)
  binned <- bin(x, y, 400, w)
  expect_equal(sum(binned$w), 1)
  expect_equal(max(binned$x), max(x))
  expect_equal(min(binned$x), min(x))
  expect_equal(sum(binned$y * binned$weight), sum(w*y)) 
})

test_that("Check bin split is correct for 2D", {
  n <- 1000
  x1 <- runif(n, 0 ,10)
  x2 <- runif(n, 0, 10)
  x3 <- rnorm(n, 0, 0.3)
  y = sin(x1) + cos(x2) + x3
  x <- cbind(x1, x2)
  w <- rep(1/n, n)
  binned <- bin(x, y, 400 , w)
  expect_equal(sum(binned$w), 1)
  expect_equal(max(binned$x), max(x))
  expect_equal(min(binned$x), min(x))
  expect_equal(sum(binned$y * binned$weight), sum(w*y)) 
})

test_that("Check bin split is correct for 1D under dif weights", {
  n <- 300
  x <- runif(n,0,10)
  x1 <- rnorm(n,0,0.2)
  y <- sin(x) + x1
  rand <- runif(n, 0, 1)
  w <- rand/sum(rand)
  binned <- bin(x, y, 400, w)
  expect_equal(sum(binned$w), 1)
  expect_equal(max(binned$x), max(x))
  expect_equal(min(binned$x), min(x))
  expect_equal(sum(binned$y * binned$weight), sum(w*y)) 
})


test_that("Check bin split is correct for 2D under rand weight", {
  n <- 1000
  x1 <- runif(n, 0 ,10)
  x2 <- runif(n, 0, 10)
  x3 <- rnorm(n, 0, 0.3)
  y = sin(x1) + cos(x2) + x3
  x <- cbind(x1, x2)
  rand <- runif(n, 0, 1)
  w <- rand/sum(rand)
  binned <- bin(x, y, 400 , w)
  expect_equal(sum(binned$w), 1)
  expect_equal(max(binned$x), max(x))
  expect_equal(min(binned$x), min(x))
  expect_equal(sum(binned$y * binned$weight), sum(w*y)) 
})
