test_that("Check that the estimates for cv and gcv are reasonable", {
  x <- mcycle$times
  y <- mcycle$accel
  d <- data.frame(x,y)
  bw <- seq(0.01, 0.2, by = 0.05)
  w <- rep(1/length(x), length(x))
  loocv <- loocv.llr(x, y, w, bandwidth = bw)
  cv <- cv.llr(x, y, w, k = 5, bandwidth = bw)
  expect_equal(loocv$SSE[1], Inf)
  expect_equal(cv$SSE[1], Inf)
})

test_that("Test that small bw_opt will give a warning",{
  x <- mcycle$times
  y <- mcycle$accel
  d <- data.frame(x,y)
  bw <- seq(0.001, 0.002)
  w <- rep(1/length(x), length(x))
  expect_error(loocv <- loocv.llr(x, y, w, bandwidth = bw), "Choose a larger range of bandwidth")
  expect_error(cv <- cv.llr(x, y, w, k = 5, bandwidth = bw), "Choose a larger range of bandwidth")
  expect_error(autoloocv.llr(x, y, w, bw_range = c(0.001, 0.002), 
               control = list(numPopulation=5, maxIter=10, Vmax=2, ci=1.49445, cg=1.49445, w=0.729)),
               "Convergence failed, choose a larger bandwidth or use workers")
})

test_that("cv using kdtree and non-kdtree is the same", { 
  x <- mcycle$times
  y <- mcycle$accel
  d <- data.frame(x,y)
  bw <- seq(0.01, 0.15, by = 0.01)
  w <- rep(1/length(x), length(x))
  set.seed(1)
  cv_tree <- cv.llr(x, y, w, k =5, bandwidth = bw, kdtree = TRUE, approx = FALSE)
  set.seed(1)
  cv <- cv.llr(x, y, w, k = 5, bandwidth = bw, kdtree = TRUE)
  expect_equal(cv_tree, cv)
})

