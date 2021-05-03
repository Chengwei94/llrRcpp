test_that("Check openmp setthread working", {
  expect_equal(set_num_threads(5), 5)
})
