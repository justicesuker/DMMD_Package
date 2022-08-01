test_that("Sanity check of DMMD-i", {
  expect_error(DMMD_i(X1 = matrix(rnorm(12), nrow = 4), X2 = matrix(rnorm(12), nrow = 4), r1 = 2, r2 = 3, rc = 4), "The specified joint rank is not legal, please check.")
})
