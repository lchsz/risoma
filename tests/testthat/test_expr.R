
test_that("calc_one_tsi example code runs correctly", {
  # Exclusive expression
  exp1 <- c(0, 0, 10, 0)
  expect_equal(calc_one_tsi(exp1), 1)

  # Uniform expression
  exp2 <- c(5, 5, 5, 5)
  expect_equal(calc_one_tsi(exp2), 0)

  # Mixed expression
  exp3 <- c(2, 8, 4, 6)
  expected_tsi <- sum(1 - exp3 / max(exp3)) / (length(exp3) - 1)
  expect_equal(calc_one_tsi(exp3), expected_tsi)

  # Invalid input (negative values)
  exp4 <- c(-1, 3, 5)
  expect_error(calc_one_tsi(exp4), "Expression values must be non-negative")

  # Invalid input (all zeros)
  exp5 <- c(0, 0, 0)
  expect_error(calc_one_tsi(exp5), "At least one value must be non-negative")

  # Invalid input (single value)
  exp6 <- c(10)
  expect_error(calc_one_tsi(exp6))
})
