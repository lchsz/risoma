
test_that("TSI works as expected", {
  expr <- c(0, 8, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0)
  expect_equal(calc_tsi(expr), 0.95, tolerance = 0.01)
})
