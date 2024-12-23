
test_that("CIGAR works as expected", {
  expect_equal(get_cigar("-ACGTA", "GACGT-"), "1D4M1I")
  expect_equal(get_cigar("TCGGACCAGGCTTCATTCCCC", "ACGGACCAGGCTTCATTCCCC"), "1X20M")
})
