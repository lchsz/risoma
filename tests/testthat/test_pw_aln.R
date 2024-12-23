
test_that("Pairwise alignment works as expected", {
  sub_mat <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = 0)
  expect_equal(pw_aln("ACGTA", "GACGT", sub_mat), c("-ACGTA", "GACGT-"))
  expect_equal(pw_aln("TCGGACCAGGCTTCATTCCCC", "ACGGACCAGGCTTCATTCCCC", sub_mat),
               c("TCGGACCAGGCTTCATTCCCC", "ACGGACCAGGCTTCATTCCCC"))
})
