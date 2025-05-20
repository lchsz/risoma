data("isomirs")

# Get isomiR information ======================================================

## get_isoform_num -----
test_that("get_isoform_num calculates the number of isoforms correctly", {
  # Calculate the number of isoforms per reference miRNA
  isoform_counts <- get_isoform_num(isomirs)

  # Check if the output is a named numeric vector
  expect_true(is.numeric(isoform_counts))
  expect_true(!is.null(names(isoform_counts)))

  # Check if the counts are correct
  expect_equal(isoform_counts[["vvi-miR172c"]], 2)
  expect_equal(isoform_counts[["vvi-miR156f"]], 3)

  # Invalid input with NULL IsomirDataSet
  expect_error(get_isoform_num(NULL))

  # Invalid input with empty ref2isomir slot
  isomir_dataset <- new("IsomirDataSet", ref2isomir = list())
  expect_equal(length(get_isoform_num(isomir_dataset)), 0)
})

## get_isomir_by_ref -------------
test_that("get_isomir_by_ref extracts isomiR data correctly", {
  # Extract isomiR data for an existing reference miRNA
  result <- get_isomir_by_ref(isomirs, ref = "vvi-miR172c")
  expect_true(inherits(result, "Isomir"))
  expect_equal(result@mature_id, "vvi-miR172c")
  expect_equal(result@read_seqs[1], "GGAATCTTGATGATGCTG")
  expect_equal(result@dist, c(3, 2, 0))

  # Extract isomiR data for a non-existing reference miRNA
  result <- get_isomir_by_ref(isomirs, ref = "miR-3")
  # Output should be NULL for non-existing reference miRNAs
  expect_null(result)

  # Invalid input with NULL IsomirDataSet
  expect_error(get_isomir_by_ref(NULL, ref = "vvi-miR172c"))

  # Invalid input with NULL ref
  expect_error(get_isomir_by_ref(isomirs, ref = NULL))
})

## get_alignment_by_ref ------
test_that("get_alignment_by_ref extracts alignment data correctly", {
  # Extract alignment data for an existing reference miRNA
  result <- get_alignment_by_ref(isomirs, ref = "vvi-miR156f")
  expect_true(inherits(result, "DNAStringSet"))
  expect_equal(length(result), 5)

  # Extract alignment data for a non-existing reference miRNA
  result <- get_alignment_by_ref(isomirs, ref = "miR-2")
  expect_null(result)

  # Test Case 3: Invalid input with NULL IsomirDataSet
  expect_error(get_alignment_by_ref(NULL, ref = "vvi-miR156f"))

  # Test Case 4: Invalid input with NULL ref
  expect_error(get_alignment_by_ref(isomirs, ref = NULL))
})

# Expression Data ===========================================================

## calc_group_expr_num -----------------
test_that("calc_group_expr_num example code runs correctly", {
  # Check if calc_expr_num runs without errors
  expect_no_error({
    expr_num <- calc_group_expr_num(isomirs)
  })

  # Check if the output is a data frame
  expect_true(is.data.frame(expr_num))

  # Check if the data frame contains the expected columns
  expected_columns <- c("group", "total", "organ_specific", "common")
  expect_true(all(expected_columns %in% colnames(expr_num)))

  # Check if the total expressed isomiRs are non-negative
  expect_true(all(expr_num$total >= 0))

  # Check if the organ-specific isomiRs are non-negative
  expect_true(all(expr_num$organ_specific >= 0))

  # Check if the common isomiRs are non-negative
  expect_true(all(expr_num$common >= 0))

  # Check if the sum of organ-specific and common isomiRs is less than or equal
  # to total
  expect_true(all(expr_num$organ_specific + expr_num$common <= expr_num$total))
})

## get_ref_expr --------------
test_that("get_ref_expr extracts reference miRNA expression data correctly", {
  # Extract expression data for reference miRNAs
  result <- get_ref_expr(isomirs)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 44)
  expect_equal(rownames(result)[c(1, 2)], c("vvi-miR156b", "vvi-miR156f"))
  expect_equal(colnames(result)[1], "Berry_FS")

  # Check if the expression values are correct
  expect_equal(result["vvi-miR172c", "Flower_FB"], 7.5)
  expect_equal(result["vvi-miR156f", "Inflorescence_WD"], 252)

  # Invalid input with NULL IsomirDataSet
  expect_error(get_ref_expr(NULL))

  # Check handling of IsomirDataSet with empty ref2isomir slot
  isomir_dataset <- new("IsomirDataSet", ref2isomir = list())
  expect_error(get_ref_expr(isomir_dataset))
})

## get_expr_by_ref -------
test_that("get_expr_by_ref extracts expression data correctly", {
  # Extract expression data for an existing reference miRNA
  result <- get_expr_by_ref(isomirs, ref = "vvi-miR172c")
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 3.0)
  expect_equal(rownames(result)[1], "18M3D")
  expect_equal(colnames(result)[1], "Berry_FS")

  #  Extract expression data for a non-existing reference miRNA
  result <- get_expr_by_ref(isomirs, ref = "miR-3")
  expect_null(result)

  # Invalid input with NULL IsomirDataSet
  expect_error(get_expr_by_ref(NULL, ref = "vvi-miR172c"))

  # Invalid input with NULL ref
  expect_error(get_expr_by_ref(vvi_isomirs, ref = NULL))
})

# Calculate Tissue Specificity Index =================================

## calc_tsi -----
test_that("calc_tsi example code runs correctly", {
  # Check if calc_tsi runs without errors
  expect_no_error({
    tsi_results <- calc_tsi(isomirs)
  })

  # Check if the output is a data frame
  expect_true(is.data.frame(tsi_results))

  # Check if the data frame contains the expected columns
  expected_columns <- c("read_seq", "type", "tsi")
  expect_true(all(expected_columns %in% colnames(tsi_results)))

  # Check if TSI values are within the valid range [0, 1]
  expect_true(all(tsi_results$tsi >= 0 & tsi_results$tsi <= 1))

  # Check if the top 5 most specific isoforms are correctly ordered
  top_5 <- head(tsi_results[order(-tsi_results$tsi), ], 5)
  expect_true(all(diff(top_5$tsi) <= 0))
})

# DEG ====================================================================

test_that("get_deg_data returns correct structure and data", {
  # Test normal case: treatment = "Flower_FB", control = "Inflorescence_WD"
  result <- get_deg_data(isomirs, "Flower_FB", "Inflorescence_WD")

  # Validate output structure
  expect_type(result, "list")
  expect_named(result, c("expr", "sample_info"))
  expect_s3_class(result$expr, "data.frame")
  expect_s3_class(result$sample_info, "data.frame")

  # Verify sample_info filtering and ordering
  expect_equal(nrow(result$sample_info), 4)
  expect_equal(result$sample_info$name,
               c("SRR1528376", "SRR1528377", "SRR1528374", "SRR1528375"))
  expect_equal(levels(result$sample_info$group), c("Flower_FB", "Inflorescence_WD"))
  expect_equal(result$sample_info$group,
               factor(c("Flower_FB", "Flower_FB", "Inflorescence_WD", "Inflorescence_WD"),
                      levels = c("Flower_FB", "Inflorescence_WD")))

  # Verify expression matrix
  expect_equal(colnames(result$expr),
               c("SRR1528376", "SRR1528377", "SRR1528374", "SRR1528375"))
  expect_equal(rownames(result$expr)[1:2],
               c("AAGCTCAGGAGGGATAGCGCC", "AAGCTCAGGAGGGATAGCGCCA"))

  # Test reversed group order (treatment = "Inflorescence_WD", control = "Flower_FB")
  result_rev <- get_deg_data(isomirs, "Inflorescence_WD", "Flower_FB")
  expect_equal(result_rev$sample_info$group,
               factor(c("Inflorescence_WD", "Inflorescence_WD", "Flower_FB", "Flower_FB"),
                      levels = c("Inflorescence_WD", "Flower_FB")))
  expect_equal(result_rev$sample_info$name,
               c("SRR1528374", "SRR1528375", "SRR1528376", "SRR1528377"))
})

test_that("get_deg_data handles invalid groups", {
  # Non-existent treatment group
  expect_error(get_deg_data(isomirs, "X", "Inflorescence_WD"),
               "Groups not found in sample_info: X")
  # Non-existent control group
  expect_error(get_deg_data(isomirs, "Flower_FB", "X"),
               "Groups not found in sample_info: X")
})

test_that("expr columns match sample_info rows", {
  result <- get_deg_data(isomirs, "Flower_FB", "Inflorescence_WD")
  # Verify column names in expr match row names in sample_info
  expect_equal(colnames(result$expr), result$sample_info$name)
})
