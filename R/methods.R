
# Get isomiR information ======================================================

## get_isoform_num -----
#' Count Isoforms per Reference miRNA
#'
#' This method calculates the number of isoforms associated with each reference
#' miRNA in an \linkS4class{IsomirDataSet} object. It returns a vector of counts
#' for reference miRNAs.
#'
#' @docType methods
#' @name get_isoform_num
#' @rdname get_isoform_num
#'
#' @param x An \linkS4class{IsomirDataSet} object containing isomiR data.
#'   Requires the \strong{ref2isomir} slot, which stores a list of reference
#'   miRNAs with associated isoforms.
#'
#' @return A named numeric vector where:
#'   \itemize{
#'     \item Names: Reference miRNA IDs.
#'     \item Values: Number of isoforms associated with each reference miRNA.
#'   }
#'
#' @examples
#' # Load example IsomirDataSet
#' data(isomirs)
#'
#' # Count isoforms per reference miRNA
#' get_isoform_num(isomirs)
#'
#' @export
setMethod("get_isoform_num", "IsomirDataSet", function(x) {
  if (is.null(x))
    stop("x in null")

  stopifnot(is(x, "IsomirDataSet"))

  isoform_num <- sapply(x@ref2isomir, function(isomir) {
    read_seqs <- isomir@read_seqs[isomir@dist != 0]
    length(read_seqs)
  })

  isoform_num[isoform_num != 0]
})


## get_isomir_by_ref -------------
#' Extract IsomiR Data for a Specific Reference miRNA
#'
#' This method extracts the isomiR data associated with a specified reference
#' miRNA from an \linkS4class{IsomirDataSet} object. It returns the isomiR
#' information stored in the \strong{ref2isomir} slot for the given reference
#' miRNA.
#'
#' @docType methods
#' @name get_isomir_by_ref
#' @rdname get_isomir_by_ref
#'
#' @param x An \linkS4class{IsomirDataSet} object containing isomiR data.
#'   Requires the \strong{ref2isomir} slot, which stores a list of reference
#'   miRNAs with associated isomiR information.
#' @param ref A character string specifying the reference miRNA ID.
#'
#' @return A list containing the isomiR data for the specified reference miRNA.
#'
#' @examples
#' # Load example IsomirDataSet
#' data(isomirs)
#'
#' # Extract isomiR data for a specific reference miRNA
#' get_isomir_by_ref(isomirs, ref = "vvi-miR396c")
#'
#' @export
setMethod("get_isomir_by_ref", "IsomirDataSet", function(x, ref) {
  if (is.null(x))
    stop("x is null")

  stopifnot(is(x, "IsomirDataSet"))

  x@ref2isomir[[ref]]
})


## get_alignment_by_ref ------
#' Extract Aligned Isoforms to a Reference Sequence
#'
#' This method extract the alignment of the sequences of isomiRs and their
#' reference template sequence using multiple sequence alignment (MSA). It is
#' useful for visualizing sequence variations and identifying potential editing
#' events in isomiRs.
#'
#' @docType methods
#' @name get_alignment_by_ref
#' @rdname get_alignment_by_ref
#'
#' @param x An object of class \linkS4class{IsomirDataSet} containing isomiR data.
#' @param ref A character string specifying the reference mature microRNA ID.
#'
#' @return A \linkS4class{DNAStringSet} object containing the aligned sequences
#'   of isomiRs and the reference template sequence.
#'
#' @examples
#' # Load example data
#' data(isomirs)
#'
#' # Align isoforms for a specific reference microRNA
#' get_alignment_by_ref(isomirs, ref = "vvi-miR396c")
#'
#' @export
setMethod("get_alignment_by_ref", "IsomirDataSet", function(x, ref) {
  if (is.null(x))
    stop("x is null")

  stopifnot(is(x, "IsomirDataSet"))

  result <- NULL

  isomir <- get_isomir_by_ref(x, ref)

  if (!is.null(isomir))
    result <- isomir@aln

  return(result)
})

# Expression Data ==============================================================

## calc_group_expr_num ----------------------
#' Calculate the Number of Expressed IsomiRs
#'
#' This method calculates the number of expressed isomiRs for each group in an
#' \linkS4class{IsomirDataSet} object.
#' It categorizes isomiRs into three types:
#' 1. Total: The total number of isomiRs expressed in each group.
#' 2. Organ-Specific: The number of isomiRs expressed in only one group.
#' 3. Common: The number of isomiRs expressed in all groups.
#'
#' @docType methods
#' @name calc_group_expr_num
#' @rdname calc_group_expr_num
#'
#' @param x An object of class \linkS4class{IsomirDataSet} containing isomiR
#'   expression data.
#'
#' @return A data frame with the following columns:
#'   \itemize{
#'     \item group: The name of the group.
#'     \item total: The total number of isomiRs expressed in the group.
#'     \item organ_specific: The number of isomiRs specific to the group
#'           (expressed in only this group).
#'     \item common: The number of isomiRs expressed in all groups.
#'   }
#'
#' @examples
#' # Load example data
#' data(isomirs)
#'
#' # Calculate the number of expressed isomiRs
#' calc_group_expr_num(isomirs)
#'
#' @export
setMethod("calc_group_expr_num", "IsomirDataSet", function(x) {
  stopifnot(is(x, "IsomirDataSet"))

  mature_seqs <- sapply(x@ref2isomir, function(isomir)
    isomir@mature_seq)

  expr <- x@group_expr
  expr <- expr[!(rownames(expr) %in% mature_seqs), ]

  specific_expr <- expr[apply(expr, 1, function(x)
    sum(x > 0) == 1), ]
  common_expr <- expr[apply(expr, 1, function(x)
    sum(x > 0) == length(x)), ]
  total <- apply(expr, 2, function(x)
    sum(x > 0))
  organ_specific <- apply(specific_expr, 2, function(x)
    sum(x > 0))
  common = nrow(common_expr)
  df <- data.frame(total, organ_specific, common)
  df$group <- rownames(df)
  rownames(df) <- NULL

  df[c("group", "total", "organ_specific", "common")]
})


## get_ref_expr --------------
#' Extract Expression Data for Reference miRNAs
#'
#' This method extracts the expression data for reference miRNAs from an
#' \linkS4class{IsomirDataSet} object. It retrieves the expression values of
#' reference miRNAs across experimental groups and returns a data frame with
#' reference miRNA IDs as row names.
#'
#' @docType methods
#' @name get_ref_expr
#' @rdname get_ref_expr
#'
#' @param x An \linkS4class{IsomirDataSet} object containing isomiR data.
#'
#' @return A numeric data frame where:
#'   \itemize{
#'     \item Rownames: Reference miRNA IDs.
#'     \item Columns: Experimental groups.
#'     \item Values: Expression values (TPM).
#'   }
#'
#' @examples
#' # Load example IsomirDataSet
#' data(isomirs)
#'
#' # Extract expression data for reference miRNAs
#' get_ref_expr(isomirs)
#'
#' @export
setMethod("get_ref_expr", "IsomirDataSet", function(x) {
  stopifnot(is(x, "IsomirDataSet"))

  expr <- x@group_expr

  mature_seqs <- sapply(x@ref2isomir, function(isomir)
    isomir@mature_seq)

  mature_df <- as.data.frame(mature_seqs)
  mature_df <- unique(mature_df)
  mature_df <- subset(mature_df, mature_seqs %in% rownames(expr))

  expr <- expr[mature_df$mature_seqs, ]
  rownames(expr) <- rownames(mature_df)

  expr
})


## get_expr_by_ref -------

#' Extract Expression Data for IsomiRs of a Specific Reference
#'
#' This method extracts expression data for isomiRs associated with a specified
#' reference mature microRNA. It is useful for analyzing the expression patterns
#' of isomiRs derived from a specific microRNA across different conditions or
#' groups.
#'
#' @docType methods
#' @name get_expr_by_ref
#' @rdname get_expr_by_ref
#'
#' @param x An object of class \linkS4class{IsomirDataSet} containing isomiR
#'   data.
#' @param ref A character string specifying the reference mature microRNA ID.
#'
#' @return A data frame containing the expression data for isomiRs associated
#'   with the specified reference microRNA. The row names of the data frame are
#'   replaced with unique CIGAR strings for easier interpretation.
#'
#' @examples
#' # Load example data
#' data(isomirs)
#'
#' # Extract expression data for a specific reference microRNA
#' get_expr_by_ref(isomirs, ref = "vvi-miR396c")
#'
#' @export
setMethod("get_expr_by_ref", "IsomirDataSet", function(x, ref) {
  if (is.null(x))
    stop("x is null")

  if (is.null(ref))
    stop("ref is null")

  stopifnot(is(x, "IsomirDataSet"))

  result <- NULL
  isomir <- x@ref2isomir[[ref]]
  if (!is.null(isomir)) {
    expr <- x@group_expr
    result <- subset(expr, rownames(expr) %in% isomir@read_seqs)
    ids <- isomir@ids[match(rownames(result), isomir@read_seqs)]
    rownames(result) <- ids
  }

  result
})


# Calculate Tissue Specificity Index =================================

## calc_tsi -----
#' Calculate Tissue Specificity Index (TSI) for Isoforms and Reference miRNAs
#'
#' This method computes the Tissue Specificity Index (TSI) for both isomiRs and
#' reference miRNAs across experimental groups in an \linkS4class{IsomirDataSet}
#' object. TSI quantifies how specifically a sequence is expressed in a single
#' group (1 = group-specific, 0 = ubiquitously expressed).
#'
#' @docType methods
#' @name calc_tsi
#' @rdname calc_tsi
#'
#' @param x An \linkS4class{IsomirDataSet} object containing isoform expression
#'   data and reference miRNA annotations.
#'
#' @return A data frame with three columns:
#'   \itemize{
#'     \item read_seq: Character vector of sequence identifiers (row names of
#'           the expression matrix).
#'     \item type: Factor indicating if the sequence is a reference miRNA
#'            ("ref") or isoform ("isoform").
#'     \item tsi: Numeric vector of TSI values rounded to 2 decimal places.
#'   }
#'
#' @examples
#' # Load example IsomirDataSet
#' data(isomirs)
#'
#' # Calculate TSI values
#' tsi_results <- calc_tsi(isomirs)
#'
#' # View top 5 most specific isoforms
#' head(tsi_results[order(-tsi_results$tsi), ], 5)
#'
#' @export
setMethod("calc_tsi", "IsomirDataSet", function(x) {
  stopifnot(is(x, "IsomirDataSet"))

  mature_seqs <- sapply(x@ref2isomir, function(isomir)
    isomir@mature_seq)
  expr <- x@group_expr
  type <- rep("isoform", nrow(expr))
  type[rownames(expr) %in% mature_seqs] <- "reference"
  tau <- round(apply(expr, 1, calc_one_tsi), 2)
  names(tau) <- NULL

  data.frame(read_seq = rownames(expr),
             type = type,
             tsi = tau)
})


# DEG ====================================================================

## get_deg_data ---------------

#' Extract differential expression data from IsomirDataSet
#'
#' This function extracts expression data and sample information from an
#' IsomirDataSet object for specified treatment and control groups, preparing
#' it for differential expression analysis. The function allows selecting
#' either reference sequences or isoform sequences.
#'
#' @docType methods
#' @name get_deg_data
#' @rdname get_deg_data
#'
#' @param x An object of class \linkS4class{IsomirDataSet}
#' @param treatment Character vector specifying the treatment group name
#' @param control Character vector specifying the control group name
#'
#' @return A list containing two elements:
#' \itemize{
#'   \item \code{expr} - Expression matrix with sequences as rows and samples
#'         as columns. For type="ref", contains reference sequences;
#'         for type="isoform", contains isomiR sequences.
#'   \item \code{sample_info} - Data frame containing sample metadata with
#'         group factor ordered as c(treatment, control)
#' }
#'
#' @examples
#' data(isomirs)
#' deg_data <- get_deg_data(isomirs,
#'                         treatment = "Flower_FB",
#'                         control = "Inflorescence_WD")
#'
#' @export
setMethod("get_deg_data", "IsomirDataSet", function(x,
                                                treatment,
                                                control) {
  stopifnot(is(x, "IsomirDataSet"))

  if (!all(c(treatment, control) %in% x@sample_info$group)) {
    missing_groups <- setdiff(c(treatment, control), x@sample_info$group)
    stop(paste("Groups not found in sample_info:", paste(missing_groups, collapse = ", ")))
  }

  condiction <- data.frame(group = c(treatment, control))

  sample_info <- x@sample_info
  sample_info$fq <- NULL
  sample_info <- dplyr::left_join(condiction, sample_info, by = "group")
  sample_info$group <- factor(sample_info$group, levels = c(treatment, control))
  rownames(sample_info) <- sample_info$name

  expr <- x@sample_expr

  expr <- expr[rownames(sample_info)]
  expr <- expr[rowSums(expr) != 0, ]

  list(expr = expr, sample_info = sample_info)
})
