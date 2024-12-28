
setMethod("has_ref", "Isomir", function(x) {
  any(x@dist == 0)
})


#' =============================IsomirDataSet =================================

setMethod("get_ref", "IsomirDataSet", function(x) {
  sapply(x@isomir_list, function(x) {
    length(x@read_seqs)
  })
})


setMethod("ex_isomir", "IsomirDataSet", function(x, ref) {
  x@isomir_list[[ref]]
})


setMethod("ex_expr", "IsomirDataSet", function(x, ref) {
  result <- NULL
  isomir <- x@isomir_list[[ref]]
  if (!is.null(isomir)) {
    expr <- x@expr
    result <- subset(expr, rownames(expr) %in% isomir@read_seqs)
  }

  cigars <- isomir@cigars[match(rownames(result), isomir@read_seqs)]
  rownames(result) <- make.unique(cigars)

  return(result)
})


setMethod("aln_isoforms", "IsomirDataSet", function(x, ref) {
  isomir <- ex_isomir(x, ref)
  seqs <- isomir@read_seqs
  names(seqs) <- make.unique(isomir@cigars)
  seqs <- Biostrings::DNAStringSet(seqs, use.names = TRUE)
  seq_aln <- msa::msa(seqs, order = "input")
  Biostrings::DNAStringSet(seq_aln)
})


setMethod("get_tissues", "IsomirDataSet", function(x) {
  result <- NULL
  colnames(x@expr)
})


setMethod("ex_ref_expr", "IsomirDataSet", function(x) {
  mature_seqs <- sapply(x@isomir_list, function(x) x@mature_seq)
  expr <- x@expr[mature_seqs, ]
  rownames(expr) <- names(mature_seqs)

  return(expr)
})
