
setMethod("plotExpr", "IsomirDataSet", function(x) {
  pheatmap::pheatmap(x@expr, scale = "row")
})


setMethod("exIsomirByRef", "IsomirDataSet", function(x, ref) {
  x@isomirList[[ref]]
})


setMethod("hasRef", "IsomirDataSet", function(x) {
  sapply(x@isomirList, function(isomir) {
    any(isomir@dist == 0)
  })
})


setMethod("getRef", "IsomirDataSet", function(x) {
  names(x@isomirList)
})


setMethod("exExpr", "IsomirDataSet", function(x, ref) {
  result <- NULL
  isomir <- x@isomirList[[ref]]
  if (!is.null(isomir)) {
    expr <- x@expr
    result <- subset(expr, rownames(expr) %in% isomir@readSeqs)
  }

  cigars <- isomir@cigars[match(rownames(result), isomir@readSeqs)]
  rownames(result) <- make.unique(cigars)

  return(result)
})


setMethod("alnIsoforms", "IsomirDataSet", function(x, ref) {
  isomir <- exIsomirByRef(x, ref)
  seqs <- isomir@readSeqs
  names(seqs) <- make.unique(isomir@cigars)
  seqs <- Biostrings::DNAStringSet(seqs, use.names = TRUE)
  seqAln <- msa::msa(seqs, order = "input")
  Biostrings::DNAStringSet(seqAln)
})
