
setMethod("has_ref", "Isomir", function(x) {
  any(x@dist == 0)
})


#' ============================= IsomirDataSet ================================

setMethod("has_isomir", "IsomirDataSet", function(x) {
  sapply(x@ref2isomir, function(isomir) length(isomir@read_seqs))
})


setMethod("ex_isomir", "IsomirDataSet", function(x, ref) {
  x@ref2isomir[[ref]]
})


setMethod("ex_expr", "IsomirDataSet", function(x, ref) {
  result <- NULL
  isomir <- x@ref2isomir[[ref]]
  if (!is.null(isomir)) {
    expr <- x@group_expr
    result <- subset(expr, rownames(expr) %in% isomir@read_seqs)
  }

  cigars <- isomir@cigars[match(rownames(result), isomir@read_seqs)]
  rownames(result) <- make.unique(cigars)

  result
})


setMethod("aln_isoforms", "IsomirDataSet", function(x, ref) {
  isomir <- ex_isomir(x, ref)
  seqs <- isomir@read_seqs
  names(seqs) <- isomir@ids
  template_seq <- isomir@template_seq
  names(template_seq) <- "template"
  seqs <- c(seqs, template_seq)
  seqs <- Biostrings::DNAStringSet(seqs, use.names = TRUE)
  seq_aln <- msa::msa(seqs, order = "input")
  Biostrings::DNAStringSet(seq_aln)
})


setMethod("get_groups", "IsomirDataSet", function(x) {
  result <- NULL
  colnames(x@group_expr)
})


setMethod("ex_ref_expr", "IsomirDataSet", function(x) {
  mature_seqs <- sapply(x@ref2isomir, function(isomir) isomir@mature_seq)
  expr <- x@group_expr[mature_seqs, ]
  rownames(expr) <- names(mature_seqs)

  expr
})


setMethod("calc_tsi", "IsomirDataSet", function(x) {
  mature_seqs <- sapply(x@ref2isomir, function(isomir) isomir@mature_seq)
  expr <- x@group_expr
  type <- rep("isoform", nrow(expr))
  type[rownames(expr) %in% mature_seqs] <- "ref"
  tau <- round(apply(expr, 1, calc_one_tsi), 2)
  names(tau) <- NULL

  data.frame(read_seq = rownames(expr), type = type, tsi = tau)
})


setMethod("calc_deg", "IsomirDataSet",
          function(x, treatment, control, padj=0.05) {
  sample_info <- x@sample_info
  sample_info <- sample_info[sample_info$group %in% c(treatment, control), ]
  sample_info$group <- factor(sample_info$group, levels = c(control, treatment))
  expr <- x@sample_expr
  expr <- expr[rownames(sample_info)]
  expr <- expr[rowSums(expr) != 0, ]

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = expr,
                                colData = sample_info,
                                design = ~ group)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::lfcShrink(dds, coef = DESeq2::resultsNames(dds)[[2]],
                   type = "apeglm")
  res <- res[!(res$baseMean == 0 &
                 is.na(res$padj) & is.na(res$pvalue) &
                 is.na(res$log2FoldChange)), ]
  res <- res[!res$log2FoldChange == 0, ]
  res <- res[!(is.na(res$pvalue) & is.na(res$padj)), ]
  res <- res[!is.na(res$padj), ]
  res <- res[res$padj <= padj, ]

  res[order(res$log2FoldChange), ]
})


setMethod("has_deg", "IsomirDataSet", function(x, deg) {
  deg_num <- sapply(x@ref2isomir,
                    function(isomir) sum(isomir@read_seqs %in% rownames(deg)))

  deg_num[deg_num != 0]
})


setMethod("ex_deg", "IsomirDataSet", function(x, deg, ref) {
  result <- NULL
  isomir <- x@ref2isomir[[ref]]

  if (!is.null(isomir)) {
    result <- subset(deg, rownames(deg) %in% isomir@read_seqs)
  }

  ids <- isomir@ids[match(rownames(result), isomir@read_seqs)]
  rownames(result) <- ids

  result
})
