#' Get CIGAR string represents alignment between query and ref
#'
#' @param query query string
#' @param ref reference
#' @return CIGAR
#'
#' @export
calc_one_cigar <- function(query, ref) {
  query <- strsplit(query, "")[[1]]
  ref <- strsplit(ref, "")[[1]]
  cigar <- rep("M", length(query))
  cigar[query == "-"] <- "D"
  query[cigar == "D"] <- "D"
  ref[cigar == "D"] <- "D"
  cigar[ref == "-"] <- "I"
  query[cigar == "I"] <- "I"
  ref[cigar == "I"] <- "I"
  cigar[query != ref] <- "X"
  r <- rle(cigar)
  paste0(as.character(r$lengths), r$values, collapse = "")
}


pw_aln <- function(query, subject, sub_mat) {
  alns <- pwalign::pairwiseAlignment(
    pattern = query,
    subject = subject,
    substitutionMatrix = sub_mat,
    gapOpening = 0,
    gapExtension = 1,
    type = "global"
  )
  aln <- alns[1]
  str1 <- as.character(pwalign::alignedPattern(aln))
  str2 <- as.character(pwalign::alignedSubject(aln))

  return(c(str1, str2))
}


calc_cigar <- function(isoforms, mirnas) {
  cigars <- NULL
  if (nrow(isoforms > 0)) {
    queries <- Biostrings::DNAStringSet(isoforms$read_seq)
    subjects <- Biostrings::DNAStringSet(isoforms$mature_seq)

    query_num <- length(queries)
    cigars <- vector("character", query_num)

    sub_mat <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = 0)

    for (i in seq_len(query_num)) {
      query <- queries[i]
      subject <- subjects[i]
      ss <- pw_aln(query, subject, sub_mat)
      cigar <- calc_one_cigar(ss[1], ss[2])
      cigars[i] <- cigar
    }
  }

  return(cigars)
}
