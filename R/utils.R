#' Get CIGAR string represents alignment between query and ref
#'
#' @param query query string
#' @param ref reference
#' @return CIGAR
#'
#' @export
calcOneCigar <- function(query, ref) {
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


pwAln <- function(query, subject, subMat) {
  alns <- pwalign::pairwiseAlignment(
    pattern = query,
    subject = subject,
    substitutionMatrix = subMat,
    gapOpening = 0,
    gapExtension = 1,
    type = "global"
  )
  aln <- alns[1]
  str1 <- as.character(pwalign::alignedPattern(aln))
  str2 <- as.character(pwalign::alignedSubject(aln))

  return(c(str1, str2))
}


calcCigar <- function(isoforms) {
  cigars <- NULL
  if (nrow(isoforms > 0)) {
    queries <- Biostrings::DNAStringSet(isoforms$read_seq)
    subjects <- Biostrings::DNAStringSet(isoforms$mature_seq)

    queryNum <- length(queries)
    cigars <- vector("character", queryNum)

    subMat <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = 0)

    for (i in seq_len(queryNum)) {
      query <- queries[i]
      subject <- subjects[i]
      ss <- pwAln(query, subject, subMat)
      cigar <- calcOneCigar(ss[1], ss[2])
      cigars[i] <- cigar
    }
  }

  return(cigars)
}
