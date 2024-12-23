
#' Get CIGAR string represents alignment between query and ref
#'
#' @param query query string
#' @param ref reference
#' @return CIGAR
#'
#' @export
get_cigar <- function(query, ref) {
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
