#' Load microRNA Data from miRBase
#'
#' This function retrieves microRNA (miRNA) data from the miRBase database for a
#' specified species. It extracts mature miRNA sequences, seed sequences, and 
#' flanking regions, and returns the data in a structured data frame.
#'
#' @param spe A character string specifying the species name (e.g., "ath" for
#'   Arabidopsis thaliana).
#' @param flank_5p_len An integer specifying the length of the 5' flanking
#' region (default: 3).
#' @param flank_3p_len An integer specifying the length of the 3' flanking
#' region (default: 3).
#' @param seed_len Lenght of the seed sequence. (default: 16)
#'
#' @return A data frame containing the following columns:
#'  \itemize{
#'    \item mature_id: The ID of the mature miRNA.
#'    \item mature_seq: The sequence of the mature miRNA.
#'    \item mature_start: The start position of the mature miRNA in the precursor sequence.
#'    \item mature_end: The end position of the mature miRNA in the precursor sequence.
#'    \item seed_seq: The seed sequence of the mature miRNA.
#'    \item pre_id: The ID of the precursor miRNA.
#'    \item template_seq: The template sequence, including the mature miRNA and its flanking regions.
#'    \item flank_5p_seq: The 5' flanking sequence.
#'    \item flank_3p_seq: The 3' flanking sequence.
#'  }
#'
#' @examples
#' # Load miRNA data for Arabidopsis thaliana
#' mirnas <- load_mirnas(spe = "ath")
#'
#' # View the first few rows of the data
#' head(mirnas)
#'
#' @export
load_mirnas <- function(spe, flank_5p_len = 3, flank_3p_len = 3, seed_len = 16) {
  id2spe <- as.list(mirbase.db::mirbaseID2SPECIES)
  pre_ids <- AnnotationDbi::mappedkeys(id2spe[id2spe == spe])
  id2pre_seq <- as.list(mirbase.db::mirbaseSEQUENCE[pre_ids])
  id2matrue <- AnnotationDbi::mget(pre_ids, mirbase.db::mirbaseMATURE)

  mirnas <- vector("list", length(pre_ids))
  for (i in seq_along(pre_ids)) {
    pre_id <- pre_ids[[i]]
    pre_seq <- gsub("U", "T", toupper(id2pre_seq[[pre_id]]))
    matures <- id2matrue[[pre_id]]
    mature_starts = mirbase.db::matureFrom(matures)
    mature_ends = mirbase.db::matureTo(matures)
    mature_seqs <- substr(pre_seq, mature_starts, mature_ends)
    mature_seqs <- mature_seqs[nchar(mature_seqs) >= seed_len]

    seed_starts <- floor((nchar(mature_seqs) - seed_len) / 2) + 1
    seed_ends <- seed_starts + seed_len - 1
    seed_seqs <- substr(mature_seqs, seed_starts, seed_ends)

    template_stars <- mature_starts - flank_5p_len
    template_ends <- mature_ends + flank_3p_len
    template_stars[template_stars < 1] <- 1
    template_ends[template_ends > nchar(pre_seq)] <- nchar(pre_seq)

    flank_5p_seqs <- substr(pre_seq, template_stars, mature_starts - 1)
    flank_5p_seqs <- stringr::str_pad(flank_5p_seqs, width=flank_5p_len,
                                      pad = "N", side = "left")
    flank_3p_seqs <- substr(pre_seq, mature_ends + 1, template_ends)
    flank_3p_seqs <- stringr::str_pad(flank_3p_seqs, width=flank_3p_len,
                                      pad = "N", side = "right")
    template_seqs <- paste0(flank_5p_seqs, mature_seqs, flank_3p_seqs)

    df <- data.frame(
      mature_id = mirbase.db::matureName(matures),
      mature_seq = mature_seqs,
      mature_start = mature_starts,
      mature_end = mature_ends,
      seed_seq = seed_seqs,
      pre_id = pre_id,
      template_seq = template_seqs,
      flank_5p_seq = flank_5p_seqs,
      flank_3p_seq = flank_3p_seqs
    )

    mirnas[[i]] <- df
  }

  dplyr::bind_rows(mirnas)
}
