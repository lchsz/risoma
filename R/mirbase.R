#' Load miRNAs of a species from mirbase.db
#'
#' @param spe species
#' @param seed_size Length of seed in mature miRNA
#' @return A data.frame storing mature miRNA ID (matrue_id),
#'   mature miRNA sequence (mature_seq), mature start (mature_start),
#'   seed on miRNA (seed_seq), precursor ID (pre_id) and precursor sequence
#'   (pre_seq)
#' @examples
#' load_mirna("vvi", 13)
#'
#' @export
load_mirnas <- function(spe, seed_size = 13) {
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

    seed_starts <- floor((nchar(mature_seqs) - seed_size) / 2) + 1
    seed_ends <- seed_starts + seed_size - 1
    seed_seqs <- substr(mature_seqs, seed_starts, seed_ends)

    template_stars <- mature_starts - 3
    template_ends <- mature_ends + 3
    template_stars[template_stars < 1] <- 1
    template_ends[template_ends > nchar(pre_seq)] <- nchar(pre_seq)
    template_seqs <-substr(pre_seq, template_stars, template_ends)

    df <- data.frame(
      mature_id = mirbase.db::matureName(matures),
      mature_seq = mature_seqs,
      mature_start = mature_starts,
      mature_end = mature_ends,
      seed_seq = seed_seqs,
      pre_id = pre_id,
      pre_seq = pre_seq,
      template_seq = template_seqs
    )

    mirnas[[i]] <- df
  }

  dplyr::bind_rows(mirnas)
}
