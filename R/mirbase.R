#' Load miRNAs of a species from mirbase.db
#'
#' @param spe species
#' @param word_size Length of seed in mature miRNA
#' @return A data.frame storing mature miRNA ID (matrue_id),
#'   mature miRNA sequence (mature_seq), mature start (mature_start),
#'   seed on miRNA (seed_seq), precursor ID (pre_id) and precursor sequence
#'   (pre_seq)
#' @examples
#' load_mirna("vvi", 13)
#'
#' @export
load_mirna <- function(spe, word_size = 13) {
  id2spe <- as.list(mirbase.db::mirbaseID2SPECIES)
  pre_ids <- AnnotationDbi::mappedkeys(id2spe[id2spe == spe])

  pre_seqs <- mirbase.db::mirbaseSEQUENCE
  id2pre_seq <- as.list(pre_seqs[pre_ids])

  id2matrue <- AnnotationDbi::mget(pre_ids, mirbase.db::mirbaseMATURE)

  out <- vector("list", length(pre_ids))

  for (i in seq_along(pre_ids)) {
    pre_id <- pre_ids[[i]]
    pre_seq <- gsub("U", "T", toupper(id2pre_seq[[pre_id]]))
    mature <- id2matrue[[pre_id]]
    mature_seq <- substr(pre_seq,
                         mirbase.db::matureFrom(mature),
                         mirbase.db::matureTo(mature))

    # should check if the length of mirna_seq is shorted than word_size?
    seed_start <- floor((nchar(mature_seq) - word_size) / 2) + 1
    seed_end <- seed_start + word_size - 1
    seed_seq <- substr(mature_seq, seed_start, seed_end)

    df <- data.frame(
      mature_id = mirbase.db::matureName(mature),
      mature_seq = mature_seq,
      mature_start = mirbase.db::matureFrom(mature),
      seed_seq = seed_seq,
      pre_id = pre_id,
      pre_seq = pre_seq
    )

    out[[i]] <- df
  }

  dplyr::bind_rows(out)
}
