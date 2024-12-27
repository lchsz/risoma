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
loadMirna <- function(spe, wordSize = 13) {
  id2spe <- as.list(mirbase.db::mirbaseID2SPECIES)
  preIds <- AnnotationDbi::mappedkeys(id2spe[id2spe == spe])

  preSeqs <- mirbase.db::mirbaseSEQUENCE
  id2preSeq <- as.list(preSeqs[preIds])

  id2matrue <- AnnotationDbi::mget(preIds, mirbase.db::mirbaseMATURE)

  mirnas <- vector("list", length(preIds))

  for (i in seq_along(preIds)) {
    preId <- preIds[[i]]
    preSeq <- gsub("U", "T", toupper(id2preSeq[[preId]]))
    mature <- id2matrue[[preId]]
    matureSeq <- substr(preSeq,
                         mirbase.db::matureFrom(mature),
                         mirbase.db::matureTo(mature))

    # should check if the length of mirna_seq is shorted than word_size?
    seedStart <- floor((nchar(matureSeq) - wordSize) / 2) + 1
    seedEnd <- seedStart + wordSize - 1
    seedSeq <- substr(matureSeq, seedStart, seedEnd)

    df <- data.frame(
      mature_ID = mirbase.db::matureName(mature),
      mature_seq = matureSeq,
      mature_start = mirbase.db::matureFrom(mature),
      seed_seq = seedSeq,
      pre_ID = preId,
      pre_seq = preSeq
    )

    mirnas[[i]] <- df
  }

  dplyr::bind_rows(mirnas)
}
