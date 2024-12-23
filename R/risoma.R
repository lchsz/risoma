#' Find isoform from a sample
#'
#' @param fq_file FASTQ of a samle
#' @param mirnas A data.frame storing mature miRNA ID (matrue_id),
#'   mature miRNA sequence (mature_seq), mature start (mature_start),
#'   seed on miRNA (seed_seq), precursor ID (pre_id) and precursor sequence
#'   (pre_seq)
#' @param max_ed_5p the maximum distance between 5’ region of reads with reference miRNA
#' @param max_ed_3p the maximum distance between 3' region of reads with reference miRNA
#' @return a data.frame with columns: mature_id, read_seq, read_num and dist
#'
#' @export
detect_sample <- function(fq_file, mirnas, max_ed_5p, max_ed_3p) {
  reads <- mark_duplicates(fq_file, 1)
  isoforms <- find_isoforms(mirnas, reads, max_ed_5p, max_ed_3p)
  return(isoforms)
}


#' Detect isoforms from a tissue consist of several duplicates
#'
#' @param sample_info description
#' @param mirnas A data.frame storing mature miRNA ID, mature miRNA sequence,
#'    mature start, seed on miRNA, precursor ID and precursor sequence
#' @param max_ed_5p the maximum distance between 5’ region of reads with reference miRNA
#' @param max_ed_3p the maximum distance between 3' region of reads with reference miRNA
#' @param min_tpm a isomiR is defined as “expressed” only when the TPM of all
#'    biological replicates were greater than or equal to the threshold set
#' @return a list
#'
#' @export
detect_tissue <- function(sample_info,
                          mirnas,
                          max_ed_5p,
                          max_ed_3p,
                          min_tpm) {
  tissue_isoforms <- apply(sample_info, 1, function(x) {
    sample_name <- x["sample"]
    fq_file <- x["fq"]
    message("detecting isoform for sample: ", sample_name)
    sample_isoforms <- detect_sample(fq_file, mirnas, max_ed_5p, max_ed_3p)
    sample_isoforms$tpm <- calc_expr(sample_isoforms)
    subset(sample_isoforms, tpm >= min_tpm)
  })

  read_seqs <- tissue_isoforms[[1]]$read_seq
  for (i in seq_along(tissue_isoforms)[-1]) {
    read_seqs <- intersect(read_seqs, tissue_isoforms[[i]]$read_seq)
  }

  lapply(tissue_isoforms, function(sample_isoforms) {
    sample_isoforms <- subset(sample_isoforms, read_seq %in% read_seqs)
    cigars <- calc_cigar(sample_isoforms)
    sample_isoforms$id <- paste0(sample_isoforms$mature_id, "-", cigars)
    sample_isoforms
  })
}


#' Detect isoform based on the metadata file containing sample name, tissue and
#' FASTQ path
#'
#' @param sample_info_file description
#' @param spe description
#' @param word_size description
#' @param max_ed_5p description
#' @param max_ed_3p description
#' @param min_tpm description
#' @return description
#'
#' @export
detect_isoform <- function(sample_info_file,
                           spe,
                           word_size = 13,
                           max_ed_5p = 2,
                           max_ed_3p = 3,
                           min_tpm = 2) {
  mirnas <- load_mirna(spe, word_size)
  meta_data <- read.csv(sample_info_file, stringsAsFactors = FALSE)
  rownames(meta_data) <- meta_data$sample
  group <- split(meta_data, meta_data$tissue)
  lapply(
    group,
    detect_tissue,
    mirnas = mirnas,
    max_ed_5p = max_ed_5p,
    max_ed_3p = max_ed_3p,
    min_tpm = min_tpm
  )
}


#' Extract isomiRs from isoforms
#'
#' @param isoforms description
#' @return description
#'
#' @export
generate_isomir <- function(isoforms) {
  lapply(isoforms, function(tissu_isoform) {
    lapply(tissu_isoform, function(sample_isoform) {
      group <- dplyr::group_by(sample_isoform, read_seq)
      group <- dplyr::mutate(group, ids = paste(id, collapse = ","))
      group <- dplyr::distinct(group, read_seq, .keep_all = T)
      as.data.frame(dplyr::ungroup(dplyr::select(group, -mature_id)))
    })
  })
}

#' Extract expression profile from isomiRs
#'
#' @param isomirs description
#' @return a data.frame
#'
#' @export
generate_expr <- function(isomirs) {
  tissue_expr <- lapply(isomirs, function(tissu_isomir) {
    sampel_name <- names(tissu_isomir[1])
    sample_isomir <- tissu_isomir[[1]]
    sample_isomir <- subset(sample_isomir, dist > 0)
    sample_expr <- sample_isomir[, c("read_seq", "tpm")]
    colnames(sample_expr)[colnames(sample_expr) == "tpm"] <- sampel_name
    for (i in  seq_along(tissu_isomir)[-1]) {
      sampel_name <- names(tissu_isomir[i])
      sample_isomir <- tissu_isomir[[i]]
      sample_isomir <- subset(sample_isomir, dist > 0)
      sample_expr_temp <- sample_isomir[, c("read_seq", "tpm")]
      colnames(sample_expr_temp)[colnames(sample_expr_temp) == "tpm"] <-
        sampel_name
      sample_expr <- merge(x = sample_expr, y = sample_expr_temp, by = "read_seq")
    }
    sample_expr
  })

  expr <- tissue_expr[[1]]
  for (i in  seq_along(tissue_expr)[-1]) {
    expr <- merge(
      x = expr,
      y = tissue_expr[[i]],
      by = "read_seq",
      all = TRUE
    )
  }
  row.names(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  expr
}


#' Cluster and visualize expression profile
#'
#' @param expr description
#' @param scale description
#'
#' @export
draw_heatmap <- function(expr, scale = "row") {
  pheatmap::pheatmap(expr, scale = scale)
}


calc_cigar <- function(isoforms) {
  result <- NULL
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
      cigar <- get_cigar(ss[1], ss[2])
      cigars[i] <- cigar
    }
  }

  return(cigars)
}
