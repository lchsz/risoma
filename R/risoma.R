#' ============================== exported functions =============================


#' Find isoform from a sample
#'
#' @param fq_file FASTQ of a samle
#' @param mirnas A data.frame storing mature miRNA ID (matrue_id),
#'   mature miRNA sequence (mature_seq), mature start (mature_start),
#'   seed on miRNA (seed_seq), precursor ID (pre_id) and precursor sequence
#'   (pre_seq)
#' @param max_ed_5p the maximum distance between 5’ region of reads with
#'   reference miRNA
#' @param max_ed_3p the maximum distance between 3' region of reads with
#'   reference miRNA
#' @return a data.frame with columns: mature_id, read_seq, read_num and dist
#'
#' @export
detect_one_sample <- function(fq_file, mirnas, max_ed_5p, max_ed_3p) {
  reads <- ifelse(
    endsWith(fq_file, ".gz"),
    mark_duplicates_gz(fq_file, 1),
    mark_duplicates(fq_file, 1)
  )

  isoforms <- find_isoforms(mirnas, reads, max_ed_5p, max_ed_3p)
  return(isoforms)
}


#' Detect isoform based on the sample information file containing sample name,
#' tissue/treatment and FASTQ path
#'
#' @param sample_info_file description
#' @param fq_dir description
#' @param mirnas description
#' @param max_ed_5p description
#' @param max_ed_3p description
#' @param min_tpm description
#' @return IsomirDataSet
#'
#' @export
detect_isomirs <- function(sample_info_file,
                          fq_dir,
                          mirnas,
                          max_ed_5p = 2,
                          max_ed_3p = 3,
                          min_tpm = 2) {

  meta_data <- read.csv(sample_info_file,
                       stringsAsFactors = FALSE,
                       comment.char = "#")

  dup_name <- meta_data$name[duplicated(meta_data$name)]
  if (length(dup_name) != 0) {
    stop("Duplicated sample name: ", paste(dup_name, collapse = ", "))
  }

  meta_data$fq <- file.path(fq_dir, meta_data$fq)

  dup_fq <- meta_data$fq[duplicated(meta_data$fq)]
  if (length(dup_fq) != 0) {
    stop("Duplicated FASTQ files:", paste(dup_fq, collapse = ", "))
  }

  no_fq <- meta_data$fq[!file.exists(meta_data$fq)]
  if (length(no_fq) != 0) {
    stop("FASTQ files do not exists: ", paste(no_fq, collapse = ", "))
  }

  rownames(meta_data) <- meta_data$name
  group <- split(meta_data, meta_data$group)

  rep_num <- sapply(group, function(x)
    nrow(x))

  rep_one <- rep_num[rep_num < 2]
  if (length(rep_one) != 0) {
    stop(
      "At least two replicates are required for tissu/treatment: ",
      paste(names(rep_one), collapse = ", ")
    )
  }

  isoform_list <- lapply(
    group,
    detect_one_tissue,
    mirnas = mirnas,
    max_ed_5p = max_ed_5p,
    max_ed_3p = max_ed_3p,
    min_tpm = min_tpm
  )

  isoform_list <- isoform_list[sapply(isoform_list, function(x) nrow(x) != 0)]

  message("Generating expression profiile from isomiR list")
  expr <- generate_expr(isoform_list)

  message("Clustriing isoforms based on reference mature microRNAs")
  isomir_list <- cluster_isoforms(isoform_list)

  new(
    "IsomirDataSet",
    isoform_list = isoform_list,
    expr = expr,
    isomir_list = isomir_list
  )
}


#' ====================== private functions ====================================


#' Detect isoforms from a tissue consist of several duplicates
#'
#' @param sample_info A data.frame with sample name, fq and group information
#'   for replicates in a tissu/treatment
#' @param mirnas A data.frame storing mature miRNA ID, mature miRNA sequence,
#'   mature start, seed on miRNA, precursor ID and precursor sequence
#' @param max_ed_5p the maximum distance between 5’ region of reads with
#'   reference miRNA
#' @param max_ed_3p the maximum distance between 3' region of reads with
#'   reference miRNA
#' @param min_tpm a isomiR is defined as “expressed” only when the TPM of all
#'    biological replicates were greater than or equal to the threshold set
#' @return a data.frame
detect_one_tissue <- function(sample_info,
                            mirnas,
                            max_ed_5p,
                            max_ed_3p,
                            min_tpm) {
  message("Detecting isoforms for tissue: ", sample_info["group"][[1]])
  tissue_isoforms <- apply(sample_info, 1, function(x) {
    sample_name <- x["name"]
    fq_file <- x["fq"]
    isoforms <- detect_one_sample(fq_file, mirnas, max_ed_5p, max_ed_3p)
    isoforms$tpm <- calc_tpm(isoforms)
    isoforms <- subset(isoforms, tpm >= min_tpm)
    isoforms$cigar <- calc_cigar(isoforms)
    isoforms
  })

  isoforms <- tissue_isoforms[[1]]
  for (i in seq_along(tissue_isoforms)[-1]) {
    isoforms_temp <- tissue_isoforms[[i]][, c("read_seq", "tpm")]
    isoforms <- merge(x = isoforms, y = isoforms_temp, by = "read_seq")
    isoforms$tpm <- rowMeans(isoforms[, c("tpm.x", "tpm.y")])
    isoforms$tpm.x <- NULL
    isoforms$tpm.y <- NULL
  }

  isoforms
}


cluster_isoforms <- function(isoform_list) {
  df_list <- lapply(isoform_list, function(x) {
    x[, c("mature_id", "mature_seq", "read_seq", "dist", "cigar")]
  })

  df <- dplyr::bind_rows(df_list)
  df <- dplyr::distinct(df)

  group <- split(df, df$mature_id)
  isomir_list <- lapply(group, function(x) {
    x <- dplyr::arrange(x, desc(dist))
    new(
      "Isomir",
      mature_id = x$mature_id[1],
      mature_seq = x$mature_seq[1],
      read_seqs = x$read_seq,
      dist = x$dist,
      cigars = x$cigar
    )
  })

  isomir_list[sapply(isomir_list, function(x) has_ref(x))]
}


#' Extract expression profile of each sample replicate from isomiRs
#'
#' @param isomir_list isomiRs of all samples
#' @return a data.frame of expression
generate_expr <- function(isoform_list) {
  tissue_name <- names(isoform_list[1])
  expr <- isoform_list[[1]][, c("read_seq", "tpm")]
  expr <- dplyr::distinct(expr, read_seq, .keep_all = TRUE)
  colnames(expr)[colnames(expr) == "tpm"] <- tissue_name
  for (i in seq_along(isoform_list)[-1]) {
    tissue_name <- names(isoform_list[i])
    expr_temp <- isoform_list[[i]][, c("read_seq", "tpm")]
    expr_temp <- dplyr::distinct(expr_temp, read_seq, .keep_all = TRUE)
    colnames(expr_temp)[colnames(expr_temp) == "tpm"] <- tissue_name
    expr <- merge(x = expr, y = expr_temp, by = "read_seq", all = TRUE)
  }

  rownames(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  return(expr)
}
