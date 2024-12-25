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
detect_one <- function(fq_file, mirnas, max_ed_5p, max_ed_3p) {
  reads <- mark_duplicates(fq_file, 1)
  isoforms <- find_isoforms(mirnas, reads, max_ed_5p, max_ed_3p)
  return(isoforms)
}


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
#' @return a list
#'
#' @export
detect_tissue <- function(sample_info,
                          mirnas,
                          max_ed_5p,
                          max_ed_3p,
                          min_tpm) {
  tissue_isoforms <- apply(sample_info, 1, function(x) {
    sample_name <- x["name"]
    fq_file <- x["fq"]
    message("detecting isoform for sample: ", sample_name)
    isoforms <- detect_one(fq_file, mirnas, max_ed_5p, max_ed_3p)
    isoforms$tpm <- calc_expr(isoforms)
    subset(isoforms, tpm >= min_tpm)
  })

  read_seqs <- tissue_isoforms[[1]]$read_seq
  for (i in seq_along(tissue_isoforms)[-1]) {
    read_seqs <- intersect(read_seqs, tissue_isoforms[[i]]$read_seq)
  }

  lapply(tissue_isoforms, function(isoforms) {
    isoforms <- subset(isoforms, read_seq %in% read_seqs)
    cigars <- calc_cigar(isoforms)
    isoforms$id <- paste0(isoforms$mature_id, "-", cigars)
    isoforms
  })
}


#' Detect isoform based on the sample infomation file containing sample name,
#' tissue/treatment and FASTQ path
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
  meta_data <- read.csv(sample_info_file, stringsAsFactors = FALSE,
                        comment.char = "#")

  dup_name <- meta_data$name[duplicated(meta_data$name)]
  if (length(dup_name) != 0) {
    stop("Duplicated sample name: ", paste(dup_name, collapse = ", "))
  }

  dup_fq <- meta_data$fq[duplicated(meta_data$fq)]
  if (length(dup_fq) != 0) {
    stop("Duplicated FASTQ files:", paste(dup_fq, collapse = ", "))
  }

  no_fq <- meta_data$fq[!file.exists(meta_data$fq)]
  if (length(no_fq) != 0) {
    stop("FASTQ files do not exists: ", paste(no_fq, collapse = ", "))
  }

  row.names(meta_data) <- meta_data$name
  group <- split(meta_data, meta_data$group)

  rep_num <- sapply(group, function(x) nrow(x))
  rep_one <- rep_num[rep_num < 2]
  if (length(rep_one) != 0) {
    stop("At least two replicates are required for tissu/treatment: ",
         paste(names(rep_one), collapse = ", "))
  }

  lapply(
    group,
    detect_tissue,
    mirnas = mirnas,
    max_ed_5p = max_ed_5p,
    max_ed_3p = max_ed_3p,
    min_tpm = min_tpm
  )
}


#' Extract isomiRs from a list containing all isoforms
#'
#' @param all_isoforms a list containing all isoforms
#' @return a list of all isomiRs
#'
#' @export
generate_isomir <- function(all_isoforms) {
  lapply(all_isoforms, function(tissu_isoforms) {
    lapply(tissu_isoforms, function(isoforms) {
      group <- dplyr::group_by(isoforms, read_seq)
      group <- dplyr::mutate(group, ids = paste(id, collapse = ","))
      group <- dplyr::distinct(group, read_seq, .keep_all = T)
      as.data.frame(dplyr::ungroup(dplyr::select(group, -mature_id)))
    })
  })
}





#' Extract expression profile of each tissue/treatment from isomiRs
#'
#' @param isomirs isomiRs of all samples
#' @return a data.frame
#'
#' @export
group_expr <- function(isomirs) {
  expr_list <- lapply(isomirs, function(tissu_isomirs) {
    sample_isomir <- tissu_isomirs[[1]]
    # sample_isomir <- subset(sample_isomir, dist > 0)
    sample_expr <- sample_isomir[, c("read_seq", "tpm")]
    for (i in seq_along(tissu_isomirs)[-1]) {
      sample_isomir <- tissu_isomirs[[i]]
      # sample_isomir <- subset(sample_isomir, dist > 0)
      sample_expr_temp <- sample_isomir[, c("read_seq", "tpm")]
      sample_expr <- merge(x = sample_expr, y = sample_expr_temp, by = "read_seq")
    }
    row.names(sample_expr) <- sample_expr$read_seq
    sample_expr$read_seq <- NULL
    data.frame(read_seq = row.names(sample_expr),
               tpm = apply(sample_expr, 1, function(x) mean(x)))
  })

  expr <- expr_list[[1]]
  colnames(expr)[colnames(expr) == "tpm"] <- names(expr_list[1])

  for (i in  seq_along(expr_list)[-1]) {
    expr_temp <- expr_list[[i]]
    colnames(expr_temp)[colnames(expr_temp) == "tpm"] <- names(expr_list[i])
    expr <- merge(
      x = expr,
      y = expr_temp,
      by = "read_seq",
      all = TRUE
    )
  }
  row.names(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  return(expr)
}


#' Cluster and visualize expression profile
#'
#' @param expr description
#' @param scale description
#'
#' @export
plot_expr <- function(expr, scale = "row") {
  pheatmap::pheatmap(expr, scale = scale)
}


#' Extract expression profile of each sample replicate from isomiRs
#'
#' @param isomirs isomiRs of all samples
#' @return a data.frame
#'
#' @export
generate_expr <- function(isomirs) {
  expr_list <- lapply(isomirs, function(tissu_isomirs) {
    sample_name <- names(tissu_isomirs[1])
    sample_isomirs <- tissu_isomirs[[1]]
    sample_isomirs <- subset(sample_isomirs, dist > 0)
    sample_expr <- sample_isomirs[, c("read_seq", "tpm")]
    colnames(sample_expr)[colnames(sample_expr) == "tpm"] <- sample_name
    for (i in  seq_along(tissu_isomirs)[-1]) {
      sample_name <- names(tissu_isomirs[i])
      sample_isomirs <- tissu_isomirs[[i]]
      sample_isomirs <- subset(sample_isomirs, dist > 0)
      sample_expr_temp <- sample_isomirs[, c("read_seq", "tpm")]
      colnames(sample_expr_temp)[colnames(sample_expr_temp) == "tpm"] <- sample_name
      sample_expr <- merge(x = sample_expr, y = sample_expr_temp, by = "read_seq")
    }
    sample_expr
  })

  expr <- expr_list[[1]]
  for (i in  seq_along(expr_list)[-1]) {
    expr <- merge(
      x = expr,
      y = expr_list[[i]],
      by = "read_seq",
      all = TRUE
    )
  }
  row.names(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  return(expr)
}


#' Extract CIGAR of each sample replicate from isomiRs
#'
#' @param isomirs isomiRs of all samples
#' @return a data.frame
#'
#' @export
generate_cigar <- function(isomirs) {
  cigar_list <- lapply(isomirs, function(tissu_isomirs) {
    sample_name <- names(tissu_isomirs[1])
    sample_isomirs <- tissu_isomirs[[1]]
    sample_cigar <- sample_isomirs[, c("read_seq", "id")]
    colnames(sample_cigar)[colnames(sample_cigar) == "id"] <- sample_name
    for (i in  seq_along(tissu_isomirs)[-1]) {
      sample_name <- names(tissu_isomirs[i])
      sample_isomirs <- tissu_isomirs[[i]]
      sample_cigar_temp <- sample_isomirs[, c("read_seq", "id")]
      colnames(sample_cigar_temp)[colnames(sample_cigar_temp) == "id"] <- sample_name
      sample_cigar <- merge(x = sample_cigar, y = sample_cigar_temp, by = "read_seq")
    }
    sample_cigar
  })

  cigar <- cigar_list[[1]]
  for (i in  seq_along(cigar_list)[-1]) {
    cigar <- merge(
      x = cigar,
      y = cigar_list[[i]],
      by = "read_seq",
      all = TRUE
    )
  }
  row.names(cigar) <- cigar$read_seq
  cigar$read_seq <- NULL
  cigar[is.na(cigar)] = ""

  return(cigar)
}
