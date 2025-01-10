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
  if (max_ed_5p < 0) stop("max_ed_5p should >= 0")
  if( max_ed_3p < 0) stop("max_ed_3p should >= 0")
  reads  <- NULL
  if(endsWith(fq_file, ".gz")) {
    reads  <- mark_duplicates_gz(fq_file, 1)
  } else {
    reads  <-mark_duplicates(fq_file, 1)
  }

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
                          max_ed_5p = 3,
                          max_ed_3p = 6,
                          min_tpm = 2) {

  meta_data <- read.csv(sample_info_file,
                       stringsAsFactors = FALSE,
                       comment.char = "#")

  dup_name <- meta_data$name[duplicated(meta_data$name)]
  if (length(dup_name) != 0) {
    stop("Duplicated sample name: ", paste(dup_name, collapse = ", "))
  }

  meta_data$fq <- file.path(normalizePath(fq_dir), meta_data$fq)

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

  g2s_isoforms <- lapply(
    group,
    detect_one_tissue,
    mirnas = mirnas,
    max_ed_5p = max_ed_5p,
    max_ed_3p = max_ed_3p,
    min_tpm = min_tpm
  )

  g2isoforms <- lapply(g2s_isoforms, function(s2isoforms) {
    isoforms <- s2isoforms[[1]]
    for (i in seq_along(s2isoforms)[-1]) {
      isoforms_temp <- s2isoforms[[i]][, c("id", "tpm")]
      isoforms <- merge(x = isoforms, y = isoforms_temp, by = "id")
      isoforms$tpm <- rowMeans(isoforms[, c("tpm.x", "tpm.y")])
      isoforms$tpm.x <- NULL
      isoforms$tpm.y <- NULL
    }
    isoforms
  })

  g2isoforms <- g2isoforms[sapply(g2isoforms, function(x) nrow(x) != 0)]

  message("Generating expression profiile for samples")
  sample_expr <- generate_sample_expr(g2s_isoforms)

  message("Generating expression profiile for groups")
  group_expr <- generate_group_expr(g2isoforms)

  message("Clustriing isoforms based on reference mature microRNAs")
  ref2isomir <- cluster_isoforms(g2isoforms)

  new(
    "IsomirDataSet",
    sample_info = meta_data,
    group2isoforms = g2isoforms,
    sample_expr = sample_expr,
    group_expr = group_expr,
    ref2isomir = ref2isomir
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
  s2isoforms <- apply(sample_info, 1, function(x) {
    message("Detecting isoforms for sample: ", x["name"])
    fq_file <- x["fq"]
    isoforms <- detect_one_sample(fq_file, mirnas, max_ed_5p, max_ed_3p)
    isoforms$id <- paste(isoforms$mature_id, isoforms$read_seq, sep = ":")
    isoforms$tpm <- calc_tpm(isoforms)
    isoforms <- subset(isoforms, tpm >= min_tpm)
    isoforms$cigar <- calc_cigar(isoforms)
    isoforms
  })

  ids <- unique(unlist(lapply(s2isoforms, function(isoforms) isoforms$id)))

  lapply(s2isoforms, function(isoforms) {
    isoforms[isoforms$id %in% ids, ]
  })
}


cluster_isoforms <- function(g2isoforms) {
  g2isoforms <- lapply(g2isoforms, function(isoforms) {
    isoforms[, c("mature_id", "mature_seq", "template_seq", "read_seq", "dist",
          "dist_5p", "dist_3p", "cigar")]
  })

  isoforms <- dplyr::bind_rows(g2isoforms)
  isoforms <- dplyr::distinct(isoforms)

  group <- split(isoforms, isoforms$mature_id)
  ref2isomir <- lapply(group, function(x) {
    x <- dplyr::arrange(x, desc(dist))
    x$cigar[x$dist == 0] <- x$mature_id[[1]]
    new(
      "Isomir",
      mature_id = x$mature_id[[1]],
      mature_seq = x$mature_seq[[1]],
      template_seq = x$template_seq[[1]],
      read_seqs = x$read_seq,
      dist = x$dist,
      dist_5p = x$dist_5p,
      dist_3p = x$dist_3p,
      cigars = x$cigar,
      ids = make.unique(x$cigar)
    )
  })

  ref2isomir[sapply(ref2isomir, function(isomir) has_ref(isomir))]
}


#' Extract expression profile of each sample replicate from isomiRs
#'
#' @param g2isoforms isomiRs of all samples
#' @return a data.frame of expression
generate_group_expr <- function(g2isoforms) {
  g_name <- names(g2isoforms[1])
  expr <- g2isoforms[[1]][, c("read_seq", "tpm")]
  expr <- expr[!duplicated(expr$read_seq), ]
  colnames(expr)[colnames(expr) == "tpm"] <- g_name
  for (i in seq_along(g2isoforms)[-1]) {
    g_name <- names(g2isoforms[i])
    expr_temp <- g2isoforms[[i]][, c("read_seq", "tpm")]
    expr_temp <- expr_temp[!duplicated(expr_temp$read_seq), ]
    colnames(expr_temp)[colnames(expr_temp) == "tpm"] <- g_name
    expr <- merge(x = expr, y = expr_temp, by = "read_seq", all = TRUE)
  }

  rownames(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  return(expr)
}


generate_sample_expr <- function(g2s_isoforms) {
  g2expr <- lapply(g2s_isoforms, function(s2isoforms) {
    s_name <- names(s2isoforms[1])
    expr <- s2isoforms[[1]][, c("read_seq", "read_num")]
    expr <- unique(expr)
    colnames(expr)[colnames(expr) == "read_num"] <- s_name
    for (i in seq_along(s2isoforms)[-1]) {
      s_name <- names(s2isoforms[i])
      expr_temp <- s2isoforms[[i]][, c("read_seq", "read_num")]
      expr_temp <- unique(expr_temp)
      colnames(expr_temp)[colnames(expr_temp) == "read_num"] <- s_name
      expr <- merge(x = expr, y = expr_temp, by = "read_seq")
    }
    expr
  })

  expr <- g2expr[[1]]
  for (i in seq_along(g2expr)[-1]) {
    expr_temp <- g2expr[[i]]
    expr <- merge(x = expr, y = expr_temp, by = "read_seq", all = TRUE)
  }

  rownames(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  return(expr)
}
