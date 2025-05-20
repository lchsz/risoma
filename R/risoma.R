
#' Detect Isoforms from a Single Small RNA-Seq Sample
#'
#' This function detects miRNA isoforms (isomiRs) from a single small RNA-Seq
#' sample. It reads the sequencing data from a FASTQ file, caculate reads TPM,
#' and identifies isomiRs by comparing reads to reference miRNAs with specified
#' InDel and SNP thresholds.
#'
#' @param fq_file A character string specifying the path to the FASTQ file. The
#'   file can be compressed (`.gz`) or uncompressed.
#' @param mirnas A data frame containing reference miRNA information. Must
#' include the following columns:
#'   \itemize{
#'     \item mature_id: Unique identifier of the mature miRNA.
#'     \item mature_seq: Sequence of the mature miRNA.
#'     \item seed_seq: Seed sequence (center region of the mature miRNA).
#'     \item template_seq: Template sequence including flanking regions.
#'     \item flank_5p_seq: 5' flanking sequence of the mature miRNA.
#'     \item flank_3p_seq: 3' flanking sequence of the mature miRNA.
#'   }
#'
#' @param max_indel_5p An integer specifying the maximum allowed indels at
#'   the 5' end. Must be non-negative.
#' @param max_indel_3p An integer specifying the maximum allowed indels at the
#'   3' end. Must be non-negative.
#' @param max_snp_5p An integer specifying the maximum allowed mismatches at the
#'   5' end. Must be non-negative.
#' @param max_snp_3p An integer specifying the maximum allowed mismatches at the
#'   3' end. Must be non-negative.
#' @param max_snp_seed An integer specifying the maximum allowed mismatches in
#'   the seed region. Must be non-negative.
#' @param max_snp An integer specifying the maximum allowed mismatches. Must be
#'   non-negative.
#'
#' @return A data frame containing detected isoforms with the following columns:
#'   \itemize{
#'     \item mature_id: Reference miRNA ID.
#'     \item mature_seq: Reference miRNA sequence.
#'     \item template_seq: Template sequence including flanking regions.
#'     \item read_seq: Read sequence matched to the reference.
#'     \item tpm: Numeric vector of TPM values.
#'     \item indel_5p: Indels at the 5' end.
#'     \item indel_3p: Indels at the 3' end.
#'     \item snp_5p: Mismatches at the 5' end.
#'     \item snp_3p: Mismatches at the 3' end.
#'     \item snp_seed: Mismatches in the seed.
#'     \item dist: Total edit distance.
#'   }
#'
#' @details
#' The function performs the following steps:
#' 1. Checks if max_indel_5p, max_indel_3p, max_snp_5p, max_snp_3p and max_snp
#'    are non-negative.
#' 2. Reads the FASTQ file and caculate reads TPM using \code{load_reads} or
#'    \code{load_reads_gz}.
#' 3. Identifies isomiRs by comparing reads to reference miRNAs using the
#'    \code{find_isoforms} function.
#'
#' @export
detect_one_sample <- function(fq_file,
                              mirnas,
                              max_indel_5p = 3,
                              max_indel_3p = 3,
                              max_snp_5p = 1,
                              max_snp_3p = 1,
                              max_snp_seed = 1,
                              max_snp = 3) {
  if (is.null(fq_file))
    stop("fq_file is NULL")
  if (!file.exists(fq_file))
    stop(paste(fq_file, "does not exist"))
  if (is.null(mirnas))
    stop("mirnas is NULL")
  if (max_indel_5p < 0)
    stop("max_ed_5p should >= 0")
  if (max_indel_5p < 0)
    stop("max_ed_3p should >= 0")
  if (max_snp_5p < 0)
    stop("max_snp_5p should >= 0")
  if (max_snp_3p < 0)
    stop("max_snp_3p should >= 0")
  if (max_snp_seed < 0)
    stop("max_snp_seed should >= 0")
  if (max_snp < 0)
    stop("max_snp should >= 0")

  reads  <- NULL
  if (endsWith(fq_file, ".gz")) {
    reads <- load_reads_gz(fq_file)
  } else {
    reads <- load_reads(fq_file)
  }

  find_isoforms(
    mirnas,
    reads,
    max_indel_5p,
    max_indel_3p,
    max_snp_5p,
    max_snp_3p,
    max_snp_seed,
    max_snp
  )
}

#' Detect Isoforms for a Tissue or Treatment Group
#'
#' This function detects miRNA isoforms (isomiRs) for all samples within a
#' tissue or treatment group. It processes each sample's FASTQ file, identifies
#' isomiRs, calculates TPM, filters low-expression isoforms.
#'
#' @param sample_info A data frame containing sample information for the
#'   tissue/treatment group. Must include the following columns:
#'   \itemize{
#'      \item name: Character vector of sample names.
#'      \item fq: Character vector of paths to FASTQ files.
#'   }
#'
#' @param mirnas A data frame containing reference miRNA information. Must
#'   include the following columns:
#'   \itemize{
#'     \item mature_id: Unique identifier of the mature miRNA.
#'     \item mature_seq: Sequence of the mature miRNA.
#'     \item seed_seq: Seed sequence.
#'     \item template_seq: Template sequence including flanking regions.
#'     \item flank_5p_seq: 5' flanking sequence of the mature miRNA.
#'     \item flank_3p_seq: 3' flanking sequence of the mature miRNA.
#'   }
#'
#' @param max_indel_5p An integer specifying the maximum allowed indels at
#'   the 5' end. Must be non-negative.
#' @param max_indel_3p An integer specifying the maximum allowed indels at the
#'   3' end. Must be non-negative.
#' @param max_snp_5p An integer specifying the maximum allowed mismatches at the
#'  5' end. Must be non-negative.
#' @param max_snp_3p An integer specifying the maximum allowed mismatches at the
#'  3' end. Must be non-negative.
#' @param max_snp_seed An integer specifying the maximum allowed mismatches in
#'  the seed region. Must be non-negative.
#' @param max_snp An integer specifying the maximum allowed mismatches. Must be
#'  non-negative.
#' @param min_tpm A numeric value specifying the minimum TPM threshold for
#'   filtering low-expression isoforms. Isoforms with TPM below this threshold
#'   are excluded.
#'
#' @return A list of data frames, one for each sample. Each data frame contains
#'   the following columns:
#'   \itemize{
#'     \item mature_id: Reference miRNA ID.
#'     \item mature_seq: Reference miRNA sequence.
#'     \item template_seq: Template sequence including flanking regions.
#'     \item read_seq: Read sequence matched to the reference.
#'     \item tpm: Numeric vector of TPM values.
#'     \item indel_5p: Indels at the 5' end.
#'     \item indel_3p: Indels at the 3' end.
#'     \item snp_5p: Mismatches at the 5' end.
#'     \item snp_3p: Mismatches at the 3' end.
#'     \item snp_seed: Mismatches in the seed.
#'     \item dist: Total edit distance.
#'     \item id: A unique ID to each isoform.
#'   }
#'
#' @details
#' The function performs the following steps for each sample:
#' 1. Detects isomiRs using `detect_one_sample`.
#' 2. Assigns a unique ID to each isoform (`mature_id:read_seq`).
#' 3. Filters isoforms with TPM below `min_tpm`.
#' 4. Ensures consistency across samples by retaining only isoforms present in
#'    all samples.
#'
#' @export
detect_one_group <- function(sample_info,
                             mirnas,
                             max_indel_5p = 3,
                             max_indel_3p = 3,
                             max_snp_5p = 1,
                             max_snp_3p = 1,
                             max_snp_seed = 1,
                             max_snp = 3,
                             min_tpm = 1) {
  sample2isoforms <- apply(sample_info, 1, function(x) {
    message("Detecting isoforms for sample: ", x["name"])
    fq_file <- x["fq"]
    isoforms <- detect_one_sample(
      fq_file,
      mirnas,
      max_indel_5p,
      max_indel_3p,
      max_snp_5p,
      max_snp_3p,
      max_snp_seed,
      max_snp
    )
    isoforms$id <- paste(isoforms$mature_id, isoforms$read_seq, sep = ":")
    isoforms
  })

  expr <- sample2isoforms[[1]][, c("id", "tpm")]
  for (i in seq_along(sample2isoforms)[-1]) {
    expr_temp <- sample2isoforms[[i]][, c("id", "tpm")]
    expr <- merge(x = expr, y = expr_temp, by = "id", all = TRUE)
  }

  rownames(expr) <- expr$id
  expr$id <- NULL
  expr[is.na(expr)] = 0

  row_medians <- apply(expr, 1, median)
  expr$median <- row_medians
  expr <- subset(expr, median >= min_tpm)

  lapply(sample2isoforms, function(isoforms) {
    isoforms[isoforms$id %in% rownames(expr), ]
  })
}

#' Generate Expression Matrix for Samples
#'
#' This function generates an expression matrix for individual samples based on
#' isoform read TPM. It combines read TPM for each isoform across samples into
#' a single matrix, with isoforms as rows and samples as columns.
#'
#' @param group2sample_isoforms A list of lists, where each inner list contains
#'   isoform data frame for a specific group and sample. Each data frame must
#'   include the following columns:
#'   \itemize{
#'    \item read_seq: Character vector of isoform sequences.
#'    \item TPM: Numeric vector of TPM.
#'   }
#'
#' @return A numeric data frame where:
#'   \itemize{
#'     \item Rownames: Isoform sequences (read_seq).
#'     \item Columns: Sample names.
#'     \item Values: TPM. Missing values are replaced with 0.
#'   }
#'
#' @export
generate_sample_expr <- function(group2sample_isoforms) {
  group2expr <- lapply(group2sample_isoforms, function(sample2isoforms) {
    sample_name <- names(sample2isoforms[1])
    expr <- sample2isoforms[[1]][, c("read_seq", "tpm")]
    expr <- unique(expr)
    colnames(expr)[colnames(expr) == "tpm"] <- sample_name
    for (i in seq_along(sample2isoforms)[-1]) {
      sample_name <- names(sample2isoforms[i])
      expr_temp <- sample2isoforms[[i]][, c("read_seq", "tpm")]
      expr_temp <- unique(expr_temp)
      colnames(expr_temp)[colnames(expr_temp) == "tpm"] <- sample_name
      expr <- merge(x = expr, y = expr_temp, by = "read_seq", all = TRUE)
    }
    expr
  })

  expr <- group2expr[[1]]
  for (i in seq_along(group2expr)[-1]) {
    expr_temp <- group2expr[[i]]
    expr <- merge(x = expr,
                  y = expr_temp,
                  by = "read_seq",
                  all = TRUE)
  }

  rownames(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  expr
}

#' Generate Expression Matrix for Groups
#'
#' This function generates an expression matrix for groups (e.g., tissues or
#' treatments) based on isoform expression data. It combines TPM
#' (Transcripts Per Million) values for each isoform across groups into a single
#' matrix, with isoforms as rows and groups as columns.
#'
#' @param group2isoforms A list of data frames, where each data frame contains
#'   isoform data for a specific group. Each data frame must include the
#'   following columns:
#'   \itemize{
#'     \item read_seq: Character vector of isoform sequences.
#'     \item tpm: Numeric vector of mean TPM values across all the sampel in the
#'           a group.
#'   }
#'
#' @return A numeric data frame where:
#'  \itemize{
#'    \item Rownames: Isoform sequences (read_seq).
#'    \item Columns: Group names (names of group2isoforms).
#'    \item Values: TPM values. Missing values are replaced with 0.
#'  }
#'
#' @export
generate_group_expr <- function(group2isoforms) {
  group_name <- names(group2isoforms[1])
  expr <- group2isoforms[[1]][, c("read_seq", "tpm")]
  expr <- expr[!duplicated(expr$read_seq), ]
  colnames(expr)[colnames(expr) == "tpm"] <- group_name
  for (i in seq_along(group2isoforms)[-1]) {
    group_name <- names(group2isoforms[i])
    expr_temp <- group2isoforms[[i]][, c("read_seq", "tpm")]
    expr_temp <- expr_temp[!duplicated(expr_temp$read_seq), ]
    colnames(expr_temp)[colnames(expr_temp) == "tpm"] <- group_name
    expr <- merge(x = expr,
                  y = expr_temp,
                  by = "read_seq",
                  all = TRUE)
  }

  rownames(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  expr
}

#' Cluster Isoforms by Reference miRNAs
#'
#' This function clusters isoforms by their reference miRNAs and organizes the
#' data into a list of \linkS4class{Isomir} objects. It performs multiple
#' sequence alignment (MSA) using the provided sequences and generates CIGAR
#' strings to represent the alignment differences. It ensures that isoforms are
#' grouped by their reference miRNA and assigns unique identifiers to each
#' isoform based on its CIGAR string.
#'
#' @param group2isoforms A list of data frames, where each data frame contains
#'   isoform information for a specific group (e.g., tissue or treatment). Each
#'   data frame must include the following columns:
#'   \itemize{
#'     \item mature_id: Reference miRNA ID.
#'     \item mature_seq: Reference miRNA sequence.
#'     \item seed_seq: Seed sequence (center region of the mature miRNA).
#'     \item template_seq: Template sequence including flanking regions.
#'     \item read_seq: Read sequence matched to the reference.
#'     \item indel_5p: Number of indels at the 5' end.
#'     \item indel_3p: Number of indels at the 3' end.
#'     \item snp_5p: Number of SNPs at the 5' end.
#'     \item snp_3p: Number of SNPs at the 5' end.
#'     \item snp_seed:
#'     \item dist
#'   }
#'
#' @return A list of \linkS4class{Isomir}.
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts relevant columns (mature_id, mature_seq, template_seq, read_seq,
#'    indel_5p, indel_3p, snp_5p, snp_3p, snp_seed, dist) from each group's
#'    isoform data.
#' 2. Combines all isoforms into a single data frame and removes duplicates.
#' 3. Splits the isoforms by reference miRNA (mature_id).
#' 4. Performs multiple sequence alignment (MSA) using the isoform, reference
#'    miRNA and template sequences, and generates CIGAR strings to represent
#'    the alignment differences.
#' 5. For each reference miRNA, creates an \linkS4class{Isomir} object
#'    containing the associated isoforms.
#'
#' @export
cluster_isoforms <- function(group2isoforms) {
  if (is.null(group2isoforms))
    stop("group2isoforms is NULL")
  if (length(group2isoforms) == 0)
    stop ("group2isoforms is empty")

  group2isoforms <- lapply(group2isoforms, function(isoforms) {
    isoforms[, c(
      "mature_id",
      "mature_seq",
      "seed_seq",
      "template_seq",
      "read_seq",
      "indel_5p",
      "indel_3p",
      "snp_5p",
      "snp_3p",
      "snp_seed",
      "dist"
    )]
  })

  isoforms <- dplyr::bind_rows(group2isoforms)
  isoforms <- dplyr::distinct(isoforms)

  group <- split(isoforms, isoforms$mature_id)
  ref2isomir <- lapply(group, function(x) {
    x <- dplyr::arrange(x, desc(dist))

    mature_id = x$mature_id[[1]]
    mature_seq = x$mature_seq[[1]]
    template_seq = x$template_seq[[1]]

    seqs <- c(x$read_seq, mature_seq, template_seq)
    seqs <- Biostrings::DNAStringSet(seqs)
    msa <- msa::msa(seqs, order = "input", verbose = FALSE)
    seq_aln <- as.character(msa)

    query_num <- length(seq_aln) - 2
    cigars <- vector("character", query_num)

    seq2 <- seq_aln[length(seq_aln) - 1]
    for (i in seq_len(query_num)) {
      seq1 <- seq_aln[i]
      cigars[i] <- generate_cigar(seq1, seq2)
    }

    ids <- make.unique(cigars)
    ids[x$dist == 0] <- mature_id

    aln <- Biostrings::DNAStringSet(msa)
    aln[length(aln) - 1] <- NULL
    names(aln) <- c(ids, "template")

    new(
      "Isomir",
      mature_id = x$mature_id[[1]],
      mature_seq = x$mature_seq[[1]],
      seed_seq = x$seed_seq[[1]],
      template_seq = x$template_seq[[1]],
      read_seqs = x$read_seq,
      indel_5p = x$indel_5p,
      indel_3p = x$indel_3p,
      snp_5p = x$snp_5p,
      snp_3p = x$snp_3p,
      snp_seed = x$snp_seed,
      dist = x$dist,
      cigars = cigars,
      ids = ids,
      aln = aln
    )
  })

  ref2isomir
}

#' Detect IsomiRs from Small RNA-Seq Data
#'
#' This function detects isomiRs from small RNA-Seq data for multiple samples.
#' It processes FASTQ files, identifies isomiRs, and generates expression
#' profiles for samples and groups. The results are returned as an
#' \linkS4class{IsomirDataSet} object, which can be used for further analysis.
#'
#' @param sample_info_file A character string specifying the path to the sample
#'   information file. The file should be a CSV file with columns: name, group,
#'   and fq (FASTQ file name).
#' @param fq_dir A character string specifying the directory containing the
#'   FASTQ files.
#' @param mirnas A data frame containing microRNA data, typically generated by
#'   \code{load_mirnas}
#' @param max_indel_5p An integer specifying the maximum allowed indels at
#'   the 5' end. Must be non-negative (default: 3).
#' @param max_indel_3p An integer specifying the maximum allowed indels at the
#'   3' end. Must be non-negative (default: 3).
#' @param max_snp_5p An integer specifying the maximum allowed mismatches at the
#'   5' end. Must be non-negative (default: 1).
#' @param max_snp_3p An integer specifying the maximum allowed mismatches at the
#'   3' end. Must be non-negative (default: 1).
#' @param max_snp_seed An integer specifying the maximum allowed mismatches in
#'   the seed region. Must be non-negative (default: 1).
#' @param max_snp An integer specifying the maximum allowed mismatches. Must be
#'   non-negative (default: 3).
#' @param min_tpm A numeric value specifying the minimum TPM threshold for
#'   filtering low-expression isomiRs (default: 1).
#'
#' @return An object of class \linkS4class{IsomirDataSet}.
#'
#' @details
#' This function performs the following steps:
#' 1. Reads the sample information file and validates the input data (e.g.,
#'    checks for duplicate sample names or missing FASTQ files).
#' 2. Groups samples by their experimental conditions (e.g., tissue or treatment).
#' 3. Detects isomiRs for each sample using the \code{\link{detect_one_tissue}}
#'    function.
#' 4. Aggregates isomiR expression data across replicates within each group.
#' 5. Generates expression profiles for individual samples and groups.
#' 6. Clusters isomiRs based on reference mature microRNAs.
#' 7. Returns an \linkS4class{IsomirDataSet} object for further analysis.
#'
#' @export
detect_isomirs <- function(sample_info_file,
                           fq_dir,
                           mirnas,
                           max_indel_5p = 3,
                           max_indel_3p = 3,
                           max_snp_5p = 1,
                           max_snp_3p = 1,
                           max_snp_seed = 1,
                           max_snp = 3,
                           min_tpm = 1) {
  if (!file.exists(sample_info_file))
    stop(paste(sample_info_file, "does not exist"))

  if (!file.exists(fq_dir))
    stop(paste(fq_dir, "does not exist"))

  meta_data <- read.csv(
    sample_info_file,
    stringsAsFactors = FALSE,
    comment.char = "#",
    sep = ",",
    strip.white = TRUE
  )

  dup_name <- meta_data$name[duplicated(meta_data$name)]
  if (length(dup_name) != 0)
    stop("Duplicated sample name: ", paste(dup_name, collapse = ", "))

  meta_data$fq <- file.path(normalizePath(fq_dir), meta_data$fq)

  dup_fq <- meta_data$fq[duplicated(meta_data$fq)]
  if (length(dup_fq) != 0)
    stop("Duplicated FASTQ files:", paste(dup_fq, collapse = ", "))

  no_fq <- meta_data$fq[!file.exists(meta_data$fq)]
  if (length(no_fq) != 0)
    stop("FASTQ files do not exists: ", paste(no_fq, collapse = ", "))

  rownames(meta_data) <- meta_data$name
  group2sample_info <- split(meta_data, meta_data$group)

  # rep_num <- sapply(group2sample_info, function(x)
  #   nrow(x))

  # rep_one <- rep_num[rep_num < 2]
  # if (length(rep_one) != 0) {
  #   stop(
  #     "At least two replicates are required for tissu/treatment: ",
  #     paste(names(rep_one), collapse = ", ")
  #   )
  # }

  group2sample_isoforms <- lapply(
    group2sample_info,
    detect_one_group,
    mirnas = mirnas,
    max_indel_5p = max_indel_5p,
    max_indel_3p = max_indel_3p,
    max_snp_5p = max_snp_5p,
    max_snp_3p = max_snp_3p,
    max_snp_seed = max_snp_seed,
    max_snp = max_snp,
    min_tpm = min_tpm
  )

  group2isoforms <- lapply(group2sample_isoforms, function(sample2isoforms) {
    isoforms <- sample2isoforms[[1]]
    for (i in seq_along(sample2isoforms)[-1]) {
      isoforms_temp <- sample2isoforms[[i]]
      isoforms <- merge(x = isoforms, y = isoforms_temp,
                        by = c("id", "mature_id", "mature_seq", "seed_seq",
                               "template_seq", "read_seq", "indel_5p",
                               "indel_3p", "snp_5p", "snp_3p", "snp_seed",
                               "dist"))

      isoforms$tpm.x[is.na(isoforms$tpm.x)] = 0
      isoforms$tpm.y[is.na(isoforms$tpm.y)] = 0
      isoforms$tpm <- rowMeans(isoforms[, c("tpm.x", "tpm.y")])
      isoforms$tpm.x <- NULL
      isoforms$tpm.y <- NULL
    }
    isoforms
  })

  group2isoforms <- group2isoforms[sapply(group2isoforms, function(x)
    nrow(x) != 0)]

  message("Generating expression profiile for samples")
  sample_expr <- generate_sample_expr(group2sample_isoforms)

  message("Generating expression profiile for groups")
  group_expr <- generate_group_expr(group2isoforms)

  message("Clustriing isoforms based on reference mature microRNAs")
  ref2isomir <- cluster_isoforms(group2isoforms)

  new(
    "IsomirDataSet",
    sample_info = meta_data,
    group2isoforms = group2isoforms,
    sample_expr = sample_expr,
    group_expr = group_expr,
    ref2isomir = ref2isomir
  )
}
