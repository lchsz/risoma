#' ============================== exported functions =============================

#' Find isoform from a sample
#'
#' @param fqFile FASTQ of a samle
#' @param mirnas A data.frame storing mature miRNA ID (matrue_id),
#'   mature miRNA sequence (mature_seq), mature start (mature_start),
#'   seed on miRNA (seed_seq), precursor ID (pre_id) and precursor sequence
#'   (pre_seq)
#' @param maxEd5p the maximum distance between 5’ region of reads with reference miRNA
#' @param maxEd3p the maximum distance between 3' region of reads with reference miRNA
#' @return a data.frame with columns: mature_id, read_seq, read_num and dist
#'
#' @export
detectOneSample <- function(fqFile, mirnas, maxEd5p, maxEd3p) {
  reads <- markDuplicates(fqFile, 1)
  isoforms <- findIsoforms(mirnas, reads, maxEd5p, maxEd3p)
  return(isoforms)
}

#' Detect isoform based on the sample information file containing sample name,
#' tissue/treatment and FASTQ path
#'
#' @param sampleInfoFile description
#' @param spe description
#' @param wordSize description
#' @param maxEd5p description
#' @param maxEd3p description
#' @param minTpm description
#' @return IsomirDataSet
#'
#' @export
detectIsomirs <- function(sampleInfoFile,
                          spe,
                          wordSize = 13,
                          maxEd5p = 2,
                          maxEd3p = 3,
                          minTpm = 2) {
  mirnas <- loadMirna(spe, wordSize)
  metaData <- read.csv(sampleInfoFile,
                       stringsAsFactors = FALSE,
                       comment.char = "#")

  dupName <- metaData$name[duplicated(metaData$name)]
  if (length(dupName) != 0) {
    stop("Duplicated sample name: ", paste(dupName, collapse = ", "))
  }

  dupFq <- metaData$fq[duplicated(metaData$fq)]
  if (length(dupFq) != 0) {
    stop("Duplicated FASTQ files:", paste(dupFq, collapse = ", "))
  }

  noFq <- metaData$fq[!file.exists(metaData$fq)]
  if (length(noFq) != 0) {
    stop("FASTQ files do not exists: ", paste(noFq, collapse = ", "))
  }

  row.names(metaData) <- metaData$name
  group <- split(metaData, metaData$group)

  repNum <- sapply(group, function(x)
    nrow(x))
  repOne <- repNum[repNum < 2]
  if (length(repOne) != 0) {
    stop(
      "At least two replicates are required for tissu/treatment: ",
      paste(names(repOne), collapse = ", ")
    )
  }

  isoformList <- lapply(
    group,
    detectOneTissue,
    mirnas = mirnas,
    maxEd5p = maxEd5p,
    maxEd3p = maxEd3p,
    minTpm = minTpm
  )

  message("Generating expression profiile from isomiR list")
  expr <- generateExpr(isoformList)

  message("Clustriing isoforms based on reference mature microRNAs")
  isomirList <- clusterIsoforms(isoformList)

  new(
    "IsomirDataSet",
    isoformList = isoformList,
    expr = expr,
    isomirList = isomirList
  )
}

#' ====================== private functions ====================================

#' Detect isoforms from a tissue consist of several duplicates
#'
#' @param sampleInfo A data.frame with sample name, fq and group information
#'   for replicates in a tissu/treatment
#' @param mirnas A data.frame storing mature miRNA ID, mature miRNA sequence,
#'   mature start, seed on miRNA, precursor ID and precursor sequence
#' @param maxEd5p the maximum distance between 5’ region of reads with
#'   reference miRNA
#' @param maxEd3p the maximum distance between 3' region of reads with
#'   reference miRNA
#' @param minTpm a isomiR is defined as “expressed” only when the TPM of all
#'    biological replicates were greater than or equal to the threshold set
#' @return a data.frame
#'
detectOneTissue <- function(sampleInfo,
                            mirnas,
                            maxEd5p,
                            maxEd3p,
                            minTpm) {
  tissueIsoforms <- apply(sampleInfo, 1, function(x) {
    sampleName <- x["name"]
    fqFile <- x["fq"]
    message("detecting isoform for sample: ", sampleName)
    isoforms <- detectOneSample(fqFile, mirnas, maxEd5p, maxEd3p)
    isoforms$TPM <- calcTpm(isoforms)
    isoforms <- subset(isoforms, TPM >= minTpm)
    cigars <- calcCigar(isoforms)
    isoforms$CIGAR <- cigars
    isoforms$ID <- paste0(isoforms$mature_ID, "-", cigars)
    isoforms
  })

  # read_seqs <- tissue_isoforms[[1]]$read_seq
  # for (i in seq_along(tissue_isoforms)[-1]) {
  #   read_seqs <- intersect(read_seqs, tissue_isoforms[[i]]$read_seq)
  # }
  #
  # lapply(tissue_isoforms, function(isoforms) {
  #   isoforms <- subset(isoforms, read_seq %in% read_seqs)
  #   cigars <- calc_cigar(isoforms)
  #   isoforms$id <- paste0(isoforms$mature_id, "-", cigars)
  #   isoforms
  # })

  isoforms <- tissueIsoforms[[1]]
  for (i in seq_along(tissueIsoforms)[-1]) {
    isoformsTemp <- tissueIsoforms[[i]][, c("read_seq", "TPM")]
    isoforms <- merge(x = isoforms, y = isoformsTemp, by = "read_seq")
    isoforms$TPM <- rowMeans(isoforms[, c("TPM.x", "TPM.y")])
    isoforms$TPM.x <- NULL
    isoforms$TPM.y <- NULL
  }

  return(isoforms)
}

clusterIsoforms <- function(isoformList) {
  dfList <- lapply(isoformList, function(x) {
    x[, c("ID", "mature_ID", "read_seq", "dist", "CIGAR")]
  })

  df <- dplyr::bind_rows(dfList)
  df <- dplyr::distinct(df)

  group <- split(df, df$mature_ID)
  lapply(group, function(x) {
    x <- dplyr::arrange(x, desc(dist))
    new(
      "Isomir",
      matureId = x$mature_ID[1],
      readSeqs = x$read_seq,
      isoformIds = x$ID,
      dist = x$dist,
      cigars = x$CIGAR
    )
  })
}

#' Extract isomiRs from a list containing all isoforms
#'
#' @param isoformList a list containing all isoforms
#' @return a list of all isomiRs
generateIsomirs <- function(isoformList) {
  lapply(isoformList, function(isoforms) {
    group <- dplyr::group_by(isoforms, read_seq)
    group <- dplyr::mutate(group, IDs = paste(make.unique(sort(ID)), collapse = ","))
    # group <- dplyr::mutate(group, ID = sort(IDs)[1])
    group <- dplyr::distinct(group, read_seq, .keep_all = TRUE)
    as.data.frame(dplyr::ungroup(group))
  })
}

#' Extract expression profile of each sample replicate from isomiRs
#'
#' @param isomirList isomiRs of all samples
#' @return a data.frame of expression
generateExpr <- function(isoformList) {
  tissueName <- names(isoformList[1])
  expr <- isoformList[[1]][, c("read_seq", "TPM")]
  expr <- dplyr::distinct(expr, read_seq, .keep_all = TRUE)
  colnames(expr)[colnames(expr) == "TPM"] <- tissueName
  for (i in seq_along(isoformList)[-1]) {
    tissueName <- names(isoformList[i])
    exprTemp <- isoformList[[i]][, c("read_seq", "TPM")]
    exprTemp <- dplyr::distinct(exprTemp, read_seq, .keep_all = TRUE)
    colnames(exprTemp)[colnames(exprTemp) == "TPM"] <- tissueName
    expr <- merge(x = expr, y = exprTemp, by = "read_seq", all = TRUE)
  }

  row.names(expr) <- expr$read_seq
  expr$read_seq <- NULL
  expr[is.na(expr)] = 0

  return(expr)
}
