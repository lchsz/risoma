#' An S4 class to represent a dataset containing all the isomiR information
#'   detected from RNA-Seq samples.
#'
#' @slot isoformList A list of isoforms
#' @slot expr A data.fram of expression profile
#' @slot isomirList A list of isomiRs
#'
#' @export
setClass(
  "IsomirDataSet",
  slots = c(
    sample_info = "ANY",
    group2isoforms = "ANY",
    sample_expr = "ANY",
    group_expr = "ANY",
    ref2isomir = "ANY"
  )
)


#' An S4 class to represent a isomiR
#'
#' @slot mature_id microRNA ID
#' @slot read_seqs Vector of read sequences
#' @slot dist vector of edit distance
#' @slot cigars vector of CIGAR
#'
#' @export
setClass(
  "Isomir",
  slots = c(
    mature_id = "character",
    mature_seq = "character",
    template_seq = "character",
    read_seqs = "character",
    dist = "integer",
    dist_5p = "integer",
    dist_3p = "integer",
    cigars = "character",
    ids = "character"
  )
)

