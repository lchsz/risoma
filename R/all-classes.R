#' An S4 class to represent a dataset containing all the isomiR information
#' detected from RNA-Seq samples.
#'
#' @slot sample_info A data.frame containing sample information
#' @slot group2isoforms A list of isoforms grouped by sample groups
#' @slot sample_expr A data.frame of expression profile for samples
#' @slot group_expr A data.frame of expression profile for groups
#' @slot ref2isomir A list of isomiRs grouped by reference mature microRNAs
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


#' An S4 class to represent an isomiR
#'
#' @slot mature_id: Reference miRNA ID.
#' @slot mature_seq: Reference miRNA sequence.
#' @slot seed_seq: The seed sequence of the mature miRNA.
#' @slot template_seq: Template sequence including flanking regions.
#' @slot read_seqs: Character vector of read sequences.
#' @slot read_num: Abundance of the read.
#' @slot indel_5p: Indels at the 5' end.
#' @slot indel_3p: Indels at the 3' end.
#' @slot snp_5p: Mismatches at the 5' end.
#' @slot snp_3p: Mismatches at the 3' end.
#' @slot snp_seed: Mismatches in the seed.
#' @slot dist: Total edit distance.
#' @slot cigars: CIGAR string representing alignment details.
#' @slot ids: Character vector of unique isoform IDs.
#' @slot aln: Multiple sequence alignment (MSA) of the isoform, reference miRNA
#'            and template sequences
#'
#' @export
setClass(
  "Isomir",
  slots = c(
    mature_id = "character",
    mature_seq = "character",
    seed_seq = "character",
    template_seq = "character",
    read_seqs = "character",
    indel_5p = "integer",
    indel_3p = "integer",
    snp_5p = "integer",
    snp_3p = "integer",
    snp_seed = "integer",
    dist = "integer",
    cigars = "character",
    ids = "character",
    aln = "ANY"
  )
)

