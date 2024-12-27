#' An S4 class to represent a dataset containing all the isomiR information
#'   detected from RNA-Seq samples.
#'
#' @slot isoformList A list of isoforms
#' @slot expr A data.fram of expression profile
#' @slot isomirList A list of isomiRs
setClass(
  "IsomirDataSet",
  slots = c(
    isoformList = "ANY",
    expr = "ANY",
    isomirList = "ANY"
  )
)


setClass(
  "Isomir",
  slots = c(
    matureId = "character",
    readSeqs = "character",
    isoformIds = "character",
    dist = "integer",
    cigars = "character"
  )
)

