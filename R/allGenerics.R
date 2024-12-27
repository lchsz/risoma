#' Cluster and visualize expression profile
#'
#' @param x description
#'
#' @export
setGeneric("plotExpr", function(x)
  standardGeneric("plotExpr"))


#' @export
setGeneric("exIsomirByRef", function(x, ref)
  standardGeneric("exIsomirByRef"))


#' @export
setGeneric("hasRef", function(x)
  standardGeneric("hasRef"))


#' @export
setGeneric("getRef", function(x)
  standardGeneric("getRef"))


#' @export
setGeneric("exExpr", function(x, ref)
  standardGeneric("exExpr"))


#' @export
setGeneric("alnIsoforms", function(x, ref)
  standardGeneric("alnIsoforms"))
