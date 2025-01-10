#' @export
setGeneric("has_ref", function(x, ref)
  standardGeneric("has_ref"))


#' @export
setGeneric("has_isomir", function(x)
  standardGeneric("has_isomir"))


#' @export
setGeneric("ex_isomir", function(x, ref)
  standardGeneric("ex_isomir"))


#' @export
setGeneric("ex_expr", function(x, ref)
  standardGeneric("ex_expr"))


#' @export
setGeneric("aln_isoforms", function(x, ref)
  standardGeneric("aln_isoforms"))

#' @export
setGeneric("get_groups", function(x, ref)
  standardGeneric("get_groups"))


#' @export
setGeneric("ex_ref_expr", function(x)
  standardGeneric("ex_ref_expr"))


#' @export
setGeneric("calc_tsi", function(x)
  standardGeneric("calc_tsi"))


#' @export
setGeneric("calc_deg", function(x, treatment, control, padj=0.05)
  standardGeneric("calc_deg"))


#' @export
setGeneric("has_deg", function(x, deg)
  standardGeneric("has_deg"))


#' @export
setGeneric("ex_deg", function(x, deg, ref)
  standardGeneric("ex_deg"))

