

# Get isomiR information ======================================================

## get_isoform_num -----
#' @rdname get_isoform_num
#' @export
setGeneric("get_isoform_num", function(x)
  standardGeneric("get_isoform_num"))

## get_isomir_by_ref -------------
#' @rdname get_isomir_by_ref
#' @export
setGeneric("get_isomir_by_ref", function(x, ref)
  standardGeneric("get_isomir_by_ref"))

## get_alignment_by_ref ------
#' @rdname get_alignment_by_ref
#' @export
setGeneric("get_alignment_by_ref", function(x, ref)
  standardGeneric("get_alignment_by_ref"))

# Expression Data ===========================================================

## calc_group_expr_num -----------------
#' @rdname calc_group_expr_num
#' @export
setGeneric("calc_group_expr_num", function(x)
  standardGeneric("calc_group_expr_num"))

## get_ref_expr --------------
#' @rdname get_ref_expr
#' @export
setGeneric("get_ref_expr", function(x)
  standardGeneric("get_ref_expr"))

## get_expr_by_ref -------
#' @rdname get_expr_by_ref
#' @export
setGeneric("get_expr_by_ref", function(x, ref)
  standardGeneric("get_expr_by_ref"))


# Calculate Tissue Specificity Index =================================

## calc_tsi -----
#' @rdname calc_tsi
#' @export
setGeneric("calc_tsi", function(x)
  standardGeneric("calc_tsi"))


# DEG ====================================================================

## get_deg_data ---------------
#' @rdname get_deg_data
#' @export
setGeneric("get_deg_data", function(x,
                                treatment,
                                control)
  standardGeneric("get_deg_data"))
