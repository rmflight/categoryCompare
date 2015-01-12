setGeneric("combine_enrichments", function(...) standardGeneric("combine_enrichments"))

setGeneric("combine_annotations", function(annotation_list) standardGeneric("combine_annotations"))

#' get significant annotations
#' 
#' @export
setGeneric("get_significant_annotations", function(combined_enrichment_or_stat_results, ...) standardGeneric("get_significant_annotations"))

#' get statistics
#' 
#' @export
setGeneric("extract_statistics", function(in_results) standardGeneric("extract_statistics"))