setGeneric("combine_enrichments", function(...) standardGeneric("combine_enrichments"))

setGeneric("combine_annotations", function(annotation_list) standardGeneric("combine_annotations"))

#' get significant annotations
#' 
#' @export
setGeneric("get_significant_annotations", function(in_results, ...) standardGeneric("get_significant_annotations"))

#' get statistics
#' 
#' @export
setGeneric("extract_statistics", function(in_results) standardGeneric("extract_statistics"))

#' generate annotation graph
#' 
#' @export
setGeneric("generate_annotation_graph", function(comb_enrichment, annotation_similarity = "combined", low_cut = 5, hi_cut = 500) standardGeneric("generate_annotation_graph"))