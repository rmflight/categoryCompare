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

#' generate statistical table
#' 
#' @export
setGeneric("generate_table", function(comb_enrichment, link_type = "explicit") standardGeneric("generate_table"))

#' unique annotation combinations
#' 
#' @export
setGeneric("annotation_combinations", function(in_graph) standardGeneric("annotation_combinations"))
# we may unexport this eventually

setGeneric("remove_edges", function(edge_obj, cutoff, edge_attr = "weight", value_direction = "under") standardGeneric("remove_edges"))
