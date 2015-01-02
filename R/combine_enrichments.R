#' combine enrichments
#' 
#' This is one of the primary workhorse functions behind \pkg{categoryCompare}.
#' The primary function of \code{categoryCompare} is to enable \emph{comparisons}
#' of different enrichment analyses. To facilitate that, we must first 
#' \strong{combine} one (really, we can do this with a single) or more 
#' \linkS4class{enriched_results}
#' 
#' @param ... one or more \linkS4class{enriched_results}
#' 
#' @return \linkS4class{combined_enrichment}
#' @export
setMethod("combine_enrichments", signature = "enriched_result", function(...) .combine_enrichments(...))

.combine_enrichments <- function(...){
  enriched <- list(...)
  
  enriched_types <- lapply(enriched, function(x){x@annotation@type})
  
  # stop if there are more than one type
  n_type <- length(unique(unlist(enriched_types)))
  if (n_type != 1){
    stop("Cannot combine enriched_result's with more than one annotation type.", call.=FALSE)
  }
  
  all_annotation <- combine_annotation_features(lapply(enriched, function(x){x@annotation@annotation_features}))
  
}

#' combine annotation-features
#' 
#' For the generation of a proper annotation-annotation relationship graph, we
#' need to combine the annotation-feature relationships across multiple
#' \linkS4class{annotation} objects
#' 
#' @param annotation_features list of annotation_features to combine
#' 
#' @export
#' @return list of combined annotations
#' 
combine_annotation_features <- function(annotation_features){
  annotation_names <- lapply(annotation_features, function(x){names(x)})
  annotation_names <- unique(unlist(annotation_names))
  
  annotation_out <- vector("list", length(annotation_names))
  names(annotation_out) <- annotation_names
  
  for (i_annot in seq(1, length(annotation_features))){
    use_names <- names(annotation_features[[i_annot]])
    annotation_out[use_names] <- lapply(use_names, function(x){union(annotation_out[[x]], annotation_features[[i_annot]][[x]])})
  }
  
  return(annotation_out)
}