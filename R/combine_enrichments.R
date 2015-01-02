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
  
  annotation_graph <- generate_annotation_similarity_graph(all_annotation)
  
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

#' annotation similarity graph
#' 
#' given an annotation-feature list, generate a similarity graph between all of 
#' the annotations
#' 
#' @param annotation_features list where each entry is a set of features to that annotation
#' @param overlap_type which type of overlap coefficient to report
#' 
#' @export
#' @return graphNEL
#' 
#' @import graph
generate_annotation_similarity_graph <- function(annotation_features, overlap_type = "combined"){
  num_features <- sapply(annotation_features, length)
  
  keep_annotations <- (num_features >= 5) & (num_features <= 1000)
  annotation_features <- annotation_features[keep_annotations]
  
  use_annotations <- names(annotation_features)
  n_annotation <- length(use_annotations)
  out_graph <- new("graphNEL", nodes = use_annotations, edgemode = "directed")
  
  all_comparisons <- expand.grid(seq(1, n_annotation), seq(1, n_annotation))
  all_comparisons <- all_comparisons[(all_comparisons[,2] > all_comparisons[,1]), ]
    
  similarity <- sapply(seq(1, nrow(all_comparisons)), function(x){
    do_comparison <- all_comparisons[x, ]
    n1 <- annotation_features[[do_comparison[1]]]
    n2 <- annotation_features[[do_comparison[2]]]
    
    use_similarity <- switch(overlap_type,
                             overlap = overlap_coefficient(n1, n2),
                             jaccard = jaccard_coefficient(n1, n2),
                             combined = combined_coefficent(n1, n2))
    
    use_similarity
  })
  
  similarity_non_zero <- similarity != 0
  all_comparisons <- all_comparisons[, similarity_non_zero]
  similarity <- similarity[similarity_non_zero]
  
  from_edge <- use_annotations[all_comparisons[,1]]
  to_edge <- use_annotations[all_comparisons[,2]]
  
  out_graph <- addEdge(from_edge, to_edge, out_graph, similarity)
  return(out_graph)
}

#' overlap coefficient
#' 
#' calculates the similarity using the "overlap" coefficient, which is
#' 
#' length(intersect(n1, n2)) / length(union(n1, n2))
#' 
#' @param n1 group 1 of objects
#' @param n2 group 2 of objects
#' 
#' @return double
#' @export
overlap_coefficient <- function(n1, n2){
  length(intersect(n1, n2)) / length(union(n1, n2))
}

#' jaccard coefficient
#' 
#' calculates similarity of two groups of objects using "jaccard" coefficient,
#' defined as:
#' 
#' length(intersect(n1, n2)) / min(c(length(n1), length(n2)))
#' 
#' @param n1 group 1
#' @param n2 group 2
#' 
#' @return double
#' @export
jaccard_coefficient <- function(n1, n2){
  length(intersect(n1, n2)) / min(c(length(n1), length(n2)))
}

#' combined coefficient
#' 
#' takes an average of the \code{overlap} and \code{jaccard} coefficients
#' 
#' @param n1 group 1
#' @param n2 group 2
#' 
#' @return double
#' @export
combined_coefficient <- function(n1, n2){
  o_coef <- overlap_coefficient(n1, n2)
  j_coef <- overlap_coefficient(n1, n2)
  return((0.5 * o_coef) + (0.5 * j_coef))
}
