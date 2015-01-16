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
setMethod("combine_enrichments", signature = list(annotation_similarity = "character"), 
          function(annotation_similarity, ...) .combine_enrichments(annotation_similarity, ...))

setMethod("combine_enrichments", signature = "missing",
          function(...) .combine_enrichments(annotation_similarity = "combined", ...))

.combine_enrichments <- function(annotation_similarity = "combined", ...){
  enriched <- list(...)
  enriched_class <- unique(sapply(enriched, class))
  if (enriched_class != "enriched_result"){
    stop("objects passed are not all of class enriched_result", call. = FALSE)
  }
  
  all_annotation <- combine_annotations(lapply(enriched, function(x){x@annotation}))
  
  annotation_graph <- generate_annotation_similarity_graph(all_annotation@annotation_features, annotation_similarity)
  
  # create a dummy that may be a little smaller for passing to other functions
  out_combined <- new("combined_enrichment",
                      enriched = enriched)
  
  combined_stats <- extract_statistics(out_combined)
  
  out_combined <- new("combined_enrichment",
                      enriched = enriched,
                      annotation = all_annotation,
                      graph = annotation_graph,
                      statistics = combined_stats)
  
  out_combined
}

#' combine annotations
#' 
#' Takes multiple \linkS4class{annotation} objects and combines them so that there
#' is a consistent sole set for creating the \code{annotation_graph} and providing
#' other information about each annotation entry.
#' 
#' @param ... one or more \linkS4class{annotation}
#' 
#' @return \linkS4class{annotation}
#' @export
setMethod("combine_annotations", signature = "list", function(annotation_list) .combine_annotations(annotation_list))

.combine_annotations <- function(annotation_list){
  
  combined_type <- unique(unlist(lapply(annotation_list, function(x){x@type})))
  
  # stop if there are more than one type
  n_type <- length(combined_type)
  if (n_type != 1){
    stop("Cannot combine annotation's with more than one annotation type.", call.=FALSE)
  }
  
  combined_features <- combine_annotation_features(lapply(annotation_list, function(x){x@annotation_features}))
  
  all_annotation_names <- names(combined_features)
  
  combined_description <- combine_text(lapply(annotation_list, function(x){x@description}),
                                       all_annotation_names,
                                       "description")
  combined_links <- combine_text(lapply(annotation_list, function(x){x@links}),
                                 all_annotation_names,
                                 "links")
  
  combined_counts <- sapply(combined_features, length)
  
  out_annotation <- new("annotation",
                        annotation_features = combined_features,
                        description = combined_description,
                        counts = combined_counts,
                        links = combined_links,
                        type = combined_type)
  
  out_annotation
  
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

#' combine text
#' 
#' Given lists of named character objects, and a character vector of names to be
#' in the final object, either get the character string from the list that has
#' the names, or check that the character string is the same across all of the
#' lists.
#' 
#' @param list_characters list containing named character strings
#' @param names_out the full list of names to use
#' 
#' @export
#' @return named character vector
combine_text <- function(list_characters, names_out, text_id){
  n_char <- length(names_out)
  list_names <- names(list_characters)
  n_list <- length(list_characters)
  
  list_text <- matrix("", n_char, n_list)
  rownames(list_text) <- names_out
  colnames(list_text) <- names(list_characters)
  
  for (i_name in list_names){
    match_list <- intersect(names_out, names(list_characters[[i_name]]))
    list_text[match_list, i_name] <- list_characters[[i_name]][match_list]
  }
  
  list_has_char <- nchar(list_text) > 0
  
  out_char <- lapply(names_out, function(x){
    unique(list_text[x, list_has_char[x,]])
  })
  
  names(out_char) <- names_out
  
  n_result <- sapply(out_char, length)
  
  n_multi <- sum(n_result > 1)
  if (n_multi != 0){
    errmsg <- paste(n_multi, text_id, "are not identical across annotation objects", collapse = " ")
    stop(errmsg, call. = FALSE)
  }
  
  unlist(out_char)
  
}

#' annotation similarity graph
#' 
#' given an annotation-feature list, generate a similarity graph between all of 
#' the annotations
#' 
#' @param annotation_features list where each entry is a set of features to that annotation
#' @param similarity_type which type of overlap coefficient to report
#' 
#' @export
#' @return graphNEL
#' 
#' @import graph
generate_annotation_similarity_graph <- function(annotation_features, similarity_type = "combined"){
  num_features <- sapply(annotation_features, length)
  
  keep_annotations <- (num_features >= 5) & (num_features <= 1000)
  annotation_features <- annotation_features[keep_annotations]
  
  use_annotations <- names(annotation_features)
  n_annotation <- length(use_annotations)
  out_graph <- new("graphNEL", nodes = use_annotations, edgemode = "directed")
  
  all_comparisons <- expand.grid(seq(1, n_annotation), seq(1, n_annotation))
  all_comparisons <- all_comparisons[(all_comparisons[,2] > all_comparisons[,1]), ]
    
  similarity <- choose_mclapply(seq(1, nrow(all_comparisons)), function(x){
    do_comparison <- all_comparisons[x, ]
    n1 <- annotation_features[[do_comparison[1,1]]]
    n2 <- annotation_features[[do_comparison[1,2]]]
    
    use_similarity <- switch(similarity_type,
                             overlap = overlap_coefficient(n1, n2),
                             jaccard = jaccard_coefficient(n1, n2),
                             combined = combined_coefficient(n1, n2))
    
    use_similarity
  })
  similarity <- unlist(similarity)
  similarity_non_zero <- similarity != 0
  all_comparisons <- all_comparisons[similarity_non_zero, ]
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
  j_coef <- jaccard_coefficient(n1, n2)
  return((0.5 * o_coef) + (0.5 * j_coef))
}


#' extract statistics
#' 
#' extract all statistics for a \linkS4class{statistical_results} object. These
#' can then be combined into a \code{data.frame} that can be returned or used
#' to annotate the graph of annotations.
#' 
#' @param in_results the \linkS4class{statistical_results} object
#' @return data.frame
#' @exportMethod extract_statistics
setMethod("extract_statistics", signature = list(in_results = "statistical_results"),
          function(in_results) .extract_statistics_statistical_results(in_results))

.extract_statistics_statistical_results <- function(in_results){
  out_data <- as.data.frame(in_results@statistic_data)
  row.names(out_data) <- in_results@annotation_id
  
  out_data
}

#' extract statistics
#' 
#' extract all statistics from a \linkS4class{combined_enrichment} object and 
#' create a \linkS4class{combined_statistics} where each statistic from the underlying 
#' \linkS4class{statistical_results} object in each of the enrichments
#' is named according to which enrichment it was in and what statistic it was.
#' 
#' @param in_results the \linkS4class{combined_enrichment} object
#' @return combined_statistics
#' @exportMethod extract_statistics
setMethod("extract_statistics", signature = list(in_results = "combined_enrichment"),
          function(in_results) .extract_statistics_combined_enrichment(in_results))

.extract_statistics_combined_enrichment <- function(in_results){
  all_data <- lapply(in_results@enriched, function(x){extract_statistics(x@statistics)})
  
  all_rows <- unique(unlist(lapply(all_data, rownames)))
  
  enrichment_names <- rep(names(all_data), unlist(lapply(all_data, ncol)))
  
  statistic_names <- unlist(lapply(all_data, names), use.names = FALSE)
  
  new_named_data <- lapply(names(all_data), function(x){
    tmp_x <- all_data[[x]]
    names(tmp_x) <- paste(x, names(tmp_x), sep = ".")
    tmp_x
  })
  names(new_named_data) <- names(all_data)
  
  all_cols <- (unlist(lapply(new_named_data, colnames)))
  
  data_matrix <- matrix(NA, length(all_rows), length(all_cols))
  data_frame <- as.data.frame(data_matrix)
  
  rownames(data_frame) <- all_rows
  names(data_frame) <- all_cols
  
  for (i_data in names(new_named_data)){
    tmp_data <- new_named_data[[i_data]]
    data_frame[rownames(tmp_data), names(tmp_data)] <- tmp_data
  }
  
  new("combined_statistics",
      statistic_data = data_frame,
      annotation_id = rownames(data_frame),
      which_enrichment = enrichment_names,
      which_statistic = statistic_names)
}

#' all data extraction
#' 
#' given a \linkS4class{combined_enrichment} object, get the significant and present
#' annotations, their descriptions, links (if present), and all the statistics, and 
#' generate an object that can be used to 