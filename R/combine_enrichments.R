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
setMethod("combine_enrichments", signature = "enriched_result", 
          function(...) .combine_enrichments(...))

.combine_enrichments <- function(...){
  enriched <- list(...)
  
  all_annotation <- combine_annotations(lapply(enriched, function(x){x@annotation}))
  
  # create a dummy that may be a little smaller for passing to other functions
  out_combined <- new("combined_enrichment",
                      enriched = enriched)
  
  combined_stats <- extract_statistics(out_combined)
  
  out_combined <- new("combined_enrichment",
                      enriched = enriched,
                      annotation = all_annotation,
                      statistics = combined_stats)
  
  out_combined
}

#' generate the annotation graph
#' 
#' given a \linkS4class{combined_enrichment}, generate the annotation similarity graph
#' 
#' @param comb_enrichment the combined_enrichment object
#' @param annotation_similarity which similarity measure to use
#' @param low_cut keep only those annotations in the graph with at least this many annotated features
#' @param hi_cut keep only those annotations with less than this many annotated features
#' 
#' @return \linkS4class{cc_graph}
#' 
#' @export
setMethod("generate_annotation_graph", signature = list(comb_enrichment = "combined_enrichment"),
          function(comb_enrichment, annotation_similarity, low_cut, hi_cut) .generate_annotation_graph(comb_enrichment, annotation_similarity, low_cut, hi_cut))

.generate_annotation_graph <- function(comb_enrichment, annotation_similarity = "combined", low_cut = 5, hi_cut = 500){
  
  keep_features <- comb_enrichment@statistics@annotation_id
  annotation_features <- comb_enrichment@annotation@annotation_features[keep_features]
  n_features <- sapply(annotation_features, length)
  
  keep_annotations <- names(annotation_features)[(n_features >= low_cut) & (n_features <= hi_cut)]
  
  annotation_table <- generate_table(comb_enrichment, link_type = "explicit")
  
  in_graph_annotation <- intersect(keep_annotations, rownames(annotation_table))
  
  annotation_graph <- generate_annotation_similarity_graph(annotation_features[in_graph_annotation], annotation_similarity)
  
  annotation_graph <- add_data_to_graph(annotation_graph, annotation_table)
  
  annotation_graph@significant <- comb_enrichment@statistics@significant@significant[in_graph_annotation, , drop = FALSE]
  annotation_graph
}

#' add table data to graph
#' 
#' given the annotation_graph and a data.frame, add all of the data in the data.frame
#' to the graph so it is available elsewhere. Note that for NA integer and numerics,
#' the value is modified to -100, and for infinite values, it is modified to 1e100.
#' 
#' @param graph the graph to work on
#' @param data the data to add to it
#' @import graph
#' 
#' @return graphNEL
add_data_to_graph <- function(graph, data){
  type_convert <- c('STRING','INTEGER','FLOATING','STRING')
  type_defaults <- list(character = "NA", integer = -100, numeric = -100, logical = "NA")
  names(type_convert) <- c('character','integer','numeric','logical')
  
  data_types <- lapply(data, class)
  
  graph_entries <- nodes(graph)
  data_entries <- rownames(data)
  
  match_entries <- intersect(graph_entries, data_entries)
  
  for (i_data in names(data_types)){
    use_type <- data_types[[i_data]][1]
    nodeDataDefaults(graph, i_data) <- type_defaults[[use_type]][1]
    attr(nodeDataDefaults(graph, i_data), "class") <- type_convert[use_type]
    
    # NA is not nice in RCy, so convert to the default, which is -100, pretty non-sensical
    tmp_data <- data[match_entries, i_data]
    if (use_type %in% c("integer", "numeric")){
      tmp_data[is.na(tmp_data)] <- type_defaults[[use_type]][1]
      tmp_data[is.infinite(tmp_data)] <- 1e100
    }
    
    nodeData(graph, match_entries, i_data) <- tmp_data
  }
  
  graph
}

#' generate table
#' 
#' given a \linkS4class{combined_enrichment} object, get out the data.frame
#' either for investigation or to add data to the \linkS4class{annotation_graph}.
#' 
#' @param comb_enrichment the combined_enrichment object
#' @param link_type should their be an "explicit" link (see details)
#' @details the \code{link_type} controls whether to create an "explicit" link
#'   that is actually a column in the data.frame, or create an "implicit" html link
#'   that is part of the \code{@@name} column in the returned data.frame. Useful
#'   if you are embedding the data.frame in an html report.
#' @return data.frame
#' @export
setMethod("generate_table", signature = list(comb_enrichment = "combined_enrichment"),
          function(comb_enrichment, link_type) .generate_table(comb_enrichment, link_type))

.generate_table <- function(comb_enrichment, link_type = "explicit"){
  # get the bits we need
  # from the combined_statistcs we take the statistic_data and significant_annotations,
  # and from the annotation slot we take the description and links, and we put this
  # all together into a single data.frame
  
  base_data <- comb_enrichment@statistics@statistic_data
  which_enrichment <- names(comb_enrichment@enriched)
  n_enrichment <- length(which_enrichment)
  sig_data <- comb_enrichment@statistics@significant@significant
  meas_data <- comb_enrichment@statistics@significant@measured
  
  keep_data <- rowSums(sig_data) > 0
  
  # ideally we only keep things that were significant in at least one of the
  # enrichments we did, but if nothing is significant in any (implies not having
  # done sig cutoffs yet), then we will just return everything
  if (sum(keep_data) == 0){
    keep_data <- rep(TRUE, length(keep_data))
  }
  
  base_data <- base_data[keep_data, ]
  sig_data <- sig_data[keep_data, ]
  meas_data <- meas_data[keep_data, ]
  
  keep_annot <- comb_enrichment@statistics@annotation_id[keep_data]
  
  annot_obj <- comb_enrichment@annotation
  has_desc <- length(annot_obj@description) > 0
  has_link <- length(annot_obj@links) > 0
  
  obj_desc <- data.frame(name = keep_annot, stringsAsFactors = FALSE)
  
  if (has_desc){
    if ((link_type != "explicit") & (has_link)){
      text_desc <- paste('<a href="', annot_obj@links[keep_annot], '">', annot_obj@desc[keep_annot], '</a>', sep = "")
      obj_desc$description <- text_desc
    } else if ((link_type == "explicit") & has_link){
      obj_desc$description <- annot_obj@description[keep_annot]
      obj_desc$link <- annot_obj@links[keep_annot]
    } else {
      obj_desc$description <- annot_obj@description[keep_annot]
    }
  }
  
  sig_cols <- as.data.frame(sig_data)
  colnames(sig_cols) <- paste(colnames(sig_cols), "sig", sep = ".")
  
  meas_cols <- as.data.frame(meas_data)
  colnames(meas_cols) <- paste(colnames(meas_cols), "meas", sep = ".")
  
  if (has_desc){
    out_data <- cbind(obj_desc[, c("name", "description")], sig_cols, meas_cols, base_data)
  } else {
    out_data <- cbind(obj_desc[, "name"], sig_cols, meas_cols, base_data)
  }
  
  if (exists("link", where = -1)){
    out_data <- cbind(out_data, obj_desc$link)
  }
  out_data
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
#' @return cc_graph
#' 
#' @import graph
generate_annotation_similarity_graph <- function(annotation_features, similarity_type = "combined"){
    
  use_annotations <- names(annotation_features)
  n_annotation <- length(use_annotations)
  out_graph <- new("cc_graph", nodes = use_annotations, edgemode = "directed")
  
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
  
  edgeDataDefaults(out_graph, attr = "weight") <- 0
  attr(edgeDataDefaults(out_graph, attr = "weight"), "class") <- "FLOATING"
  
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

#' combined_statistics
#' 
#' constructor function for the combined_statistics object, makes sure that
#' empty things get initialized correctly
#' 
#' @param statistic_data the data.frame of statistics
#' @param which_enrichment which enrichment gave the results
#' @param which_statistic which statistics were calculated in each case
#' @param annotation_id the annotations for which we are returning statistics
#' @param significant the significant annotations
#' @param measured the measured annotations
#' @param use_names the order of naming
#' 
#' @export
#' @return combined_statistics
combined_statistics <- function(statistic_data, which_enrichment, which_statistic,
                                annotation_id, significant = NULL, measured = NULL, use_names = NULL){
  
  nul_sig <- is.null(significant)
  nul_meas <- is.null(measured)
  nul_names <- is.null(use_names)
  
  enrich_names <- unique(which_enrichment)
  
  if (nul_names){
    use_names <- enrich_names
  } else {
    use_names <- use_names[(use_names %in% enrich_names)]
  }
  
  if (nul_sig){
    significant <- matrix(FALSE, nrow = length(annotation_id), ncol = length(use_names))
    rownames(significant) <- annotation_id
    colnames(significant) <- use_names
  }
  
  if (nul_meas){
    measured <- matrix(FALSE, nrow = length(annotation_id), ncol = length(use_names))
    rownames(measured) <- annotation_id
    colnames(measured) <- use_names
  }
  
  sig_data <- new("significant_annotations",
                  significant = significant,
                  measured = measured)
  
  new("combined_statistics",
      statistic_data = statistic_data,
      which_enrichment = which_enrichment,
      which_statistic = which_statistic,
      annotation_id = annotation_id,
      significant = sig_data)
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
  
  base_enrich_names <- names(all_data)
  
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
  
 combined_statistics(statistic_data = data_frame,
                     annotation_id = rownames(data_frame),
                     which_enrichment = enrichment_names,
                     which_statistic = statistic_names,
                     use_names = base_enrich_names)
}

#' add data to graph
#' 
#' given a \linkS4class{combined_enrichment} object, add the data about the significant 
#' and present annotations, their descriptions, links (if present), and all the statistics, 
#' and add that data to the annotation graph object.
#' 
#' @param combined_enrichment a \linkS4class{combined_enrichment} object
#' @exportMethod