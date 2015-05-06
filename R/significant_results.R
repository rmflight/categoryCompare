#' index a list
#' 
#' Provided a list, and a condition, returns the logical indices into the named
#' part of the list provided. Uses \code{subset} like non-standard evaluation
#' so that we can define appropriate expressions.
#' 
#' @param list_to_query the list to run the query on
#' @param ... the expressions that do the queries
#' 
#' @export
#' @return logical "&" of all queries
#' 
multi_query_list <- function(list_to_query, ...){
  queries <- as.list(substitute(list(...)))[-1L]
  
  n_query <- length(queries)

  stopifnot(length(queries) > 0)
  
  
  query_results <- lapply(queries, eval, list_to_query)
  
  # how many objects do we have in each query
  # they should all be the same to allow merging
  n_objects <- unique(unlist(lapply(query_results, length)))
  
  if (length(n_objects) != 1){
    stop("cannot merge queries on objects of different length", call. = FALSE)
  }
  
  result_logical <- rep(TRUE, n_objects)
  
  for (i_query in seq(1, n_query)){
    result_logical <- result_logical & query_results[[i_query]]
  }
  
  names(result_logical) <- NULL
  result_logical
}

#' get significant annotations
#' 
#' given a \linkS4class{statistical_results} object and some conditional expressions,
#' return the significant annotations
#' 
#' @param in_results the \linkS4class{statistical_results} object
#' @param ... conditional expressions
#' 
#' @examples
#' 
#' test_stat <- new("statistical_results",
#'                  annotation_id = c("a1", "a2", "a3"),
#'                  statistics = list(pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
#'                    counts = c(a1 = 5, a2 = 10, a3 = 1),
#'                    odds = c(a1 = 20, a2 = 100, a3 = 0)))
#' get_significant_annotations(test_stat, pvalues < 0.05)
#' get_significant_annotations(test_stat, odds > 10)
#' get_significant_annotations(test_stat, pvalues < 0.05, counts >= 1)
#' 
#' @return vector of significant annotation_id's
#' @exportMethod get_significant_annotations
setMethod("get_significant_annotations", 
          signature = list(in_results = "statistical_results"),
          function(in_results, ...) .get_significant_stat_results(in_results, ...))

.get_significant_stat_results <- function(in_results, ...){
  queries <- as.list(substitute(list(...)))[-1L]
  stopifnot(length(queries) > 0)
  
  out_ids <- in_results@annotation_id
  
  sig_entries <- multi_query_list(in_results@statistic_data, ...)
  
  out_ids[sig_entries]
}

#' get significant annotations
#' 
#' In the case where we have a \linkS4class{combined_enrichment} and we want
#' to get all of the significant annotations from each of them, and put them
#' together so we can start doing real meta-analysis.
#' 
#' Note that this function returns the original \linkS4class{combined_enrichment} object with a modified
#' \linkS4class{combined_statistics} slot where the significant annotations have been added in. 
#' 
#' @param in_results a \linkS4class{combined_enrichment} object
#' @param ... conditional expressions
#' 
#' @return \linkS4class{combined_enrichment} object
#' @exportMethod get_significant_annotations
setMethod("get_significant_annotations",
          signature = list(in_results = "combined_enrichment"),
          function(in_results, ...) .get_significant_combined_enrichment(in_results, ...))

.get_significant_combined_enrichment <- function(in_results, ...){
  queries <- as.list(substitute(list(...)))[-1L]
  stopifnot(length(queries) > 0)
  
  all_measured <- lapply(in_results@enriched,
                         function(x){x@statistics@annotation_id})
  
  all_significant <- lapply(in_results@enriched,
                            function(x){get_significant_annotations(x@statistics, ...)})
  
  annotation_measured <- unique(unlist(all_measured))
  n_measured <- length(annotation_measured)
  n_enriched <- length(in_results@enriched)
  
  out_measured <- matrix(FALSE, n_measured, n_enriched)
  rownames(out_measured) <- annotation_measured
  colnames(out_measured) <- names(all_measured)
  
  out_significant <- out_measured
  
  for (i_meas in names(all_measured)){
    out_measured[all_measured[[i_meas]], i_meas] <- TRUE
  }
  
  for (i_meas in names(all_significant)){
    out_significant[all_significant[[i_meas]], i_meas] <- TRUE
  }
  
  sig_annotation <- new("significant_annotations",
                        significant = out_significant,
                        measured = out_measured,
                        sig_calls = sapply(queries, deparse))
  
  in_results@statistics@significant <- sig_annotation
  
  in_results
}

#' filter graph by significant entries
#' 
#' If a graph has already been generated, it may be faster to filter a previously
#' generated one than generate a new one from significant data.
#' 
#' @param in_graph the \linkS4Class{cc_graph} previously generated
#' @param comb_enrich the \linkS4Class{combined_enrichment} that you want to use to filter with
#' 
#' @export
#' @importFrom graph subGraph nodes
#' @return cc_graph
filter_annotation_graph <- function(in_graph, comb_enrich){
  sig_matrix <- comb_enrich@statistics@significant@significant
  
  annotation_list <- rownames(sig_matrix)
  
  sig_entries <- rowSums(sig_matrix) > 0
  
  keep_annotation <- annotation_list[sig_entries]
  
  # use intersect in case there is something odd of the graph and significant entries
  keep_intersect <- intersect(keep_annotation, nodes(in_graph))
  
  if (length(keep_intersect) > 0){
    in_graph <- subGraph(keep_intersect, in_graph)
    out_graph <- as(in_graph, "cc_graph")
    out_graph@significant <- sig_matrix[keep_intersect, ]
    
    return(out_graph)
  } else {
    warning("No matching nodes and annotations found", call. = TRUE)
  }
  
}
