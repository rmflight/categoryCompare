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

#' get signficant annotations
#' 
#' given a \linkS4class{statistical_results} object and some conditional expressions,
#' return the significant annotations
#' 
#' @param stat_results the statistical_results object
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
#' @export
#' @return vector of significant annotation_id's
get_significant_annotations <- function(stat_results, ...){
  out_ids <- stat_results@annotation_id
  
  sig_entries <- multi_query_list(stat_results@statistics, ...)
  
  out_ids[sig_entries]
}