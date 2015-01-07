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
  
  result_logical
}