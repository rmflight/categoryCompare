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
  queries <- as.list(substitute(list(...)))
  
  browser(expr = TRUE)
}