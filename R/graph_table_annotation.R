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

.extract_statistics_statistical_results(in_results){
  
}

#' extract statistics
#' 
#' extract all statistics from a \linkS4class{combined_enrichment} object and 
#' create a \code{data.frame} where each statistic from the underlying 
#' \linkS4class{statistical_results} object in each of the enrichments
#' is named according to which enrichment it was in and what statistic it was.
#' 
#' @param in_results the \linkS4class{combined_enrichment} object
#' @return data.frame
#' @exportMethod extract_statistics
setMethod("extract_statistics", signature = list(in_results = "combined_enrichment"),
          function(in_results) .extract_statistics_combined_enrichment(in_results))

.extract_statistics_combined_enrichment(in_results){
  
}