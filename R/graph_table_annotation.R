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
  out_data <- as.data.frame(in_results@statistics)
  rownames(out_data) <- in_results@annotation_id
  
  out_data
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

.extract_statistics_combined_enrichment <- function(in_results){
  all_data <- lapply(in_results@enriched, function(x){extract_statistics(x@statistics)})
  
  all_rows <- unique(unlist(lapply(all_data, rownames)))
  
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
  
  data_frame
}