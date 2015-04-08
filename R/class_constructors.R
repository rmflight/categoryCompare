#' creates enriched result
#' 
#' given all the slots for an \linkS4Class{enriched_result}, checks that all
#' the data is self-consistent, and creates the \code{enriched_result} object.
#' 
#' @param features the features that were differentially expressed (see details)
#' @param universe all of the features that were measured
#' @param annotation an \linkS4Class{annotation} object
#' @param statistics a \linkS4Class{statistical_results} object
#' 
#' @export
#' @return enriched_result
#' 
enriched_result <- function(features, universe, annotation, statistics){
  stopifnot(class(statistics) == "statistical_results")
  stopifnot(class(annotation) == "annotation")
  
  annotation_names <- names(annotation@annotation_features)
  stat_names <- statistics@annotation_id
  
  annotation_stat <- intersect(annotation_names, stat_names)
  
  if (sum(stat_names %in% annotation_stat) != length(stat_names)){
    stop("There are missing annotations!", call. = TRUE)
  }
  
  annotation@annotation_features <- annotation@annotation_features[stat_names]
  if (length(annotation@description) != 0){
    annotation@description <- annotation@description[stat_names]
  }
  
  if (length(annotation@counts) != 0){
    annotation@counts <- annotation@counts[stat_names]
  }
  
  if (length(annotation@links) != 0){
    annotation@links <- annotation@links[stat_names]
  }
  
  return(new("enriched_result",
             features = features,
             universe = universe,
             annotation = annotation,
             statistics = statistics))
}