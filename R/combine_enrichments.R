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
  
}