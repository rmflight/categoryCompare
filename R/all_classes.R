#' annotation class
#' 
#' This class holds an annotation object that defines how annotations relate to
#' features, as well as various pieces about each annotation
#' 
#' These objects may be created by hand, or may result from specific functions
#' to create them. Most notably, this package provides functions for creating
#' a Gene Ontology annotation.
#' 
#' @slot annotation_features list of annotation to feature relationships
#' @slot description character vector providing descriptive text about the annotation
#' @slot counts numeric vector of how many features are in each annotation
#' @slot links character vector defining html links for each annotation (may be empty)
#' @slot type a one word short description of the "type" of annotation
#' 
#' @export
annotation <- setClass("annotation",
                       slots = list(annotation_features = "list",
                                    description = "character",
                                    counts = "numeric",
                                    links = "character",
                                    type = "character"))

#' statistical results class
#' 
#' This class holds the part of an enrichment that is the statistical results.
#' It has two pieces, a \code{list} of \code{statistics} that is a named list
#' with the actual numerical results of applying the statistics. The other piece
#' is the \code{annotation_id} vector defining which entry in each vector of the
#' \code{statistics} is.
#' 
#' @slot statistics list of numerical statistics
#' @slot annotation_id vector of ids
#' 
#' @export
statistical_results <- setClass("statistical_results",
                                slots = list(statistics = "list",
                                             annotation_id = "ANY"))

#' the enriched results class
#' 
#' @slot features the "features" of interest, a vector of class ANY
#' @slot universe all of the "features" in the background
#' @slot annotation list giving the annotation to feature relationship
#' @slot statistics named list with various statistics of the enrichment results
#' 
#' @export
enriched_result <- setClass("enriched_result",
                            slots = list(features = "ANY",
                                         universe = "ANY",
                                         annotation = "annotation",
                                         statistics = "statistical_results"))

#' combined enrichments
#' 
#' The \code{combined_enrichment} class holds the results of combining several 
#' \linkS4class{enriched_result}s together, which includes the original 
#' \linkS4class{enriched_result}s, as well as the \code{annotation_graph}
#' and combined \linkS4class{annotation} objects.
#' 
#' @slot enriched list of enriched objects
#' @slot enriched_type character describing the enrichment annotation
#' @slot annotation \linkS4class{annotation} where the annotation_features
#' have been combined across the \linkS4class{enriched_results}
#' @slot graph the annotation graph, where nodes are annotations, and edges are
#' weighted by similarity between annotation_features.
#' 
#' @export
combined_enrichment <- setClass("combined_enrichment",
                                slots = list(enriched = "list",
                                             annotation = "annotation",
                                             graph = "graphNEL"))