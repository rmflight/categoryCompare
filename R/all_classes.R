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
                                         statistics = "list"))

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