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
#' @slot statistic_data list of numerical statistics
#' @slot annotation_id vector of ids
#' 
#' @export
statistical_results <- setClass("statistical_results",
                                slots = list(statistic_data = "list",
                                             annotation_id = "ANY"))

#' the enriched results class
#' 
#' @slot features the "features" of interest, a vector of class ANY
#' @slot universe all of the "features" in the background
#' @slot annotation list giving the annotation to feature relationship
#' @slot statistics a \linkS4class{statistical_results} object
#' 
#' @export
enriched_result <- setClass("enriched_result",
                            slots = list(features = "ANY",
                                         universe = "ANY",
                                         annotation = "annotation",
                                         statistics = "statistical_results"))

#' significant annotations
#' 
#' The \code{significant_annotations} class holds which annotations from which
#' enrichment were both \strong{measured} and \strong{significant}. Each of these
#' slots is a \emph{logical matrix} with rows named by \emph{annotation_id} and 
#' columns named by the names of the \linkS4class{enriched_result} that was combined.
#' 
#' @slot significant logical matrix
#' @slot measured logical matrix
#' @slot sig_calls character representations of calls used to filter the data
#' 
#' @export
significant_annotations <- setClass("significant_annotations",
                                    slots = list(significant = "matrix",
                                                 measured = "matrix",
                                                 sig_calls = "character"))

#' combined statistics
#' 
#' holds the results of extracting a bunch of statistics from a \linkS4class{combined_enrichment}
#' into one entity. This is useful because we want to enable multiple data representations and
#' simple filtering on the actual \code{data.frame} of statistics, and this provides flexibility
#' to enable that.
#' 
#' @slot statistic_data a \code{data.frame} of all of the statistics from all of the enrichments
#' @slot significant a \linkS4class{significant_annotations} object, that may be empty
#' @slot which_enrichment a \code{vector} giving which enrichment each column of the statistics came from
#' @slot which_statistic a \code{vector} providing which statistic each column contains
#' 
#' @export
combined_statistics <- setClass("combined_statistics",
                                contains = "statistical_results",
                                slots = list(statistic_data = "data.frame",
                                             annotation_id = "character",
                                             significant = "significant_annotations",
                                             which_enrichment = "character",
                                             which_statistic = "character"))

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
                                             graph = "graphNEL",
                                             statistics = "combined_statistics"))

#' cc_graph
#' 
#' A \code{cc_graph} class is a \code{graphNEL} with the added slot of
#' \code{significant}, a matrix of rows (nodes / annotations) and whether
#' they were found to be significant in a given enrichment (columns). This
#' matrix is used for classifying the annotations into different groups, and
#' generating either pie-charts or coloring the nodes in a visualization.
#' 
#' @slot significant numeric matrix of ones and zeros
#' 
#' @export
#' @importFrom graph graphNEL
cc_graph <- setClass("cc_graph",
                     contains = "graphNEL",
                     slots = list(significant = "matrix"))

#' node_assign
#' 
#' The \code{node_assign} class holds the unique annotation combinations and the
#' assignment of the nodes to those combinations for use in visualization.
#' 
#' @slot groups the unique groups, as a logical matrix
#' @slot assignments named character vector providing association with groups
#' @slot colors named character vector of hex colors for groups or experiments
#' @slot color_type whether doing group or experiment based colors
#' @slot pie_locs if doing experiment colors, then pie graphs were generated here
#' 
#' @export
node_assign <- setClass("node_assign",
                        slots = list(groups = "matrix",
                                     assignments = "character",
                                     colors = "character",
                                     color_type = "character",
                                     pie_locs = "character"))
