#' hypergeom feature class
#' 
#' class to hold features undergoing hypergeometric enrichment
#' 
#' @slot significant the significant features
#' @slot universe all of the features measured
#' @slot annotation annotation object
#' 
#' @export
setClass("hypergeom_features",
         slots = list(significant = "ANY",
                      universe = "ANY",
                      annotation = "annotation"))

#' do hypergeometric enrichment
#' 
#' @param hypergeom_features a hypergeometric_features object
#' @param direction which direction to do the enrichment (over or under)
#' @export
#' @return enriched_result
#' 
hypergeometric_feature_enrichment <- function(hypergeom_features, direction = "over"){
  
  # cleanup the features and annotations (should be in separate function)
  hypergeom_features@universe <- unique(hypergeom_features@universe)
  
  tmp_annot_feature <- hypergeom_features@annotation@annotation_features
  annotation_universe <- unique(unlist(tmp_annot_feature))
  
  hypergeom_features@universe <- intersect(hypergeom_features@universe, annotation_universe)
  tmp_annot_feature <- lapply(tmp_annot_feature, intersect, hypergeom_features@universe)
  
  n_feature <- sapply(tmp_annot_feature, length)
  keep_annot <- n_feature > 0
  
  tmp_annot_feature <- tmp_annot_feature[keep_annot]
  
  hypergeom_features@significant <- intersect(hypergeom_features@significant, hypergeom_features@universe)
  hypergeom_features@annotation@annotation_features <- tmp_annot_feature
  
  
  # this probably needs its own function eventually
  if (length(hypergeom_features@annotation@description) != 0){
    hypergeom_features@annotation@description <- hypergeom_features@annotation@description[names(tmp_annot_feature)]
  }
  
  if (length(hypergeom_features@annotation@links) != 0){
    hypergeom_features@annotation@links <- hypergeom_features@annotation@links[names(tmp_annot_feature)]
  }
  
  # now get the counts annotated
  num_white_drawn <- sapply(hypergeom_features@annotation@annotation_features, function(x) sum(hypergeom_features@significant %in% x))
  
  if (length(num_white_drawn) == 0){
    num_white_drawn <- 0
  }
  
  num_white <- Biobase:::listLen(hypergeom_features@annotation@annotation_features)
  
  if (length(num_white) == 0){
    num_white <- 0
  }
  
  num_black <- length(hypergeom_features@universe) - num_white
  num_drawn <- length(hypergeom_features@significant)
  
  hyper_stats <- hypergeometric_basic(num_white, num_black, num_drawn, num_white_drawn, direction)
  
  hyper_stats$counts <- num_white_drawn[names(hyper_stats$p)]
  
  out_stats <- new("statistical_results",
                   statistic_data = hyper_stats,
                   annotation_id = names(hypergeom_features@annotation@annotation_features))
  
  out_enrich <- new("enriched_result",
                    features = hypergeom_features@significant,
                    universe = hypergeom_features@universe,
                    statistics = out_stats,
                    annotation = hypergeom_features@annotation)
  
  out_enrich
  
}

#' generate link text
#' 
#' given a named vector of links, generate an actual html link formatted for output in html documents
#' 
#' @param links the vector of links
#' 
#' @export
#' @return character
generate_link <- function(links){
  link_names <- names(links)
  
  if (max(nchar(link_names)) == 0){
    link_names <- "link"
  }
  
  paste('<a href="', links, '">', link_names, '</a>', sep = "")
  
}

#' do hypergeometric test
#' 
#' does a hypergeometric enrichment test
#' 
#' @param num_white number of white balls in urn
#' @param num_black number of black balls in urn
#' @param num_drawn number of balls taken from urn
#' @param num_white_drawn number of white balls taken from urn
#' @param direction which direction is the test
#' 
#' @export
#' @return list
hypergeometric_basic <- function(num_white, num_black, num_drawn, num_white_drawn, direction = "over"){
  n_2_1 <- num_white - num_white_drawn
  n_1_2 <- num_drawn - num_white_drawn
  n_2_2 <- num_black - n_1_2
  
  odds_ratio <- (num_white_drawn * n_2_2) / (n_1_2 * n_2_1)
  
  expected <- ((num_white_drawn + n_1_2) * (num_white_drawn + n_2_1)) / (num_white_drawn + n_1_2 + n_2_1 + n_2_2)
  
  p_values <- switch(direction,
                     over =  phyper(num_white_drawn - 1L, num_white, num_black, num_drawn, lower.tail = FALSE),
                     under = phyper(num_white_drawn, num_white, num_black, num_drawn, lower.tail = TRUE)
  )
  
  list(p = p_values, odds = odds_ratio, expected = expected)
}
