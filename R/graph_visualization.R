#' unique combinations
#' 
#' determine the unique combinations of annotations that exist in the 
#' significant matrix of the \S4Class{cc_graph} and assign each node in the graph
#' to a group.
#' 
#' @param in_graph the \S4Class{cc_graph} to work on
#' 
#' @export
#' @return node_assignment
annotation_combinations <- function(in_graph){
  sig_matrix <- in_graph@significant
  
  unique_combinations <- unique(sig_matrix)
  name_combinations <- paste("G", seq(1, nrow(unique_combinations)), sep = "")
  
  rownames(unique_combinations) <- name_combinations
  
  combination_assign <- rep("G", nrow(sig_matrix))
  names(combination_assign) <- rownames(sig_matrix)
  
  for (in_comb in name_combinations){
    has_match <- apply(sig_matrix, 1, function(in_sig){
      identical(in_sig, unique_combinations[in_comb, ])
    })
    combination_assign[has_match] <- in_comb
  }
  new("node_assign", groups = unique_combinations, assignments = combination_assign)
}
