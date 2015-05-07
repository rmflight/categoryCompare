#' unique combinations
#' 
#' determine the unique combinations of annotations that exist in the 
#' significant matrix of the \S4Class{cc_graph} and assign each node in the graph
#' to a group.
#' 
#' @param in_graph the \S4Class{cc_graph} to work on
#' 
#' @return node_assignment
#' @exportMethod annotation_combinations
setMethod("annotation_combinations", 
          signature = list(in_graph = "cc_graph"),
          function(in_graph) .annotation_combinations(in_graph))

.annotation_combinations <- function(in_graph){
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

#' generate colors
#' 
#' given a bunch of items, generate a set of colors for either single node colorings
#' or pie-chart annotations. Colors are generated using the \emph{hcl} colorspace,
#' and for \code{n_color >= 5}, the colors are re-ordered in an attempt to create
#' the largest contrasts between colors, as they result from being picked on a
#' circle in \emph{hcl} space.
#' 
#' @param n_color
#' 
#' @export
#' @importFrom colorspace hcl_rainbow
generate_colors <- function(n_color){
  out_color <- rainbow_hcl(n_color, c = 100)
  
  if (n_color <= 4){
    return(out_color)
  }
  
  out_index <- seq(1, n_color)
  tmp_index <- out_index
  
  if ((n_color %% 2) == 1){
    swap_index_1 <- seq(2, n_color - 1, 2)
    swap_index_2 <- swap_index_1 + 1
    
    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }
    
    swap_index_1 <- seq(3, n_color - 1, 2)
    swap_index_2 <- swap_index_1 + 1
    
    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }
    
  } else {
    swap_index_1 <- seq(1, n_color, 2)
    swap_index_2 <- swap_index_1 + 1
    
    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }
    
    swap_index_1 <- seq(2, n_color - 2, 3)
    swap_index_2 <- seq(n_color - 1, 3, -3)
    
    for (i_swap in seq(1, length(swap_index_1))){
      out_index[swap_index_1[i_swap]] <- tmp_index[swap_index_2[i_swap]]
      out_index[swap_index_2[i_swap]] <- tmp_index[swap_index_1[i_swap]]
      tmp_index <- out_index
    }
    
  }
  
  out_color <- out_color[out_index]
  return(out_color)
}
