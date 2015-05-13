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
#' @importFrom colorspace rainbow_hcl
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

#' assign colors
#' 
#' given a \linkS4Class{node_assign}, assign colors to either the independent groups
#' of unique annotations, or to each of the experiments independently.
#' 
#' @param in_assign the \linkS4Class{node_assign} object generated from a \linkS4Class{cc_graph}
#' @param type either "group" or "experiment"
#' 
#' @export
#' @return node_assign with colors
assign_colors <- function(in_assign, type = "experiment"){
  grp_matrix <- in_assign@groups
  
  if (type == "experiment"){
    n_color <- ncol(grp_matrix)
    use_color <- generate_colors(n_color)
    names(use_color) <- colnames(grp_matrix)
    
    in_assign@colors <- use_color
    in_assign@color_type <- "pie"
    in_assign@pie_locs <- generate_piecharts(grp_matrix, use_color)
  } else {
    n_color <- nrow(grp_matrix)
    use_color <- generate_colors(n_color)
    names(use_color) <- rownames(grp_matrix)
    in_assign@colors <- use_color
    in_assign@color_type <- "solid"
  }
  
  return(in_assign)
}

#' create piecharts for visualization
#' 
#' given a group matrix and the colors for each experiment, generate the pie graphs
#' that will be used as glyphs in Cytoscape
#' 
#' this should \emph{not be exported in the final version}
#' 
#' @param grp_matrix the group matrix
#' @param use_color the colors for each experiment
#' 
#' @export 
#' @return list of png files that are pie graphs
#' @importFrom colorspace desaturate
#' @import Cairo
generate_piecharts <- function(grp_matrix, use_color){
  n_grp <- nrow(grp_matrix)
  n_color <- length(use_color)
  
  # defines how many pie segments are needed, common to all the pie-charts
  pie_area <- rep(1 / n_color, n_color)
  names(pie_area) <- rep("", n_color) # add blank names so nothing gets printed
  
  # use desaturated version of colors when there is non-significance
  desat_color <- desaturate(use_color)
  names(desat_color) <- names(use_color)
  piecharts <- sapply(rownames(grp_matrix), function(i_grp){
    tmp_logical <- grp_matrix[i_grp, ]
    tmp_color <- use_color
    
    # add the proper desaturated versions of the colors
    tmp_color[!tmp_logical] <- desat_color[!tmp_logical]
    
    # use a tempfile so that multiple runs should generate their own files
    out_file <- tempfile(i_grp, fileext = ".png")
    Cairo(width = 640, height = 640, file = out_file, type = "png", bg = "transparent")
    par(mai = c(0, 0, 0, 0))
    pie(pie_area, col = tmp_color, clockwise = TRUE)
    dev.off()
    out_file
  })
  return(piecharts)
}

#' visualize in cytoscape
#' 
#' given a graph, and the node assignments, visualize the graph in cytoscape
#' for manipulation
#' 
#' @param in_graph the cc_graph to visualize
#' @param in_assign the node_assign generated
#' @param description something descriptive about the vis (useful when lots of different visualizations)
#' @param ... other parameters for \code{CytoscapeWindow}
#' 
#' @import RCytoscape
#' @export
#' @return something
vis_in_cytoscape <- function(in_graph, in_assign, description = "", ...){
  
  # initialize and add the visual attribute so we can color according to the
  # data that lives in in_assign
  nodeDataDefaults(in_graph, "visattr") <- ""
  attr(nodeDataDefaults(in_graph, "visattr"), "class") <- 'STRING'
  nodeData(in_graph, names(in_assign@assignments), "visattr") <- in_assign@assignments
  
  cyt_window <- CytoscapeWindow(description, graph = in_graph, ...)
  displayGraph(cyt_window)
  setLayoutProperties(cyt_window, 'force-directed', list(edge_attribute='weight'))
  layoutNetwork(cyt_window, 'force-directed')
  
  if (in_assign@color_type == "solid"){
    setNodeColorRule(cyt_window, "visattr", names(in_assign@colors), in_assign@colors, mode = "lookup")
    redraw(cyt_window)
  } else if (in_assign@color_type == "pie"){
    pie_images <- in_assign@pie_locs[in_assign@assignments]
    names(pie_images) <- NULL
    pie_images <- paste("file://localhost", pie_images, sep = "")
    setNodeImageDirect(cyt_window, names(in_assign@assignments), pie_images)
    setDefaultNodeColor(cyt_window, 'transparent')
    setNodeOpacityDirect(cyt_window, names(in_assign@assignments), 0)
    setDefaultNodeShape(cyt_window, "diamond")
    redraw(cyt_window)
  }
  return(cyt_window)
}

#' remove edges
#' 
#' given a \linkS4Class{CytoscapeWindowClass}, remove edges according to provided
#' values.
#' 
#' @param edge_obj a CytoscapeWindowClass
#' @param cutoff what cutoff to use to remove edges
#' @param edge_attr what attribute has the values
#' @param value_direction remove those edges "under" or "over" the value
#' 
#' @export
#' @return nothing
setMethod("remove_edges", signature=list(edge_obj="CytoscapeWindowClass", cutoff="numeric"), function(edge_obj, cutoff, edge_attr, value_direction) 
  .remove_edges_cw(edge_obj, cutoff, edge_attr, value_direction))

.remove_edges_cw <-	function(cyt_window, cutoff, edge_attr = "weight", value_direction = "under"){
  edge_data <- getAllEdgeAttributes(cyt_window)
  
  switch(value_direction,
         under = edge_data <- edge_data[(as.numeric(edge_data[, edge_attr]) < cutoff),],
         over = edge_data <- edge_data[(as.numeric(edge_data[, edge_attr]) > cutoff),]
  )
  
  attr_names <- names(edge_data)
  if (!('edgeType' %in% attr_names)){
    edge_names <- paste(edge_data$source,' (unspecified) ',edge_data$target,sep='')
  } else {
    edge_names <- paste(edge_data$source,' (',edge_data$edgeType,') ',edge_data$target, sep='')
  }
  selectEdges(cyt_window,edge_names)
  deleteSelectedEdges(cyt_window)
  
  layoutNetwork(cyt_window, 'force-directed')
  
  message("Removed ", length(edge_names), " edges from graph\n")
}

#' remove graph edges
#' 
#' @param edge_obj cc_graph
#' @param cutoff the cutoff to use
#' @param edge_attr which attribute to use
#' @param value_direction remove edges with value under or over
#' 
#' @export
#' @return cc_graph
setMethod("remove_edges", signature=list(edge_obj="cc_graph", cutoff="numeric"), function(edge_obj, cutoff, edge_attr, value_direction) 
  .remove_edges_ccgraph(edge_obj, cutoff, edge_attr, value_direction))

.remove_edges_ccgraph <-	function(in_graph, cutoff, edge_attr = "weight", value_direction = "under"){
  edge_data <- unlist(edgeData(in_graph, , , edge_attr))
  
  switch(value_direction,
         under = del_edges <- names(edge_data)[edge_data < cutoff],
         over = del_edges <- names(edge_data)[edge_data > cutoff]
  )
  
  if (length(del_edges) > 0){
    del_edges <- strsplit(del_edges, "|", fixed = TRUE)
    from_node <- sapply(del_edges, function(x){x[1]})
    to_node <- sapply(del_edges, function(x){x[2]})
    in_graph <- removeEdge(from_node, to_node, in_graph)
  }
  
  message("Removed ", length(del_edges), " edges from graph\n")
  return(in_graph)
}
