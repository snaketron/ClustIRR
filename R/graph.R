get_graph <- function(clust_irr) {
    
    edges <- get_edges(clust_irr)
    
    # if we don't compute this, igraph will break some from/to connections
    # and effectively make the graph useless by assigning wrong indices
    v <- unique(c(unique(edges$from), unique(edges$to)))
    
    ig <- igraph::graph_from_data_frame(edges, 
                                        directed = FALSE,
                                        vertices = v)
    return(ig)
}


plot_graph <- function(ig) {
  
  # prepare nodes & edges for visNetwork
  nodes <- igraph::as_data_frame(ig, what = "vertices")
  edges <- igraph::as_data_frame(ig, what = "edges")
  
  names(nodes) <- "id"
  nodes$color.border <- "black"
  nodes$color.highlight <- "red"
  nodes$size <- 20
  

  
  return(NULL)
}