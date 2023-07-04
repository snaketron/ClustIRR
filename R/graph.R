get_graph <- function(clust_irr) {
    
    edges <- get_edges(clust_irr)
    
    igraph_object <- igraph::graph_from_data_frame(edges, 
                                                   directed = FALSE)
    return(igraph_object)
}
