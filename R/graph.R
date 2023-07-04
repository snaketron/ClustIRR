get_graph <- function(clustirr) {
    
    edges <- get_edges(clustirr)
    
    igraph_object <- igraph::graph_from_data_frame(edges, 
                                                   directed = FALSE)
    return(igraph_object)
}