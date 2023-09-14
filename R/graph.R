
get_graph <- function(clust_irr) {
    check_clustirr(clust_irr = clust_irr)

    edges <- get_edges(clust_irr = clust_irr)

    if (is.null(edges)) {
        warning("No local or global edges to build igraph from \n")
        return(NULL)
    }

    ig <- graph_from_data_frame(edges, directed = FALSE)
    
    # these are the nodes that are part of s but do not have connections
    # in the graph, so we have to add them as individual vertices
    n <- as.numeric(V(ig)$name)
    n <- setdiff(x = slot(clust_irr, "inputs")$s$id, y = n)
    ig <- igraph::add_vertices(graph = ig,
                               nv = length(n),
                               name = as.character(n))
    
    return(ig)
}


plot_graph <- function(clust_irr, expand_clones = FALSE) {
    check_clustirr(clust_irr = clust_irr)
    ig <- get_graph(clust_irr = clust_irr)
    if (is.null(ig)) {
        warning("No graph to plot \n")
        return(NULL)
    }
    edges <- as_data_frame(ig, what = "edges")
    nodes <- as_data_frame(ig, what = "vertices")
    chains <- names(slot(clust_irr, "clust"))
    types <- unique(edges$type)
    edges <- configure_edges(edges = edges, chains = chains, types = types)
    nodes <- configure_nodes(
        nodes = nodes, edges = edges, chains = chains,
        types = types, s = slot(clust_irr, "inputs")$s
    )
    if (!expand_clones) {
        nodes <- nodes[!duplicated(nodes$label), ]
        edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, ]
    }
    if(nrow(edges) != 0){
        ledges <- data.frame(
            color = unique(edges$color),
            label = unique(edges$type),
            arrows = "", width = 4
        )
    } else {
        ledges <- data.frame()
    }
    if (length(chains) > 1) {
        chains <- append(chains, "Both chains")
    }
    lnodes <- data.frame(
        label = chains,
        color = "", shape = "dot", size = 8
    )
    lnodes$color <- ifelse(test = lnodes$label == "CDR3b" |
        lnodes$label == "CDR3g" |
        lnodes$label == "CDR3h",
    yes = "blue",
    no = "yellow"
    )
    lnodes$color <- ifelse(test = lnodes$label == "Both chains",
        yes = "green",
        no = lnodes$color
    )
    shapes <- c("dot", "diamond")
    labels <- c("Expanded", "Singleton")
    for (i in seq_along(shapes)) {
        if (shapes[i] %in% unique(nodes$shape)) {
            lnodes <- rbind(
                lnodes,
                data.frame(
                    label = labels[i],
                    color = "black",
                    shape = shapes[i],
                    size = 8
                )
            )
        }
    }
    return(configure_network(
        nodes = nodes, edges = edges,
        lnodes = lnodes, ledges = ledges
    ))
}


join_graphs <- function(clust_irr_1, clust_irr_2) {
    if(slot(clust_irr_1, "inputs")$control$global_max_dist!=
       slot(clust_irr_2, "inputs")$control$global_max_dist) {
       stop("variable global_max_dist used in clust_irr_1 and clust_irr_2") 
    }
    
    ig1 <- get_graph(clust_irr = clust_irr_1)
    ig2 <- get_graph(clust_irr = clust_irr_2)
    
    d1 <- igraph::get.data.frame(ig1, what = "both")
    d1$edges$sample <- "s1"
    d1$edges <- d1$edges[, c("from", "to", "chain", "sample")]
    d1$edges$from <- paste0("s1|", d1$edges$from)
    d1$edges$to <- paste0("s1|", d1$edges$to)
    
    d1$vertices$id <- d1$vertices$name
    d1$vertices$name <- paste0("s1|", d1$vertices$name)
    
    d2 <- igraph::get.data.frame(ig2, what = "both")
    d2$edges$sample <- "s2"
    d2$edges <- d2$edges[, c("from", "to", "chain", "sample")]
    d2$edges$from <- paste0("s2|", d2$edges$from)
    d2$edges$to <- paste0("s2|", d2$edges$to)
    
    d2$vertices$id <- d2$vertices$name
    d2$vertices$name <- paste0("s2|", d2$vertices$name)
    
    d1$vertices <- rbind(d1$vertices, d2$vertices)
    d1$edges <- rbind(d1$edges, d2$edges)
    
    global_max_dist <- slot(clust_irr_1, "inputs")$control$global_max_dist
    ige <- get_intergraph_edges(s1=slot(clust_irr_1, "inputs")$s,
                                s2=slot(clust_irr_2, "inputs")$s,
                                global_max_dist = global_max_dist)
    ige$from <- paste0("s1|", ige$from)
    ige$to <- paste0("s2|", ige$to)
    d1$edges <- rbind(d1$edges, ige)
    d <- d1
    
    # build graph
    return(graph_from_data_frame(
        d$edges, 
        directed = FALSE,
        vertices = d$vertices))
}
