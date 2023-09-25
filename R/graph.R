
get_graph <- function(clust_irr) {
    
    get_clone_edges <- function(x, edges, cs) {
        e <- edges[edges$chain == x, ]
        if(nrow(e)==0) {
            return(NULL)
        }
        
        e <- merge(x = e, y = cs[, c(x, "clone_id", "clone_size")], 
                   by.x = "from_cdr3", by.y = x, all.x = TRUE)
        e$from <- e$clone_id
        e$clone_id <- NULL
        e$clone_size <- NULL
        
        e <- merge(x = e, y = cs[, c(x, "clone_id", "clone_size")], 
                   by.x = "to_cdr3", by.y = x, all.x = TRUE)
        e$to <- e$clone_id
        e$clone_id <- NULL
        e$clone_size <- NULL
        
        return(e)
    }
    
    check_clustirr(clust_irr = clust_irr)

    edges <- get_edges(clust_irr = clust_irr)
    
    if(is.null(edges)) {
        warning("No local or global edges to build igraph from \n")
        return(NULL)
    }
    
    # cells
    s <- slot(clust_irr, "inputs")$s
    
    # clones
    cs <- s
    cs$clone_size <- 1
    cs$id <- NULL
    cs <- aggregate(clone_size~., data = cs, FUN = sum)
    cs$clone_id <- 1:nrow(cs)
    
    # get chains to be analyzed
    clone_edges <- lapply(X = get_chains(colnames(s)), 
                          FUN = get_clone_edges,
                          edges = edges, cs = cs)
    clone_edges <- do.call(rbind, clone_edges)
    cs$name <- cs$clone_id
    cs <- cs[, rev(colnames(cs))]
    
    # remove self-edges (apply only to clones)   
    i <- which(clone_edges$from==clone_edges$to)
    if(length(i) != 0) {
        clone_edges <- clone_edges[-i,]
    }
    
    # build graph
    ig <- igraph::graph_from_data_frame(
        clone_edges[, c("from", "to", "chain", "type", 
                        "motif", "from_cdr3", "to_cdr3")],
        directed = FALSE,
        vertices = cs)

    return(ig)
}



join_graphs <- function(clust_irr_1, 
                        clust_irr_2) {
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
    
    global_max_dist <- slot(clust_irr_1, "inputs")$control$global_max_dist
    chains <- colnames(slot(clust_irr_1, "inputs")$s)
    ige <- get_intergraph_edges(s1=d1$vertices,
                                s2=d2$vertices,
                                global_max_dist = global_max_dist,
                                chains = chains)
    # ige$from <- paste0("s1|", ige$from)
    # ige$to <- paste0("s2|", ige$to)
    
    d1$vertices <- rbind(d1$vertices, d2$vertices)
    d1$edges <- rbind(d1$edges, d2$edges, ige)
    d <- d1
    
    # build graph
    return(graph_from_data_frame(
        d$edges, 
        directed = FALSE,
        vertices = d$vertices))
}


plot_graph <- function(clust_irr) {
    
    check_clustirr(clust_irr = clust_irr)
    
    ig <- get_graph(clust_irr = clust_irr)
    if (is.null(ig)) {
        warning("No graph to plot \n")
        return(NULL)
    }
    plot(ig)
}

plot_graph_old <- function(clust_irr) {
    
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
    nodes <- nodes[!duplicated(nodes$label), ]
    edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, ]
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