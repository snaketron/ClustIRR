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

join_graphs <- function(clust_irr_1, clust_irr_2) {
    
    if(slot(clust_irr_1, "inputs")$control$global_max_dist!=
       slot(clust_irr_2, "inputs")$control$global_max_dist) {
       stop("variable global_max_dist used in clust_irr_1 and clust_irr_2") 
    }
    
    ig1 <- get_graph(clust_irr = clust_irr_1)
    ig2 <- get_graph(clust_irr = clust_irr_2)
    
    # get the vertices/edges of graph 1
    d1 <- igraph::get.data.frame(ig1, what = "both")
    if(nrow(d1$edges)!=0) {
        d1$edges$sample <- "s1"
        d1$edges <- d1$edges[, c("from", "to", "chain", "sample", "type")]
        d1$edges$from <- paste0("s1|", d1$edges$from)
        d1$edges$to <- paste0("s1|", d1$edges$to)
    }
    d1$vertices$id <- d1$vertices$name
    d1$vertices$name <- paste0("s1|", d1$vertices$name)
    d1$vertices$sample <- "s1"
    
    # get the vertices/edges of graph 2
    d2 <- igraph::get.data.frame(ig2, what = "both")
    if(nrow(d2$edges)!=0) {
        d2$edges$sample <- "s2"
        d2$edges <- d2$edges[, c("from", "to", "chain", "sample", "type")]
        d2$edges$from <- paste0("s2|", d2$edges$from)
        d2$edges$to <- paste0("s2|", d2$edges$to)
    }
    d2$vertices$id <- d2$vertices$name
    d2$vertices$name <- paste0("s2|", d2$vertices$name)
    d2$vertices$sample <- "s2"
    
    global_max_dist <- slot(clust_irr_1, "inputs")$control$global_max_dist
    chains <- colnames(slot(clust_irr_1, "inputs")$s)
    chains <- chains[chains!="id"]
    ige <- get_intergraph_edges(s1=d1$vertices,
                                s2=d2$vertices,
                                global_max_dist = global_max_dist,
                                chains = chains)
    
    d1$vertices <- rbind(d1$vertices, d2$vertices)
    d1$edges <- rbind(d1$edges, d2$edges, ige)
    d <- d1
    
    d1$edges <- configure_edges(es = d1$edges)
    
    # build joint graph
    g <- graph_from_data_frame(d1$edges, directed=FALSE, vertices=d$vertices)
    
    # make graph look visually better
    g <- configure_vertices_plot(g = g, is_jg = TRUE)
    g <- configure_edges_plot(g = g, is_jg = TRUE)
    
    return(g)
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

configure_edges <- function(es) {
    
    # edge color by type
    get_type_edge <- function(x) {
        x <- sort(unique(x))
        if(length(x)==2) {
            return("black")
        }
        if(x == "local") {
            return("red")
        }
        if(x == "global") {
            return("purple")
        }
        return("black")
    }
    
    # edge shape by chain
    get_chain_edge <- function(x) {
        x <- sort(unique(x))
        if(length(x)==2) {
            return("solid")
        }
        if(x %in% c("CDR3b", "CDR3g", "CDR3h")) {
            return("dashed")
        }
        if(x %in% c("CDR3a", "CDR3d", "CDR3l")) {
            return("dotted")
        }
    }
    
    if(nrow(es) == 0) {
        return(NULL)
    }
    
    attr_type <- aggregate(type~from+to, data = es, FUN = get_type_edge)
    attr_type$type_color <- attr_type$type
    attr_type$type <- NULL
    
    attr_chain <- aggregate(chain~from+to, data = es, FUN = get_chain_edge)
    attr_chain$chain_shape <- attr_chain$chain
    attr_chain$chain <- NULL
    
    return(merge(x = attr_chain, y = attr_type, by = c("from", "to")))
}

configure_edges_plot <- function(g, is_jg) {
    n_e <- length(g)
    if(n_e != 0) {
        E(g)$color <- E(g)$type_color
        E(g)$lty <- E(g)$chain_shape
    }
    return(g)
}

configure_vertices_plot <- function(g, is_jg) {
    # default features
    V(g)$size <- 1.5+log10(V(g)$clone_size)
    
    if(is_jg==TRUE) {
        V(g)$color <- ifelse(test = V(g)$sample == "s1", 
                             yes = "steelblue", no = "darkorange")
        V(g)$frame.color <- V(g)$color
    } 
    else {
        V(g)$color <- "white"
        V(g)$frame.color <- "black"
    }
    
    return(g)
}



get_visnet <- function(nodes, edges, ledges, lnodes) {
    V(g)$size <- V(g)$size*10
    visIgraph(igraph = g,
              idToLabel = TRUE,
              layout = "layout_components",
              randomSeed = 1234,
              physics = FALSE,
              smooth = FALSE,
              type = "square")
}








