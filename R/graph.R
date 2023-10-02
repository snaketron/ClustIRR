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
    
    get_local_edges <- function(clust_irr) {
        
        get_lp <- function(x, lp, chain) {
            if (sum(lp$motif == x) > 1) {
                p <- t(combn(x = lp$cdr3[lp$motif == x], m = 2))
                return(data.frame(
                    from_cdr3 = p[, 1], to_cdr3 = p[, 2],
                    motif = x, type = "local",
                    chain = chain
                ))
            }
            return(NULL)
        }
        
        get_from_to <- function(x, lp) {
            y <- lp$cdr3[lp$motif == x]
            if(length(y)>1) {
                y <- t(utils::combn(x = y, m = 2))
                y <- data.frame(y)
                colnames(y) <- c("from_cdr3", "to_cdr3")
                y$motif <- x
                return(y)
            }
        }
        
        el<-vector(mode="list", length = length(get_clustirr_clust(clust_irr)))
        names(el) <- names(get_clustirr_clust(clust_irr))
        for(chain in names(get_clustirr_clust(clust_irr))) {
            lp <- get_clustirr_clust(clust_irr)[[chain]]$local$lp
            if(is.null(lp) == FALSE && nrow(lp) != 0) {
                lp <- lapply(X = unique(lp$motif), FUN = get_from_to, lp = lp)
                lp <- do.call(rbind, lp)
                if(is.null(lp)==FALSE) {
                    lp$chain <- chain
                    lp$type <- "local"
                    el[[chain]] <- lp
                }
            }
        }
        el <- do.call(rbind, el)
        return(unique(el))
    }
    
    get_global_edges <- function(clust_irr) {
        eg<-vector(mode="list", length = length(get_clustirr_clust(clust_irr)))
        names(eg) <- names(get_clustirr_clust(clust_irr))
        for(chain in names(get_clustirr_clust(clust_irr))) {
            g <- get_clustirr_clust(clust_irr)[[chain]]$global
            if(is.null(g) == FALSE && nrow(g) != 0) {
                eg[[chain]] <- data.frame(from_cdr3 = g[,1], 
                                          to_cdr3 = g[,2], 
                                          motif = NA,
                                          chain = chain,
                                          type = "global")
            }
        }
        eg <- do.call(rbind, eg)
        return(eg)
    }
    
    get_edge_order <- function(x) {
        return(paste0(sort(c(x[1], x[2])), collapse = '-'))
    }
    
    check_clustirr(clust_irr = clust_irr)
    
    # cells
    s <- get_clustirr_inputs(clust_irr)$s
    
    # clones
    cs <- s
    cs$clone_size <- 1
    cs$id <- NULL
    cs <- aggregate(clone_size~., data = cs, FUN = sum)
    cs$clone_id <- seq_len(nrow(cs))
    
    # get edges
    edges <- rbind(get_local_edges(clust_irr = clust_irr),
                   get_global_edges(clust_irr = clust_irr))
    
    # build graph with only vertices
    if(is.null(edges)) {
        cs$name <- cs$clone_id
        cs <- cs[, rev(colnames(cs))]
        
        ig <- graph_from_data_frame(d = data.frame(from = 1, to = 1),
                                    directed = FALSE,
                                    vertices = cs)
        ig <- delete_edges(ig, edges = 1)
        return(ig)
    }
    
    # get chains to be analyzed
    clone_edges <- lapply(X = get_chains(colnames(s)), 
                          FUN = get_clone_edges,
                          edges = edges, 
                          cs = cs)
    clone_edges <- do.call(rbind, clone_edges)
    clone_edges$edge_id <- apply(X = clone_edges[, c("from", "to")],
                                 MARGIN = 1,
                                 FUN = get_edge_order)
    
    cs$name <- cs$clone_id
    cs <- cs[, rev(colnames(cs))]
    
    # remove self-edges (apply only to clones)   
    i <- which(duplicated(clone_edges$edge_id))
    if(length(i) != 0) {
        clone_edges <- clone_edges[-i,]
    }
    
    # build graph
    ig <- graph_from_data_frame(clone_edges[, c("from","to","chain","type", 
                                                "motif","from_cdr3","to_cdr3")],
                                directed = FALSE,
                                vertices = cs)
    
    return(ig)
}

join_graphs <- function(clust_irr_1, clust_irr_2) {
    
    if(get_clustirr_inputs(clust_irr_1)$control$global_max_dist!=
       get_clustirr_inputs(clust_irr_2)$control$global_max_dist) {
       stop("variable global_max_dist used in clust_irr_1 and clust_irr_2") 
    }
    
    ig1 <- get_graph(clust_irr = clust_irr_1)
    ig2 <- get_graph(clust_irr = clust_irr_2)
    
    # get the vertices/edges of graph 1
    d1 <- get.data.frame(ig1, what = "both")
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
    d2 <- get.data.frame(ig2, what = "both")
    if(nrow(d2$edges)!=0) {
        d2$edges$sample <- "s2"
        d2$edges <- d2$edges[, c("from", "to", "chain", "sample", "type")]
        d2$edges$from <- paste0("s2|", d2$edges$from)
        d2$edges$to <- paste0("s2|", d2$edges$to)
    }
    d2$vertices$id <- d2$vertices$name
    d2$vertices$name <- paste0("s2|", d2$vertices$name)
    d2$vertices$sample <- "s2"
    
    global_max_dist <- get_clustirr_inputs(clust_irr_1)$control$global_max_dist
    chains <- colnames(get_clustirr_inputs(clust_irr_1)$s)
    chains <- chains[chains!="id"]
    ige <- get_intergraph_edges(s1=d1$vertices,
                                s2=d2$vertices,
                                global_max_dist = global_max_dist,
                                chains = chains)
    
    d1$vertices <- rbind(d1$vertices, d2$vertices)
    d1$edges <- rbind(d1$edges, d2$edges, ige)
    d <- d1
    
    d1$edges <- config_edges(es = d1$edges)
    
    # build joint graph
    g <- graph_from_data_frame(d1$edges, directed=FALSE, vertices=d$vertices)
    
    # make graph look visually better
    g <- config_vertices_plot(g = g, is_jg = TRUE)
    g <- config_edges_plot(g = g, is_jg = TRUE)
    
    return(g)
}

plot_graph <- function(clust_irr, as_visnet = FALSE) {
    
    check_clustirr(clust_irr = clust_irr)
    
    ig <- get_graph(clust_irr = clust_irr)
    if(is.null(ig)) {
        warning("No graph to plot \n")
        return(NULL)
    }
    
    # get the vertices/edges of the graph
    d <- get.data.frame(ig, what = "both")
    d$vertices$id <- d$vertices$name
    if(nrow(d$edges)!=0) {
        d$edges <- d$edges[, c("from", "to", "chain", "type")]
        d$edges <- config_edges(es = d$edges)
        ig <- graph_from_data_frame(d$edges, directed=FALSE,vertices=d$vertices)
    } 
    else {
        d$edges <- data.frame(from = 1, to = 1)
        ig <- graph_from_data_frame(d$edges, directed=FALSE,vertices=d$vertices)
        ig <- delete_edges(ig, edges = 1)
    }
    
    # make graph look visually better
    ig <- config_vertices_plot(g = ig, is_jg = FALSE)
    ig <- config_edges_plot(g = ig, is_jg = FALSE)
    
    if(as_visnet == FALSE) {
        return(plot(ig, vertex.label = NA))
    }
    
    if(as_visnet == TRUE) {
        V(ig)$size <- V(ig)$size*10
        return(visIgraph(igraph = ig,
                         idToLabel = TRUE,
                         layout = "layout_components",
                         randomSeed = 1234,
                         physics = FALSE,
                         smooth = FALSE,
                         type = "square"))
    }
}
