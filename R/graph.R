

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
        return(list(graph = ig, clones = cs))
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
    
    
    return(list(graph = ig, clones = cs))
}


get_joint_graph <- function(clust_irrs) {
    
    get_vertices <- function(x) {
        return(x$vertices)
    }
    
    get_edges <- function(x) {
        return(x$edges)
    }
    
    check_input <- function(clust_irrs) {
        if(missing(clust_irrs)==TRUE) {
            stop("clust_irrs input missing")
        }
        if(is.list(clust_irrs)==FALSE) {
            stop("clust_irrs must be a list of clust_irr objects")
        }
        if(length(clust_irrs)<=1) {
            stop("get_joint_graph needs >= 2 clust_irr outputs")
        }
        gmd <- numeric(length = length(clust_irrs))
        for(i in 1:length(clust_irrs)) {
            check_clustirr(clust_irr = clust_irrs[[i]])
            gmd[i]<-get_clustirr_inputs(clust_irrs[[i]])$control$global_max_dist
        }
        if(length(unique(gmd))!=1) {
            stop("all global_max_dist should be equal")
        }
        
        # check if same chain names
        cs <- colnames(get_clustirr_inputs(clust_irrs[[1]])$s)
        for(i in 2:length(clust_irrs)) {
            if(any(colnames(get_clustirr_inputs(clust_irrs[[i]])$s)!=cs)) {
                stop("different chains in graphs")
            }
        }
    }
    
    # check input
    check_input(clust_irrs = clust_irrs)
    
    # get chains
    igs <- lapply(X = clust_irrs, FUN = get_graph)
    chains <- colnames(get_clustirr_inputs(clust_irrs[[1]])$s)
    chains <- chains[chains!="id"]
    
    # get global_max_dist
    gmd <- get_clustirr_inputs(clust_irrs[[1]])$control$global_max_dist
    
    ige <- get_intergraph_edges(igs=igs, global_max_dist=gmd, chains=chains)
    
    df_v <- do.call(rbind, lapply(X = ige$igs_df, FUN = get_vertices))
    df_e <- do.call(rbind, lapply(X = ige$igs_df, FUN = get_edges))
    df_e <- rbind(df_e, ige$ige)
    
    df_e <- config_edges(es = df_e)
    
    # build joint graph
    g <- graph_from_data_frame(df_e, directed=FALSE, vertices=df_v)
    
    # make graph look visually better
    g <- config_vertices_plot(g = g, is_jg = TRUE)
    g <- config_edges_plot(g = g, is_jg = TRUE)
    
    # prepare clones
    # cs_1 <- ig1$clones
    # cs_1$sample <- "s1"
    # cs_2 <- ig2$clones
    # cs_2$sample <- "s2"
    
    return(list(graph = g, clones = NA))
}


plot_graph <- function(clust_irr, 
                       as_visnet = FALSE) {
    
    check_clustirr(clust_irr = clust_irr)
    
    ig <- get_graph(clust_irr = clust_irr)
    clones <- ig$clones
    if(is.null(ig$graph)) {
        warning("No graph to plot \n")
        return(list(graph = NA, clones = clones))
    }
    
    # get the vertices/edges of the graph
    d <- get.data.frame(ig$graph, what = "both")
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
    
    # plot
    if(as_visnet == FALSE) {
        plot(ig, vertex.label = NA)
    }
    if(as_visnet == TRUE) {
        V(ig)$size <- V(ig)$size*10
        if(length(E(ig))==0) {
            # apparently if no edges, visnetwork can't plot
            ig <- add_edges(graph = ig, edges = c(1,1))
        }
        visIgraph(igraph = ig,
                  idToLabel = TRUE,
                  layout = "layout_components",
                  randomSeed = 1234,
                  physics = FALSE,
                  smooth = FALSE,
                  type = "square")
        
    }
}


plot_joint_graph <- function(clust_irrs, 
                             as_visnet = FALSE) {
    
    
    check_input <- function(clust_irrs) {
        if(missing(clust_irrs)==TRUE) {
            stop("clust_irrs input missing")
        }
        if(is.list(clust_irrs)==FALSE) {
            stop("clust_irrs must be a list of clust_irr objects")
        }
        if(length(clust_irrs)<=1) {
            stop("get_joint_graph needs >= 2 clust_irr outputs")
        }
        gmd <- numeric(length = length(clust_irrs))
        for(i in 1:length(clust_irrs)) {
            check_clustirr(clust_irr = clust_irrs[[i]])
            gmd[i]<-get_clustirr_inputs(clust_irrs[[i]])$control$global_max_dist
        }
        if(length(unique(gmd))!=1) {
            stop("all global_max_dist should be equal")
        }
        
        # check if same chain names
        cs <- colnames(get_clustirr_inputs(clust_irrs[[1]])$s)
        for(i in 2:length(clust_irrs)) {
            if(any(colnames(get_clustirr_inputs(clust_irrs[[i]])$s)!=cs)) {
                stop("different chains in graphs")
            }
        }
    }
    
    # check input
    check_input(clust_irrs = clust_irrs)
    
    jg <- get_joint_graph(clust_irrs) 
    
    if(is.null(jg$graph)) {
        warning("No graph to plot \n")
        return(jg)
    }
    
    # make graph look visually better
    jg$graph <- config_vertices_plot(g = jg$graph, is_jg = TRUE)
    jg$graph <- config_edges_plot(g = jg$graph, is_jg = TRUE)
    
    # plot
    if(as_visnet == FALSE) {
        plot(jg$graph, vertex.label = NA)
    }
    if(as_visnet == TRUE) {
        V(jg$graph)$size <- V(jg$graph)$size*10
        visIgraph(igraph = jg$graph,
                  idToLabel = TRUE,
                  layout = "layout_components",
                  randomSeed = 1234,
                  physics = FALSE,
                  smooth = FALSE,
                  type = "square")
    }
}

