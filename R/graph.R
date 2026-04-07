
get_graph <- function(clust_irr) {
    
    get_edges <- function(clust_irr) {
        eg <- vector(mode="list", length =length(get_clustirr_clust(clust_irr)))
        names(eg) <- names(get_clustirr_clust(clust_irr))
        for(chain in names(get_clustirr_clust(clust_irr))) {
            g <- get_clustirr_clust(clust_irr)[[chain]]
            if(is.null(g) == FALSE && nrow(g) != 0) {
                eg[[chain]] <- g
            }
        }
        eg <- do.call(rbind, eg)
        return(eg)
    }
    
    get_edge_order <- function(x) {
        return(paste0(sort(c(x[1], x[2])), collapse = '-'))
    }
    
    build_graph <- function(e, cs, sample_id, chains) {
        
        add_edges <- function(e, ig) {
            es <- as.vector(apply(X = e[, c("from", "to")], MARGIN = 1,
                                     FUN = function(x) {return(x)}))
            
            ig <- igraph::add_edges(graph = ig, 
                                    edges = es,
                                    weight = e$weight,
                                    cweight = e$cweight,
                                    nweight = e$nweight,
                                    ncweight = e$ncweight,
                                    type = "within-repertoire",
                                    chain = e$chain)
            
            return(ig)
        }

        ig <- graph_from_data_frame(d = data.frame(from = cs$name[1], 
                                                   to = cs$name[1]),
                                    directed = FALSE,
                                    vertices = cs)
        ig <- delete_edges(ig, edges = 1)
        
        # add global edges
        if(is.null(e)==FALSE && nrow(e)!=0) {
            for(chain in chains) {
                chain_ge <- e[e$chain == chain, ]
                if(nrow(chain_ge)!=0) {
                    ig <- add_edges(e = chain_ge, ig = ig)
                }
            }
        }
        
        return(ig)
    }
    
    check_clustirr(clust_irr = clust_irr)
    
    # get chains
    chains <- get_chains(x = colnames(get_clustirr_inputs(clust_irr)$s))
    
    # cells
    s <- get_clustirr_inputs(clust_irr)$s
    # setting up the sample id
    sample_id <- s$sample[1]
    
    # meta
    meta <- get_clustirr_inputs(clust_irr)$meta
    
    # get clones
    cs <- get_clones(sample_id = sample_id, s = s, meta = meta)
    
    # get edges between clones
    e <- get_edges(clust_irr = clust_irr)
    
    # build graph with only vertices
    if(is.null(e)) {
        ig <- graph_from_data_frame(d = data.frame(from = cs$name[1], 
                                                   to = cs$name[1]), 
                                    directed = FALSE, vertices = cs)
        ig <- delete_edges(ig, edges = 1)
        
        return(list(graph = ig, clust_irrs = clust_irr, 
                    clones = cs, multigraph = FALSE))
    }
    
    # build graph
    message("creating graph (", sample_id, ") \n")
    ig <- build_graph(e = e, cs = cs, sample_id = sample_id, chains = chains)
    
    return(list(graph = ig, clust_irrs = clust_irr, 
                clones = cs, multigraph = FALSE))
}


get_joint_graph <- function(clust_irrs) {
    
    check_input <- function(clust_irrs) {
        if(missing(clust_irrs)==TRUE) {
            stop("clust_irrs input missing")
        }
        if(is.list(clust_irrs)==FALSE) {
            stop("clust_irrs must be a list of clust_irr objects")
        }
        for(i in 1:length(clust_irrs)) {
            check_clustirr(clust_irr = clust_irrs[[i]])
        }
        
        # check if same chain names
        cs <- colnames(get_clustirr_inputs(clust_irrs[[1]])$s)
        for(i in 2:length(clust_irrs)) {
            if(any(colnames(get_clustirr_inputs(clust_irrs[[i]])$s)!=cs)) {
                stop("different chains in graphs")
            }
        }
    }
    
    get_clust_irrs_names <- function(clust_irrs) {
        clust_irrs_names <- names(clust_irrs)
        if(is.null(clust_irrs_names)) {
            names(clust_irrs) <- paste0("s", 1:length(clust_irrs))
        }
        
        for(i in 1:length(clust_irrs)) {
            clust_irrs[[i]]@inputs[["sample_id"]] <- names(clust_irrs)[i]
        }
        
        return(clust_irrs)
    }
    
    get_v_e <- function(x, what) {
        return(igraph::as_data_frame(x$graph, what = what))
    }
    
    # check input
    check_input(clust_irrs = clust_irrs)
    
    # get clust_irrs names
    clust_irrs <- get_clust_irrs_names(clust_irrs = clust_irrs)
    
    # get joint controls
    control <- get_joint_controls(clust_irrs = clust_irrs)
    
    # build graphs
    message("[1/3] generating individual graphs... \n")
    igs <- lapply(X = clust_irrs, FUN = get_graph)
    names(igs) <- names(clust_irrs)
    
    # get chains
    chains <- get_chains(x = colnames(get_clustirr_inputs(clust_irrs[[1]])$s))
    
    message("[2/3] joining graphs... \n")
    ige <- get_intergraph_edges(igs = igs, chains = chains, control = control)
    
    # get the vertices/edges of the graph
    df_v <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "vertices"))
    df_e <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "edges"))
    
    # these are the cols we want to keep in this order
    cols <- c("from", "to", "weight", "cweight", "nweight", "ncweight",
              "type", "chain")
    
    if(nrow(df_e)!=0) {
        if(is.null(ige)==FALSE && nrow(ige)!=0) {
            df_e <- rbind(df_e[, cols], ige[, cols])
        }
    } else {
        if(is.null(ige)==FALSE && nrow(ige)!=0) {
            df_e <- ige[, cols]
        }
    }
    
    message("[3/3] optimizing joint graph... \n")
    # delete duplicate edges (if any)
    key <- data.frame(pmin(df_e$from, df_e$to), 
                      pmax(df_e$from, df_e$to), df_e$chain)
    df_e <- df_e[!duplicated(key), ]
    
    # clean unused vars
    rm(igs, ige)
    
    # build joint graph
    g <- graph_from_data_frame(df_e, directed = FALSE, vertices = df_v)
    
    return(list(graph = g, clust_irrs = clust_irrs, multigraph = TRUE))
}


plot_graph <- function(g,
                       select_by = "Ag_species",
                       as_visnet = FALSE, 
                       show_singletons = TRUE,
                       node_opacity = 1) {
    
    check_graph <- function(g) {
        if(missing(g)) {
            stop("missing input g")
        }
        if(is.list(g)==FALSE) {
            stop("missing input g")
        }
        if(is_igraph(g$graph)==FALSE) {
            stop("wrong input g")
        }
        if(is.null(g$graph)) {
            stop("wrong input g")
        }
        if(is.logical(g$multigraph)==FALSE) {
            stop("wrong input g")
        }
    }
    
    check_graph(g = g)
    check_as_visnet(as_visnet = as_visnet)
    check_show_singletons(show_singletons = show_singletons)
    check_select_by(select_by = select_by)
    check_node_opacity(node_opacity = node_opacity)
    
    # unpack g
    ig <- g$graph
    cs <- igraph::as_data_frame(ig, what = "vertices")
    is_jg <- g$multigraph
    
    if(!show_singletons){
        k <- which(degree(ig) == 0 & V(ig)$clone_size <= 1)
        if(length(k)!=0) {
            ig <- delete_vertices(ig, k)
        }
    }
    
    ig <- config_vertices_plot(g = ig, is_jg = is_jg, 
                               node_opacity = node_opacity)
    if(as_visnet == FALSE) {
        plot(ig, vertex.label = NA)
    }
    if(as_visnet == TRUE) {
        V(ig)$size <- V(ig)$size*5
        if(length(E(ig))==0) {
            message("no edges found in the graph, dummy self-edge added")
            # apparently if no edges, visnetwork can't plot
            ig <- add_edges(graph = ig, edges = c(V(ig)[1],V(ig)[1]))
        } 
        else {
            if(all(E(ig)$width == 1)==FALSE) {
                E(ig)$width <- 2*(1/(1+exp(-(-6-0.1*E(ig)$weight))))
            }
        }
        
        vi <- visIgraph(igraph = ig,
                        idToLabel = TRUE,
                        layout = "layout_components",
                        randomSeed = 1234,
                        physics = FALSE,
                        smooth = FALSE,
                        type = "square")
        
        # set vertex titles -> CDR sequences
        chains <- get_chains(x = colnames(vi$x$nodes))
        vi$x$nodes$title <- apply(X = vi$x$nodes[, chains, drop = FALSE], 
                                  y = chains, 
                                  MARGIN = 1, FUN = function(x, y) {
                                      paste0(paste0("<b>", y, "</b>", ':', x), 
                                             collapse = '<br>')})
        
        u <- unique(unlist(strsplit(x=unique(cs[, select_by]),split = ',')))
        vi <- vi %>% visOptions(selectedBy = list(variable = select_by,
                                                  values = u,
                                                  multiple = TRUE))
        
        return(vi)
    }
}

get_intergraph_edges <- function(igs,
                                 chains,
                                 control) {

    get_igg <- function(x, ix, igs, control) {
        # prepare pair-rep data
        s1_name <- ix$name_i[x]
        s2_name <- ix$name_j[x]
        chain <- ix$chain[x]
        s1 <- igs[[ix$index_i[x]]]$clones
        s2 <- igs[[ix$index_j[x]]]$clones
        rm(igs)
        
        message("joining (", x, ") ", s1_name, ' and ', s2_name, "\n")
        
        s1$chain <- chain
        s1$cdr3 <- s1[, chain]
        s2$chain <- chain
        s2$cdr3 <- s2[, chain]
        
        b <- get_score_pair(s_from = s1, s_to = s2, control = control)
        
        if(is.null(b)==FALSE && nrow(b)!=0) {
            b$chain <- chain
            b$sample <- paste0(s1_name, "|", s2_name)
            b$type <- "between-repertoire"
            return(b)
        }
        
        return(NULL)
    }
    
    get_ix <- function(xs, ns, chains, control) {
        ixs <- c()
        for(c in chains) {
            if(control$knn==TRUE) {
                ix <- expand.grid(1:xs, 1:xs, KEEP.OUT.ATTRS = FALSE)
                colnames(ix) <- c("index_i", "index_j")
                ix <- ix[ix$index_i!=ix$index_j,]
                ix$name_i <- ns[ix$index_i]
                ix$name_j <- ns[ix$index_j]
                ix$chain <- c
                ixs <- rbind(ixs, ix)
            } else {
                ix <- subset(expand.grid(1:xs, 1:xs, KEEP.OUT.ATTRS = FALSE), 
                             Var1 < Var2)
                colnames(ix) <- c("index_i", "index_j")
                ix$name_i <- ns[ix$index_i]
                ix$name_j <- ns[ix$index_j]
                ix$chain <- c
                ixs <- rbind(ixs, ix)
            }
        }
        return(ixs)
    }
    
    # indices 
    ix <- get_ix(xs = length(igs), ns = names(igs), 
                 chains = chains, control = control)
    
    # find global similarities between pairs of clone tables
    message("merging clust_irrs: ", nrow(ix), "\n")
    ige <- lapply(X = seq_len(nrow(ix)),
                  ix = ix,
                  FUN = get_igg,
                  igs = igs,
                  control = control)
    ige <- do.call(rbind, ige)

    return(ige)
}

get_clones <- function(sample_id, s, meta) {
    s$name <- s$id
    s <- s[, rev(colnames(s))]
    
    # append s and meta
    if(is.null(meta)==FALSE & missing(meta)==FALSE) {
        cols <- intersect(colnames(s), colnames(meta))
        if(length(cols)!=0) {
            meta <- meta[, (colnames(meta) %in% cols)==FALSE, drop = FALSE]
        }
        if(ncol(meta)!=0) {
            s <- cbind(s, meta)
        }
    }
    
    return(s)
}

# Description:
# Compare the control lists of individual clust_irr objects. If they match,
# then extract the control list from one list.
get_joint_controls <- function(clust_irrs) {
    # if these tests pass -> extract controls of any of the elements
    for(i in 1:(length(clust_irrs)-1)) {
        for(j in (i+1):length(clust_irrs)) {
            c_i <- get_clustirr_inputs(clust_irrs[[i]])$control
            c_j <- get_clustirr_inputs(clust_irrs[[j]])$control
            
            if(length(c_i)!=length(c_j)) {
                stop("different controls used in individual clust_irrs")
            }
            
            c_ij <- vapply(X = names(c_i), c_i = c_i, c_j = c_j,
                           FUN = function(x, c_i, c_j) {
                               if(is.null(c_i[[x]]) & is.null(c_j[[x]])) {
                                   return(TRUE)
                               }
                               if(is.na(c_i[[x]]) & is.na(c_j[[x]])) {
                                   return(TRUE)
                               }
                               return(c_i[[x]]==c_j[[x]])}, 
                           FUN.VALUE = logical(1))
            
            if(any(c_ij==FALSE)) {
                stop("different controls used in individual clust_irrs")
            }
        } 
    }
    
    return(get_clustirr_inputs(x = clust_irrs[[1]])$control)
}

