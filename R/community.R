detect_communities <- function(graph, 
                               algorithm = "leiden", 
                               resolution = 1,
                               weight_type = "ncweight",
                               chains) {
    
    message("1/5 formatting graph (g)...")
    cg <- get_formatted_graph(g = graph, 
                              weight_type = weight_type, 
                              chains = chains) 
    
    message("2/5 community detection...")
    cg <- get_community_detection(g = cg, 
                                  algorithm = algorithm, 
                                  resolution = resolution)
    
    message("3/5 community summary (cs)...")
    V(graph)$community <- V(cg)$community
    cs <- get_community_summary(g = graph, 
                                chains = chains,
                                weight_type = weight_type)
    
    message("4/5 extracting community matrix (cm)...")
    cm <- get_community_matrix(g = cg)
    
    message("5/5 extracting vertices...")
    vs <- as_data_frame(x = cg, what = "vertices")
    
    # save configs
    config <- list(input_g = graph, 
                   algorithm = algorithm, 
                   resolution = resolution,
                   weight_type = weight_type, 
                   chains = chains)
    
    return(list(community_matrix = cm, 
                community_summary = cs, 
                node_summary = vs, 
                graph = cg, 
                input_config = config))
}

get_formatted_graph <- function(g, 
                                weight_type, 
                                chains) {
    
    set_weight <- function(g, weight_type) {
        if(weight_type == "cweight") {
            E(g)$weight <- E(g)$cweight
        }
        if(weight_type == "ncweight") {
            E(g)$weight <- E(g)$ncweight
        }
        if(weight_type == "weight") {
            E(g)$weight <- E(g)$ncweight
        }
        if(weight_type == "nweight") {
            E(g)$weight <- E(g)$cweight
        }
        return(g)
    }
    
    set_chain <- function(g, chains) {
        i <- which(!E(g)$chain %in% chains)
        if(length(i)!=0) {
            g <- delete_edges(graph = g, i)
        }
        return(g)
    }
    
    set_concat <- function(x) {
        paste0(unique(sort(x)), collapse = '+')
    }
    
    g <- set_weight(g = g, weight_type = weight_type)
    g <- set_chain(g = g, chains = chains)
    
    # global graph
    i <- which(E(g)$clustering == "local")
    if(length(i)!=0) {
        g <- delete_edges(graph = g, i)
    }
    
    if(is_simple(g)==FALSE) {
        g <- simplify(g, edge.attr.comb = list(weight = "sum",
                                               type = "first",
                                               chain = "concat",
                                               clustering = "concat",
                                               "ignore"))
        # average weight
        E(g)$weight <- E(g)$weight/2
        E(g)$clustering <- unlist(lapply(X = E(g)$clustering, FUN = set_concat))
        E(g)$chain <- unlist(lapply(X = E(g)$chain, FUN = set_concat))
    }
    
    g <- delete_edges(graph = g, which(E(g)$weight <= 0))
    return(g)
}

get_community_detection <- function(g, 
                                    algorithm, 
                                    resolution) {
    
    if(algorithm == "louvain") {
        c <- cluster_louvain(graph = g, 
                             weights = E(g)$weight, 
                             resolution = resolution)
        V(g)$community <- c$membership
    }
    if(algorithm == "leiden") {
        c <- cluster_leiden(graph = g, 
                            weights = E(g)$weight, 
                            resolution_parameter = resolution)
        V(g)$community <- c$membership
    }
    return(g)
}

get_community_summary <- function(g, chains, weight_type) {
    
    get_vs_stats <- function(vs) {
        
        vs$f <- 1
        
        # number of cells
        vs_cells <- aggregate(clone_size~community+sample, data = vs, FUN = sum)
        vs_cells <- acast(data = vs_cells, formula = community~sample, 
                          value.var = "clone_size", fill = 0)
        vs_cells <- data.frame(vs_cells)
        vs_cells$n <- apply(X = vs_cells, MARGIN = 1, FUN = sum)
        colnames(vs_cells) <- paste0("cells_", colnames(vs_cells))
        vs_cells$community <- rownames(vs_cells)
        vs_cells$community <- as.numeric(as.character(vs_cells$community))
        vs_cells <- vs_cells[order(vs_cells$community, decreasing = FALSE), ]
        
        
        # number of clones
        vs_clones <- aggregate(f~community+sample, data = vs, FUN = sum)
        vs_clones <- acast(data = vs_clones, formula = community~sample, 
                           value.var = "f", fill = 0)
        vs_clones <- data.frame(vs_clones)
        vs_clones$n <- apply(X = vs_clones, MARGIN = 1, FUN = sum)
        colnames(vs_clones) <- paste0("clones_", colnames(vs_clones))
        vs_clones$community <- rownames(vs_clones)
        vs_clones$community <- as.numeric(as.character(vs_clones$community))
        vs_clones <- vs_clones[order(vs_clones$community, decreasing = FALSE), ]
        
        # merge clones and cells
        vs_stats <- merge(x = vs_clones, y = vs_cells, by = "community")
        return(vs_stats)
    }
    
    get_es_stats <- function(es, vs, chains, weight_type) {
        
        # keep only global
        es <- es[es$clustering == "global", ]
        
        # add community id to 'from node'
        es <- merge(x = es, y = vs[, c("name", "community")], 
                    by.x = "from", by.y = "name", all.x = TRUE)
        es$from_community <- es$community
        es$community <- NULL
        
        # add community id to 'to node'
        es <- merge(x = es, y = vs[, c("name", "community")], 
                    by.x = "to", by.y = "name", all.x = TRUE)
        es$to_community <- es$community
        es$community <- NULL
        
        # within-community stats is what we care about
        es <- es[which(es$from_community == es$to_community), ]
        if(nrow(es)!=0) {
            es$n_edges <- 1
            es$community <- es$from_community
        } else {
            stop("no edges to summarize")
        }
        
        es$key <- apply(X = es[,c("from", "to")], MARGIN = 1, FUN=function(x) {
            return(paste0(sort(x), collapse = '-'))
        })
        
        # set weight
        es$weight <- es[, weight_type]
        
        if(length(chains)==2) {
            # keys
            k1 <- paste0("w_", chains[1])
            k2 <- paste0("w_", chains[2])
            k12 <- paste0("w_", chains[1], '_', chains[2])
            n1 <- paste0("n_", chains[1])
            n2 <- paste0("n_", chains[2])
            
            # chain 1
            es_1 <- es[es$chain == chains[1], c("key", "weight", "community")]
            es_1[, k1] <- es_1$weight
            
            # chain 2
            es_2 <- es[es$chain == chains[2], c("key", "weight", "community")]
            es_2[, k2] <- es_2$weight
            
            # merge chains
            es_w <- merge(x = es_1[c("key", k1, "community")], 
                          y = es_2[c("key", k2, "community")], 
                          by = c("key", "community"),
                          all = TRUE)
            
            es_w[, n1] <- ifelse(test = is.na(es_w[, k1])|es_w[, k1]<0, 
                                 yes = 0, no = 1)
            es_w[, n2] <- ifelse(test = is.na(es_w[, k2])|es_w[, k2]<0, 
                                 yes = 0, no = 1)
            
            es_w[, k1] <- ifelse(test = is.na(es_w[, k1])|es_w[, k1]<0, 
                                 yes = 0, no = es_w[, k1])
            es_w[, k2] <- ifelse(test = is.na(es_w[, k2])|es_w[, k2]<0, 
                                 yes = 0, no = es_w[, k2])
            es_w[, k12] <- (es_w[, k1]+es_w[, k2])/2
            
            
            a <- merge(
                x = merge(
                    x = aggregate(es_w[[k1]]~community, data = es_w, FUN=mean),
                    y = aggregate(es_w[[k1]]~community, data = es_w, FUN=var),
                    by = "community"),
                y = aggregate(es_w[[n1]]~community, data = es_w, FUN = sum),
                by = "community")
            colnames(a) <- c("community", paste0(k1, "_mean"), 
                             paste0(k1, "_var"), n1)
            
            b <- merge(
                x = merge(
                    x = aggregate(es_w[[k2]]~community, data = es_w, FUN=mean),
                    y = aggregate(es_w[[k2]]~community, data = es_w, FUN=var),
                    by = "community"),
                y = aggregate(es_w[[n2]]~community, data = es_w, FUN = sum),
                by = "community")
            colnames(b) <- c("community", paste0(k2, "_mean"), 
                             paste0(k2, "_var"), n2)
            
            c <- merge(
                x = aggregate(es_w[[k12]]~community, data = es_w, FUN=mean),
                y = aggregate(es_w[[k12]]~community, data = es_w, FUN=var),
                by = "community")
            colnames(c) <- c("community", paste0(k12, "_mean"), 
                             paste0(k12, "_var"))
            
            edges_stats <- merge(x = merge(x = a, y = b, by = "community"),
                                 y = c, by = "community")
        }
        if(length(chains)==1) {
            # keys
            k1 <- paste0("w_",chains[1])
            n1 <- paste0("n_",chains[1])
            
            es_w <- es[es$chain == chains[1], c("key", "weight", "community")]
            es_w[, k1] <- es_w$weight
            es_w <- es_w[c("key", k1, "community")]
            
            es_w[, n1] <- ifelse(test = is.na(es_w[, k1])|es_w[, k1]<0, 
                                 yes = 0, no = 1)
            
            edges_stats <- es_w
        }
        
        return(edges_stats)
    }
    
    es <- as_data_frame(x = g, what = "edges")
    vs <- as_data_frame(x = g, what = "vertices")
    
    vs_stats <- get_vs_stats(vs = vs)
    es_stats <- get_es_stats(es = es, 
                             vs = vs, 
                             chains = chains,
                             weight_type = weight_type)
    
    o <- merge(x = vs_stats, y = es_stats, by = "community", all.x = TRUE)
    o <- o[order(o$community, decreasing = FALSE), ]
    return(o)
}

get_community_matrix <- function(g) {
    vs <- as_data_frame(x = g, what = "vertices")
    
    cm <- acast(data = vs, formula = community~sample, 
                value.var = "clone_size", 
                fun.aggregate = sum, fill = 0)
    
    return(cm)
}
