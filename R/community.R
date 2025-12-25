detect_communities <- function(graph, 
                               weight = "nweight",
                               algorithm = "leiden", 
                               metric = "average",
                               resolution = 1,
                               iterations = 100,
                               chains) {
    
    check_inputs(graph = graph,
                 algorithm = algorithm, 
                 metric = metric,
                 resolution = resolution,
                 iterations = iterations,
                 weight = weight,
                 chains = chains)
    
    message("[1/5] formatting graph...")
    cg <- get_formatted_graph(graph = graph, 
                              metric = metric,
                              weight = weight, 
                              chains = chains) 
    
    message("[2/5] community detection...")
    cg <- get_community_detection(g = cg, 
                                  algorithm = algorithm, 
                                  resolution = resolution,
                                  iterations = iterations)
    
    message("[3/5] community summary...")
    cs <- get_community_summary(g = cg$graph, chains = chains)
    
    message("[4/5] extracting community occupancy matrix...")
    cm <- get_community_matrix(g = cg$graph)
    
    message("[5/5] extracting nodes")
    vs <- igraph::as_data_frame(x = cg$graph, what = "vertices")
    
    # save configs
    config <- list(input_g = graph, 
                   algorithm = algorithm, 
                   resolution = resolution,
                   weight = weight,
                   metric = metric,
                   chains = chains)
    
    return(list(community_occupancy_matrix = cm, 
                community_summary = cs, 
                node_summary = vs, 
                graph = cg$graph, 
                graph_structure_quality = list(modularity = cg$modularity,
                                               quality = cg$quality),
                input_config = config))
}

get_formatted_graph <- function(graph,
                                metric,
                                weight,
                                chains) {
    
    set_chain <- function(graph, chains) {
        i <- which(!E(graph)$chain %in% chains)
        if(length(i)!=0) {
            graph <- delete_edges(graph = graph, i)
        }
        return(graph)
    }
    
    graph <- set_chain(graph = graph, chains = chains)
    
    graph <- simplify(graph, edge.attr.comb = list(ncweight = "concat",
                                                   nweight = "concat",
                                                   chain = "concat",
                                                   "ignore"))
    
    if(metric == "average") {
        E(graph)$nweight <- vapply(X = E(graph)$nweight, 
                                   FUN.VALUE = numeric(1),
                                   FUN = function(x) {return(sum(x)/2)})
        E(graph)$ncweight <- vapply(X = E(graph)$ncweight, 
                                    FUN.VALUE = numeric(1),
                                    FUN = function(x) {return(sum(x)/2)})
    }
    if(metric == "max") {
        E(graph)$nweight <- vapply(X = E(graph)$nweight, 
                                   FUN.VALUE = numeric(1),
                                   FUN = function(x) {return(max(x))})
        E(graph)$ncweight <- vapply(X = E(graph)$ncweight, 
                                    FUN.VALUE = numeric(1),
                                    FUN = function(x) {return(max(x))})
     }
    
    if(weight == "nweight") {
        E(graph)$w <- E(graph)$nweight
    } else {
        E(graph)$w <- E(graph)$ncweight
    }
    
    
    # if trim*2 > CDR3 lengths -> NA
    i <- which(E(graph)$w <= 0 | is.na(E(graph)$w))
    if(length(i)!=0) {
        graph <- delete_edges(graph = graph, i)
    }
    return(graph)
}

get_community_detection <- function(g, 
                                    algorithm, 
                                    resolution,
                                    iterations) {
    
    m <- NA
    q <- NA
    if(algorithm == "louvain") {
        c <- cluster_louvain(graph = g, 
                             weights = E(g)$w, 
                             resolution = resolution)
        V(g)$community <- c$membership
        m <- c$modularity[length(c$modularity)]
    }
    if(algorithm == "leiden") {
        c <- cluster_leiden(graph = g, 
                            weights = E(g)$w, 
                            resolution = resolution,
                            n_iterations = iterations)
        V(g)$community <- c$membership
        m <- modularity(g, membership = c$membership, weights = E(g)$w, 
                        resolution = resolution, directed = FALSE)
        q <- c$quality
    }
    if(algorithm == "infomap") {
        c <- cluster_infomap(graph = g, 
                             e.weights = E(g)$w, 
                             nb.trials = iterations)
        V(g)$community <- c$membership
        m <- c$modularity[length(c$modularity)]
    }
    
    return(list(graph = g, quality = q, modularity = m))
}

get_community_summary <- function(g, chains) {
    
    get_estats <- function(e, v, chains) {
        e <- merge(x = e, y = v[,c("name", "community")], 
                   all.x = TRUE, by.x = "from", by.y = "name")
        e <- merge(x = e, y = v[,c("name", "community", "sample")], 
                   all.x = TRUE, by.x = "to", by.y = "name")
        e <- e[e$community.x==e$community.y,]
        if(nrow(e)==0) {
            return(NULL)
        }
        
        e$community <- e$community.x
        e <- e[,c("community", "chain", "nweight", "ncweight", "w")]
        e$n_edges <- 1
        es <- vector(mode = "list", length = length(chains))
        names(es) <- chains
        
        for(c in chains) {
            x <- e[e$chain == c, ] %>%
                group_by(community) %>%
                summarise(across(c(nweight, ncweight, w), mean, na.rm = TRUE),
                          n_edges = sum(n_edges, na.rm = TRUE), 
                          .groups = "drop") %>%
                ungroup()
            x$chain <- NULL  
            colnames(x) <- c("community",paste0(c("nweight", "ncweight", 
                                                  "w", "n_edges"), "_", c))
            es[[c]] <- x
        }
        
        if(length(chains)>1) {
            es <- merge(x = es[[chains[1]]], y = es[[chains[2]]], 
                        by = "community", all = TRUE)
            es[is.na(es)] <- 0
            es <- es[,sort(colnames(es))]
        } else {
            es <- es[[chains[1]]]
        }
        return(es)
    }
    
    get_vstats <- function(v, wide) {
        
        v$cells <- v$clone_size
        v$clones <- 1
        
        if(wide) {
            # number of cells
            vcells <- aggregate(cells~community+sample, data = v, FUN = sum)
            vcells <- acast(data = vcells, formula = community~sample, 
                            value.var = "cells", fill = 0)
            vcells <- data.frame(vcells)
            vcells$n <- apply(X = vcells, MARGIN = 1, FUN = sum)
            colnames(vcells) <- paste0("cells_", colnames(vcells))
            vcells$community <- rownames(vcells)
            vcells$community <- as.numeric(as.character(vcells$community))
            vcells <- vcells[order(vcells$community, decreasing = FALSE), ]
            
            # number of clones
            vclones <- aggregate(clones~community+sample, data = v, FUN = sum)
            vclones <- acast(data = vclones, formula = community~sample, 
                             value.var = "clones", fill = 0)
            vclones <- data.frame(vclones)
            vclones$n <- apply(X = vclones, MARGIN = 1, FUN = sum)
            colnames(vclones) <- paste0("clones_", colnames(vclones))
            vclones$community <- rownames(vclones)
            vclones$community <- as.numeric(as.character(vclones$community))
            vclones <- vclones[order(vclones$community, decreasing = FALSE), ]
            
            # merge clones and cells
            vstats <- merge(x = vclones, y = vcells, by = "community")
        } 
        else {
            # number of cells
            vcells <- aggregate(cells~community+sample, data = v, 
                                FUN = sum, drop = FALSE)
            
            # number of clones
            vclones <- aggregate(clones~community+sample, data = v, 
                                 FUN = sum, drop = FALSE)
            
            # merge clones and cells
            vstats <- merge(x = vclones, y = vcells, 
                            by = c("community", "sample"))
            vstats$cells[is.na(vstats$cells)] <- 0
            vstats$clones[is.na(vstats$clones)] <- 0
        }
        
        return(vstats)
    }
    
    b <- igraph::as_data_frame(x = g, what = "both")
    es <- get_estats(e = b$edges, v = b$vertices, chains = chains)
    vs <- get_vstats(v = b$vertices, wide = TRUE)
    
    # merge results
    cs <- merge(x = vs, y = es, by = "community", all.x = TRUE)
    cs <- cs[order(cs$community, decreasing = FALSE), ]
    cs[is.na(cs)] <- 0
    return(cs)
}

get_community_matrix <- function(g) {
    vs <- igraph::as_data_frame(x = g, what = "vertices")
    
    cm <- acast(data = vs, formula = community~sample, 
                value.var = "clone_size", 
                fun.aggregate = sum, fill = 0)
    
    return(cm)
}

check_inputs <- function(graph, 
                         algorithm, 
                         metric,
                         resolution,
                         iterations,
                         weight, 
                         chains) {
    
    
    # check graph
    if(missing(graph)) {
        stop("graph must be an igraph object")
    }
    if(is_igraph(graph)==FALSE) {
        stop("graph must be an igraph object")
    }
    
    # check algorithm
    if(missing(algorithm)) {
        stop("algorithm must be louvain, leiden or infomap")
    }
    if(length(algorithm)!=1) {
        stop("algorithm must be louvain, leiden or infomap")
    }
    if(is.character(algorithm)==FALSE) {
        stop("algorithm must be character")
    }
    if(!algorithm %in% c("louvain", "leiden", "infomap")) {
        stop("algorithm must be louvain, leiden or infomap")
    }
    
    # check metric
    if(missing(metric)) {
        stop("metric must be average or max")
    }
    if(length(metric)!=1) {
        stop("metric must be average or max")
    }
    if(is.character(metric)==FALSE) {
        stop("metric must be character")
    }
    if(!metric %in% c("average", "max")) {
        stop("metric must be average or max")
    }
    
    # check resolution
    if(missing(resolution)) {
        stop("resolution must be a number > 0")
    }
    if(length(resolution)!=1) {
        stop("resolution must be a number > 0")
    }
    if(is.numeric(resolution)==FALSE) {
        stop("resolution must be a number > 0")
    }
    if(is.finite(resolution)==FALSE) {
        stop("resolution must be a number > 0")
    }
    if(resolution<=0) {
        stop("resolution must be a number > 0")
    }
    
    
    # check iterations
    if(missing(iterations)) {
        stop("iterations must be a number > 0")
    }
    if(length(iterations)!=1) {
        stop("iterations must be a number > 0")
    }
    if(is.numeric(iterations)==FALSE) {
        stop("iterations must be a number > 0")
    }
    if(is.finite(iterations)==FALSE) {
        stop("iterations must be a number > 0")
    }
    if(iterations<=0) {
        stop("iterations must be a number > 0")
    }
    
    
    # check weight
    if(missing(weight)) {
        stop("weight must be ncweight or nweight")
    }
    if(length(weight)!=1) {
        stop("weight must be ncweight or nweight")
    }
    if(is.character(weight)==FALSE) {
        stop("weight must be character")
    }
    if(!weight %in% c("ncweight", "nweight")) {
        stop("weight must be ncweight or nweight")
    }
    
    
    # check chains
    if(missing(chains)) {
        stop("chains must be a character vector")
    }
    if(length(chains)<1 | length(chains)>2) {
        stop("chains must be a character vector")
    }
    if(is.character(chains)==FALSE) {
        stop("chains must be a character vector")
    }
    if(any(chains %in% c("CDR3a", "CDR3b", "CDR3g", 
                         "CDR3d", "CDR3h", "CDR3l"))==FALSE) {
        stop("chains must be a character vector")
    }
}

set_chain <- function(graph, chains) {
    i <- which(!E(graph)$chain %in% chains)
    if(length(i)!=0) {
        graph <- delete_edges(graph = graph, i)
    }
    return(graph)
}
