detect_communities <- function(graph, 
                               algorithm = "leiden", 
                               resolution = 1,
                               weight = "ncweight",
                               metric = "average",
                               chains) {
    
    check_inputs(graph = graph,
                 algorithm = algorithm, 
                 resolution = resolution,
                 weight = weight, 
                 metric = metric, 
                 chains = chains)
    
    message("1/5 formatting graph (g)...")
    cg <- get_formatted_graph(graph = graph, 
                              weight = weight, 
                              metric = metric,
                              chains = chains) 
    
    message("2/5 community detection...")
    cg <- get_community_detection(g = cg, 
                                  algorithm = algorithm, 
                                  resolution = resolution)
    
    message("3/5 community summary (cs)...")
    cs <- get_community_summary(g = cg, chains = chains)
    
    message("4/5 extracting community matrix (cm)...")
    cm <- get_community_matrix(g = cg)
    
    message("5/5 extracting vertices...")
    vs <- as_data_frame(x = cg, what = "vertices")
    
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
                graph = cg, 
                input_config = config))
}

get_formatted_graph <- function(graph, 
                                weight,
                                metric,
                                chains) {
    
    set_weight <- function(graph, weight) {
        if(weight == "ncweight") {
            E(graph)$weight <- E(graph)$ncweight
        }
        if(weight == "nweight") {
            E(graph)$weight <- E(graph)$nweight
        }
        return(graph)
    }
    
    set_chain <- function(graph, chains) {
        i <- which(!E(graph)$chain %in% chains)
        if(length(i)!=0) {
            graph <- delete_edges(graph = graph, i)
        }
        return(graph)
    }
    
    graph <- set_weight(graph = graph, weight = weight)
    graph <- set_chain(graph = graph, chains = chains)
    
    if(is_simple(graph)==FALSE) {
        graph <- simplify(graph, edge.attr.comb = list(weight = "concat",
                                                       chain = "concat",
                                                       "ignore"))
        
        if(metric == "average") {
            E(graph)$w <- vapply(X = E(graph)$weight, FUN.VALUE = numeric(1),
                                      FUN = function(x) {return(sum(x)/2)})
        }
        if(metric == "strict") {
            E(graph)$w <- vapply(X = E(graph)$weight, FUN.VALUE = numeric(1),
                                      FUN = function(x) {
                                          if(length(x)==2) {
                                              return(min(x))
                                          }
                                          return(min(x,0))})
        }
        if(metric == "loose") {
            E(graph)$w <- vapply(X = E(graph)$weight, FUN.VALUE = numeric(1),
                                      FUN = function(x) {
                                          if(length(x)==2) {
                                              return(max(x))
                                          }
                                          return(max(x,0))})
        }
    }
    
    graph <- delete_edges(graph = graph, which(E(graph)$w <= 0))
    return(graph)
}

get_community_detection <- function(g, 
                                    algorithm, 
                                    resolution) {
    
    if(algorithm == "louvain") {
        c <- cluster_louvain(graph = g, weights = E(g)$w, 
                             resolution = resolution)
        V(g)$community <- c$membership
    }
    if(algorithm == "leiden") {
        c <- cluster_leiden(graph = g, weights = E(g)$w, 
                            resolution_parameter = resolution)
        V(g)$community <- c$membership
    }
    return(g)
}

get_community_summary <- function(g, 
                                  chains, 
                                  metric) {
    
    get_vstats <- function(vs) {
        
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
    
    get_estats <- function(x, g, chains) {
        sg <- subgraph(graph = g, vids = which(V(g)$community==x))
        
        if(length(sg)==1) {
            v <- numeric(length = length(chains)*2+1)
            names(v) <- c("w", paste0("w_", chains), paste0("n_", chains))
            v <- c(x, v)
            names(v)[1] <- "community"
            return(v)
        }
        
        es <- as_data_frame(x = sg, what = "edges")
        l <- lapply(X = 1:nrow(es), es = es, chains, 
                    FUN = function(x, es, chains) {
                        v <- numeric(length = length(chains)*2+1)
                        names(v) <- c("w", paste0("w_", chains), 
                                      paste0("n_", chains))
                        for(c in chains) {
                            i <- which(es$chain[[x]]==c)
                            if(length(i)==0) {
                                v[paste0("w_", c)] <- 0
                                v[paste0("n_", c)] <- 0
                            } else {
                                v["w"] <- es$w[[x]][i]
                                v[paste0("w_", c)] <- es$weight[[x]][i]
                                v[paste0("n_", c)] <- 1
                            }
                        }
                        return(v)
                    })
        l <- data.frame(do.call(rbind, l))
        l <- c(x, colMeans(l[,which(regexpr(pattern = "w", 
                                            text = colnames(l))!=-1)]),
               colSums(l[,which(regexpr(pattern = "n", 
                                        text = colnames(l))!=-1)]))
        names(l)[1] <- "community"
        return(l)
    } 
        
    # get community statistics on edges
    es <- lapply(X = unique(V(g)$community), g = g, 
                 chains = chains, FUN = get_estats)
    es <- data.frame(do.call(rbind, es))
    
    # get community statistics on vertices
    vs <- get_vstats(vs = as_data_frame(x = g, what = "vertices"))
    
    # merge results
    o <- merge(x = vs, y = es, by = "community", all.x = TRUE)
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

check_inputs <- function(graph, 
                         algorithm, 
                         resolution,
                         weight, 
                         metric, 
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
        stop("algorithm must be louvain or leiden")
    }
    if(length(algorithm)!=1) {
        stop("algorithm must be louvain or leiden")
    }
    if(is.character(algorithm)==FALSE) {
        stop("algorithm must be character")
    }
    if(!algorithm %in% c("louvain", "leiden")) {
        stop("algorithm must be louvain or leiden")
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
    
    
    
    # check metric
    if(missing(metric)) {
        stop("metric must be average, strict or loose")
    }
    if(length(metric)!=1) {
        stop("metric must be average, strict or loose")
    }
    if(is.character(metric)==FALSE) {
        stop("metric must be character")
    }
    if(!metric %in% c("average", "strict", "loose")) {
        stop("metric must be average, strict or loose")
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
