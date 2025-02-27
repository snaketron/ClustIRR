detect_communities <- function(graph, 
                               weight = "nweight",
                               algorithm = "leiden", 
                               resolution = 1,
                               iterations = 100,
                               chains) {
    
    check_inputs(graph = graph,
                 algorithm = algorithm, 
                 resolution = resolution,
                 iterations = iterations,
                 weight = weight,
                 chains = chains)
    
    message("[1/5] formatting graph...")
    cg <- get_formatted_graph(graph = graph, weight = weight, chains = chains) 
    
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
    vs <- as_data_frame(x = cg$graph, what = "vertices")
    
    # save configs
    config <- list(input_g = graph, 
                   algorithm = algorithm, 
                   resolution = resolution,
                   weight = weight,
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
    
    E(graph)$nweight <- vapply(X = E(graph)$nweight, FUN.VALUE = numeric(1),
                               FUN = function(x) {return(sum(x)/2)})
    E(graph)$ncweight <- vapply(X = E(graph)$ncweight, FUN.VALUE = numeric(1),
                                FUN = function(x) {return(sum(x)/2)})
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

get_community_summary <- function(g, 
                                  chains) {
    
    get_vstats <- function(vs, wide) {
        
        vs$cells <- vs$clone_size
        vs$clones <- 1
        
        if(wide) {
            # number of cells
            vcells <- aggregate(cells~community+sample, data = vs, FUN = sum)
            vcells <- acast(data = vcells, formula = community~sample, 
                            value.var = "cells", fill = 0)
            vcells <- data.frame(vcells)
            vcells$n <- apply(X = vcells, MARGIN = 1, FUN = sum)
            colnames(vcells) <- paste0("cells_", colnames(vcells))
            vcells$community <- rownames(vcells)
            vcells$community <- as.numeric(as.character(vcells$community))
            vcells <- vcells[order(vcells$community, decreasing = FALSE), ]
            
            # number of clones
            vclones <- aggregate(clones~community+sample, data = vs, FUN = sum)
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
            vcells <- aggregate(cells~community+sample, data = vs, 
                                FUN = sum, drop = FALSE)
            
            # number of clones
            vclones <- aggregate(clones~community+sample, data = vs, 
                                 FUN = sum, drop = FALSE)
            
            # merge clones and cells
            vstats <- merge(x = vclones, y = vcells, 
                            by = c("community", "sample"))
            vstats$cells[is.na(vstats$cells)] <- 0
            vstats$clones[is.na(vstats$clones)] <- 0
        }
        
        return(vstats)
    }
    
    get_estats <- function(x, g, chains) {
        sg <- subgraph(graph = g, vids = which(V(g)$community==x))
        
        if(length(sg)==1) {
            v <- numeric(length = length(chains)*3+1)
            names(v) <- c("w", 
                          paste0("ncweight_", chains),
                          paste0("nweight_", chains),
                          paste0("n_", chains))
            v <- c(x, v)
            names(v)[1] <- "community"
            return(v)
        }
        
        es <- as_data_frame(x = sg, what = "edges")
        l <- lapply(X = 1:nrow(es), es = es, chains, 
                    FUN = function(x, es, chains) {
                        v <- numeric(length = length(chains)*3+1)
                        names(v) <- c("w", 
                                      paste0("ncweight_", chains),
                                      paste0("nweight_", chains),
                                      paste0("n_", chains))
                        for(c in chains) {
                            i <- which(es$chain[[x]]==c)
                            if(length(i)==0) {
                                v[paste0("ncweight_", c)] <- 0
                                v[paste0("nweight_", c)] <- 0
                                v[paste0("n_", c)] <- 0
                            } else {
                                v["w"] <- es$w[[x]]#[i]
                                v[paste0("ncweight_", c)] <- es$ncweight[[x]][i]
                                v[paste0("nweight_", c)] <- es$nweight[[x]][i]
                                v[paste0("n_", c)] <- 1
                            }
                        }
                        return(v)
                    })
        l <- data.frame(do.call(rbind, l))
        l <- c(x, colMeans(l[,which(regexpr(pattern = "w", 
                                            text = colnames(l))!=-1)],
                           na.rm = TRUE),
               colSums(l[,which(regexpr(pattern = "n\\_", 
                                        text = colnames(l))!=-1), drop=FALSE]))
        
        names(l)[1] <- "community"
        return(l)
    } 
    
    # get community statistics on edges
    es <- lapply(X = unique(V(g)$community), g = g, 
                 chains = chains, FUN = get_estats)
    es <- data.frame(do.call(rbind, es))
    
    # get community statistics on vertices (wide and tall format)
    vs_wide <- get_vstats(vs = as_data_frame(x = g, what = "vertices"), 
                          wide = TRUE)
    vs_tall <- get_vstats(vs = as_data_frame(x = g, what = "vertices"), 
                          wide = FALSE)
    
    # merge results
    cs_wide <- merge(x = vs_wide, y = es, by = "community", all.x = TRUE)
    cs_wide <- cs_wide[order(cs_wide$community, decreasing = FALSE), ]
    
    cs_tall <- merge(x = vs_tall, y = es, by = "community", all.x = TRUE)
    cs_tall <- cs_tall[order(vs_tall$community, decreasing = FALSE), ]
    
    
    return(list(wide = cs_wide, 
                tall = cs_tall))
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

# Description:
# find components, cliques, subgraphs in a community
decode_communities <- function(community_id,
                               graph, 
                               edge_filter,
                               node_filter) {
    
    apply_op <- function(vec, op, val) {
        ops <- list("==" = `==`,
                    "!=" = `!=`,
                    "<"  = `<`,
                    ">"  = `>`,
                    "<=" = `<=`,
                    ">=" = `>=`)
        
        if (!op %in% names(ops)) {
            stop("Invalid operator, choose: '==', '!=', '<', '>', '<=', '>='")
        }
        
        return(ops[[op]](vec, val))
    }
    
    if(any(vertex_attr_names(graph)=="community")==FALSE) {
        stop("no community ID as node attribute")
    } else {
        graph <- subgraph(graph = graph, 
                          vids = V(graph)$community == community_id)
    }
    if(length(graph)==1 | length(E(graph)) == 0) {
        warning("community has only one vertex")
    }

    # consider edges
    if(nrow(edge_filter)!=0) {
        # this is where the edge filter results will be kept
        etm <- matrix(data = 0, 
                      nrow = nrow(edge_filter), 
                      ncol = length(E(graph)))
        for(i in seq_len(nrow(edge_filter))) {
            a_name <- edge_filter$name[i]
            a_value <- edge_filter$value[i]
            a_operation <- edge_filter$operation[i]
            
            j <- which(edge_attr_names(graph) == a_name)
            if(length(j) != 0) {
                v <- edge_attr(graph = graph, name = a_name)
                etm[i,] <- apply_op(vec = v, val = a_value, op = a_operation)
            }
        }
        etm <- apply(X = etm, MARGIN = 2, FUN = prod)
        i <- which(etm == FALSE)
        if(length(i) != 0) {
            graph <- delete_edges(graph = graph, edges = i)
        }
    }
    
    # now partition based on node-attributes
    vs <- as_data_frame(x = graph, what = "vertices")
    V(graph)$key <- apply(X = vs[, node_filter$name, drop=FALSE], 
                          MARGIN = 1, FUN = paste, collapse = '|')
    sgs <- lapply(
        X = unique(V(graph)$key), g = graph, 
        FUN = function(x, g) {
            # get a subgraph with shared node attributes
            sg <- subgraph(graph = g, vids = which(V(g)$key == x))
            
            # find connected components
            V(sg)$components <- components(graph = sg)$membership
            V(sg)$component_id <- paste0(V(sg)$key, '|', 
                                         V(sg)$components)
            V(sg)$component_id <- as.numeric(as.factor(V(sg)$component_id))
            return(disjoint_union(lapply(
                X = unique(V(sg)$component_id), g = sg, 
                FUN = function(x, g) {
                    vids <- which(V(g)$component_id == x)
                    return(subgraph(graph = g, vids = vids))
                })))
        })
    sgs <- disjoint_union(sgs)
    
    return(sgs)
}

