
# Description:
# find components, cliques, subgraphs in a community
decode_community <- function(community_id,
                             graph, 
                             edge_filter = NULL,
                             node_filter = NULL) {
    
    check_graph(graph = graph)
    check_community_id(graph = graph, community_id = community_id)
    check_edge_filter(edge_filter = edge_filter, graph = graph)
    check_node_filter(node_filter = node_filter, graph = graph)
    check_node_and_edges_filter(node_filter = node_filter,
                                edge_filter = edge_filter)
    
    # subset graph
    graph <- subgraph(graph = graph, vids = V(graph)$community == community_id)
    
    graph <- filter_edges(g = graph, edge_filter = edge_filter)
    graph <- filter_nodes(g = graph, node_filter = node_filter)
    return(graph)
}

decode_all_communities <- function(graph, 
                                   edge_filter = NULL, 
                                   node_filter = NULL) {
    check_graph(graph = graph)
    check_edge_filter(edge_filter = edge_filter, graph = graph)
    check_node_filter(node_filter = node_filter, graph = graph)
    check_node_and_edges_filter(node_filter = node_filter,
                                edge_filter = edge_filter)
    
    cs <- unique(V(graph)$community)
    o <- lapply(X = cs, 
                FUN = decode_community, 
                graph = graph, 
                edge_filter = edge_filter,
                node_filter = node_filter)
    names(o) <- cs
    return(o)
}

get_component_stats <- function(x) {
    
    # what is he
    community_id <- V(x)$community[1]
    
    # number of nodes per component
    count_n <- table(V(x)$component_id)
    
    # components has more than one node
    i <- which(count_n>1)
    stats_components <- c()
    if(length(i)>0) {
        stats_components <- lapply(X = names(i), g = x, FUN = function(x, g) {
            o <- subgraph(g, vids = which(V(g)$component_id==x))
            
            ncweight <- E(o)$ncweight
            if(length(ncweight)==0) {
                ncweight <- NA
            }
            nweight <- E(o)$nweight
            if(length(nweight)==0) {
                nweight <- NA
            }
            
            return(data.frame(component_id = x,
                              community = V(o)$community[1],
                              mean_ncweight = mean(ncweight),
                              mean_nweight = mean(nweight),
                              n_nodes = length(o),
                              n_edges = length(E(o)),
                              n_clique_edges = length(o)*(length(o)-1)/2,
                              diameter = diameter(o)))
        })
        stats_components <- do.call(rbind, stats_components)
    }
    
    i <- which(count_n==1)
    stats_singletons <- c()
    if(length(i)>0) {
        stats_singletons <- data.frame(component_id = names(i),
                                       community = community_id,
                                       mean_ncweight = NA,
                                       mean_nweight = NA,
                                       n_nodes = 1,
                                       n_edges = 0,
                                       n_clique_edges = 0,
                                       diameter = 0)
    }
    
    return(rbind(stats_components, stats_singletons))
}

check_graph <- function(graph) {
    if(missing(graph)) {
        stop("missing graph")
    }
    if(is_igraph(graph)==FALSE) {
        stop("graph is not igraph")
    }
    if(any(vertex_attr_names(graph)=="community")==FALSE) {
        stop("no community as node attribute in graph")
    }
    if(is.null(vertex_attr(graph = graph, name = "community"))) {
        stop("no community as node attribute in graph")
    }
}

check_community_id <- function(community_id, 
                               graph) {
    if(missing(community_id)) {
        stop("no community_id")
    }
    if(length(community_id)!=1) {
        stop("community_id must have length one")
    }
    if(is.numeric(community_id)==FALSE) {
        stop("community_id is not numeric")
    }
    if(community_id<0) {
        stop("community_id is negative")
    }
    if(is.infinite(community_id)) {
        stop("community_id is infinite")
    }
    j <- which(vertex_attr(graph = graph, name = "community")==community_id)
    if(length(j)==0) {
        stop("input community_id not found in graph")
    }
}

check_edge_filter <- function(edge_filter,
                              graph) {
    if(missing(edge_filter)|is.null(edge_filter)) {
        warning("no edge_filter")
    } 
    else {
        if(is.data.frame(edge_filter)==FALSE) {
            stop("edge_filter provided but not data frame")  
        }
        if(nrow(edge_filter)==0) {
            stop("edge_filter provided but no rows")  
        }
        if(!all(colnames(edge_filter) %in% c("name", "value", "operation"))) {
            stop("edge_filter columns must be: name, value and operation")
        }
        if(!all(edge_filter$operation %in% c("<",">",">=","<=","==","!="))) {
            stop("allowed edge_filter operation: <, >, >=, <=, ==, !=")
        }
        if(any(is.character(edge_filter$name)==FALSE)) {
            stop("edge_filter name must be character")
        }
        if(any(is.character(edge_filter$name)==FALSE)) {
            stop("edge_filter name must be character")
        }
        if(any(is.numeric(edge_filter$value)==FALSE)) {
            stop("edge_filter value must be numeric")
        }
        if(any(is.infinite(edge_filter$value)==TRUE)) {
            stop("edge_filter value must be numeric")
        }
        
        # graph
        if(!all(edge_filter$name %in% edge_attr_names(graph = graph))) {
            stop("some edge_filter names not found as edge attributes in graph")
        }
    }
}

check_node_filter <- function(node_filter, 
                              graph) {
    if(missing(node_filter)|is.null(node_filter)) {
        warning("no node_filter")
    } 
    else {
        if(is.data.frame(node_filter)==FALSE) {
            stop("node_filter provided but not data frame")  
        }
        if(nrow(node_filter)==0) {
            stop("node_filter provided but no rows")  
        }
        if(!all(colnames(node_filter) %in% c("name"))) {
            stop("node_filter columns must be: name")
        }
        if(any(is.character(node_filter$name)==FALSE)) {
            stop("node_filter name must be character")
        }
        if(any(is.character(node_filter$name)==FALSE)) {
            stop("node_filter name must be character")
        }
        
        # graph
        if(!all(node_filter$name %in% vertex_attr_names(graph = graph))) {
            stop("some node_filter names not found as node attributes in graph")
        }
    }
}

check_node_and_edges_filter <- function(node_filter,
                                        edge_filter) {
    
    if(missing(node_filter) & missing(edge_filter)) {
        stop("no nodes AND edge filter provided")
    }
}

apply_filter_op <- function(vec, op, val) {
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

filter_edges <- function(g, 
                         edge_filter) {
    if(is.null(edge_filter)==FALSE) {
        etm <- matrix(data = 0, nrow = nrow(edge_filter), ncol = length(E(g)))
        for(i in seq_len(nrow(edge_filter))) {
            a_n <- edge_filter$name[i]
            a_v <- edge_filter$value[i]
            a_o <- edge_filter$operation[i]
            
            j <- which(edge_attr_names(g) == a_n)
            if(length(j) != 0) {
                v <- edge_attr(graph = g, name = a_n)
                etm[i,] <- apply_filter_op(vec = v, val = a_v, op = a_o)
            }
        }
        etm <- apply(X = etm, MARGIN = 2, FUN = prod)
        i <- which(etm == FALSE)
        if(length(i) != 0) {
            g <- delete_edges(graph = g, edges = i)
        }
    }
    return(g)
}

filter_nodes <- function(g, 
                         node_filter) {
    if(is.null(node_filter)==FALSE) {
        vs <- as_data_frame(x = g, what = "vertices")
        V(g)$key <- apply(X = vs[, node_filter$name, drop=FALSE], 
                          MARGIN = 1, FUN = paste, collapse = '|')
        
        sgs <- lapply(X = unique(V(g)$key), g = g, 
                      FUN = function(x, g) {
                          sg <- subgraph(graph = g, vids = which(V(g)$key == x))
                          
                          V(sg)$components <- components(graph = sg)$membership
                          V(sg)$component_id <- paste0(V(sg)$key, '|', 
                                                       V(sg)$components)
                          return(disjoint_union(lapply(
                              X = unique(V(sg)$component_id), g = sg, 
                              FUN = function(x, g) {
                                  vids <- which(V(g)$component_id == x)
                                  return(subgraph(graph = g, vids = vids))
                              })))
                      })
        sgs <- disjoint_union(sgs)
        V(sgs)$component_id <- as.numeric(as.factor(V(sgs)$component_id))
        
        sgs_stat <- get_component_stats(x = sgs)
        vs <- as_data_frame(x = sgs, what = "vertices")
        
        return(list(community_graph = sgs, 
                    component_stats = sgs_stat,
                    node_summary = vs))
    } 
    else {
        sgs <- lapply(
            X = unique(V(g)$community), g = g, 
            FUN = function(x, g) {
                # get a subgraph with shared node attributes
                sg <- subgraph(graph = g, vids = which(V(g)$community == x))
                
                # find connected components
                V(sg)$components <- components(graph = sg)$membership
                V(sg)$component_id <- paste0(V(sg)$community, '|', 
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
        
        sgs_stat <- get_component_stats(x = sgs)
        
        vs <- as_data_frame(x = sgs, what = "vertices")
        
        return(list(community_graph = sgs, 
                    component_stats = sgs_stat,
                    node_summary = vs))
    }
}
