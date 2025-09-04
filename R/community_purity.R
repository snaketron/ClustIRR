# Description:
# compute purity of communitie with respect to a node feature: 
# numeric (e.g., gene expression) via coefficient of variation 
# or categorical (e.g., tissue type) via gini impurity
get_community_purity <- function(graph, 
                                 node_feature) {
    
    if(missing(graph) || is.null(graph)) {
        stop("missing graph")
    }
    if(is_igraph(graph) == FALSE) {
        stop("graph must be igraph")
    }
    
    if(missing(node_feature) || is.null(node_feature)) {
        stop("missing node_feature")
    }
    if(length(node_feature) != 1) {
        stop("node_feature length != 1")
    }
    if(is.character(node_feature)==FALSE) {
        stop("node_feature must be character")
    }
    
    if(!any(vertex_attr_names(graph = graph) %in% node_feature)) {
        stop("node_feature not found among node attributes")
    }
    c <- vertex_attr(graph = graph, name = "community")
    if(is.null(c)) {
        stop("communities not found in graph")
    }
    f <- vertex_attr(graph = graph, name = node_feature)
    if(is.null(f)) {
        stop("feature not found in graph")
    }
    
    feature_type <- NA
    if(is.numeric(f)) {
        feature_type <- "num"
    }
    if(is.character(f) | is.factor(f) | is.logical(f)) {
        feature_type <- "cat"
    }
    if(is.na(feature_type)) {
        stop("node_feature must be categorical or numeric")  
    }
    
    if(feature_type == "cat") {
        return(get_gini_cat(ls = f, cs = c))
    } 
    return(get_gini_num(fs = f, cs = c))
}


get_gini_cat <- function(ls, cs) {
    
    get_gi <- function(l) {
        ls <- unique(l)
        l_len <- length(l)
        s <- 0
        for(i in seq_len(length.out = length(ls))) {
            s <- s + (sum(l == ls[i])/l_len)^2
        }
        return(s)
    }
    
    ucs <- unique(cs)
    uls <- unique(ls)
    n <- length(ucs)
    
    gi <- numeric(length = n)
    sz <- numeric(length = n)
    for(i in seq_len(length.out = n)) {
        j <- which(cs == ucs[i])
        gi[i] <- 1-get_gi(l = ls[j])
        sz[i] <- length(j)
    }
    
    gi <- data.frame(GI = gi, n = sz)
    gi$community <- ucs
    
    return(gi)
}


get_gini_num <- function(fs, cs) {
    
    get_gi <- function(f) {
        o <- numeric(length = 4)
        o[1] <- mean(f)
        o[2] <- sd(f)
        o[3] <- o[2]/o[1]
        o[4] <- length(f)
        return(o)
    }
    
    ucs <- unique(cs)
    n <- length(ucs)
    gi <- matrix(data = 0, nrow = n, ncol = 4)
    for(i in seq_len(length.out = n)) {
        j <- which(cs == ucs[i])
        gi[i,] <- get_gi(f = fs[j])
    }
    
    gi <- data.frame(gi)
    colnames(gi) <- c("mean", "sd", "cv", "n")
    gi$community <- ucs
    
    return(gi)
}
