# Description:
# Compute the statistics of a given community feature
get_community_feature_stats <- function(node_summary, 
                                        node_feature,
                                        community_id) {
    
    if(missing(node_summary) || is.null(node_summary)) {
        stop("missing node_summary")
    }
    if(is.data.frame(node_summary) == FALSE) {
        stop("node_summary must be a data.frame")
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
    
    if(missing(community_id) || is.null(community_id)) {
        stop("missing community_id")
    }
    if(length(community_id) != 1) {
        stop("community_id length != 1")
    }
    if(is.numeric(community_id)==FALSE) {
        stop("community_id must be character")
    }
    
    
    if(!any(colnames(node_summary) %in% node_feature)) {
        stop("node_feature not found among node attributes")
    }
    if(!any(colnames(node_summary) %in% "community")) {
        stop("community not found in node_summary")
    }
    c <- node_summary[, "community"]
    f <- node_summary[, node_feature]
    
    i <- which(c == community_id)
    if(length(i)==0) {
        stop("community not found")
    }
    c <- c[i]
    f <- f[i]
    
    feature_type <- NA
    if(is.numeric(f)) {
        feature_type <- "num"
        return(data.frame(community = community_id, 
                          feature = node_feature,
                          feature_mean = mean(f), 
                          feature_median = median(f), 
                          feature_sum = sum(f),
                          feature_type = feature_type,
                          n = length(f)))
    }
    if(is.character(f) | is.factor(f) | is.logical(f)) {
        feature_type <- "cat"
        f <- data.frame(table(f))
        colnames(f) <- c("feature", "feature_count")
        f$feature_prop <- f$feature_count/sum(f$feature_count)
        f$community <- community_id
        f$feature_type <- feature_type
        f$n <- length(c)
        f <- f[,c("community", "feature", "feature_count", 
                  "feature_prop", "feature_type", "n")]
        f <- f[order(f$feature_prop, decreasing = TRUE),]
        return(f)
    }
}




# Description:
# compute purity of communities with respect to a node feature: 
# numeric (e.g., gene expression) via coefficient of variation 
# or categorical (e.g., tissue type) via gini impurity
get_community_feature_purity <- function(node_summary, 
                                         node_feature) {
    
    if(missing(node_summary) || is.null(node_summary)) {
        stop("missing node_summary")
    }
    if(is.data.frame(node_summary) == FALSE) {
        stop("node_summary must be a data.frame")
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
    
    
    if(!any(colnames(node_summary) %in% node_feature)) {
        stop("node_feature not found among node attributes")
    }
    if(!any(colnames(node_summary) %in% "community")) {
        stop("community not found in node_summary")
    }
    c <- node_summary[, "community"]
    f <- node_summary[, node_feature]
    
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
    # Vectorized GI function
    get_gi <- function(l) {
        p <- tabulate(match(l, unique(l))) / length(l)
        return(sum(p^2))
    }
    
    get_h <- function(l) {
        p <- tabulate(match(l, unique(l))) / length(l)
        return(-sum(p*log(p)))
    }
    
    
    # split ls by cs
    groups <- split(ls, cs)
    
    # Compute GI and size per community
    return(data.frame(community = as.numeric(names(groups)),
                      GI = 1 - vapply(groups, get_gi, numeric(1)),
                      H = vapply(groups, get_h, numeric(1)),
                      n = lengths(groups),
                      row.names = NULL))
}

get_gini_num <- function(fs, cs) {
    
    groups <- split(fs, cs)
    
    gi <- data.frame(community = as.numeric(names(groups)),
                     mean = vapply(groups, mean, numeric(1)),
                     sd   = vapply(groups, sd,   numeric(1)),
                     n    = lengths(groups),
                     row.names = NULL)
    
    gi$cv <- gi$sd/gi$mean
    return(gi)
}
