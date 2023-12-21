
config_edges <- function(es) {
  
  # edge color by type
  get_type_edge <- function(x) {
    x <- sort(unique(x))
    if(length(x)==2) {
      return("black")
    }
    if(x == "local") {
      return("purple")
    }
    if(x == "global") {
      return("#27AE60")
    }
    return("black")
  }
  
  # edge shape by chain
  get_chain_edge <- function(x) {
    x <- sort(unique(x))
    if(length(x)==2) {
      return("solid")
    }
    if(x %in% c("CDR3b", "CDR3g", "CDR3h")) {
      return("dashed")
    }
    if(x %in% c("CDR3a", "CDR3d", "CDR3l")) {
      return("dotted")
    }
  }
  
  if(nrow(es) == 0) {
    return(NULL)
  }
  
  attr_type <- aggregate(type~from+to, data = es, FUN = get_type_edge)
  attr_type$type_color <- attr_type$type
  attr_type$type <- NULL
  
  attr_chain <- aggregate(chain~from+to, data = es, FUN = get_chain_edge)
  attr_chain$chain_shape <- attr_chain$chain
  attr_chain$chain <- NULL
  
  return(merge(x = attr_chain, y = attr_type, by = c("from", "to")))
}

config_edges_plot <- function(g, is_jg) {
  n_e <- length(E(g))
  if(n_e != 0) {
    E(g)$color <- E(g)$type_color
    E(g)$lty <- E(g)$chain_shape
  }
  return(g)
}

config_vertices_plot <- function(g, is_jg) {
  # default features
  V(g)$size <- 1.5+log2(V(g)$clone_size)
  
  if(is_jg==TRUE) {
    V(g)$color_num <- as.numeric(as.factor(V(g)$sample))
    max_n <- max(V(g)$color_num)
    V(g)$color <- hcl.colors(n=max_n, palette = "Roma")[V(g)$color_num]
    V(g)$frame.color <- V(g)$color
  } 
  else {
    V(g)$color <- "black"
    V(g)$frame.color <- "black"
  }
  
  return(g)
}

get_intergraph_edges <- function(igs, global_max_dist, chains) {
  
  get_igg <- function(x, i, igs, global_max_dist, chain) {
    
    get_hdist <- function(x, 
                          id_x, 
                          id_y, 
                          seq_x, 
                          seq_y,
                          sample_x,
                          sample_y,
                          global_max_dist) {
      d <- stringdist(a = seq_x[x], b = seq_y, method = "hamming")
      js <- which(d <= global_max_dist)
      if(length(js) == 0) {
        return(NULL)
      }
      return(data.frame(from = id_x[x], 
                        to = id_y[js], 
                        sample_x = sample_x, 
                        sample_y = sample_y))
    }
    
    get_hamming_dist <- function(x, 
                                 id_x, 
                                 id_y, 
                                 seq_x, 
                                 seq_y, 
                                 len_x, 
                                 len_y, 
                                 sample_x,
                                 sample_y,
                                 global_max_dist) {
      
      is_x <- which(len_x == x)
      is_y <- which(len_y == x)
      
      if(length(is_x)==0|length(is_y)==0) {
        return(NULL)
      }
      
      hd <- lapply(X = seq_along(is_x),
                   FUN = get_hdist,
                   id_x = id_x[is_x], 
                   id_y = id_y[is_y], 
                   seq_x = seq_x[is_x], 
                   seq_y = seq_y[is_y],
                   sample_x = sample_x,
                   sample_y = sample_y,
                   global_max_dist = global_max_dist)
      hd <- do.call(rbind, hd)
      return(hd)
    }
    
    
    s1_name <- names(igs)[i]
    s2_name <- names(igs)[x]
    
    s1 <- igs[[i]]$clones
    s2 <- igs[[x]]$clones
    
    seq_x <- s1[,chain]
    seq_y <- s2[,chain]
    id_x <- s1[,"name"] 
    id_y <- s2[,"name"]
    len_x <- nchar(seq_x)
    len_y <- nchar(seq_y)
    # browser()
    hd <- lapply(X = unique(c(len_x, len_y)),
                 FUN = get_hamming_dist,
                 id_x = id_x,
                 id_y = id_y,
                 seq_x = seq_x,
                 seq_y = seq_y,
                 len_x = len_x,
                 len_y = len_y,
                 sample_x = s1_name,
                 sample_y = s2_name,
                 global_max_dist = global_max_dist)
    # browser()
    hd <- do.call(rbind, hd)
    if(is.null(hd)==FALSE && nrow(hd)!=0) {
      hd$chain <- chain
      hd$sample <- paste0(hd$sample_x, "|", hd$sample_y)
      hd$sample_x <- NULL
      hd$sample_y <- NULL
      hd$type <- "inter-sample"
      return(hd)
    }
    return(NULL)
  }
  
  # find global similarities between pairs of clone tables
  ige <- vector(mode = "list", length = length(chains)*(length(igs)-1))
  
  count <- 1
  for(i in 1:(length(igs)-1)) {
    message("merging clust_irr index: ", i, "\n")
    for(chain in chains) {
      ige[[count]] <- do.call(rbind,
                              lapply(X = (i+1):length(igs), 
                                     i = i,
                                     FUN = get_igg,
                                     igs = igs,
                                     chain = chain,
                                     global_max_dist = global_max_dist))
      count <- count + 1
    }
  }
  ige <- do.call(rbind, ige)
  return(ige)
}
