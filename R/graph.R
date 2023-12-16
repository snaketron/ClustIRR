

get_graph <- function(x, clust_irrs) {
  
  get_local_edges <- function(clust_irr, cs) {
    
    get_motif_to_id <- function(x, cs, lp, chain) {
      y <- lp$cdr3[lp$motif == x]
      if(length(y)>0) {
        return(unique(cs$clone_id[cs[,chain] %in% y]))
      }
      return(NA)
    }
    
    el <- vector(mode="list", length = length(get_clustirr_clust(clust_irr)))
    names(el) <- names(get_clustirr_clust(clust_irr))
    for(chain in names(get_clustirr_clust(clust_irr))) {
      lp <- get_clustirr_clust(clust_irr)[[chain]]$local$lp
      cs <- cs[cs[, chain] %in% lp$cdr3,]
      
      if(is.null(lp) == FALSE && nrow(lp) != 0  && nrow(s) != 0) {
        
        lm <- lapply(X = unique(lp$motif), 
                     FUN = get_motif_to_id, 
                     lp = lp, 
                     cs = cs, 
                     chain = chain)
        names(lm) <-  unique(lp$motif)
        
        if(is.null(lm)==FALSE) {
          el[[chain]] <- lm
        }
      }
    }
    return(el)
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
  
  build_graph <- function(le, ge, cs, sample_id) {
    
    add_edges_ig <- function(le, ig, sample_id) {
      
      add_motif_edges <- function(x, ig, sample_id) {
        xp <- utils::combn(x = paste0(sample_id, '|', x), m = 2)
        xp <- as.vector(xp)
        return(igraph::add.edges(graph = ig, edges = xp))
      }
      
      for(i in 1:length(le)) {
        # if only one element in le then do not add edge, else add
        if(length(le[[i]])!=1) {
          ig <- add_motif_edges(x = le[[i]], ig = ig, sample_id = sample_id)
        }
      }
      return(ig)
    }
    
    ig <- graph_from_data_frame(d = data.frame(from = cs$name[1], 
                                               to = cs$name[1]),
                                directed = FALSE,
                                vertices = cs)
    ig <- delete_edges(ig, edges = 1)
    
    # browser()
    
    # add local edges
    if(length(le)!=0) {
      for(chain in names(le)) {
        if(length(le[[chain]])!=0) {
          ig <- add_edges_ig(le = le[[chain]], ig = ig, sample_id = sample_id)
        }
      }
    }
    
    return(ig)
  }
  
  clust_irr <- clust_irrs[[x]]
  
  check_clustirr(clust_irr = clust_irr)
  
  # cells
  s <- get_clustirr_inputs(clust_irr)$s
  
  # clones
  cs <- s
  cs$clone_size <- 1
  cs$id <- NULL
  cs <- aggregate(clone_size~., data = cs, FUN = sum)
  cs$clone_id <- seq_len(nrow(cs))
  cs$sample <- x
  cs$name <- paste0(x, '|', cs$clone_id)
  cs <- cs[, rev(colnames(cs))]
  
  # get local and global edges between clones
  le <- get_local_edges(clust_irr = clust_irr, cs = cs)
  ge <- get_global_edges(clust_irr = clust_irr)
  
  # TODO: what does it mean is.null(ge) or legnth(le)==0
  
  # build graph with only vertices
  if(length(le)==0 & is.null(ge)) {
    ig <- graph_from_data_frame(d = data.frame(from = cs$name[1], 
                                               to = cs$name[1]),
                                directed = FALSE,
                                vertices = cs)
    ig <- delete_edges(ig, edges = 1)
    
    return(list(graph = ig, clones = cs))
  }
  
  # build graph
  ig <- build_graph(le = le, ge = ge, cs = cs, sample_id = x)
  
  return(list(graph = ig, clones = cs))
}


get_joint_graph <- function(clust_irrs) {
  
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
  
  get_clust_irrs_names <- function(clust_irrs) {
    clust_irrs_names <- names(clust_irrs)
    if(is.null(clust_irrs_names)) {
      names(clust_irrs) <- paste0("s", 1:length(clust_irrs))
    }
    return(clust_irrs)
  }
  
  get_v_e <- function(x, what) {
    return(get.data.frame(x$graph, what = what))
  }
  
  # check input
  check_input(clust_irrs = clust_irrs)
  
  # get clust_irrs names
  clust_irrs <- get_clust_irrs_names(clust_irrs = clust_irrs)
  
  # get igs
  igs <- lapply(X = names(clust_irrs), clust_irrs=clust_irrs, FUN = get_graph)
  
  # rename igs
  names(igs) <- names(clust_irrs)
  
  # get chains
  chains <- get_chains(x = colnames(get_clustirr_inputs(clust_irrs[[1]])$s))
  
  # get global_max_dist
  gmd <- get_clustirr_inputs(clust_irrs[[1]])$control$global_max_dist
  
  # get intergraph edges (global)
  ige <- get_intergraph_edges(igs=igs, global_max_dist=gmd, chains=chains)
  
  # get the vertices/edges of the graph
  df_v <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "vertices"))
  df_e <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "edges"))
  df_e$type <- "intra-sample"
  df_e <- rbind(df_e, ige$ige[, c("from", "to", "type")])
  # df_e <- config_edges(es = df_e)
  
  # build joint graph
  g <- graph_from_data_frame(df_e, directed=FALSE, vertices=df_v)
  
  # make graph look visually better
  # g <- config_vertices_plot(g = g, is_jg = TRUE)
  # g <- config_edges_plot(g = g, is_jg = TRUE)
  
  return(list(graph = g, clones = df_v))
}


plot_graph <- function(clust_irr, 
                       as_visnet = FALSE) {
  
  check_clustirr(clust_irr = clust_irr)
  
  clust_irr <- list(clust_irr)
  names(clust_irr) <- "S"
  ig <- get_graph(x = "S", clust_irrs = clust_irr)
  clones <- ig$clones
  if(is.null(ig$graph)) {
    warning("No graph to plot \n")
    return(list(graph = NA, clones = clones))
  }
  ig <- ig$graph
  ig <- config_vertices_plot(g = ig, is_jg = FALSE)
  # plot
  if(as_visnet == FALSE) {
    plot(ig, vertex.label = NA)
  }
  if(as_visnet == TRUE) {
    V(ig)$size <- V(ig)$size*10
    if(length(E(ig))==0) {
      # apparently if no edges, visnetwork can't plot
      ig <- add_edges(graph = ig, edges = c("S1","S1"))
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
  # jg$graph <- config_edges_plot(g = jg$graph, is_jg = TRUE)
  
  # plot
  if(as_visnet == FALSE) {
    plot(jg$graph, vertex.label = NA)
  }
  if(as_visnet == TRUE) {
    browser()
    V(jg$graph)$size <- V(jg$graph)$size*5
    visIgraph(igraph = jg$graph,
              idToLabel = TRUE,
              layout = "layout_components",
              randomSeed = 1234,
              physics = FALSE,
              smooth = FALSE,
              type = "square")
  }
}

