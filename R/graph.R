get_edges <- function(clust_irr) {
  
  if(!methods::is(clust_irr, 'clust_irr')){
    base::stop("Input has to be object of class clust_irr")
  }
  
  le <- get_local_edges(clust_irr = clust_irr)
  ge <- get_global_edges(clust_irr = clust_irr)
  e <- base::rbind(le, ge)
  
  if(base::is.null(e) || base::nrow(e)==0) {
    return(NULL)
  }
  
  me <- base::vector(mode = "list", length = base::length(clust_irr$clust))
  base::names(me) <- base::names(clust_irr$clust)
  for(chain in base::names(clust_irr$clust)) {
    tmp_d <- clust_irr$inputs$s[, c(chain, "ID")]
    tmp_e <- e[e$chain == chain, ]
    
    if(base::nrow(tmp_e)!=0 & base::nrow(tmp_d)!=0) {
      base::colnames(tmp_d) <- c("CDR3", "ID")
      me[[chain]] <- base::merge(x = base::merge(
        x = tmp_e, y = tmp_d, by.x = "from_cdr3", by.y = "CDR3"),
        y = tmp_d, by.x = "to_cdr3", by.y = "CDR3")
    }
  }
  me <- base::do.call(base::rbind, me)
  if(base::is.null(me)==FALSE && base::nrow(me)!=0) {
    me$from <- me$ID.x
    me$to <- me$ID.y
    me <- me[,c("from", "to", "from_cdr3","to_cdr3", 
                "motif", "type", "chain")]
    return(me)
  }
  return(NULL)
}

get_local_edges <- function(clust_irr) {
  
  get_lp <- function(x, lp, chain) {
    if(base::sum(lp$motif == x)>1) {
      p <- base::t(utils::combn(x = lp$cdr3[lp$motif == x], m = 2))
      return(base::data.frame(from_cdr3 = p[,1], to_cdr3 = p[,2],
                              motif = x, type = "local",
                              chain = chain))
    }
    return(NULL)
  }
  
  edges_local <- vector(mode = "list", length = length(clust_irr$clust))
  names(edges_local) <- names(clust_irr$clust)
  for(chain in names(clust_irr$clust)) {
    lp <- clust_irr$clust[[chain]]$local$lp
    if(base::is.null(lp)==FALSE && base::nrow(lp)!=0) {
      edges_local[[chain]] <- base::do.call(base::rbind, base::lapply(
        X = base::unique(lp$motif), FUN = get_lp, lp = lp,
        chain = chain))
    }
  }
  edges_local <- base::do.call(base::rbind, edges_local)
  return(edges_local)
}

get_global_edges <- function(clust_irr) {
  
  get_diff_str <- function(m) {
    c1 <- strsplit(m[1], "")[[1]]
    c2 <- strsplit(m[2], "")[[1]]
    c1[c1 != c2] <- "-"
    return(paste(c1, collapse = ""))
  }
  
  get_gp <- function(gp, chain) {
    return(base::data.frame(from_cdr3 = gp[,1], to_cdr3 = gp[,2],
                            motif = apply(gp, 1, get_diff_str), 
                            type = "global",
                            chain = chain))
  }
  
  edges_global <- vector(mode = "list", length = length(clust_irr$clust))
  names(edges_global) <- names(clust_irr$clust)
  for(chain in names(clust_irr$clust)) {
    g <- clust_irr$clust[[chain]]$global
    if(base::is.null(g)==FALSE && base::nrow(g)!=0) {
      edges_global[[chain]] <- get_gp(gp = g, chain = chain)
    }
  }
  edges_global <- base::do.call(base::rbind, edges_global)
  return(edges_global)
}


get_graph <- function(clust_irr) {
    
    edges <- get_edges(clust_irr = clust_irr)
    
    ig <- igraph::graph_from_data_frame(edges, 
                                        directed = TRUE)
    return(ig)
}


plot_graph <- function(clust_irr) {
  
  ig <- get_graph(clust_irr = clust_irr)
  nodes <- igraph::as_data_frame(ig, what = "vertices")
  edges <- igraph::as_data_frame(ig, what = "edges")
  
  nodes <- configure_nodes(nodes = nodes, edges = edges)
  edges <- configure_edges(edges = edges)
  
  ledges <- base::data.frame(color = base::unique(edges$color), 
                             label = base::unique(edges$type),
                             arrows = "", width = 4)
  lenodes <- base::data.frame(label = base::unique(edges$chain), 
                              color = "", shape = "dot", size = 10)
  lenodes$color <- base::ifelse(test = lenodes$label == "CDR3b",
                                yes = "green",
                                no = "blue")
  
  return(
    visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
    visNetwork::visIgraphLayout(layout = "layout_components", 
                                randomSeed = 1234) %>%
    visNetwork::visOptions(highlightNearest = 
                             base::list(enabled = TRUE, 
                                        degree = 1,
                                        algorithm = "hierarchical"),
                           selectedBy = base::list(variable = "group", 
                                                   multiple = TRUE), 
                           manipulation = FALSE) %>%
    visNetwork::visLegend(addEdges = ledges, addNodes = lenodes, 
                          useGroups = FALSE, position = "right", 
                          width=0.15, zoom = FALSE)
    )
}


configure_nodes <- function(nodes, edges) {
  
  names(nodes) <- "id"
  nodes$color.background <- "green"
  nodes$color.border <- "black"
  nodes$color.highlight <- "red"
  nodes$size <- 20
  id_cdr3 <- base::data.frame(id = base::c(edges$from, edges$to),
                              label = base::c(edges$from_cdr3, edges$to_cdr3))
  nodes <- base::merge(nodes, base::unique(id_cdr3), by = "id")
  nodes$title <- base::paste("<p><b>", nodes$label, "</b></p>")
  nodes$group <- nodes$label
  nodes$shape <- "dot"
  nodes$shadow <- FALSE
  return(nodes)
}


configure_edges <- function(edges) {
  
  edges$length <- 15
  edges$width <- 10
  edges$color <- base::ifelse(test = (edges$type == "local"),
                              yes = "gray",
                              no = "orange")
  for(i in 1:base::nrow(edges)){
    t <- base::unique(edges$type[((edges$from == edges$from[i]) & 
                                    (edges$to == edges$to[i])) |
                                   ((edges$to == edges$from[i]) &
                                      (edges$from == edges$to[i]))])
    if(base::length(t) == 2){
      edges$color[i] <- "#9A0000"
      edges$type[i] <- "local & global"
    }
  }
  edges$arrows <- ""
  edges$dashes <- FALSE
  edges$title <- edges$motif
  edges$smooth <- FALSE
  edges$shadow <- FALSE
  return(edges)
}