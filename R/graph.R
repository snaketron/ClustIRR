get_edges <- function(clust_irr) {
  
  check_clustirr(clust_irr = clust_irr)
  
  le <- get_local_edges(clust_irr = clust_irr)
  ge <- get_global_edges(clust_irr = clust_irr)
  e <- base::rbind(le, ge)
  
  if(base::is.null(e) || base::nrow(e)==0) {
    base::warning("No local or global edges found \n")
    return(NULL)
  }
  
  me <- base::vector(mode = "list", length = base::length(clust_irr$clust))
  base::names(me) <- base::names(clust_irr$clust)
  for(chain in base::names(clust_irr$clust)) {
    tmp_d <- clust_irr$inputs$s[, c(chain, "id")]
    # remove clonal expanded cdr3s for version 1 and 2
    if(clust_irr$inputs$version != 3) {
      tmp_d <- tmp_d[!base::duplicated(tmp_d[[chain]]), ]
    }
    tmp_e <- e[e$chain == chain, ]
    
    if(base::nrow(tmp_e)!=0 & base::nrow(tmp_d)!=0) {
      base::colnames(tmp_d) <- c("CDR3", "id")
      me[[chain]] <- base::merge(x = base::merge(
        x = tmp_e, y = tmp_d, by.x = "from_cdr3", by.y = "CDR3"),
        y = tmp_d, by.x = "to_cdr3", by.y = "CDR3")
    }
  }
  me <- base::do.call(base::rbind, me)
  if(base::is.null(me)==FALSE && base::nrow(me)!=0) {
    me$from <- me$id.x
    me$to <- me$id.y
    me <- me[,c("from", "to", "from_cdr3","to_cdr3", 
                "motif", "type", "chain")]
    # remove self-reference connections to prevent circles (only v3)
    me <- me[!(me$from == me$to),]
    return(me)
  }
  base::warning("No local or global edges found \n")
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
  
  return(base::unique(edges_local))
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
  
  check_clustirr(clust_irr = clust_irr)
  
  edges <- get_edges(clust_irr = clust_irr)

  if(base::is.null(edges)) {
    base::warning("No local or global edges to build igraph from \n")
    return(NULL)
  }
  
  ig <- igraph::graph_from_data_frame(edges, 
                                      directed = TRUE)
  return(ig)
}


plot_graph <- function(clust_irr) {
  
  check_clustirr(clust_irr = clust_irr)
  
  ig <- get_graph(clust_irr = clust_irr)
  
  if(base::is.null(ig)) {
    base::warning("No graph to plot \n")
    return(NULL)
  }
  
  edges <- igraph::as_data_frame(ig, what = "edges")
  nodes <- igraph::as_data_frame(ig, what = "vertices")
  chains <- base::unique(edges$chain)
  
  edges <- configure_edges(edges = edges)
  nodes <- configure_nodes(nodes = nodes, edges = edges, chains = chains,
                           s = clust_irr$inputs$s)
  
  ledges <- base::data.frame(color = base::unique(edges$color), 
                             label = base::unique(edges$type),
                             arrows = "", width = 4)
  
  
  if(base::length(chains)>1) { chains <- base::append(chains, "CDR3a+b") }
  lnodes <- base::data.frame(label = chains, 
                              color = "", shape = "dot", size = 8)
  lnodes$color <- base::ifelse(test = lnodes$label == "CDR3b",
                                yes = "blue",
                                no = "yellow")
  lnodes$color <- base::ifelse(test = lnodes$label == "CDR3a+b",
                               yes = "green",
                               no = lnodes$color)
  lnodes <- base::rbind(lnodes, 
                        base::data.frame(label = c("Expanded", "Singleton"), 
                                         color = "black", 
                                         shape = c("dot", "diamond"), 
                                         size = 8))
  
  return(configure_network(nodes = nodes, edges = edges,
                           lnodes = lnodes, ledges = ledges))
}


configure_nodes <- function(nodes, edges, chains, s) {
  
  base::names(nodes) <- "id"
  nodes <- base::merge(nodes, s, by = "id")
  min_size <- 20
  clone_boost <- 10
  if(base::length(chains)>1){
    nodes$label <- base::paste(nodes$CDR3b, "(b) -", nodes$CDR3a, "(a)")
    nodes$clone_count <- base::apply(X = nodes, MARGIN = 1, function(x) 
      sum(s$CDR3a == x['CDR3a'] & s$CDR3b == x['CDR3b']))
  } else {
    if (chains == "CDR3a"){
      nodes$label <- nodes$CDR3a
      nodes$clone_count <- apply(X = nodes, MARGIN = 1, function(x) 
        sum(s$CDR3a == x['CDR3a']))
    }
    if (chains == "CDR3b"){
      nodes$label <- nodes$CDR3b
      nodes$clone_count <- apply(X = nodes, MARGIN = 1, function(x) 
        sum(s$CDR3b == x['CDR3b']))
    }
  }
  nodes$size <- base::log2(nodes$clone_count)*clone_boost+min_size
  nodes$color.border <- "black"
  nodes$color.highlight <- "red"
  
  get_u_motifs <- function(x){
    
    get_unique_str <- function(x){
      res <- base::unique(e[x][e[x] != "-"])
      return(base::paste(base::unique(base::unlist(
        base::strsplit(res, ", "))), collapse = ", "))
    }
    
    id <- base::as.numeric(x["id"])
    e <- edges[edges$from == id | edges$to == id,]
    l <- list("CDR3b_global", "CDR3b_local", "CDR3a_global", "CDR3a_local")
    return(lapply(X=l, FUN = get_unique_str))
  }
  
  set_node_title <- function(x) {
    
    ma <- "<mark style=\"background-color: white; color: black;\">"
    
    bma <- base::paste(
    " </mark><mark style=\"background-color: blue; color: white;\">",
    "<i>(&beta;)</i><br>", ma)
    
    ama <- base::paste(
      " </mark><mark style=\"background-color: yellow; color: black;\">",
      "<i>(&alpha;)</i>", ma)

    l <- x[["label"]]
    l <- base::gsub(pattern = " \\(b\\) - ", replacement = bma, x = l) 
    l <- base::gsub(pattern = " \\(a\\)", replacement = ama, x = l)
    l <- base::paste(ma, "<b>", l, "</b></mark><br><br>")
    
    cc <- x[["clone_count"]]
    cc <- base::paste("<b> Clone count:", cc, "</b><br>")
    
    m <- get_motif_list(x)
    
    return(base::paste(l, cc, "<br>", m))
  }
  
  get_color <- function(x){
    b <- x['CDR3b_global'] != "" | x['CDR3b_local'] != ""
    a <- x['CDR3a_global'] != "" | x['CDR3a_local'] != ""
    if(a&b) { return("green") }
    if(b) { return("blue") }
    return("yellow")
  }
  
  nodes$group <- base::apply(X = nodes, MARGIN = 1, FUN = get_u_motifs)
  nodes$CDR3b_global <- base::unlist(lapply(nodes$group, `[`, 1))
  nodes$CDR3b_local <- base::unlist(lapply(nodes$group, `[`, 2))
  nodes$CDR3a_global <- base::unlist(lapply(nodes$group, `[`, 3))
  nodes$CDR3a_local <- base::unlist(lapply(nodes$group, `[`, 4))
  nodes$group <- base::apply(X = nodes, MARGIN = 1, function(x) {
    u <- base::unlist(x["group"])
    base::paste(u[base::nzchar(u)], collapse = ", ")
  })
  nodes$shape <- base::ifelse(test = nodes$clone_count > 1,
                              yes = "dot", no = "diamond")
  
  nodes$color.background <- base::apply(X = nodes, MARGIN = 1, FUN = get_color)
  nodes$title <- base::apply(X = nodes, MARGIN = 1, FUN = set_node_title)
  
  nodes$shadow <- FALSE
  nodes <- nodes[order(nodes$label),]
  return(nodes)
}


configure_edges <- function(edges) {
  
  get_edge_id <- function(x){
    t <- base::c(base::as.numeric(x["from"]), base::as.numeric(x["to"]))
    return(base::paste(base::min(t), base::max(t)))
  }
  
  edges$edge_id <- base::apply(X = edges, MARGIN = 1, FUN = get_edge_id)
  
  edges <- base::unique(edges[, c("edge_id", "motif", "type", "chain")])
  
  alpha_g <- edges[(edges$chain == "CDR3a") & edges$type == "global",
                     c("edge_id", "motif")]
  alpha_l <- edges[(edges$chain == "CDR3a") & edges$type == "local",
                     c("edge_id", "motif")]
  beta_g <- edges[(edges$chain == "CDR3b") & edges$type == "global",
                   c("edge_id", "motif")]
  beta_l <- edges[(edges$chain == "CDR3b") & edges$type == "local",
                   c("edge_id", "motif")]

  # hotfix, to be fixed in a better way once rewritten rest of the code
  if(!nrow(alpha_g) == 0) {
    alpha_g <- do.call(cbind.data.frame, 
                       aggregate(motif ~ edge_id, data = alpha_g, 
                                 FUN = function(x) { 
                                   count <- length(x)
                                   c(motif = paste(x, collapse = ", "), 
                                     count = count)})
    )
    names(alpha_g) <- c("edge_id", "CDR3a_global", "CDR3a_global_count")
  }
  if(!nrow(alpha_l) == 0) {
    alpha_l <- do.call(cbind.data.frame, 
                       aggregate(motif ~ edge_id, data = alpha_l, 
                                 FUN = function(x) { 
                                   count <- length(x)
                                   c(motif = paste(x, collapse = ", "), 
                                     count = count)})
    )
    names(alpha_l) <- c("edge_id", "CDR3a_local", "CDR3a_local_count")
  }
  if(!nrow(beta_g) == 0) {
    beta_g <- do.call(cbind.data.frame, 
                      aggregate(motif ~ edge_id, data = beta_g, 
                                FUN = function(x) { 
                                  count <- length(x)
                                  c(motif = paste(x, collapse = ", "), 
                                    count = count)})
    )
    names(beta_g) <- c("edge_id", "CDR3b_global", "CDR3b_global_count")
    
  }
  if(!nrow(beta_l) == 0) {
    beta_l <- do.call(cbind.data.frame, 
                      aggregate(motif ~ edge_id, data = beta_l, 
                                FUN = function(x) { 
                                  count <- length(x)
                                  c(motif = paste(x, collapse = ", "), 
                                    count = count)})
    )
    names(beta_l) <- c("edge_id", "CDR3b_local", "CDR3b_local_count")
    
  }
  
  a <- merge(alpha_g, alpha_l, by = "edge_id", all = TRUE)
  b <- merge(beta_g, beta_l, by = "edge_id", all = TRUE)
  ab <- merge(a, b, by = "edge_id", all = TRUE)
  
  # another hotfix, to be rewritten into a function
  if(!("CDR3a_global" %in% colnames(ab))) ab$CDR3a_global = NA
  if(!("CDR3b_global" %in% colnames(ab))) ab$CDR3b_global = NA
  if(!("CDR3a_local" %in% colnames(ab))) ab$CDR3a_local = NA
  if(!("CDR3b_local" %in% colnames(ab))) ab$CDR3b_local = NA
  if(!("CDR3a_global_count" %in% colnames(ab))) ab$CDR3a_global_count = NA
  if(!("CDR3b_global_count" %in% colnames(ab))) ab$CDR3b_global_count = NA
  if(!("CDR3a_local_count" %in% colnames(ab))) ab$CDR3a_local_count = NA
  if(!("CDR3b_local_count" %in% colnames(ab))) ab$CDR3b_local_count = NA
  
  ab$CDR3a_global[is.na(ab$CDR3a_global)] <- "-"
  ab$CDR3a_global_count[is.na(ab$CDR3a_global_count)] <- 0
  ab$CDR3a_local[is.na(ab$CDR3a_local)] <- "-"
  ab$CDR3a_local_count[is.na(ab$CDR3a_local_count)] <- 0
  ab$CDR3b_global[is.na(ab$CDR3b_global)] <- "-"
  ab$CDR3b_global_count[is.na(ab$CDR3b_global_count)] <- 0
  ab$CDR3b_local[is.na(ab$CDR3b_local)] <- "-"
  ab$CDR3b_local_count[is.na(ab$CDR3b_local_count)] <- 0
  
  ab$CDR3a_global_count <- as.numeric(ab$CDR3a_global_count)
  ab$CDR3a_local_count <- as.numeric(ab$CDR3a_local_count)
  ab$CDR3b_global_count <- as.numeric(ab$CDR3b_global_count)
  ab$CDR3b_local_count <- as.numeric(ab$CDR3b_local_count)
  
  ab$global_count <- ab$CDR3a_global_count + ab$CDR3b_global_count
  ab$local_count <- ab$CDR3a_local_count + ab$CDR3b_local_count
  
  ab$total_count <- ab$global_count + ab$local_count
  
  edges <- ab
  
  edges$from <- base::as.numeric(base::apply(X = edges, MARGIN = 1, function(x)
    base::strsplit(x["edge_id"], " ")[[1]][1]))
  edges$to <- base::as.numeric(base::apply(X = edges, MARGIN = 1, function(x)
    base::strsplit(x["edge_id"], " ")[[1]][2]))
    
  edges$length <- 15
  edges$width <- edges$total_count
  edges$type <- base::ifelse(test = (edges$global_count > 0),
                             yes = "global",
                             no = "local")
  edges$type <- base::ifelse(test = (edges$global_count > 0 & 
                                       edges$local_count > 0),
                             yes = "local & global",
                             no = edges$type)

  edges$color <- base::ifelse(test = (edges$type == "local"),
                              yes = "gray",
                              no = "orange")
  edges$color <- base::ifelse(test = (edges$type == "local & global"),
                              yes = "#9A0000",
                              no =  edges$color )
  
  edges$title <- base::apply(X = edges, MARGIN = 1, FUN = get_motif_list)
  
  edges$arrows <- ""
  edges$dashes <- FALSE
  edges$smooth <- FALSE
  edges$shadow <- FALSE
  return(edges)
}


get_motif_list <- function(x) {
  t_b   <- ""
  t_b_g <- ""
  t_b_l <- ""
  t_a   <- ""
  t_a_g <- ""
  t_a_l <- ""
  b   <- FALSE
  b_g <- FALSE
  b_l <- FALSE
  a   <- FALSE
  a_g <- FALSE
  a_l <- FALSE
  m_y <- "<mark style=\"background-color: yellow; color: black;\">"
  m_b <- "<mark style=\"background-color: blue; color: white;\">"
  
  if(x["CDR3b_global"] != "-" & x["CDR3b_global"] != "") {  
    t_b_g <- base::paste(
      "<b>Global:</b><br>", m_b, x["CDR3b_global"], "</mark><br>")
    b_g <- TRUE }
  if(x["CDR3b_local"] != "-" & x["CDR3b_local"] != "") {
    t_b_l <- base::paste(
      "<b>Local:</b><br>", m_b, x["CDR3b_local"], "</mark><br>")
    b_l <- TRUE }
  if(b_g | b_l) { 
    t_b <- "<b>CDR3-&beta;</b><br>"
    b <- TRUE }
  if(x["CDR3a_global"] != "-" & x["CDR3a_global"] != "") {
    t_a_g <- base::paste(
      "<b>Global:</b><br>", m_y, x["CDR3a_global"], "</mark><br>")
    a_g <- TRUE }
  if(x["CDR3a_local"] != "-" & x["CDR3a_local"] != "") { 
    t_a_l <- base::paste(
      "<b>Local:</b><br>", m_y, x["CDR3a_local"], "</mark><br>")
    a_l <- TRUE }
  if(a_g | a_l) { 
    t_a <- "<b>CDR3-&alpha;</b><br>"
    a <- TRUE }
  if(b & a) { return(base::paste(t_b, t_b_g, t_b_l, "<br>", 
                                 t_a, t_a_g, t_a_l)
  )}
  if(b & !a) { return(base::paste(t_b, t_b_g, t_b_l)
  )}
  if(!b & a) { return(base::paste(t_a, t_a_g, t_a_l)
  )}
  return("")
}


configure_network <- function(nodes, edges, ledges, lnodes){
  id_style <- base::paste("'width: ", 
                          base::max(base::nchar(nodes$label))*7.5,
                          "px;'")
  return(
    visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
    visNetwork::visIgraphLayout(layout = "layout_components", 
                                randomSeed = 1234) %>%
    visNetwork::visOptions(highlightNearest = 
                             base::list(enabled = TRUE, 
                                        degree = 1,
                                        algorithm = "hierarchical"),
                           selectedBy = base::list(variable = "group", 
                                                   multiple = TRUE,
                                                   sort = TRUE), 
                           manipulation = FALSE,
                           nodesIdSelection = base::list(
                             enabled = TRUE,
                             style = id_style)) %>%
    visNetwork::visLegend(addEdges = ledges, addNodes = lnodes, 
                          useGroups = FALSE, position = "right", 
                          width=0.15, zoom = FALSE)
  )
}


check_clustirr <- function(clust_irr){
  if(!methods::is(clust_irr, 'clust_irr')){
    base::stop("Input has to be object of class clust_irr")
  }
}
