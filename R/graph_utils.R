
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

config_vertices_plot <- function(g, is_jg, node_opacity) {
  # default features
  V(g)$size <- 1.5+log2(V(g)$clone_size)
  
  if(is_jg==TRUE) {
    V(g)$color_num <- as.numeric(as.factor(V(g)$sample))
    max_n <- max(V(g)$color_num)
    V(g)$color <- hcl.colors(n=max(5, max_n), palette = "Roma",
                             alpha = node_opacity)[V(g)$color_num]
    V(g)$frame.color <- V(g)$color
  } 
  else {
    V(g)$color <- adjustcolor("black", alpha.f = node_opacity)
    V(g)$frame.color <- V(g)$color
  }
  
  return(g)
}


save_graph <- function(graph, file_name, output_folder){
    
    graph$width <- "100%"
    graph$height <- "100vh"
    
    saveWidget(
        widget = graph, 
        file = paste0(file_name, ".html"), 
        selfcontained = TRUE, 
        background = "white",
        title = file_name
    )
    
    base_file <- paste0(file_name, ".html")
    target_path <- file.path(output_folder, base_file)
    
    # this step is necessary to ensure self-contained files
    # if we put the full path in saveWidget, it will not delete the files folder
    if(dir.exists(output_folder)){
        file.rename(from = base_file, to = target_path)
        message(file_name, " exported to: ", target_path)
    } else {
        message("Directory '", output_folder, "' does not exist. \n",
                "File saved to current working directory as: ", base_file)
    }
    
}
