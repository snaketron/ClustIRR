cluster_irr <- function(s,
                        control = list(global_min_identity = 0.7,
                                       trim_flank_aa = 3,
                                       global_hamming = FALSE,
                                       global_max_hdist = 1,
                                       low_mem = FALSE)) {
  # control check
  control <- get_control(control_in = control)
  
  # input check
  input_check(s = s, control = control)

  # get chains to be analyzed
  chains <- get_chains(colnames(s))
  
  # add clone size info, sample and clone id
  s <- get_clone_size(s)
  s <- get_sample(s)
  s$id <- seq_len(length.out = nrow(s))
  
  # do clustering
  clust <- lapply(X = chains, FUN = get_clust, s = s, control = control)
  names(clust) <- chains
  
  # setup clustirr object
  return(get_clustirr_output_obj(clust = clust, 
                                 s = s, 
                                 control = control))
}


# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust <- function(x, s, control) {
  
  get_cdr3s_for_global <- function(x, chain) {
    # multiple cdr3 x clone_count -> cells
    # in the global analysis this is reverted
    x <- x[, chain]
    x <- x[is.na(x)==FALSE]
    if(length(x)==0) {
      stop("no CDR3s found")
    }
    return(x)
  }
  
  g <- get_global_clust(cdr3 = get_cdr3s_for_global(x = s, chain = x),
                        control = control)
  
  return(g)
}
