cluster_irr <- function(s,
                        r,
                        ks = 4,
                        cores = 1,
                        control = list(global_smart = FALSE,
                                       global_max_hdist = 1,
                                       local_max_fdr = 0.05,
                                       local_min_ove = 2,
                                       local_min_o = 1,
                                       trim_flank_aa = 3,
                                       global_pairs = NULL,
                                       low_mem = FALSE)) {
  # control check
  control <- get_control(control_in = control)
  
  # input check
  input_check(s = s, r = r, ks = ks, cores = cores, control = control)
  
  # if r is null or missing -> only global analysis of s
  global_only <- ifelse(test = missing(r)||is.null(r), yes = TRUE, no = FALSE)
  if(global_only) {
    r <- NULL
  }
  
  # get chains to be analyzed
  chains <- get_chains(colnames(s))
  
  # add clone size info
  s <- get_clone_size(s)
  r <- get_clone_size(r)
  
  # add ID to data
  s$id <- seq_len(length.out = nrow(s))
  
  # run analysis for each chain
  clust <- lapply(X = chains, 
                  FUN = run_chain_analysis,
                  s = s, 
                  r = r,
                  ks = ks, 
                  cores = cores, 
                  control = control,
                  global_only = global_only)
  names(clust) <- chains
  
  # setup clustirr object
  return(get_clustirr_output_obj(clust = clust, 
                                 s = s, 
                                 r = r,
                                 ks = ks, 
                                 cores = cores, 
                                 control = control))
}


# Description:
# Wrapper to be run within the main call and will perform chain-specific
# analysis
run_chain_analysis <- function(x, 
                               s, 
                               r,
                               ks, 
                               cores, 
                               control,
                               global_only) {
  
  
  get_cdr3s <- function(x, chain) {
    # multiple cdr3 x clone_count -> cells
    # in the global analysis this is reverted
    x <- rep(x[, chain], times = x$clone_size)
    x <- x[is.na(x)==FALSE]
    if(length(x)==0) {
      stop("no CDR3s found")
    }
    return(x)
  }
  
  chain <- x
  if(global_only) {
    cdr3 <- get_cdr3s(x = s, chain = chain)
    cdr3_ref <- NULL
  } 
  else {
    cdr3 <- get_cdr3s(x = s, chain = chain)
    cdr3_ref <- get_cdr3s(x = r, chain = chain)
  }
  return(get_clust(cdr3 = cdr3,
                   cdr3_ref = cdr3_ref,
                   ks = ks,
                   cores = cores,
                   control = control,
                   global_only = global_only))
}



# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust <- function(cdr3,
                      cdr3_ref,
                      ks,
                      cores,
                      control,
                      global_only) {
  # 1. local
  l <- NULL
  if(global_only==FALSE) {
    l <- get_localclust(cdr3 = cdr3, 
                        cdr3_ref = cdr3_ref, 
                        ks = ks, 
                        control = control)
  }
  
  # 2. global
  g <- get_global_clust(cdr3 = cdr3,
                        control = control)
  
  return(list(local = l, global = g))
}
