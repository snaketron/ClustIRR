cluster_irr <- function(s,
                        r,
                        ks = 4,
                        cores = 1,
                        control = list(global_smart = FALSE,
                                       global_max_hdist = 1,
                                       global_min_identity = 0.7,
                                       local_max_fdr = 0.05,
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
  
  clust <- bplapply(X = chains, 
                    FUN = get_clust,
                    s = s, 
                    r = r,
                    ks = ks,
                    control = control,
                    global_only = global_only,
                    BPPARAM = MulticoreParam(workers = cores))
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
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust <- function(x,
                      s,
                      r,
                      ks,
                      control,
                      global_only) {
  
  get_cdr3s_for_local <- function(x, chain) {
    # multiple cdr3 x clone_count -> cells
    # in the global analysis this is reverted
    x <- rep(x[, chain], times = x$clone_size)
    x <- x[is.na(x)==FALSE]
    if(length(x)==0) {
      stop("no CDR3s found")
    }
    return(x)
  }
  
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
  
  # 1. local
  l <- NULL
  if(global_only==FALSE) {
    l <- get_localclust(cdr3 = get_cdr3s_for_local(x = s, chain = x), 
                        cdr3_ref = get_cdr3s_for_local(x = r, chain = x), 
                        ks = ks, 
                        control = control)
  }
  
  # 2. global
  g <- get_global_clust(cdr3 = get_cdr3s_for_global(x = s, chain = x),
                        control = control)
  
  return(list(local = l, global = g))
}

