cluster_irr <- function(s,
                        r,
                        version = 2,
                        ks = 4,
                        cores = 1,
                        control = list(global_max_dist = 1,
                                       local_max_fdr = 0.05,
                                       local_min_ove = 2,
                                       local_min_o = 1,
                                       trim_flank_aa = 0,
                                       global_pairs = NULL,
                                       low_mem = FALSE)) {
    
    # control check
    control <- get_control(control_in = control)
    
    # input check
    input_check(s = s, r = r, version = version, ks = ks,
                cores = cores, control = control)
    
    # get chains to be analyzed
    chains <- get_chains(colnames(s))
    
    # add ID to data
    s$id <- seq_len(length.out = nrow(s))
    
    # run analysis for each chain
    clust <- lapply(X = chains, 
                    FUN = run_chain_analysis,
                    s = s, r = r, version = version,
                    ks = ks, cores = cores, control = control)
    names(clust) <- chains
    
    # setup clustirr object
    return(get_clustirr_output_obj(clust = clust, 
                                   s = s, 
                                   r = r, 
                                   version = version, 
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
                               version, 
                               ks, 
                               cores, 
                               control) {
    
    chain <- x
    s_t <- s[!is.na(s[, chain]),chain]
    r_t <- r[!is.na(r[, chain]),chain]
    if(version == 2) {
        cdr3 <- s_t
        cdr3_ref <- r_t
    }
    else {
        cdr3 <- unique(s_t)
        cdr3_ref <- unique(r_t)
    }
    return(get_clust(cdr3 = cdr3,
                     cdr3_ref = cdr3_ref,
                     version = version,
                     ks = ks,
                     cores = cores,
                     control = control))
}



# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust <- function(cdr3,
                      cdr3_ref,
                      version,
                      ks,
                      cores,
                      control) {
    # 1. local
    l <- get_localclust(cdr3 = cdr3,
                        cdr3_ref = cdr3_ref,
                        ks = ks,
                        control = control)
    
    # 2. global
    # if global_pairs are provided as input use them, else compute them
    if(!is.null(control$global_pairs)) {
        g <- control$global_pairs
    }
    else {
      if(control$global_smart) {
        g <- get_global_clust_smart(cdr3 = unique(cdr3))
      } 
      else {
        g <- get_global_clust(cdr3 = unique(cdr3),
                              global_max_dist = control$global_max_dist,
                              low_mem = control$low_mem)
      }
    }
    return(list(local = l, global = g))
}
