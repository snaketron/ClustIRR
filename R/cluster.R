cluster_irr <- function(s,
                        r,
                        version = 3,
                        ks = 4,
                        cores = 1,
                        control = list(
                            B = 1000,
                            global_max_dist = 1,
                            local_max_fdr = 0.05,
                            local_min_ove = 2,
                            local_min_o = 1,
                            trim_flank_aa = 0,
                            global_pairs = NULL,
                            low_mem = FALSE)) {
  
    # 0. control check
    control <- get_control(control_in = control)
    # 1. input check
    input_check(s = s, r = r, version = version, ks = ks, 
                cores = cores, control = control)
    # get chains to be analyzed
    chains <- get_chains(base::colnames(s))
    # add ID to data
    s$id <- base::seq_len(length.out = base::nrow(s))

    # run analysis for each chain (if available)
    clust <- base::vector(mode = "list", length = base::length(chains))
    base::names(clust) <- chains

    for(chain in chains) {
        if(version==3) {
            cdr3 <- s[, chain]
            cdr3_ref <- r[, chain]
        }
        else {
            cdr3 <- base::unique(s[, chain])
            cdr3_ref <- base::unique(r[, chain])
        }
        if(version==1) {
            clust[[chain]] <- get_clust_v1(cdr3 = cdr3,
                                           cdr3_ref = cdr3_ref,
                                           ks = ks,
                                           cores = cores,
                                           control = control)
        }
        if(version==2|version==3) {
            clust[[chain]] <- get_clust_v23(cdr3 = cdr3,
                                            cdr3_ref = cdr3_ref,
                                            ks = ks,
                                            cores = cores,
                                            control = control)
        }
    }
    return(base::structure(class = "clust_irr",
                           base::list(clust = clust,
                                      inputs = base::list(
                                          s = s,
                                          r = r,
                                          version = version,
                                          ks = ks,
                                          cores = cores,
                                          control = control))))
}




# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust_v1 <- function(cdr3,
                         cdr3_ref,
                         ks,
                         cores,
                         control) {

    # 1. local
    l <- get_localclust_v1(cdr3 = cdr3,
                           cdr3_ref = cdr3_ref,
                           ks = ks,
                           cores = cores,
                           control = control)

    # 2. global
    # if global_pairs are provided as input use them, else compute them
    if(!base::is.null(control$global_pairs)) {
        g <- control$global_pairs
    }
    else {
        if(control$low_mem) {
            g <- get_global_clust_mem(cdr3 = base::unique(cdr3),
                                      global_max_dist = control$global_max_dist)
        }
        else {
            g <- get_global_clust(cdr3 = base::unique(cdr3),
                                  global_max_dist = control$global_max_dist)
        }
    }
    return(base::list(local = l, global = g))
}



# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust_v23 <- function(cdr3,
                          cdr3_ref,
                          ks,
                          cores,
                          control) {

    # 1. local
    l <- get_localclust_v23(cdr3 = cdr3,
                            cdr3_ref = cdr3_ref,
                            ks = ks,
                            cores = cores,
                            control = control)

    # 2. global
    # if global_pairs are provided as input use them, else compute them
    if(!base::is.null(control$global_pairs)) {
        g <- control$global_pairs
    }
    else {
        if(control$low_mem) {
            g <- get_global_clust_mem(cdr3 = cdr3,
                                      global_max_dist = control$global_max_dist)

        }
        else {
            g <- get_global_clust(cdr3 = cdr3,
                                  global_max_dist = control$global_max_dist)
        }
    }
    return(base::list(local = l, global = base::unique(g)))
}
