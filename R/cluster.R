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
                            low_mem = FALSE
                        )) {
    # 0. control check
    control <- get_control(control_in = control)
    # 1. input check
    input_check(
        s = s, r = r, version = version, ks = ks,
        cores = cores, control = control
    )
    # get chains to be analyzed
    chains <- get_chains(colnames(s))
    # add ID to data
    s$id <- seq_len(length.out = nrow(s))
    # run analysis for each chain (if available)
    clust <- vector(mode = "list", length = length(chains))
    names(clust) <- chains

    for (chain in chains) {
        s_t <- s[!is.na(s[, chain]),chain]
        r_t <- r[!is.na(r[, chain]),chain]
        if (version == 3) {
            cdr3 <- s_t
            cdr3_ref <- r_t
        } else {
            cdr3 <- unique(s_t)
            cdr3_ref <- unique(r_t)
        }
        if (version == 1) {
            clust[[chain]] <- get_clust_v1(
                cdr3 = cdr3,
                cdr3_ref = cdr3_ref,
                ks = ks,
                cores = cores,
                control = control
            )
        }
        if (version == 2 | version == 3) {
            clust[[chain]] <- get_clust_v23(
                cdr3 = cdr3,
                cdr3_ref = cdr3_ref,
                version = version,
                ks = ks,
                cores = cores,
                control = control
            )
        }
    }
    return(structure(
        class = "clust_irr",
        list(
            clust = clust,
            inputs = list(
                s = s,
                r = r,
                version = version,
                ks = ks,
                cores = cores,
                control = control
            )
        )
    ))
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
    l <- get_localclust_v1(
        cdr3 = cdr3,
        cdr3_ref = cdr3_ref,
        ks = ks,
        cores = cores,
        control = control
    )

    # 2. global
    # if global_pairs are provided as input use them, else compute them
    if (!is.null(control$global_pairs)) {
        g <- control$global_pairs
    } else {
        if (control$low_mem) {
            g <- get_global_clust_mem(
                cdr3 = unique(cdr3),
                global_max_dist = control$global_max_dist
            )
        } else {
            g <- get_global_clust(
                cdr3 = unique(cdr3),
                global_max_dist = control$global_max_dist
            )
        }
    }
    return(list(local = l, global = g))
}



# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust_v23 <- function(cdr3,
                            cdr3_ref,
                            version,
                            ks,
                            cores,
                            control) {
    # 1. local
    l <- get_localclust_v23(
        cdr3 = cdr3,
        cdr3_ref = cdr3_ref,
        ks = ks,
        cores = cores,
        control = control
    )

    # 2. global
    # if global_pairs are provided as input use them, else compute them
    if (!is.null(control$global_pairs)) {
        g <- control$global_pairs
    } else {
        if (control$low_mem) {
            g <- get_global_clust_mem(
                cdr3 = unique(cdr3),
                global_max_dist = control$global_max_dist
            )
        } else {
            g <- get_global_clust(
                cdr3 = unique(cdr3),
                global_max_dist = control$global_max_dist
            )
        }
        if (version == 3) {
            clones <- unique(cdr3[duplicated(cdr3)])
            clones <- matrix(rep(clones, 2),
                nrow = length(clones), ncol = 2
            )
            if (is.null(g)) {
                g <- clones
            } else {
                g <- rbind(clones, g)
            }
        }
    }
    return(list(local = l, global = g))
}
