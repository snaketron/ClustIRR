# Description:
# This algorithm finds groups of TCRs that likely have similar peptide:MHC
# specificity. This is done with two strategies: local motif clustering and
# global CDR3 sequence clustering.
#
# SK: here we need two paragraphs of how each algorithm (local vs. global)
# works, as well as *key* differences between v=1, 2 and 3
#'
#' @param data_sample
#' @param data_ref
#' @param version
#' @param ks
#' @param cores
#' @param control
#'
#' @import stringdist
#' @import future
#' @import future.apply
#' @import methods
#' @import stats
#' @import utils
#'
#' @return list(clust = clust, edges = edges, data_sample = data_sample)
#' @export
cluster_irr <- function(data_sample,
                        data_ref,
                        version = 3,
                        ks = c(2, 3, 4),
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
    input_check(data_sample = data_sample,
                data_ref = data_ref,
                version = version,
                ks = ks,
                cores = cores,
                control = control)

    # get chains to be analyzed
    chains <- get_chains(base::colnames(data_sample))

    # add ID to data
    data_sample$ID <- base::seq_len(length.out = base::nrow(data_sample))

    # run analysis for each chain (if available)
    clust <- base::vector(mode = "list", length = base::length(chains))
    base::names(clust) <- chains
    edges <- base::vector(mode = "list", length = base::length(chains))
    base::names(edges) <- chains

    for(chain in chains) {
        if(version==3) {
            cdr3 <- data_sample[, chain]
            cdr3_ref <- data_ref[, chain]
        }
        else {
            cdr3 <- base::unique(data_sample[, chain])
            cdr3_ref <- base::unique(data_ref[, chain])
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
        edges[[chain]] <- NA # TODO
    }
    return(base::list(clust = clust,
                      edges = edges,
                      inputs = base::list(
                          data_sample = data_sample,
                          version = version,
                          ks = ks,
                          cores = cores,
                          control = control)))
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
    return(base::list(local = l, gobal = g))
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
