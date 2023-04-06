# Description:
# This is rewritten code for GLIPH v1 and v2
# Details are missing, TODO
# data_sample: must be a data.frame that has the following columns
# (CDR3b, TRBV, TRBJ, CDR3a, TRAV, TRAJ, sample_id), bare minimum
# (CDR3b, sample_id)
#' Find groups of TCRs with locally or globally similar CDR3s
#'
#' @param data_sample, data.frame: TCR sample
#' @param data_ref, data.frame: reference database
#' @param version, integer: version = 1, 2 or 3, gliph version to use
#' @param ks, vector of integers: motif lengths to use (default ks=(2,3,4)),
#' @param cores, integer: number of CPU cores to use
#' @param control, list: auxiliarry input parameters (described below)
#'
#' @return list(clust = clust, edges = edges, data_sample = data_sample)
#' @export
#'
#' @examples
gliph <- function(data_sample,
                  data_ref,
                  version = 2,
                  ks = c(2, 3, 4),
                  cores = 1,
                  control = list(
                      B = 1000,
                      global_max_dist = 1,
                      local_min_p = 0.05,
                      local_min_ove = 2,
                      local_min_o = 3,
                      trim_flanks = FALSE,
                      flank_size = 3,
                      global_pairs = NULL)) {


    # a. control check
    control <- get_control(control_in = control)

    # 1. parameter check
    # parameter_check(data_sample = data_sample,
    #                 data_ref = data_ref,
    #                 version = version,
    #                 ks = ks,
    #                 cores = cores,
    #                 control = control)

    # get chains to be analyzed
    chains <- get_chains(colnames(data_sample))


    # add ID to gliph
    data_sample$ID <- 1:nrow(data_sample)


    # run analysis for each chain (if available)
    # clust and edges will store *all* clustering results and *significant*
    # cell-cell edges (connections in graph) based on the specified filtering
    # criteria
    clust <- vector(mode = "list", length = length(chains))
    names(clust) <- chains
    edges <- vector(mode = "list", length = length(chains))
    names(edges) <- chains

    for(chain in chains) {
        if(control$trim_flanks) {
            # NAs ignored by qqgram, how about global dist? TODO
            data_sample[, chain] <- get_trimmed_flanks(
                x = data_sample[, chain],
                flank_size = control$flank_size)
            data_ref[, chain] <- get_trimmed_flanks(
                x = data_ref[, chain],
                flank_size = control$flank_size)
        }

        # create edges of a graph
        if(version==3) {
            cdr3 <- data_sample[, chain]
            cdr3_ref <- data_ref[, chain]
        }
        else {
            cdr3 <- unique(data_sample[, chain])
            cdr3_ref <- unique(data_ref[, chain])
        }


        # run local + global clustering
        if(version==1) {
            clust[[chain]] <- get_chain_run_v1(
                cdr3 = cdr3,
                cdr3_ref = cdr3_ref,
                ks = ks,
                cores = cores,
                control = control)
        }
        if(version==2|version==3) {
            clust[[chain]] <- get_chain_run_v2(
                cdr3 = cdr3,
                cdr3_ref = cdr3_ref,
                ks = ks,
                cores = cores,
                control = control)
        }

        edges[[chain]] <- NA # large clonal expansions are present->memory/disk
        # edges[[chain]] <- get_edges(local_pairs=clust[[chain]]$local_pairs,
        #                             global_pairs=clust[[chain]]$global_pairs,
        #                             cdr3 = cdr3,
        #                             chain = chain)
    }

    return(list(clust = clust,
                edges = edges,
                data_sample = data_sample,
                control = control))
}
