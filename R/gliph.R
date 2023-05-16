#' TCR clustering for specificity analysis
#'
#' @description Find groups of TCRs with locally or globally similar CDR3s.
#' Based on Gliph and Gliph2 clustering method.
#'
#' @param data_sample data.frame: TCR sample. Must be a data.frame that has
#' the following columns: CDR3b, TRBV, TRBJ, CDR3a, TRAV, TRAJ, sample_id.
#' Bare minimum: CDR3b, sample_id
#' @param data_ref data.frame: reference database
#' @param version integer: version = 1, 2 or 3, gliph version to use
#' @param ks vector of integers: motif lengths to use (default ks=(2,3,4))
#' @param cores integer: number of CPU cores to use
#' @param control list: auxiliary input parameters (described below)
#' \itemize{
#'   \item B - integer: simulation depth
#'   \item global_max_dist - integer: maximum hamming distance for global
#'   clustering
#'   \item local_max_fdr - numeric: maximum cutoff p-value for random
#' generation
#'   \item local_min_ove - numeric: minimum fold enrichment
#'   \item local_min_o - numeric: minimum motif observations
#'   \item trim_flank_aa - integer: cut off value for trimming aa flanks
#'   \item low_mem - logical: low memory mode. Slower looping, lower memory
#'   footprint
#'   \item global_pairs - matrix: optional pre-computed global pairs
#' }
#'
#' @return gliphR returns a list of the following elements:
#' \itemize{
#'    \item clust - list: local + global clusters
#'    \item edges - list: local + global edges
#'    \item data_sample - data.frame: examined data sample
#'    \item version - integer: used gliph version
#'    \item ks - vector of integers: used motif lengths
#'    \item cores - integer: number of used CPU cores
#'    \item control - list: used auxiliary input parameters
#' }
#'
#' @examples
#' ## this example shows how to test gliphR with minimal input
#'
#' # load package input data
#' data("hs_CD8_ref")
#'
#' # create a minimal reference set from the first 1000 rows
#' data_ref <- hs_CD8_ref[1:1000, 1:3]
#'
#' # draw 50 samples from the reference dataset
#' data_sample <- hs_CD8_ref[
#' sample(
#' x = 1:nrow(data_ref),
#' size = 50,
#' replace = FALSE
#' ),
#' 1:3]
#'
#' # set parameters
#' ks <- c(2, 3, 4)
#' cores <- parallel::detectCores()
#' version <- 3
#'
#' control_input <- list(
#' B = 100,
#' global_max_dist = 1,
#' local_max_fdr = 0.05,
#' local_min_ove = 2,
#' local_min_o = 3,
#' trim_flank_aa = 3,
#' low_mem = FALSE,
#' global_pairs = NULL
#' )
#'
#' # run gliph and save output to gliph_output
#' gliph_output <- gliph(
#' data_sample = data_sample,
#' data_ref = data_ref,
#' ks = ks,
#' cores = cores,
#' version = version,
#' control = control_input
#' )
#'
#' @import stringdist
#' @import future
#' @import future.apply
#' @import methods
#' @import stats
#' @import utils
#' @importFrom igraph components graph_from_data_frame
#'
#' @export
gliph <- function(data_sample,
                  data_ref,
                  version = 2,
                  ks = c(2, 3, 4),
                  cores = 1,
                  control = list(
                      B = 1000,
                      global_max_dist = 1,
                      local_max_fdr = 0.05,
                      local_min_ove = 2,
                      local_min_o = 3,
                      trim_flank_aa = 0,
                      global_pairs = NULL,
                      low_mem = FALSE
                  )) {
    # a. control check
    control <- get_control(control_in = control)

    # 1. input check
    input_check(
        data_sample = data_sample,
        data_ref = data_ref,
        version = version,
        ks = ks,
        cores = cores,
        control = control
    )

    # get chains to be analyzed
    chains <- get_chains(colnames(data_sample))

    # add ID to gliph
    data_sample$ID <- seq_len(nrow(data_sample)) # 1:nrow(data_sample)

    # run analysis for each chain (if available)
    clust <- vector(mode = "list", length = length(chains))
    names(clust) <- chains
    edges <- vector(mode = "list", length = length(chains))
    names(edges) <- chains

    for (chain in chains) {
        if (version == 3) {
            cdr3 <- data_sample[, chain]
            cdr3_ref <- data_ref[, chain]
        } else {
            cdr3 <- unique(data_sample[, chain])
            cdr3_ref <- unique(data_ref[, chain])
        }
        # run local + global clustering
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
                ks = ks,
                cores = cores,
                control = control
            )
        }
        # TODO
        edges[[chain]] <- NA # large clonal expansions are present->memory/disk
        # edges[[chain]] <- get_edges(local_pairs=clust[[chain]]$local_pairs,
        #                             global_pairs=clust[[chain]]$global_pairs,
        #                             cdr3 = cdr3,
        #                             chain = chain)
    }
    return(list(
        clust = clust,
        edges = edges,
        data_sample = data_sample,
        version = version,
        ks = ks,
        cores = cores,
        control = control
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
                cdr3 = base::unique(cdr3),
                global_max_dist = control$global_max_dist
            )
        } else {
            g <- get_global_clust(
                cdr3 = base::unique(cdr3),
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
                cdr3 = base::unique(cdr3),
                global_max_dist = control$global_max_dist
            )
        } else {
            g <- get_global_clust(
                cdr3 = base::unique(cdr3),
                global_max_dist = control$global_max_dist
            )
        }
    }
    return(list(local = l, global = g))
}
