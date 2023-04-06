# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_chain_run_v1 <- function(cdr3,
                             cdr3_ref,
                             ks,
                             cores,
                             control) {

    # 1. local clustering
    # a. get local motifs
    motifs <- lapply(X = ks,
                     FUN = get_motifs_v1,
                     cdr3 = cdr3,
                     cdr3_ref = cdr3_ref,
                     B = control$B,
                     min_o = control$local_min_o,
                     cores = cores)
    names(motifs) <- as.character(ks)

    # b. compute local enrichment scores
    future::plan(future::multisession, workers = cores)
    motif_enrichment <- do.call(rbind, future.apply::future_lapply(
        X = ks,
        FUN = get_motif_enrichment_v1,
        m = motifs,
        B = control$B,
        cores = cores,
        future.seed = TRUE))
    future::plan(future::sequential())


    # c. add filter flag
    motif_enrichment <- get_motif_filter_v1(
        m = motif_enrichment,
        min_p = control$local_min_p,
        min_ove = control$local_min_ove,
        min_o = control$local_min_o)

    # d. find motifs in CDR3
    local_pairs <- get_local_pair(
        cdr3 = cdr3,
        motif = motif_enrichment$motif[motif_enrichment$filter==TRUE])

    # 2. global
    # if global_pairs are provided as input use them, else compute them
    if(!is.null(control$global_pairs)) {
        global_pairs <- control$global_pairs
    }
    else {
        if(control$low_mem) {
            global_pairs <- get_global_pairs_mem(
                cdr3 = cdr3,
                global_max_dist = control$global_max_dist)
        }
        else {
            global_pairs <- get_global_pairs(
                cdr3 = cdr3,
                global_max_dist = control$global_max_dist)
        }
    }


    # 3. return TODO: format properly
    return(list(local_pairs = local_pairs,
                global_pairs = global_pairs,
                motif_enrichment = motif_enrichment))
}


# Description:
# Computes motif frequencies for a sample and reference
# uses functions:
# * get_kmers_freq_ref
# * get_kmers_freq_sample
get_motifs_v1 <- function(cdr3,
                          cdr3_ref,
                          B,
                          ks,
                          cores,
                          min_o) {

    # Description:
    # x = k in k-mer
    get_kmers_freq_ref <- function(x,
                                   cdr3,
                                   N,
                                   B,
                                   relevant_motifs,
                                   cores){

        get_qgrams <- function(x, q, cdr3, N, relevant_motifs) {
            draw_cdr3 <- sample(x = cdr3, size = N, replace = TRUE)
            o <- stringdist::qgrams(draw_cdr3, q = q)
            if(ncol(o)==0) {
                return(NA)
            }
            o <- o[1,]
            o <- o[names(o) %in% relevant_motifs]
            if(length(o)==0) {
                return(NA)
            }
            return(o)
        }

        future::plan(future::multisession, workers = cores)
        o <- future.apply::future_lapply(
            X = 1:B,
            q = x,
            N = N,
            relevant_motifs = relevant_motifs,
            cdr3 = cdr3,
            FUN = get_qgrams,
            future.seed = TRUE)
        future::plan(future::sequential())
        return(o)
    }


    # Description:
    # x = k in k-mer
    get_kmers_freq_sample <- function(x,
                                      cdr3,
                                      min_o){

        o <- stringdist::qgrams(cdr3, q = x)
        if(ncol(o)==0) {
            return(NA)
        }
        o <- o[1, o[1,]>=min_o]
        if(length(o)==0) {
            return(NA)
        }
        return(o)
    }

    # find motifs in sample
    motif_sample <- lapply(
        X = ks, # edit here
        FUN = get_kmers_freq_sample,
        cdr3 = cdr3, # edit here
        min_o = min_o)
    names(motif_sample) <- ks


    # do sampling & find motifs
    found_kmers <- as.vector(unlist(
        lapply(X = motif_sample, FUN = names)))
    motif_ref <- lapply(
        X = ks,
        FUN = get_kmers_freq_ref,
        cdr3 = cdr3_ref,
        B = B,
        N = length(cdr3),
        relevant_motifs = found_kmers,
        cores = cores)
    names(motif_ref) <- ks

    return(list(motif_sample = motif_sample,
                motif_ref = motif_ref))
}


# Description:
# Computes motif enrichment with data collected by function:
# * get_motifs
get_motif_enrichment_v1 <- function(x,
                                    m,
                                    B,
                                    cores) {

    get_e <- function(x, B) {
        # o
        o <- x[1]

        # e
        e <- x[-1]

        # mean, max
        e_mean <- base::mean(e)
        e_max <- base::max(e)

        # OvE = /e
        ove <- o/e_mean

        # prob e>=o
        p <- sum(e>=o)/B

        # return
        return(c(e_mean, e_max, o, ove, p))
    }

    motif_sample <- m[[as.character(x)]]$motif_sample[[1]]
    motif_ref <- m[[as.character(x)]]$motif_ref[[1]]

    # matrix of k-mer counts
    f_m <- matrix(data = 0, ncol = B+1, nrow = length(motif_sample))
    rownames(f_m) <- names(motif_sample)
    f_m[names(motif_sample), 1] <- motif_sample
    for(i in 1:length(motif_ref)) {
        s <- motif_ref[[i]]
        f_m[names(s),i+1] <- s
    }

    e <- t(apply(X = f_m, MARGIN = 1, FUN = get_e, B = B))

    # format output
    e <- data.frame(e)
    colnames(e) <- c("sim_mean", "sim_max",
                     "obs", "ove", "p")
    e$motif <- rownames(e)
    e$k <- x
    return(e)
}


# Description:
# Given data from function get_motif_enrichment, this function adds a
# filter column with filter = T if motif is enriched (given the input
# standards) and filter = F if not-enriched
get_motif_filter_v1 <- function(m,
                                min_p,
                                min_ove,
                                min_o) {
    m$filter <- FALSE
    m$filter[ m$p   <= min_p   &
              m$ove >= min_ove &
              m$obs >= min_o     ] <- TRUE
    return(m)
}
