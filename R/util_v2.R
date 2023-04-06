# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_chain_run_v2 <- function(cdr3,
                             cdr3_ref,
                             ks,
                             cores,
                             control) {

    # 1. local clustering
    # a. get local motifs
    motifs <- lapply(X = ks,
                     FUN = get_motifs_v2,
                     cdr3 = cdr3,
                     cdr3_ref = cdr3_ref,
                     min_o = control$local_min_o)
    motifs <- do.call(rbind, motifs)

    # b. compute enrichment with fisher's exact test
    motif_stats <- t(apply(
        X = motifs[, c("f_sample", "f_ref", "n_sample", "n_ref")],
        MARGIN = 1, FUN = get_motif_enrichment_fet_v2))
    motifs$ove <- motif_stats[,1]
    motifs$p_value <- motif_stats[,2]
    motifs$fdr <- stats::p.adjust(p = motifs$p_value, method = "fdr")
    rm(motif_stats)

    # c. add filter flag
    motif_enrichment <- get_motif_filter_v2(
        m = motifs,
        min_fdr = control$local_min_fdr,
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
get_motifs_v2 <- function(x,
                          cdr3,
                          cdr3_ref,
                          min_o) {

    # find kmers in sample
    kmers_s <- stringdist::qgrams(cdr3, q = x)
    if(ncol(kmers_s)==0) {
        # in real applications this should not happen -> input checks should
        # catch such errors
        stop("no kmers found in sample")
    }
    kmers_s <- kmers_s[1, kmers_s[1,]>=min_o]
    if(length(kmers_s)==0) {
        # in real applications this should not happen -> input checks should
        # catch such errors
        stop("no kmers found in sample")
    }


    # find kmers in reference
    kmers_r <- stringdist::qgrams(cdr3_ref, q = x)
    if(ncol(kmers_r)==0) {
        # in real applications this should not happen -> input checks should
        # catch such errors
        stop("no kmers found in reference")
    }
    kmers_r <- kmers_r[1,]

    # we are only interested in enrichment of motifs in sample relative to
    # reference. Remove all motifs from reference not found in sample.
    kmers_r <- kmers_r[names(kmers_r) %in% names(kmers_s)]


    # convert table to data.frame
    kmers_s <- data.frame(motif = names(kmers_s),
                          f_sample = as.numeric(kmers_s))
    kmers_s$n_sample <- sum(kmers_s$f_sample)
    kmers_r <- data.frame(motif = names(kmers_r),
                          f_ref = as.numeric(kmers_r))
    kmers_r$n_ref <- sum(kmers_r$f_ref)

    m <- merge(x = kmers_s, y = kmers_r, by = "motif", all = TRUE)
    m[is.na(m[,"f_sample"]), "f_sample"] <- 0
    m[is.na(m[,"f_ref"]), "f_ref"] <- 0
    m[is.na(m[,"n_ref"]), "n_ref"] <- kmers_r$n_ref[1]
    m[is.na(m[,"n_sample"]), "n_sample"] <- kmers_s$n_sample[1]

    m$k <- x

    # return for statistical test
    return(motifs = m)
}


# Description:
# Compute motif enrichment with Fisher's exact test (fet) based on data
# collected by function get_motifs_v2
get_motif_enrichment_fet_v2 <- function(x) {
    # Description of parameters used in hypergeometric test (below)
    # f_sample = x[1]
    # f_ref = x[2]
    # n_sample = x[3]
    # n_ref = x[4]
    #
    # q = f_sample,
    # m = f_ref+f_sample,
    # n = n_ref+n_sample-(f_ref+f_sample),
    # k = n_ref+n_sample

    # ove TODO: check how is this done in gliph2
    ove <- (x[1]/x[3])/((x[2]/x[4]))

    q <- x[1]
    m <- x[1]+x[2]
    k <- x[3]+x[4]
    n <- k-m

    # TODO: do some testing on demo data
    # u <- matrix(data = NA_integer_, nrow = 2, ncol = 2)
    u <- matrix(data = 0, nrow = 2, ncol = 2)

    u[1,1] <- x[1]
    u[1,2] <- x[2]
    u[2,1] <- x[3]-x[1]
    u[2,2] <- x[4]-x[2]

    fet <- fisher.test(u, alternative = "greater")
    ova <- fet$estimate
    p <- fet$p.value
    # ova_ci <- fet$conf.int[1:2]
    return(c(ove, p))
}


# Description:
# Given data from function get_motif_enrichment, this function adds a
# filter column with filter = T if motif is enriched (given the input
# standards) and filter = F if not-enriched
get_motif_filter_v2 <- function(m,
                                min_fdr,
                                min_ove,
                                min_o) {
    m$filter <- FALSE
    m$filter[m$fdr<=min_fdr&m$ove>=min_ove&m$f_sample>=min_o] <- TRUE
    return(m)
}
