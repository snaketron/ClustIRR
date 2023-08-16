get_localclust_v1 <- function(  cdr3,
                                cdr3_ref,
                                ks,
                                cores,
                                control) {
    # 1. trim flanks only relevant for local motifs
    if (control$trim_flank_aa != 0) {
        cdr3_core <- get_trimmed_flanks(
            x = cdr3,
            flank_size = control$trim_flank_aa
        )
        cdr3_ref_core <- get_trimmed_flanks(
            x = cdr3_ref,
            flank_size = control$trim_flank_aa
        )
    } else {
        cdr3_core <- cdr3
        cdr3_ref_core <- cdr3_ref
    }

    # 2. local clustering: get local motifs
    motifs <- lapply(
        X = ks,
        FUN = get_motifs_v1,
        cdr3 = cdr3_core,
        cdr3_ref = cdr3_ref_core,
        B = control$B,
        min_o = control$local_min_o,
        cores = cores
    )
    names(motifs) <- as.character(ks)

    # 3. compute local enrichment scores
    future::plan(future::multisession(), workers = cores)
    me <- do.call(rbind, future_lapply(
        X = ks,
        FUN = get_motif_enrichment_boot,
        m = motifs,
        B = control$B,
        cores = cores,
        future.seed = TRUE
    ))
    future::plan(future::sequential())

    # 4. add pass flag
    me$pass <- FALSE
    me$pass[me$fdr <= control$local_max_fdr &
        me$ove >= control$local_min_ove &
        me$f_s >= control$local_min_o] <- TRUE

    # 5. find motifs in input CDR3
    lp <- get_motif_in_seq(
        cdr3 = cdr3,
        cdr3_core = cdr3_core,
        motif = me$motif[me$pass == TRUE]
    )
    return(list(m = me, lp = lp))
}



get_localclust_v23 <- function( cdr3,
                                cdr3_ref,
                                ks,
                                cores,
                                control) {
    # 1. trim flanks only relevant for local motifs
    if (control$trim_flank_aa != 0) {
        cdr3_core <- get_trimmed_flanks(
            x = cdr3,
            flank_size = control$trim_flank_aa
        )
        cdr3_ref_core <- get_trimmed_flanks(
            x = cdr3_ref,
            flank_size = control$trim_flank_aa
        )
    } else {
        cdr3_core <- cdr3
        cdr3_ref_core <- cdr3_ref
    }

    # 2. local clustering: get local motifs
    m <- lapply(
        X = ks,
        FUN = get_motifs_v23,
        cdr3 = cdr3_core,
        cdr3_ref = cdr3_ref_core,
        min_o = control$local_min_o
    )
    m <- do.call(rbind, m)

    # 3. compute enrichment by fisher's exact test
    ms <- t(apply(
        X = m[, c("f_s", "f_r", "n_s", "n_r")],
        MARGIN = 1, FUN = get_motif_enrichment_fet
    ))
    m$ove <- ms[, 1]
    m$ove_ci_l95 <- ms[, 2]
    m$ove_ci_h95 <- ms[, 3]
    m$p_value <- ms[, 4]
    m$fdr <- p.adjust(p = m$p_value, method = "fdr")
    rm(ms)


    # 4. add pass flag -> TRUE if motif passes tests
    m$pass <- FALSE
    m$pass[m$fdr <= control$local_max_fdr &
        m$ove >= control$local_min_ove &
        m$f_s >= control$local_min_o] <- TRUE

    # 5. find motifs in input CDR3
    lp <- get_motif_in_seq(
        cdr3 = cdr3,
        cdr3_core = cdr3_core,
        motif = m$motif[m$pass == TRUE]
    )

    return(list(m = m, lp = lp))
}



# Description:
# Given a set of sequences and motifs (shorter sequences),
# create a seq->motif map
get_motif_in_seq <- function(   cdr3_core,
                                cdr3,
                                motif) {
    # if no enriched motifs
    if (length(motif) == 0) {
        return(NULL)
    }

    find_motif <- function(x, motif, cdr3_core, cdr3) {
        j <-
            which(regexpr(
                pattern = motif[x], text = cdr3_core) != -1)
        if (length(j) != 0) {
            return(data.frame(
                cdr3 = cdr3[j],
                cdr3_core = cdr3_core[j],
                motif = motif[x],
                stringsAsFactors = FALSE
            ))
        }
        return(NULL)
    }

    return(do.call(
        rbind,
        lapply(
            X = seq_len(length(motif)),
            motif = motif,
            FUN = find_motif,
            cdr3 = cdr3,
            cdr3_core = cdr3_core
        )
    ))
}



# Description:
# Computes motif frequencies for a sample and reference
get_motifs_v23 <- function( x,
                            cdr3,
                            cdr3_ref,
                            min_o) {
    # find kmers in sample
    kmers_s <- qgrams(cdr3, q = x)
    if (ncol(kmers_s) == 0) {
        # this should not happen -> input checks should catch such errors
        stop("no kmers found in sample")
    }
    if (length(kmers_s) == 0) {
        # this should not happen -> input checks should catch such errors
        stop("no kmers found in sample")
    }
    kmers_s <- kmers_s[1, ]

    # find kmers in reference
    kmers_r <- qgrams(cdr3_ref, q = x)
    if (ncol(kmers_r) == 0) {
        # this should not happen -> input checks should catch such errors
        stop("no kmers found in reference")
    }
    kmers_r <- kmers_r[1, ]

    # convert table to data.frame
    kmers_s <- data.frame(
        motif = names(kmers_s),
        f_s = as.numeric(kmers_s)
    )
    kmers_s$n_s <- sum(kmers_s$f_s)
    kmers_r <- data.frame(
        motif = names(kmers_r),
        f_r = as.numeric(kmers_r)
    )
    kmers_r$n_r <- sum(kmers_r$f_r)

    # we are only interested in enrichment of motifs in sample relative to
    # reference. Remove all motifs from reference not found in sample.
    # Important: n_r must be sum of all motifs in pop!
    kmers_r <- kmers_r[kmers_r$motif %in% kmers_s$motif, ]

    m <- merge(x = kmers_s, y = kmers_r, by = "motif", all = TRUE)
    m[is.na(m[, "f_s"]), "f_s"] <- 0
    m[is.na(m[, "f_r"]), "f_r"] <- 0
    m[is.na(m[, "n_r"]), "n_r"] <- kmers_r$n_r[1]
    m[is.na(m[, "n_s"]), "n_s"] <- kmers_s$n_s[1]

    m$k <- x

    # return for statistical test
    return(motifs = m)
}



# Description:
# Compute motif enrichment with Fisher's exact test (fet) based on data
# collected by function get_motifs_v2
get_motif_enrichment_fet <- function(x) {
    # Description of parameters used in hypergeometric test (below)
    # f_s = x[1]
    # f_r = x[2]
    # n_s = x[3]
    # n_r = x[4]
    #
    # q = f_s,
    # m = f_r+f_s,
    # n = n_r+n_s-(f_r+f_s),
    # k = n_r+n_s

    # ove <- (x[1] / x[3]) / ((x[2] / x[4]))
    # or use OvE provided by fisher.test

    # q <- x[1]
    # m <- x[1] + x[2]
    # k <- x[3] + x[4]
    # n <- k - m

    u <- matrix(data = 0, nrow = 2, ncol = 2)
    u[1, 1] <- x[1]
    u[2, 1] <- x[2]
    u[1, 2] <- x[3] - x[1]
    u[2, 2] <- x[4] - x[2]

    fet <- fisher.test(u, alternative = "greater", conf.level = 0.95)
    ove <- fet$estimate
    p <- fet$p.value
    ove_ci <- fet$conf.int[c(1, 2)]
    return(c(ove, ove_ci, p))
}




# Description:
# Computes motif frequencies for a sample (cdr3) and reference (cdr3_ref) sets
# of CDR3s. This is done using bootstrapping with B iterations. In each
# iteration, the algorithm draws |cdr3| CDR3s without replacement from
# cdr3_ref, finds the motifs, and computes motif frequencies. This is repeated
# B times to generate a vector of motif frequencies. The mean (expected) motif
# frequency is compared against the observed motif frequencies in sample cdr3.
get_motifs_v1 <- function(cdr3, cdr3_ref, B, ks, cores, min_o) {
    # x = k in k-mer
    get_kmers_freq_ref <- function(x, cdr3, N, B, relevant_motifs, cores) {
        get_qgrams <- function(x, q, cdr3, N, relevant_motifs) {
            draw_cdr3 <- sample(x = cdr3, size = N, replace = FALSE)
            o <- qgrams(draw_cdr3, q = q)
            if (ncol(o) == 0) {
                return(NA)
            }
            o <- o[1, ]
            o <- o[names(o) %in% relevant_motifs]
            if (length(o) == 0) {
                return(NA)
            }
            return(o)
        }

        future::plan(future::multisession(), workers = cores)
        o <- future_lapply(
            X = seq_len(B),
            q = x,
            N = N,
            relevant_motifs = relevant_motifs,
            cdr3 = cdr3,
            FUN = get_qgrams,
            future.seed = TRUE
        )
        future::plan(future::sequential())
        return(o)
    }

    # x = k in k-mer
    get_kmers_freq_sample <- function(x, cdr3, min_o) {
        o <- qgrams(cdr3, q = x)
        if (ncol(o) == 0) {
            return(NA)
        }
        o <- o[1, ]
        # o <- o[1, o[1, ] >= min_o]
        if (length(o) == 0) {
            return(NA)
        }
        return(o)
    }

    # find motifs in sample
    motif_s <- lapply(
        X = ks,
        FUN = get_kmers_freq_sample,
        cdr3 = cdr3,
        min_o = min_o
    )
    names(motif_s) <- ks
    found_kmers <- as.vector(unlist(
        lapply(X = motif_s, FUN = names)
    ))
    motif_ref <- lapply(
        X = ks,
        FUN = get_kmers_freq_ref,
        cdr3 = cdr3_ref,
        B = B,
        N = length(cdr3),
        relevant_motifs = found_kmers,
        cores = cores
    )
    names(motif_ref) <- ks
    return(list(motif_s = motif_s, motif_ref = motif_ref))
}


# Description:
# Computes motif enrichment using bootstrapping (cluster v1)
get_motif_enrichment_boot <- function(  x,
                                        m,
                                        B,
                                        cores) {
    get_e <- function(x, B) {
        # o
        o <- x[1]
        # e
        e <- x[-1]
        # mean, max
        e_mean <- mean(e)
        e_min <- min(e)
        e_max <- max(e)

        # OvE = o/e
        ove <- o / e_mean
        # prob e>=o
        p_value <- sum(e >= o) / B
        # return
        return(c(e_mean, e_min, e_max, o, ove, p_value))
    }

    motif_s <- m[[as.character(x)]]$motif_s[[1]]
    motif_ref <- m[[as.character(x)]]$motif_ref[[1]]

    # matrix of k-mer counts
    f_m <- matrix(data = 0, ncol = B + 1, nrow = length(motif_s))
    rownames(f_m) <- names(motif_s)
    f_m[names(motif_s), 1] <- motif_s
    for (i in seq_len(length(motif_ref))) {
        s <- motif_ref[[i]]
        f_m[names(s), i + 1] <- s
    }

    e <- t(apply(X = f_m, MARGIN = 1, FUN = get_e, B = B))

    # format output
    e <- data.frame(e)
    colnames(e) <- c(
        "mean_f_r", "min_f_r", "max_f_r",
        "f_s", "ove", "p_value"
    )
    e$motif <- rownames(e)
    e$fdr <- p.adjust(p = e$p_value, method = "fdr")
    e$k <- x
    e[, c(
        "motif", "k", "f_s", "mean_f_r", "min_f_r", "max_f_r",
        "ove", "p_value", "fdr"
    )]
    return(e)
}
