get_localclust <- function(cdr3,
                           cdr3_ref,
                           ks,
                           control) {
    # 1. trim flanks only relevant for local motifs
    if(control$trim_flank_aa != 0) {
        cdr3_core <- get_trimmed_flanks(x = cdr3,
                                        flank_size = control$trim_flank_aa)
        cdr3_ref_core <- get_trimmed_flanks(x = cdr3_ref,
                                            flank_size = control$trim_flank_aa)
    }
    else{
        cdr3_core <- cdr3
        cdr3_ref_core <- cdr3_ref
    }

    # 2. local clustering: get local motifs
    m <- lapply(X = ks,
                FUN = get_motifs,
                cdr3 = cdr3_core,
                cdr3_ref = cdr3_ref_core,
                min_o = control$local_min_o)
    m <- do.call(rbind, m)
    
    # 3. compute enrichment by fisher's exact test
    ms <- t(apply(X = m[, c("f_s", "f_r", "n_s", "n_r")],
                  MARGIN = 1, FUN = get_motif_enrichment_fet))
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

    # get the unique cdr3s and cdr3_cores
    us <- which(duplicated(cdr3)==FALSE)
    if(length(us) != 0) {
        cdr3 <- cdr3[us]
        cdr3_core <- cdr3_core[us]
    }
    
    # 5. find motifs in input CDR3
    lp <- get_motif_in_seq(cdr3 = cdr3,
                           cdr3_core = cdr3_core,
                           motif = m$motif[m$pass == TRUE])
    
    return(list(m = m, lp = lp))
}



# Description:
# Given a set of sequences and motifs (shorter sequences),
# create a seq->motif map
get_motif_in_seq <- function(cdr3_core,
                             cdr3,
                             motif) {
    # if no enriched motifs
    if(length(motif) == 0) {
        return(NULL)
    }

    find_motif <- function(x, motif, cdr3_core, cdr3) {
        j <- which(regexpr(pattern = motif[x], text = cdr3_core) != -1)
        if(length(j) != 0) {
            return(data.frame(cdr3 = cdr3[j],
                              cdr3_core = cdr3_core[j],
                              motif = motif[x],
                              stringsAsFactors = FALSE))
        }
        return(NULL)
    }

    return(do.call(rbind, lapply(X = seq_len(length(motif)),
                                 motif = motif,
                                 FUN = find_motif,
                                 cdr3 = cdr3,
                                 cdr3_core = cdr3_core)))
}



# Description:
# Computes motif frequencies for a sample and reference
get_motifs <- function(x,
                       cdr3,
                       cdr3_ref,
                       min_o) {
    # find kmers in sample
    kmers_s <- qgrams(cdr3, q = x)
    if(ncol(kmers_s) == 0) {
        # this should not happen -> input checks should catch such errors
        stop("no kmers found in sample")
    }
    if(length(kmers_s) == 0) {
        # this should not happen -> input checks should catch such errors
        stop("no kmers found in sample")
    }
    kmers_s <- kmers_s[1, ]
    
    # find kmers in reference
    kmers_r <- qgrams(cdr3_ref, q = x)
    if(ncol(kmers_r) == 0) {
        # this should not happen -> input checks should catch such errors
        stop("no kmers found in reference")
    }
    kmers_r <- kmers_r[1, ]
    
    # convert table to data.frame
    kmers_s <- data.frame(motif = names(kmers_s),
                          f_s = as.numeric(kmers_s))
    kmers_s$n_s <- sum(kmers_s$f_s)
    kmers_r <- data.frame(motif = names(kmers_r),
                          f_r = as.numeric(kmers_r)
    )
    kmers_r$n_r <- sum(kmers_r$f_r)
    
    # we are only interested in enrichment of motifs in sample relative to
    # reference. Remove all motifs from reference not found in sample.
    # Important: n_r must be sum of all motifs in pop!
    kmers_r <- kmers_r[kmers_r$motif %in% kmers_s$motif,,drop=FALSE]
    
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


