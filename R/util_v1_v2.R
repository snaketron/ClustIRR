# Description:
# Setup control list.
# control_in: user generated list (if missing -> use default)
get_control <- function(control_in) {

    control <- list(
        B = 1000,
        global_max_dist = 1,
        local_min_fdr = 0.05,
        local_min_ove = 2,
        local_min_o = 3,
        trim_flanks = FALSE,
        flank_size = 3,
        global_pairs = NULL,
        low_mem = FALSE)

    # if missing control_in -> use default values
    if(missing(control_in)|is.null(control_in)) {
        return(control)
    }

    # edit control by user-defined control_in
    # There are packages to update list given another list, but here we can
    # live with the following "inefficiency" as the list is generally small
    # (~5 elements)
    ns <- names(control_in)
    for(i in 1:length(control_in)) {
        control[[ns[i]]] <- control_in[[ns[i]]]
    }
    return(control)
}


# Description:
# Look for enriched motifs in CDR3s. CDR3 sequence pairs of local
# clusters are returned. Used by gliph_v1 and gliph_v2.
get_local_pair <- function(cdr3,
                           motif) {
    # if no enriched motifs
    if(length(motif)==0) {
        return(NULL)
    }

    get_motif_in_cdr <- function(x, motif, cdr3) {
        j <- which(regexpr(pattern = motif[x], text = cdr3)!=-1)
        if(length(j)==1) {
            return(data.frame(from = j, to = j, motif = motif[x]))
        }
        j <- unique(j)
        u <- t(utils::combn(x = j, m = 2))
        u <- rbind(u, cbind(j,j))
        return(data.frame(from = u[,1], to = u[,2], motif = motif[x]))
    }

    return(do.call(rbind, lapply(X = 1:length(motif),
                                 motif = motif,
                                 FUN = get_motif_in_cdr,
                                 cdr3 = cdr3)))
}


# Description:
# Look for pairs of global connections. Used by gliph_v1 and gliph_v2.
get_global_pairs <- function(cdr3,
                             global_max_dist) {
    # Jan is looping over all sequences and computing distances with the rest.
    # two problems:
    # a) slower than passing vector to stringdist (but lower max. memory footprint)
    # b) hamming distance between sequences with unequal lengths computed nontheless

    # at the very least loop over cdr3.lengths -> faster, but somewhat more
    # memory is data passed together to stringdist

    cdr3_len <- base::nchar(cdr3)
    cdr3_lens <- unique(cdr3_len)

    get_hamming_dist <- function(x, cdr3, cdr3_len, global_max_dist) {
        is <- which(cdr3_len == x)
        if(length(is)==1) {
            return(NULL)
        }
        if(length(is)==2) {
            d <- stringdist::stringdist(
                a = cdr3[is[1]],
                b = cdr3[is[2]],
                method = "hamming")
            if(d>global_max_dist) {
                return(NULL)
            }
            return(c(cdr3[is[1]], cdr3[is[2]]))
        }

        d <- stringdist::stringdistmatrix(
            a = cdr3[is],
            b = cdr3[is],
            method = "hamming")
        d[upper.tri(x = d, diag = TRUE)] <- NA
        # d[1:nrow(d), 1:nrow(d)] <- NA
        js <- which(d<=global_max_dist, arr.ind = TRUE)
        if(nrow(js)==0) {
            return(NULL)
        }
        return(cbind(is[js[,1]], is[js[,2]]))
    }

    # Description:
    # same as get_hamming_dist but slower (x3), however has much
    # smaller memory footprint -> appropriate for large input
    # Similar but faster implementation than that in Jan's code
    get_hamming_dist_memory <- function(x,
                                        cdr3,
                                        cdr3_len,
                                        global_max_dist) {
        is <- which(cdr3_len == x)
        if(length(is)==1) {
            return(NULL)
        }

        get_pairdist  <- function(x, a, len_a, global_max_dist) {
            d <- stringdist::stringdist(a = a[x],
                                        b = a[(x+1):len_a],
                                        method = "hamming")
            js <- which(d<=global_max_dist)
            if(length(js)==0) {
                return(NULL)
            }
            js <- x+js
            return(cbind(rep(x = x, times = length(js)), js))
        }

        hd <- lapply(X = 1:(length(is)-1),
                     FUN = get_pairdist,
                     a = cdr3[is],
                     len_a = length(is),
                     global_max_dist = global_max_dist)
        hd <- do.call(rbind, hd)
        if(is.null(hd)) {
            return(hd)
        }
        # map to original indices
        return(cbind(is[hd[,1]], is[hd[,2]]))
    }


    hd <- lapply(X = cdr3_lens,
                 FUN = get_hamming_dist,
                 cdr3 = cdr3,
                 cdr3_len = cdr3_len,
                 global_max_dist = global_max_dist)
    hd <- do.call(rbind, hd)
    return(hd)
}



# Description:
# Look for global connections. Low-memory mode. Used by gliph_v1 and gliph_v2
get_global_pairs_mem <- function(cdr3,
                                 global_max_dist) {

    cdr3_len <- base::nchar(cdr3)
    cdr3_lens <- unique(cdr3_len)

    get_hamming_dist <- function(x, cdr3, cdr3_len, global_max_dist) {
        is <- which(cdr3_len == x)
        if(length(is)==1) {
            return(NULL)
        }
        if(length(is)==2) {
            d <- stringdist::stringdist(
                a = cdr3[is[1]],
                b = cdr3[is[2]],
                method = "hamming")
            if(d>global_max_dist) {
                return(NULL)
            }
            return(c(cdr3[is[1]], cdr3[is[2]]))
        }

        d <- stringdist::stringdistmatrix(
            a = cdr3[is],
            b = cdr3[is],
            method = "hamming")
        d[upper.tri(x = d, diag = TRUE)] <- NA
        # d[1:nrow(d), 1:nrow(d)] <- NA
        js <- which(d<=global_max_dist, arr.ind = TRUE)
        if(nrow(js)==0) {
            return(NULL)
        }
        return(cbind(is[js[,1]], is[js[,2]]))
    }

    # Description:
    # same as get_hamming_dist from get_global_pairs, but slower. However it
    # has much smaller memory footprint -> appropriate for large input sets.
    # Similar but faster implementation than that in Jan's code
    get_hamming_dist <- function(x,
                                 cdr3,
                                 cdr3_len,
                                 global_max_dist) {
        is <- which(cdr3_len == x)
        if(length(is)==1) {
            return(NULL)
        }

        get_pairdist  <- function(x, a, len_a, global_max_dist) {
            d <- stringdist::stringdist(a = a[x],
                                        b = a[(x+1):len_a],
                                        method = "hamming")
            js <- which(d<=global_max_dist)
            if(length(js)==0) {
                return(NULL)
            }
            js <- x+js
            return(cbind(rep(x = x, times = length(js)), js))
        }

        hd <- lapply(X = 1:(length(is)-1),
                     FUN = get_pairdist,
                     a = cdr3[is],
                     len_a = length(is),
                     global_max_dist = global_max_dist)
        hd <- do.call(rbind, hd)
        if(is.null(hd)) {
            return(hd)
        }
        # map to original indices
        return(cbind(is[hd[,1]], is[hd[,2]]))
    }


    hd <- lapply(X = cdr3_lens,
                 FUN = get_hamming_dist,
                 cdr3 = cdr3,
                 cdr3_len = cdr3_len,
                 global_max_dist = global_max_dist)
    hd <- do.call(rbind, hd)
    return(hd)
}



# Description:
# Given input data, get the number of chains that ought to be analyzed
# one by one.
# x = columns(data_sample)
get_chains <- function(x) {
    js <- c(which(regexpr(pattern = "CDR3b", text = x)!=-1),
           which(regexpr(pattern = "CDR3a", text = x)!=-1))
    return(x[js])
}



# Description:
# cut the left/right flanks of each CDR3 sequence by flank_size amino
# acids
get_trimmed_flanks <- function(x,
                               flank_size) {
    x <- base::substr(x=x,
                      start=flank_size+1 ,
                      stop=base::nchar(x)-flank_size)
    x[x==""] <- NA
    return(x)
}


# Description:
# Given two sets of input:
# 1) from clustering: local_pairs and global_pairs
# 2) from data_sample: data.frame + chain (TCRa==CDR3a or TCRb==CDR3b)
get_edges_todo <- function(local_pairs,
                      global_pairs,
                      data_sample,
                      chain) {

    edges <- c()
    if(nrow(local_pairs)!=0) {
        local_pairs$edge_type <- "local"
        edges <- rbind(edges, local_pairs)
    }
    if(nrow(global_pairs)!=0) {
        global_pairs <- data.frame(
            from = global_pairs[, 1],
            to = global_pairs[, 2],
            motif = NA,
            edge_type = "global")
        edges <- rbind(edges, global_pairs)
    }
    if(nrow(local_pairs)==0 &
       nrow(global_pairs)==0) {
        return(NULL)
    }

    if(chain == "CDR3b") {
        info_from <- data_sample[edges$from, c("CDR3b", "TRBV", "TRBJ")]
        colnames(info_from) <- paste0("from_", gsub(pattern = "TRB",
                                                    replacement = "TR",
                                                    colnames(info_from)))
        info_to <- data_sample[edges$to, c("CDR3b", "TRBV", "TRBJ")]
        colnames(info_to) <- paste0("to_", gsub(pattern = "TRB",
                                                replacement = "TR",
                                                colnames(info_to)))

        info_from <- cbind(info_from, info_to)
        edges$chain <- "CDR3b"
        edges <- cbind(edges, info_from)
    }
    if(chain == "CDR3a") {
        info_from <- data_sample[edges$from, c("CDR3a", "TRAV", "TRAJ")]
        colnames(info_from) <- paste0("from_", gsub(pattern = "TRA",
                                                    replacement = "TR",
                                                    colnames(info_from)))
        info_to <- data_sample[edges$to, c("CDR3", "TRAV", "TRAJ")]
        colnames(info_to) <- paste0("to_", gsub(pattern = "TRA",
                                                replacement = "TR",
                                                colnames(info_to)))

        info_from <- cbind(info_from, info_to)
        edges$chain <- "CDR3a"
        edges <- cbind(edges, info_from)
    }
    return(edges)
}


get_edges <- function(local_pairs,
                      global_pairs,
                      cdr3,
                      chain) {

    edges <- c()
    if(is.null(local_pairs)==FALSE) {
        local_pairs$edge_type <- "local"
        edges <- rbind(edges, local_pairs)
    }
    if(is.null(global_pairs)==FALSE) {
        global_pairs <- data.frame(
            from = global_pairs[,1],
            to = global_pairs[,2],
            motif = NA,
            edge_type = "global")
        edges <- rbind(edges, global_pairs)
    }
    if(is.null(local_pairs)&is.null(global_pairs)) {
        return(NULL)
    }

    if(chain == "CDR3b") {
        edges$from_cdr3 <- cdr3[edges$from]
        edges$to_cdr3 <- cdr3[edges$to]
        edges$chain <- "CDR3b"
    }
    if(chain == "CDR3a") {
        edges$from_cdr3 <- cdr3[edges$from]
        edges$to_cdr3 <- cdr3[edges$to]
        edges$chain <- "CDR3a"
    }
    return(edges)
}

