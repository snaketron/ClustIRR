
# Description:
# Look for pairs of global connections.
# Used by gliph_v1 and gliph_v2.
get_global_clust <- function(cdr3,
                             global_max_dist) {
    # Jan is looping over all sequences and computing distances with the rest.
    # two problems:
    # a) slower than passing vector to stringdist,
    #    but lower max. memory footprint
    # b) hamming distance between sequences
    #    with unequal lengths computed nonetheless

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

    hd <- lapply(X = cdr3_lens,
                 FUN = get_hamming_dist,
                 cdr3 = cdr3,
                 cdr3_len = cdr3_len,
                 global_max_dist = global_max_dist)
    hd <- do.call(rbind, hd)
    return(hd)
}



# Description:
# Look for global connections. Low-memory mode.
# Used by gliph_v1 and gliph_v2
get_global_clust_mem <- function(cdr3,
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
