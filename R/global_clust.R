get_global_clust <- function(cdr3, global_max_dist) {
    cdr3_len <- base::nchar(cdr3)
    cdr3_lens <- base::unique(cdr3_len)

    get_hamming_dist <- function(x, cdr3, cdr3_len, global_max_dist) {
        is <- base::which(cdr3_len == x)
        if(base::length(is) == 1) {
            return(NULL)
        }
        if(base::length(is) == 2) {
            d <- stringdist::stringdist(a = cdr3[is[1]],
                                        b = cdr3[is[2]],
                                        method = "hamming")
            if(d > global_max_dist) {
                return(NULL)
            }
            return(c(cdr3[is[1]], cdr3[is[2]]))
        }

        d <- stringdist::stringdistmatrix(
            a = cdr3[is],
            b = cdr3[is],
            method = "hamming"
        )
        d[base::upper.tri(x = d, diag = TRUE)] <- NA
        js <- base::which(d <= global_max_dist, arr.ind = TRUE)
        if(base::nrow(js) == 0) {
            return(NULL)
        }
        return(base::cbind(cdr3[is[js[,1]]], cdr3[is[js[,2]]]))
    }

    hd <- base::lapply(X = cdr3_lens,
                       FUN = get_hamming_dist,
                       cdr3 = cdr3,
                       cdr3_len = cdr3_len,
                       global_max_dist = global_max_dist)
    hd <- base::do.call(base::rbind, hd)
    return(hd)
}



get_global_clust_mem <- function(cdr3, global_max_dist) {

    get_hamming_dist <- function(x, cdr3, cdr3_len, global_max_dist) {

        is <- base::which(cdr3_len == x)
        if(base::length(is) == 1) {
            return(NULL)
        }

        get_pairdist <- function(x, a, len_a, global_max_dist) {
            d <- stringdist::stringdist(a = a[x],
                                        b = a[(x + 1):len_a],
                                        method = "hamming")
            js <- base::which(d <= global_max_dist)
            if(base::length(js) == 0) {
                return(NULL)
            }
            js <- x + js
            return(base::cbind(rep(x = x, times = base::length(js)), js))
        }

        hd <- base::lapply(X = seq_len(length(is) - 1),
                           FUN = get_pairdist,
                           a = cdr3[is],
                           len_a = base::length(is),
                           global_max_dist = global_max_dist)
        hd <- base::do.call(base::rbind, hd)
        if(base::is.null(hd)) {
            return(hd)
        }
        return(cbind(cdr3[is[hd[,1]]], cdr3[is[hd[,2]]]))
    }

    cdr3_len <- base::nchar(cdr3)
    cdr3_lens <- base::unique(cdr3_len)
    hd <- base::lapply(X = cdr3_lens,
                       FUN = get_hamming_dist,
                       cdr3 = cdr3,
                       cdr3_len = cdr3_len,
                       global_max_dist = global_max_dist)
    hd <- base::do.call(base::rbind, hd)
    return(hd)
}
