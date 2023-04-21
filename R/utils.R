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

