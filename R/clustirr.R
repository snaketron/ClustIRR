clustirr <- function(s,
                     meta = NULL,
                     control = list(blast_gmi = 0.8,
                                    blast_cores = 1,
                                    trim_flank_aa = 3,
                                    db_dist = 0,
                                    db_custom = NULL,
                                    knn = FALSE,
                                    k = 30)) {
    
    # control check
    control <- get_control(control_in = control)
    
    # input check
    input_check(s = s, meta = meta, control = control)
    
    # get chains to be analyzed
    chains <- get_chains(colnames(s))
    
    # add clone size info, sample and clone id
    s <- get_clone_size(s)
    s <- get_sample(s)
    s$id <- paste0(s$sample, "|", seq_len(length.out = nrow(s)))
    
    if(is.null(meta)) {
        meta <- data.frame(id = seq_len(length.out = nrow(s)), 
                           sample = s$sample)
    } else {
        meta$id <- paste0(s$sample, "|", seq_len(length.out = nrow(meta)))
        meta$sample <- s$sample
    }
    
    # run clustirr on individual samples
    message("Step 1: Repertoire analysis \n")
    co <- lapply(X = unique(s$sample), FUN = run_clustirr, s = s, 
                 meta = meta, chains = chains, control = control)
    names(co) <- unique(s$sample)
    
    
    message("Step 2: Building graphs \n")
    if(length(co) == 1) {
        g <- get_graph(clust_irr = co[[1]])
        g$clones <- NULL
        return(g)
    }
    return(get_joint_graph(clust_irrs = co))
}

# Description:
# Run ClustIRR on individual samples
run_clustirr <- function(x, s, control, chains, meta) {
    # annotate s with known CDR3-antigen data 
    message("analysis: ", x, "\n")
    
    # subset by sample
    r <- s[s$sample == x,]
    m <- meta[meta$sample == x,]
    
    r <- match_db(cs = r, control = control)
    
    # do clustering
    clust <- lapply(X = chains, FUN = get_clust, s = r, control = control)
    names(clust) <- chains
    
    # setup clustirr object
    return(get_clustirr_obj(clust = clust, s = r, meta = m, control = control))
}

# Description:
# Wrapper of the main functions performed separately for CDR3b and CDR3a
# (if available)
get_clust <- function(x, s, control) {
    
    r <- s[, c("id", x)]
    colnames(r) <- c("id", "cdr3")
    r <- r[is.na(r$cdr3)==FALSE,]
    if(nrow(r)==0) {
        return(NULL)
    }
    if(is.null(r)) {
        return(NULL)
    }
    
    r$chain <- x
    return(get_score(s = r, control = control))
}
