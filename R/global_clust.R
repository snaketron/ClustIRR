
get_global_clust_hamming <- function(cdr3, cdr3_dup, control) {
    
    get_pairdist <- function(x, a, len_a, global_max_hdist) {
        d <- stringdist(a = a[x], b = a[(x + 1):len_a], method = "hamming")
        js <- which(d <= global_max_hdist)
        if(length(js) == 0) {
            return(NULL)
        }
        js <- x + js
        return(cbind(rep(x = x, times = length(js)), js))
    }
    
    get_hd <- function(x, cdr3, cdr3_len, global_max_hdist, low_mem) {
        
        is <- which(cdr3_len == x)
        if(length(is) == 1) {
            return(NULL)
        }
        if(length(is) == 2) {
            d <- stringdist(a = cdr3[is[1]], b = cdr3[is[2]], method="hamming")
            if(d > global_max_hdist) {
                return(NULL)
            }
            
            return(data.frame(from_cdr3 = cdr3[is[1]],
                              to_cdr3 = cdr3[is[2]],
                              weight = 1,
                              cweight = 1,
                              nweight = 1,
                              ncweight = 1,
                              max_len = NA))
        }
        
        if(low_mem) {
            hd <- lapply(X = seq_len(length(is) - 1),
                         FUN = get_pairdist,
                         a = cdr3[is],
                         len_a = length(is),
                         global_max_hdist = global_max_hdist)
            hd <- do.call(rbind, hd)
            if(is.null(hd)) {
                return(hd)
            }
            
            return(data.frame(from_cdr3 = cdr3[is[hd[, 1]]],
                              to_cdr3 = cdr3[is[hd[, 2]]],
                              weight = 1,
                              cweight = 1,
                              nweight = 1,
                              ncweight = 1,
                              max_len = NA))
        }
        else {
            d <- stringdistmatrix(a = cdr3[is], b = cdr3[is], method="hamming")
            d[upper.tri(x = d, diag = TRUE)] <- NA
            js <- which(d <= global_max_hdist, arr.ind = TRUE)
            if(nrow(js) == 0) {
                return(NULL)
            }
            
            return(data.frame(from_cdr3 = cdr3[is[js[, 1]]],
                              to_cdr3 = cdr3[is[js[, 2]]],
                              weight = 1,
                              cweight = 1,
                              nweight = 1,
                              ncweight = 1,
                              max_len = NA))
        }
    }
    
    cdr3_len <- nchar(cdr3)
    cdr3_lens <- unique(cdr3_len)
    
    hd <- lapply(X = cdr3_lens,
                 FUN = get_hd,
                 cdr3 = cdr3,
                 cdr3_len = cdr3_len,
                 global_max_hdist = control$global_max_hdist,
                 low_mem = control$low_mem)
    hd <- do.call(rbind, hd)
    
    # if there are duplicated CDR3s add one entry
    if(any(cdr3_dup==1)) {
        q <- cdr3[which(cdr3_dup==1)]
        hd_dup <- data.frame(from_cdr3 = q, 
                             to_cdr3 = q,
                             weight = 1,
                             cweight = 1,
                             nweight = 1,
                             ncweight = 1,
                             max_len = NA)
        
        hd <- rbind(hd, hd_dup)
    }
    return(hd)
}


get_global_clust_blosum <- function(cdr3, cdr3_dup, control) {
    
    # computes BLOSUM62 score for pairs of sequences returned by blaster 
    get_bscore <- function(x, d, bm, db) {
        a <- db$Seq[d$QueryId[x]]
        b <- db$Seq[d$TargetId[x]]
        
        if(is.na(a)|is.na(b)) {
            return(NA)
        }
        
        return(stringDist(x = c(a, b),
                          method = "substitutionMatrix", 
                          type = "global", 
                          substitutionMatrix = bm, 
                          gapOpening = 10,
                          gapExtension = 4))
    }
    
    get_bscore_trim <- function(x, d, bm, db, trim) {
        
        a <- db$Seq[d$QueryId[x]]
        b <- db$Seq[d$TargetId[x]]
        na <- nchar(a)
        nb <- nchar(b)
        if((na-2*trim)<=0 | (nb-2*trim)<=0) {
            return(NA)
        }
        
        a <- substr(x = a, start = trim+1, stop = nchar(a)-trim)
        b <- substr(x = b, start = trim+1, stop = nchar(b)-trim)
        
        if(is.na(a)|is.na(b)) {
            return(NA)
        }
        
        return(stringDist(x = c(a, b),
                          method = "substitutionMatrix", 
                          type = "global", 
                          substitutionMatrix = bm, 
                          gapOpening = 10,
                          gapExtension = 4))
    }
    
    get_bscore_dup <- function(x, cdr3, bm, trim) {
        
        if(is.na(cdr3[x])) {
            return(data.frame(from_cdr3 = cdr3[x],
                              to_cdr3 = cdr3[x],
                              weight = NA,
                              cweight = NA))
        } 
        
        bs <- stringDist(x = c(cdr3[x], cdr3[x]),
                         method = "substitutionMatrix", 
                         type = "global", 
                         substitutionMatrix = bm, 
                         gapOpening = 10,
                         gapExtension = 4)
        
        a <- cdr3[x]
        na <- nchar(a)
        if((na-2*trim)<=0) {
            bs_core <- NA
        } 
        else {
            cdr3_core <- substr(x = a, start = trim+1, stop = na-trim)
            if(is.na(cdr3_core)) {
                bs_core <- NA
            } 
            else {
                bs_core <- stringDist(x = c(cdr3_core, cdr3_core),
                                      method = "substitutionMatrix", 
                                      type = "global", 
                                      substitutionMatrix = bm, 
                                      gapOpening = 10,
                                      gapExtension = 4)
            }
        }
        return(data.frame(from_cdr3 = cdr3[x],
                          to_cdr3 = cdr3[x],
                          weight = -bs,
                          cweight = -bs_core))
    }
    
    db <- data.frame(Id = 1:length(cdr3), Seq = cdr3)
    
    # blast
    o <- blast(query = db, 
               db = db,
               maxAccepts = 10^4,
               minIdentity = control$global_min_identity,
               alphabet = "protein", 
               output_to_file = FALSE)
    
    # remove self-hits
    o <- o[o$QueryId!=o$TargetId,]
    
    # if empty stop
    if(nrow(o)==0) {
        return(NULL)
    }
    
    key <- apply(X = o[,c("QueryId", "TargetId")], MARGIN = 1, 
                 FUN = function(x) {paste0(sort(x), collapse = '-')})
    key_js <- which(duplicated(key)==FALSE)
    if(length(key_js)!=0) {
        o <- o[key_js,,drop=FALSE]
    }
    
    # get blosum matrix from pwalign
    data_env <- new.env(parent = emptyenv())
    # data("BLOSUM62", envir = data_env, package = "pwalign")
    data("BLOSUM62", envir = data_env, package = "Biostrings")
    
    # compute BLSOUM62 score for matches
    o$bs <- vapply(X = 1:nrow(o), 
                   FUN = get_bscore, 
                   d = o, 
                   db = db, 
                   bm = data_env[["BLOSUM62"]],
                   FUN.VALUE = numeric(1))
    
    # compute BLSOUM62 score for matches
    if(control$trim_flank_aa == 0) {
        o$core_bs <- o$bs
    } 
    else {
        o$core_bs <- vapply(X = 1:nrow(o), 
                            FUN = get_bscore_trim, 
                            d = o, 
                            db = db, 
                            bm = data_env[["BLOSUM62"]],
                            trim = control$trim_flank_aa,
                            FUN.VALUE = numeric(1))
    }
    
    out <- data.frame(from_cdr3 = db$Seq[o$QueryId],
                      to_cdr3 = db$Seq[o$TargetId],
                      weight = -o$bs,
                      cweight = -o$core_bs)
    
    # if there are duplicated CDR3s add one entry
    if(any(cdr3_dup==1)) {
        q <- cdr3[which(cdr3_dup==1)]
        out_dup <- do.call(rbind, lapply(X = 1:length(q),
                                         FUN = get_bscore_dup,
                                         cdr3 = q, 
                                         bm = data_env[["BLOSUM62"]], 
                                         trim = control$trim_flank_aa))
        
        out <- rbind(out, out_dup)
    }
    
    out$max_len <- apply(X = out[, c("from_cdr3", "to_cdr3")], MARGIN = 1,
                         FUN = function(x) {return(max(nchar(x)))})
    out$nweight <- out$weight/out$max_len
    out$ncweight <- out$cweight/out$max_len
    return(out)
}


get_global_clust <- function(cdr3, control) {
    
    cdr3 <- table(cdr3)
    cdr3_dup <- ifelse(test = as.numeric(cdr3)>1, yes = 1, no = 0)
    cdr3 <- names(cdr3)
    
    # if global_pairs are provided as input use them, else compute them
    if(!is.null(control$global_pairs)) {
        g <- control$global_pairs
    }
    else {
        if(control$global_smart==TRUE) {
            g <- get_global_clust_blosum(cdr3 = cdr3, 
                                         cdr3_dup = cdr3_dup, 
                                         control = control)
        }
        if(control$global_smart==FALSE) {
            g <- get_global_clust_hamming(cdr3 = cdr3, 
                                          cdr3_dup = cdr3_dup,
                                          control = control)
        } 
    }
    return(g)
}
