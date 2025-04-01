clustirr <- function(s,
                     meta = NULL,
                     cores = 1,
                     control = list(gmi = 0.7,
                                    trim_flank_aa = 3,
                                    db_dist = 0,
                                    db_custom = NULL)) {
    # control check
    control <- get_control(control_in = control)
    
    # input check
    input_check(s = s, meta = meta, cores = cores, control = control)
    
    # get chains to be analyzed
    chains <- get_chains(colnames(s))
    
    # add clone size info, sample and clone id
    s <- get_clone_size(s)
    s <- get_sample(s)
    s$id <- seq_len(length.out = nrow(s))
    
    if(is.null(meta)) {
        meta <- data.frame(id = seq_len(length.out = nrow(s)), 
                           sample = s$sample)
    } else {
        meta$id <- seq_len(length.out = nrow(meta))
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
    return(get_joint_graph(clust_irrs = co, cores = cores))
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
    
    get_cdr3s <- function(x, chain) {
        # multiple cdr3 x clone_count -> cells
        # in the global analysis this is reverted
        x <- x[, chain]
        x <- x[is.na(x)==FALSE]
        if(length(x)==0) {
            stop("no CDR3s found")
        }
        return(x)
    }
    
    cdr3 <- get_cdr3s(x = s, chain = x)
    
    return(get_blosum(cdr3 = cdr3, control = control))
}

get_blosum <- function(cdr3, control) {
    
    get_bscore <- function(x, o, b, gap_o, gap_e, trim) {
        
        parse_cigar <- function(cigar) {
            c <- as.numeric(gregexpr(text = cigar, pattern = "X|D|I|\\=")[[1]])
            s <- vapply(X = c, cigar = cigar, FUN.VALUE = character(1), 
                        FUN = function(x, cigar) {
                            return(substr(x = cigar, start = x, stop = x))
                        })
            n <- vapply(X = 2:(length(c)+1), c = c(0,c), cigar = cigar, 
                        FUN.VALUE = character(1), 
                        FUN = function(x, c, cigar) {
                            return(substr(x = cigar, start = c[x-1]+1, stop = c[x]-1))
                        })
            n <- as.numeric(n)
            
            return(rep(x = s, times = n))
        }
        
        q <- o$QueryMatchSeq[x] 
        t <- o$TargetMatchSeq[x] 
        cigar <- o$Alignment[x] 
        
        s <- parse_cigar(cigar = cigar)
        len_s <- length(s)
        q <- unlist(strsplit(x = q, split = NULL, fixed = TRUE))
        t <- unlist(strsplit(x = t, split = NULL, fixed = TRUE))
        
        # trim region = 1, else = 0
        tr <- numeric(length = len_s)
        if(trim != 0) {
            if(len_s-2*trim <= 0) {
                tr[1:length(tr)] <- 1
            } else {
                tr[c(1:trim, (len_s-trim+1):len_s)] <- 1
            }
        }
        
        cscore <- 0
        score <- 0
        gap_o_on <- 0
        iq <- 1
        it <- 1
        for(i in seq_len(len_s)) {
            if(s[i]=="="|s[i]=="X") {
                bi <- b[q[iq],t[it]]
                score <- score + bi
                cscore <- cscore + bi * (tr[i]==0)
                gap_o_on <- 0
                iq <- iq + 1
                it <- it + 1
            }
            else {
                if(gap_o_on==1) {
                    score <- score + gap_e
                    cscore <- cscore + gap_e * (tr[i]==0)
                } 
                else {
                    score <- score + gap_o + gap_e
                    cscore <- cscore + (gap_o + gap_e) * (tr[i]==0)
                    gap_o_on <- 1
                }
                
                iq <- iq + (s[i]=="I")
                it <- it + (s[i]=="D")
            }
        }
        
        res <- numeric(length = 4)
        res[1] <- score
        res[2] <- len_s
        res[3] <- cscore
        res[4] <- sum(tr==0)
        return(res)
    }
    
    get_bscore_dup <- function(x, b, trim) {
        x <- unlist(strsplit(x = x, split = NULL, fixed = TRUE))
        len_x <- length(x)
        len_cx <- len_x
        
        if(len_x==1) {
            score <- sum(b[x,x]) 
        } else {
            score <- sum(diag(b[x,x]))
        }
        
        if(len_x <= trim*2) {
            cscore <- 0
        } else {
            if(trim==0) {
                cscore <- score
            } 
            else {
                cx <- x[(trim+1):(len_x-trim)]
                len_cx <- length(cx)
                if(len_cx==1) {
                    cscore <- sum(b[cx,cx]) 
                } 
                else {
                    cscore <- sum(diag(b[cx,cx]))
                }
            }
        }
        
        res <- numeric(length = 4)
        res[1] <- score
        res[2] <- len_x
        res[3] <- cscore
        res[4] <- len_cx
        return(res)
    }
    
    # get duplicates
    cdr3_dup <- unique(cdr3[duplicated(cdr3)==TRUE])
    
    # remove duplicates
    cdr3 <- cdr3[duplicated(cdr3)==FALSE]
    
    # create db which is also query
    db <- data.frame(Id = 1:length(cdr3), Seq = cdr3, len = nchar(cdr3))
    
    # blast
    o <- blast(query = db, 
               db = db,
               maxAccepts = 10^4,
               maxRejects = 10^3,
               minIdentity = control$gmi,
               alphabet = "protein", 
               output_to_file = FALSE)
    
    # remove self-hits
    o <- o[o$QueryId!=o$TargetId,]
    # if empty stop
    if(nrow(o)==0) {
        return(NULL)
    }
    
    # remove partial hits
    o$TargetLen <- db$len[o$TargetId]
    o$QueryLen <- db$len[o$QueryId]
    j <- which(o$QueryMatchStart!=1 | 
                   o$TargetMatchStart!=1 | 
                   o$QueryMatchEnd != o$QueryLen | 
                   o$TargetMatchEnd != o$TargetLen)
    if(length(j)!=0) {
        o <- o[-j,]
    }
    # if empty stop
    if(nrow(o)==0) {
        return(NULL)
    }
    
    # remove duplicates (there should be none by this point)
    key <- apply(X = o[,c("QueryId", "TargetId")], MARGIN = 1, 
                 FUN = function(x) {paste0(sort(x), collapse = '-')})
    j <- which(duplicated(key)==FALSE)
    if(length(j)!=0) {
        o <- o[j,,drop=FALSE]
    }
    
    # get blosum matrix
    data_env <- get_blosum62()
    
    # compute BLSOUM62 scores
    bs <- t(vapply(X = 1:nrow(o), 
                   FUN.VALUE = numeric(4),
                   FUN = get_bscore, 
                   o = o, 
                   gap_o = -10, 
                   gap_e = -4, 
                   trim = control$trim_flank_aa,
                   b = data_env[["BLOSUM62"]]))
    
    out <- data.frame(from_cdr3 = db$Seq[o$QueryId],
                      to_cdr3 = db$Seq[o$TargetId],
                      weight = bs[,1],
                      cweight = bs[,3],
                      nweight = bs[,1]/bs[,2],
                      ncweight = bs[,3]/bs[,4],
                      max_len = bs[,2],
                      max_clen = bs[,4])
    
    # compute BLSOUM62 scores between duplicates -> more efficient
    if(length(cdr3_dup)!=0) {
        bs_dup <- t(vapply(X = cdr3_dup, 
                           FUN = get_bscore_dup, 
                           trim = control$trim_flank_aa,
                           b = data_env[["BLOSUM62"]],
                           FUN.VALUE = numeric(4)))
        
        o_dup <- data.frame(from_cdr3 = cdr3_dup,
                            to_cdr3 = cdr3_dup,
                            weight = bs_dup[,1],
                            cweight = bs_dup[,3],
                            nweight = bs_dup[,1]/bs_dup[,2],
                            ncweight = bs_dup[,3]/bs_dup[,4],
                            max_len = bs_dup[,2],
                            max_clen = bs_dup[,4])
        out <- rbind(o_dup, out)
    }
    
    return(out)
}
