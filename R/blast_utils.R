
get_bscore <- function(x, o, b, gap_o, gap_e, trim) {
    q <- o$qseq[x] 
    t <- o$sseq[x]
    len_s <- o$length[x]
    q <- unlist(strsplit(x = q, split = NULL, fixed = TRUE))
    t <- unlist(strsplit(x = t, split = NULL, fixed = TRUE))
    
    # trim region = 1, else = 0
    tr <- numeric(length = len_s)
    if(trim != 0) {
        if(len_s-2*trim <= 0) {
            tr[seq_len(length(tr))] <- 1
        } else {
            tr[c(seq_len(trim), (len_s-trim+1):len_s)] <- 1
        }
    }
    
    cscore <- 0
    score <- 0
    gap_o_on <- 0
    iq <- 1
    it <- 1
    for(i in seq_len(len_s)) {
        if(q[i]!="-" & t[i]!= "-") {
            bi <- b[q[i],t[i]]
            score <- score + bi
            cscore <- cscore + bi * (tr[i]==0)
            gap_o_on <- 0
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
        }
    }
    
    res <- numeric(length = 3)
    res[1] <- score
    res[2] <- cscore
    res[3] <- sum(tr==0)
    return(res)
}

get_score <- function(s, control) {
    
    # rblast db
    db <- AAStringSet(s$cdr3)
    names(db) <- s$id
    writeXStringSet(db, filepath = "tmp.fasta")
    makeblastdb(db_name = "db/tmp", dbtype = "prot", 
                file = "tmp.fasta", verbose = FALSE)
    bdb <- blast("db/tmp", type = "blastp")
    blast_fmt <- paste("qseqid sseqid pident length qstart qend",
                       "sstart send evalue score qseq sseq qlen slen")
    blast_args <- paste0("-num_threads ", control$blast_cores, 
                         " -evalue 10 -gapopen 11 -gapextend 2", 
                         " -word_size 3 -matrix BLOSUM62 -comp_based_stats 0")
    o <- predict(bdb, db, BLAST_args = blast_args, custom_format = blast_fmt)
    
    # clean up
    unlink("tmp.fasta")
    unlink("db/tmp", recursive = TRUE)
    
    # remove self-hits
    o <- o[o$qseqid!=o$sseqid,]
    # if empty and no duplicates found stop
    if(nrow(o)==0) {
        return(NULL)
    }
    
    # remove low gmi matches
    o <- o[o$pident >= control$blast_gmi*100,]
    if(nrow(o)==0) {
        return(NULL)
    }
    
    # remove partial hits
    j <- which(o$qstart!=1 | o$sstart!=1 | o$qend != o$qlen | o$send != o$slen)
    if(length(j)!=0) {
        o <- o[-j,]
    }
    # if empty and no duplicates found stop
    if(nrow(o)==0) {
        return(NULL)
    }
    
    # remove self-hits
    key <- paste(pmin(o$qseqid, o$sseqid), pmax(o$qseqid, o$sseqid), 
                 sep = "-")
    o <- o[!duplicated(key), , drop = FALSE]
    if(nrow(o)>0) {
        out <- data.frame(from = o$qseqid,
                          to = o$sseqid,
                          len = o$length,
                          weight = o$score,
                          nweight = o$score/o$length)
        
        data_env <- get_blosum62()
        bs <- t(vapply(X = seq_len(nrow(o)), 
                       FUN.VALUE = numeric(3),
                       FUN = get_bscore, 
                       o = o, 
                       gap_o = -11, 
                       gap_e = -2, 
                       trim = control$trim_flank_aa,
                       b = data_env[["BLOSUM62"]]))
        out$cweight <- bs[,2]
        out$clen <- bs[,3]
        
        out$nweight <- out$weight/out$len
        out$ncweight <- out$cweight/out$clen
        out$chain <- s$chain[1]
    }
    
    # KNN 
    if(control$knn == TRUE) {
        out <- out[order(out$nweight, decreasing = TRUE), ]
        out <- out %>%
            group_by(from) %>%
            arrange(desc(nweight), .by_group = TRUE) %>%
            mutate(rank = row_number()) %>%
            ungroup()
        out <- out[out$rank <= control$k,]
        out$rank <- NULL
    }
    
    return(out)
}

get_score_pair <- function(s_from, s_to, control) {
    # rblast query
    db_from <- AAStringSet(s_from$cdr3)
    names(db_from) <- s_from$id
    
    # rblast db
    db_to <- AAStringSet(s_to$cdr3)
    names(db_to) <- s_to$id
    writeXStringSet(db_to, filepath = "tmp.fasta")
    makeblastdb(db_name = "db/tmp", dbtype = "prot", 
                file = "tmp.fasta", verbose = FALSE)
    bdb <- blast("db/tmp", type = "blastp")
    blast_fmt <- paste("qseqid sseqid pident length qstart qend",
                       "sstart send evalue score qseq sseq qlen slen")
    blast_args <- paste0("-num_threads ", control$blast_cores, 
                         " -evalue 10 -gapopen 11 -gapextend 2", 
                         " -word_size 3 -matrix BLOSUM62 -comp_based_stats 0")
    o <- predict(bdb, db_from, BLAST_args = blast_args, 
                 custom_format = blast_fmt)
    
    # clean up
    unlink("tmp.fasta")
    unlink("db/tmp", recursive = TRUE)
    
    # remove low gmi matches
    o <- o[o$pident >= control$blast_gmi*100,]
    if(nrow(o)==0) {
        return(NULL)
    }
    
    # remove partial hits
    j <- which(o$qstart!=1 | o$sstart!=1 | o$qend != o$qlen | o$send != o$slen)
    if(length(j)!=0) {
        o <- o[-j,]
    }
    # if empty and no duplicates found stop
    if(nrow(o)==0) {
        return(NULL)
    }
    
    if(nrow(o)>0) {
        out <- data.frame(from = o$qseqid,
                          to = o$sseqid,
                          len = o$length,
                          weight = o$score,
                          nweight = o$score/o$length)
        
        data_env <- get_blosum62()
        bs <- t(vapply(X = seq_len(nrow(o)), 
                       FUN.VALUE = numeric(3),
                       FUN = get_bscore, 
                       o = o, 
                       gap_o = -11, 
                       gap_e = -2, 
                       trim = control$trim_flank_aa,
                       b = data_env[["BLOSUM62"]]))
        out$cweight <- bs[,2]
        out$clen <- bs[,3]
        
        out$nweight <- out$weight/out$len
        out$ncweight <- out$cweight/out$clen
        out$chain <- s_from$chain[1]
    }
    
    # KNN 
    if(control$knn == TRUE) {
        out <- out[order(out$nweight, decreasing = TRUE), ]
        out <- out %>%
            group_by(from) %>%
            arrange(desc(nweight), .by_group = TRUE) %>%
            mutate(rank = row_number()) %>%
            ungroup()
        out <- out[out$rank <= control$k,]
        out$rank <- NULL
    }
    
    return(out)
}
