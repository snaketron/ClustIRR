
get_global_clust <- function(cdr3, global_max_dist, low_mem, control) {

    get_pairdist <- function(x, a, len_a, global_max_dist) {
        d <- stringdist(a = a[x], b = a[(x + 1):len_a], method = "hamming")
        js <- which(d <= global_max_dist)
        if(length(js) == 0) {
            return(NULL)
        }
        js <- x + js
        return(cbind(rep(x = x, times = length(js)), js))
    }
    
    get_hamming_dist <- function(x, cdr3, cdr3_len, global_max_dist, low_mem) {
        is <- which(cdr3_len == x)
        if(length(is) == 1) {
            return(NULL)
        }
        if(length(is) == 2) {
            d <- stringdist(a = cdr3[is[1]], b = cdr3[is[2]], method="hamming")
            if(d > global_max_dist) {
                return(NULL)
            }
            return(data.frame(from_cdr3 = cdr3[is[1]],
                              to_cdr3 = cdr3[is[2]],
                              weight = 1))
        }
        
        if(low_mem) {
            hd <- lapply(X = seq_len(length(is) - 1),
                         FUN = get_pairdist,
                         a = cdr3[is],
                         len_a = length(is),
                         global_max_dist = global_max_dist)
            hd <- do.call(rbind, hd)
            if(is.null(hd)) {
                return(hd)
            }
            return(data.frame(from_cdr3 = cdr3[is[hd[, 1]]],
                              to_cdr3 = cdr3[is[hd[, 2]]],
                              weight = 1))
        }
        else {
            d <- stringdistmatrix(a = cdr3[is], b = cdr3[is], method="hamming")
            d[upper.tri(x = d, diag = TRUE)] <- NA
            js <- which(d <= global_max_dist, arr.ind = TRUE)
            if(nrow(js) == 0) {
                return(NULL)
            }
            return(data.frame(from_cdr3 = cdr3[is[js[, 1]]],
                              to_cdr3 = cdr3[is[js[, 2]]],
                              weight = 1))
        }
    }
    
    cdr3_len <- nchar(cdr3)
    cdr3_lens <- unique(cdr3_len)
    
    hd <- lapply(X = cdr3_lens,
                 FUN = get_hamming_dist,
                 cdr3 = cdr3,
                 cdr3_len = cdr3_len,
                 global_max_dist = global_max_dist,
                 low_mem = low_mem)
    hd <- do.call(rbind, hd)
    return(hd)
}

get_global_clust_smart <- function(cdr3) {
  
  # computes blosum62 score for pairs of sequences returned by blaster 
  get_bscore <- function(x, d, bm, db) {
    return(stringDist(x = c(db$Seq[d$QueryId[x]], db$Seq[d$TargetId[x]]),
                      method = "substitutionMatrix", 
                      type = "global", 
                      substitutionMatrix = bm, 
                      gapOpening = 10,
                      gapExtension = 4))
  }
  
  db <- data.frame(Id = 1:length(cdr3), Seq = cdr3)
  
  # blast
  o <- blast(query = db, 
             db = db, 
             maxAccepts = 1000, 
             minIdentity = 0.80,
             alphabet = "protein", 
             output_to_file = FALSE)
  
  # remove self-hits
  o <- o[o$QueryId!=o$TargetId,]
  
  # if empty stop
  if(nrow(o)==0) {
    return(NULL)
  }
  
  # compute BLSOUM62 score for matches
  data("BLOSUM62", package = "Biostrings")
  o$bs <- sapply(X = 1:nrow(o), FUN = get_bscore, d = o, db = db, bm = BLOSUM62)
  o$nbs <- ifelse(test = o$bs < 0, yes = o$bs*-1, no = 0)
  # browser()
  # normalize score between 1 (best metch) and -1 (worst match)
  # o$nbs <- o$bs/(-150)
  # o$nbs <- 1/(1+exp(-(-6-0.1*o$bs)))
  
  return(data.frame(from_cdr3 = db$Seq[o$QueryId],
                    to_cdr3 = db$Seq[o$TargetId],
                    weight = o$nbs))
}
