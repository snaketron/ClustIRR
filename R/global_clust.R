
get_global_clust_hamming <- function(cdr3, control) {
  
  get_cweight <- function(x, a, b, trim, len) {
    return(stringdist(a = substr(x=a[x], start=trim+1, stop=len-trim+1), 
                      b = substr(x=b[x], start=trim+1, stop=len-trim+1), 
                      method="hamming"))
  }
  
  get_pairdist <- function(x, a, len_a, global_max_dist) {
    d <- stringdist(a = a[x], b = a[(x + 1):len_a], method = "hamming")
    js <- which(d <= global_max_dist)
    if(length(js) == 0) {
      return(NULL)
    }
    js <- x + js
    return(cbind(rep(x = x, times = length(js)), js))
  }
  
  # Computes Hamming distances
  get_hd <- function(x, cdr3, cdr3_len, global_max_dist, 
                     low_mem, trim_flank_aa) {
    
    is <- which(cdr3_len == x)
    if(length(is) == 1) {
      return(NULL)
    }
    if(length(is) == 2) {
      d <- stringdist(a = cdr3[is[1]], 
                      b = cdr3[is[2]], 
                      method = "hamming")
      if(d > global_max_dist) {
        return(NULL)
      }
      
      # compute hamming distance of core CDR3 region
      if(trim_flank_aa == 0) {
        cweight <- 1
      } 
      else {
        if(trim_flank_aa*2 <= x) {
          cweight <- 0
        } 
        else {
          cweight <- get_cweight(x = 1, 
                                 a = cdr3[is[1]], 
                                 b = cdr3[is[2]], 
                                 trim = trim_flank_aa, 
                                 len = x)
        }
      }
      
      return(data.frame(from_cdr3 = cdr3[is[1]],
                        to_cdr3 = cdr3[is[2]],
                        weight = 1,
                        cweight = cweight))
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
      
      # compute hamming distance of core CDR3 region
      if(trim_flank_aa*2 <= x | trim_flank_aa == 0) {
        cweight <- rep(x = 0, times = nrow(hd))
      }
      else {
        cweight <- vapply(X = 1:nrow(hd), 
                          a = cdr3[is[hd[, 1]]], 
                          b = cdr3[is[hd[, 2]]],
                          trim = trim_flank_aa, 
                          len = x,
                          FUN.VALUE = numeric(1))
      }
      
      return(data.frame(from_cdr3 = cdr3[is[hd[, 1]]],
                        to_cdr3 = cdr3[is[hd[, 2]]],
                        weight = 1,
                        cweight = cweight))
    }
    else {
      d <- stringdistmatrix(a = cdr3[is], b = cdr3[is], method="hamming")
      d[upper.tri(x = d, diag = TRUE)] <- NA
      js <- which(d <= global_max_dist, arr.ind = TRUE)
      if(nrow(js) == 0) {
        return(NULL)
      }
      
      # compute hamming distance of core CDR3 region
      if(trim_flank_aa == 0) {
        cweight <- rep(x = 1, times = nrow(js))
      } 
      else {
        if(trim_flank_aa*2 <= x) {
          cweight <- rep(x = 0, times = nrow(js))
        }
        else {
          cweight <- vapply(X = 1:length(js[, 1]), 
                            a = cdr3[is[js[, 1]]], 
                            b = cdr3[is[js[, 2]]],
                            trim = trim_flank_aa, 
                            len = x,
                            FUN.VALUE = numeric(1))
        }
      }
      
      return(data.frame(from_cdr3 = cdr3[is[js[, 1]]],
                        to_cdr3 = cdr3[is[js[, 2]]],
                        weight = 1,
                        cweight = cweight))
    }
  }
  
  cdr3_len <- nchar(cdr3)
  cdr3_lens <- unique(cdr3_len)
  
  hd <- lapply(X = cdr3_lens,
               FUN = get_hd,
               cdr3 = cdr3,
               cdr3_len = cdr3_len,
               global_max_dist = control$global_max_dist,
               low_mem = control$low_mem,
               trim_flank_aa = control$trim_flank_aa)
  
  hd <- do.call(rbind, hd)
  return(hd)
}


get_global_clust_blosum <- function(cdr3, control) {
  
  # computes BLOSUM62 score for pairs of sequences returned by blaster 
  get_bscore <- function(x, d, bm, db) {
    return(stringDist(x = c(db$Seq[d$QueryId[x]], db$Seq[d$TargetId[x]]),
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
    if((na-2*trim)<=0 | (nb-2*trim)) {
      return(NA)
    }
    
    a <- substr(x = a, start = trim+1, stop = nchar(a)-trim)
    b <- substr(x = b, start = trim+1, stop = nchar(b)-trim)
    
    return(stringDist(x = c(a, b),
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
             minIdentity = 0.90,
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
    o <- o[key_js, ]
  }
  
  # get blosum matrix from Biostrings
  data_env <- new.env(parent = emptyenv())
  data("BLOSUM62", envir = data_env, package = "Biostrings")
  
  # compute BLSOUM62 score for matches
  o$bs <- vapply(X = 1:nrow(o), 
                 FUN = get_bscore, 
                 d = o, 
                 db = db, 
                 bm = data_env[["BLOSUM62"]],
                 FUN.VALUE = numeric(1))
  # o$nbs <- o$bs
  # o$nbs <- ifelse(test = o$nbs < 0, yes = 0, no = o$nbs)
  # o$nbs <- ifelse(test = o$nbs > 1, yes = 1, no = o$nbs)
  
  
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
    o$core_nbs <- o$core_bs
    # o$core_nbs <- ifelse(test = o$core_nbs < 0, yes = 0, no = o$core_nbs)
    # o$core_nbs <- ifelse(test = o$core_nbs > 1, yes = 1, no = o$core_nbs)
  }
  
  
  return(data.frame(from_cdr3 = db$Seq[o$QueryId],
                    to_cdr3 = db$Seq[o$TargetId],
                    weight = o$bs,
                    cweight = o$core_bs))
}


get_global_clust <- function(cdr3, control) {
  # if global_pairs are provided as input use them, else compute them
  if(!is.null(control$global_pairs)) {
    g <- control$global_pairs
  }
  else {
    if(control$global_smart==TRUE) {
      g <- get_global_clust_blosum(cdr3 = cdr3, control = control)
    }
    if(control$global_smart==FALSE) {
      g <- get_global_clust_hamming(cdr3 = cdr3, control = control)
    } 
  }
  return(g)
}
