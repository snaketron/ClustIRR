cluster_irr <- function(s,
                        control = list(gmi = 0.7,
                                       trim_flank_aa = 3,
                                       db_dist = 0,
                                       db_custom = NULL)) {
  # control check
  control <- get_control(control_in = control)
  
  # input check
  input_check(s = s, control = control)
  
  # get chains to be analyzed
  chains <- get_chains(colnames(s))
  
  # add clone size info, sample and clone id
  s <- get_clone_size(s)
  s <- get_sample(s)
  s$id <- seq_len(length.out = nrow(s))
  
  # annotate s with known CDR3-antigen data 
  message("[1/2] annotating CDR3s... \n")
  s <- match_db(cs = s, control = control)
  
  # do clustering
  message("[2/2] clust_irr... \n")
  clust <- lapply(X = chains, FUN = get_clust, s = s, control = control)
  names(clust) <- chains
  
  # setup clustirr object
  return(get_clustirr_output_obj(clust = clust, s = s, control = control))
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
    
    res <- numeric(length = 6)
    res[1] <- score
    res[2] <- length(s)
    res[3] <- res[1]/res[2]
    res[4] <- cscore
    res[5] <- sum(tr==0)
    res[6] <- res[4]/res[5]
    return(res)
  }
  
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
  data_env <- new.env(parent = emptyenv())
  data("BLOSUM62", envir = data_env, package = "ClustIRR")
  # ensure: min(BLOSUM62)=0
  data_env[["BLOSUM62"]] <- data_env[["BLOSUM62"]] + 4
  
  # compute BLSOUM62 scores
  bs <- do.call(rbind, lapply(X = 1:nrow(o), 
                              FUN = get_bscore, 
                              o = o, 
                              gap_o = -10, 
                              gap_e = -4, 
                              trim = control$trim_flank_aa,
                              b = data_env[["BLOSUM62"]]))
  
  o$b <- bs[,1]
  o$max_len <- bs[,2]
  o$nb <- bs[,3]
  o$cb <- bs[,4]
  o$max_clen <- bs[,5]
  o$ncb <- bs[,6]
  
  out <- data.frame(from_cdr3 = db$Seq[o$QueryId],
                    to_cdr3 = db$Seq[o$TargetId],
                    weight = o$b,
                    cweight = o$cb,
                    nweight = o$nb,
                    ncweight = o$ncb,
                    max_len = o$max_len,
                    max_clen = o$max_clen)
  return(out)
}


# Description:
# integrate clustirr with data from databases: VDJdb, tcr3d, mcpas-tcr
match_db <- function(cs, control) {
  
  get_db_info <- function(cs, db, db_type, chain, index) {
    
    get_vdjdb_info <- function(x, cs, db, chain, index) {
      if(x == "") {
        return("")
      }
      xs <- as.numeric(unlist(strsplit(index[x], split = "\\|")))
      
      return(paste0("<db:VDJdb|chain:", chain, "|",
                    "Antigen_species:", 
                    paste0(db[xs, "Antigen_species"], collapse = ';'),"|",
                    "Antigen_gene:", 
                    paste0(db[xs, "Antigen_gene"], collapse = ';'), "|",
                    "CDR3_species:", 
                    paste0(db[xs, "CDR3_species"], collapse = ';'), "|",
                    "Reference:", 
                    paste0(db[xs, "PMID"], collapse = ';'), ">"))
    }
    
    get_tcr3d_info <- function(x, cs, db, chain, index) {
      if(x == "") {
        return("")
      }
      xs <- as.numeric(unlist(strsplit(index[x], split = "\\|")))
      
      return(paste0("<db:tcr3d|chain:", chain, "|",
                    "Antigen_species:", 
                    paste0(db[xs, "Antigen_species"], collapse = ';'),"|",
                    "Antigen_gene:", 
                    paste0(db[xs, "Antigen_gene"], collapse = ';'), "|",
                    "Reference:", 
                    paste0(db[xs, "Reference"], collapse = ';'), ">"))
    }
    
    get_mcpas_info <- function(x, cs, db, chain, index) {
      if(x == "") {
        return("")
      }
      xs <- as.numeric(unlist(strsplit(index[x], split = "\\|")))
      
      return(paste0("<db:mcpas|chain:", chain, "|",
                    "Antigen_species:", 
                    paste0(db[xs, "Antigen_species"], collapse = ';'),"|",
                    "Antigen_gene:", 
                    paste0(db[xs, "Antigen_gene"], collapse = ';'), "|",
                    "CDR3_species:", 
                    paste0(db[xs, "CDR3_species"], collapse = ';'), "|",
                    "Reference:", 
                    paste0(db[xs, "Reference"], collapse = ';'), ">"))
    }
    
    if(db_type == "vdjdb") {
      return(vapply(X = 1:nrow(cs), 
                    FUN = get_vdjdb_info, 
                    cs = cs,
                    db = db[,c(chain, "CDR3_species", "Antigen_species", 
                               "Antigen_gene", "Reference")],
                    chain = chain,
                    index = index,
                    FUN.VALUE = character(1)))
      
      return(unlist(lapply(X=which(cs[,paste0("db_vdjdb_", chain)]==1),
                           cs=cs,
                           db=db[,c(chain, "CDR3_species", 
                                    "Antigen_species", 
                                    "Antigen_gene", 
                                    "Reference")],
                           chain = chain,
                           index = index,
                           FUN=get_vdjdb_info)))
    }
    if(db_type == "tcr3d") {
      return(unlist(lapply(X = 1:nrow(cs),
                           cs = cs,
                           db = db[,c(chain, "Antigen_species", 
                                      "Antigen_gene", 
                                      "Reference")],
                           chain = chain,
                           index = index,
                           FUN = get_tcr3d_info)))
    }
    if(db_type == "mcpas") {
      return(unlist(lapply(X = 1:nrow(cs),
                           cs = cs,
                           db = db[,c(chain, "CDR3_species",
                                      "Antigen_species", 
                                      "Antigen_gene", 
                                      "Reference")],
                           chain = chain,
                           index = index,
                           FUN = get_mcpas_info)))
    }
  }
  
  get_db_index <- function(x, a, b, d) {
    if(d==0) {
      z <- which(a==b[x])
    } else {
      z <- which(stringdist(a=a,b=b[x],method="lv")<=d)
    }
    if(length(z)!=0) {
      return(paste0(z, collapse = '|'))
    }
    return('')
  }
  
  load_data <- function(d) {
    
    e <- new.env()
    if(d=="vdjdb") {
      name <- data("vdjdb", package = "ClustIRR", envir = e)[1]
    }
    if(d=="mcpas") {
      name <- data("mcpas", package = "ClustIRR", envir = e)[1]
    }
    if(d=="tcr3d") {
      name <- data("tcr3d", package = "ClustIRR", envir = e)[1]
    }
    return(e[[name]])
  }
  
  # what type of chaisn are there in the data -> use them to match DBs
  chains <- get_chains(x = colnames(cs))
  
  # load DBs and pack into one list -> db
  db <- list(vdjdb = load_data(d = "vdjdb"), 
             mcpas = load_data(d = "mcpas"), 
             tcr3d = load_data(d = "tcr3d"))
  
  if(is.null(control$db_custom)==FALSE) {
    db[["custom"]] <- control$db_custom
  }
  
  for(db_name in names(db)) {
    for(chain in chains) {
      key_db <- paste0("db_", db_name, "_", chain)
      key_index <- paste0("index_", db_name, "_", chain)
      key_info <- paste0("info_", db_name, "_", chain)
      
      cs[, key_db] <- 0
      cs[, key_index] <- ''
      cs[, key_info] <- ''
      
      # insert indices
      cs[, key_index] <- vapply(X = 1:nrow(cs), 
                                a = db[[db_name]][, chain],
                                b = cs[, chain],
                                d = control$db_dist,
                                FUN.VALUE = character(1), 
                                FUN = get_db_index)
      
      cs[cs[, key_index]!='', key_db] <- 1
      
      if(any(cs[,key_db]==1)) {
        cs[, key_info] <- get_db_info(cs = cs,
                                      db_type = db_name,
                                      db = db[[db_name]],
                                      index = cs[, key_index],
                                      chain = chain)
      }
      
      cs[, key_index] <- NULL
    }
  }
  
  # get aggregate infos
  x <- cs[, which(regexpr(pattern = "info_", text = colnames(cs))!=-1)]
  a <- apply(X = x, MARGIN = 1, FUN = function(x, key) {
    if(all(x=="")) {
      return(list(ag_species = '', ag_gene = ''))
    }
    
    ag_species <- c()
    ag_gene <- c()
    x <- x[x!=""]
    for(i in 1:length(x)) {
      y <- unlist(strsplit(x = x[i], split = "\\|"))
      ag_species <- c(ag_species, unlist(strsplit(x = y[3], split = '\\;')))
      ag_gene <- c(ag_gene, unlist(strsplit(x = y[4], split = '\\;')))
    }
    ag_species <- paste0(unique(gsub(pattern = "Antigen_species\\:", 
                                     replacement = '', x = ag_species)), 
                         collapse = ',')
    ag_gene <- paste0(unique(gsub(pattern = "Antigen_gene\\:", 
                                  replacement = '', x = ag_gene)),
                      collapse = ',')
    
    return(list(ag_species = ag_species, ag_gene = ag_gene))
  })
  
  cs$Ag_species <- unlist(lapply(X = a, FUN = function(x) {x[["ag_species"]]}))
  cs$Ag_gene <- unlist(lapply(X = a, FUN = function(x) {x[["ag_gene"]]}))
  
  return(cs)
}



