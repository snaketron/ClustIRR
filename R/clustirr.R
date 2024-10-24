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
  s <- match_db(cs = s, control = control)
  
  # do clustering
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
  
  cdr3 <- table(cdr3)
  cdr3_dup <- ifelse(test = as.numeric(cdr3)>1, yes = 1, no = 0)
  cdr3 <- names(cdr3)
  
  return(get_blosum(cdr3 = cdr3, cdr3_dup = cdr3_dup, control = control))
}



get_blosum <- function(cdr3, cdr3_dup, control) {
  
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
             minIdentity = control$gmi,
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
  data("BLOSUM62", envir = data_env, package = "pwalign")
  
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
  out$max_clen <- apply(X = out[, c("from_cdr3", "to_cdr3")], MARGIN = 1,
                        trim = control$trim_flank_aa, 
                        FUN = function(x, trim) {return(max(nchar(x)-trim*2))})
  out$max_clen <- ifelse(test=out$max_clen<0, yes = 0, no = out$max_clen)
  
  
  out$nweight <- out$weight/out$max_len
  out$ncweight <- out$cweight/out$max_clen
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
    z <- which(stringdist(a=a,b=b[x],method="lv")<=d)
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
