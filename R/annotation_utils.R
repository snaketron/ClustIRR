# Description:
# annotate nodes with data from databases: VDJdb, tcr3d, mcpas-tcr
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
                                          "Antigen_gene", "Reference")],
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
                                            "Antigen_species", "Antigen_gene", 
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
    
    # make sure cs is a data.frame
    cs <- data.frame(cs)
    
    # what type of chaisn are there in the data -> use them to match DBs
    chains <- get_chains(x = colnames(cs))
    
    # load DBs and pack into one list -> db
    db_names <- get_annotation_dbs()
    db <- vector(mode = "list", length = length(db_names))
    names(db) <- db_names
    for(i in seq_len(length(db))) {
        db[[i]] <- load_data(d = names(db)[i])
    }
    
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
            cs[, key_info] <- get_db_info(cs = cs,
                                          db_type = db_name,
                                          db = db[[db_name]],
                                          index = cs[, key_index],
                                          chain = chain)
            cs[cs[, key_db]==0, key_info] <- NA
            cs[, key_index] <- NULL
        }
    }
    
    # get aggregate infos
    x <- cs[, which(regexpr(pattern = "info_", text = colnames(cs))!=-1)]
    a <- apply(X = x, MARGIN = 1, FUN = function(x, key) {
        if(all(is.na(x))) {
            return(list(ag_species = '', ag_gene = ''))
        }
        
        ag_species <- c()
        ag_gene <- c()
        x <- x[is.na(x)==FALSE]
        for(i in 1:length(x)) {
            y <- unlist(strsplit(x = x[i], split = "\\|"))
            ag_species <- c(ag_species, unlist(strsplit(x = y[3], 
                                                        split = '\\;')))
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
    
    cs$Ag_species <- unlist(lapply(X = a, FUN = function(x){x[["ag_species"]]}))
    cs$Ag_gene <- unlist(lapply(X = a, FUN = function(x) {x[["ag_gene"]]}))
    
    return(cs)
}

get_annotation_dbs <- function() {
    return(c("vdjdb", "mcpas", "tcr3d"))
}

# Description:
# integrate nodes with data from databases: VDJdb, tcr3d, mcpas-tcr
get_node_ann <- function(ns, ag_key, db) {
    
    get_ann <- function(ns, ag, type, id, db) {
        
        # created column
        agk <- paste0(id, "_", db, "_summary")
        
        # clone size (multiplier)
        clone_size <- ns$clone_size
        
        # columns of interest
        i <- which(regexec(pattern = paste0(
            paste0("info_", db, "_"), collapse = "|"), 
            text = colnames(ns))!=-1)
        is <- do.call(rbind, strsplit(x = colnames(ns)[i], split = "\\_"))
        is <- data.frame(db = is[,2], chain = is[,3], 
                         info_key = colnames(ns)[i])
        
        # results
        m <- matrix(data = 0, nrow = nrow(ns), ncol = nrow(is)+length(agk))
        colnames(m) <- c(gsub(pattern = "info\\_", 
                              replacement = paste0(id, "_"), 
                              x = is$info_key), agk)
        
        for(i in 1:nrow(is)) {
            k <- do.call(rbind, strsplit(x = ns[,is$info_key[i]], 
                                         split = "\\|"))
            if(type == "species") {
                g <- gsub(pattern = "Antigen\\_species\\:", 
                          replacement = '', x = k[,3])
                m[,i] <- vapply(X = g, 
                                ag = ag, 
                                clone_size = clone_size,
                                FUN.VALUE = numeric(1),
                                FUN = function(x, ag, clone_size) {
                                    if(ag=="") {
                                        return(FALSE)
                                    }
                                    return(regexpr(pattern = ag, 
                                                   text = x)!=-1)
                                })*clone_size
            }
            else if(type == "gene") {
                g <- gsub(pattern = "Antigen\\_gene\\:", 
                          replacement = '', x = k[,4])
                m[,i] <- vapply(X = g, 
                                ag = ag, 
                                clone_size = clone_size,
                                FUN.VALUE = numeric(1),
                                FUN = function(x, ag, clone_size) {
                                    return(regexpr(pattern = ag, 
                                                   text = x)!=-1)
                                })*clone_size
            }
        }
        
        # DB-specific cells -> summary for both chains A|B cells
        ags <- paste0(id, "_", db, "_summary")
        mm <- m[, which(regexpr(pattern = db, text = colnames(m))!=-1 & 
                            regexpr(pattern = ags, text = colnames(m))==-1),
                drop = FALSE]
        m[,ags] <- apply(X = mm, MARGIN = 1, FUN = max)
        m <- data.frame(m)
        return(m)
    }
    
    m <- ns
    vars <- c()
    seq_len(nrow(ag_key))
    for(i in seq_len(nrow(ag_key))) {
        a <- get_ann(ns = ns, ag = ag_key$ag[i], type = ag_key$type[i], 
                     id = ag_key$id[i], db = db)
        vars <- c(vars, colnames(a))
        m <- cbind(m, a)
    }
    
    return(list(ann = m, vars = vars))
}

# Get annotation database matches for visualization utilizies
get_db_match <- function(ns, db, db_dist) {
    
    get_db_info <- function(ns, db, db_type, chain, index) {
        
        get_vdjdb_info <- function(x, ns, db, chain, index) {
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
        
        get_tcr3d_info <- function(x, ns, db, chain, index) {
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
        
        get_mcpas_info <- function(x, ns, db, chain, index) {
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
            return(vapply(X = 1:nrow(ns), 
                          FUN = get_vdjdb_info, 
                          ns = ns,
                          db = db[,c(chain, "CDR3_species", "Antigen_species", 
                                     "Antigen_gene", "Reference")],
                          chain = chain,
                          index = index,
                          FUN.VALUE = character(1)))
            
            return(unlist(lapply(X=which(ns[,paste0("db_vdjdb_", chain)]==1),
                                 ns=ns,
                                 db=db[,c(chain, "CDR3_species", 
                                          "Antigen_species", "Antigen_gene", 
                                          "Reference")],
                                 chain = chain,
                                 index = index,
                                 FUN=get_vdjdb_info)))
        }
        if(db_type == "tcr3d") {
            return(unlist(lapply(X = 1:nrow(ns),
                                 ns = ns,
                                 db = db[,c(chain, "Antigen_species", 
                                            "Antigen_gene", "Reference")],
                                 chain = chain,
                                 index = index,
                                 FUN = get_tcr3d_info)))
        }
        if(db_type == "mcpas") {
            return(unlist(lapply(X = 1:nrow(ns),
                                 ns = ns,
                                 db = db[,c(chain, "CDR3_species",
                                            "Antigen_species", "Antigen_gene", 
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
    chains <- get_chains(x = colnames(ns))
    
    # load DBs and pack into one list -> db
    db_name <- db
    db <- vector(mode = "list", length = 1)
    names(db) <- db_name
    db[[1]] <- load_data(d = names(db)[1])
    
    for(chain in chains) {
        key_db <- paste0("db_", db_name, "_", chain)
        key_index <- paste0("index_", db_name, "_", chain)
        key_info <- paste0("info_", db_name, "_", chain)
        
        ns[, key_db] <- 0
        ns[, key_index] <- ''
        ns[, key_info] <- ''
        
        # insert indices
        ns[, key_index] <- vapply(X = 1:nrow(ns), 
                                  a = db[[db_name]][, chain],
                                  b = ns[, chain],
                                  d = db_dist,
                                  FUN.VALUE = character(1), 
                                  FUN = get_db_index)
        
        ns[ns[, key_index]!='', key_db] <- 1
        ns[, key_info] <- get_db_info(ns = ns,
                                      db_type = db_name,
                                      db = db[[db_name]],
                                      index = ns[, key_index],
                                      chain = chain)
        ns[, key_index] <- NULL
    }
    
    # get aggregate infos
    x <- ns[, which(regexpr(pattern = "info_", text = colnames(ns))!=-1), 
            drop=FALSE]
    a <- apply(X = x, MARGIN = 1, FUN = function(x, key) {
        if(all(x=="")) {
            return(list(ag_species = '', ag_gene = ''))
        }
        
        ag_species <- c()
        ag_gene <- c()
        x <- x[x!=""]
        for(i in 1:length(x)) {
            y <- unlist(strsplit(x = x[i], split = "\\|"))
            ag_species <- c(ag_species, unlist(strsplit(x = y[3],split ='\\;')))
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
    
    ns$Ag_species <- unlist(lapply(X = a, FUN = function(x){x[["ag_species"]]}))
    ns$Ag_gene <- unlist(lapply(X = a, FUN = function(x) {x[["ag_gene"]]}))
    
    return(ns)
}

