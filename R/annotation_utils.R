

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
                                          "Antigen_species", 
                                          "Antigen_gene", 
                                          "Reference")],
                                 chain = chain,
                                 index = index,
                                 FUN=get_vdjdb_info)))
        }
        if(db_type == "tcr3d") {
            return(unlist(lapply(X = 1:nrow(ns),
                                 ns = ns,
                                 db = db[,c(chain, "Antigen_species", 
                                            "Antigen_gene", 
                                            "Reference")],
                                 chain = chain,
                                 index = index,
                                 FUN = get_tcr3d_info)))
        }
        if(db_type == "mcpas") {
            return(unlist(lapply(X = 1:nrow(ns),
                                 ns = ns,
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
    
    ns$Ag_species <- unlist(lapply(X = a, FUN = function(x){x[["ag_species"]]}))
    ns$Ag_gene <- unlist(lapply(X = a, FUN = function(x) {x[["ag_gene"]]}))
    
    return(ns)
}



# Description:
get_beta_violins <- function(beta,
                             node_summary,
                             ag_species,
                             ag_genes,
                             db = "vdjdb",
                             db_dist = 0,
                             chain = "both") {
    
    get_violins <- function(x, d, db, ag_key, chain) {
        ag_id <- ag_key$id[x]
        ag_name <- ag_key$ag[x]
        ag_type <- ag_key$type[x]
        
        if(chain=="both") {
            ags <- paste0(ag_id, "_", db, "_summary")
            d$size <- d[, ags]
            d$specificity <- ifelse(d[, ags]>0, yes = "+", no = "-")
        } 
        if(chain == "CDR3a" | chain == "CDR3b") {
            ags <- paste0(ag_id, "_", db, "_", chain)
            d$size <- d[, ags]
            d$specificity <- ifelse(d[, ags]>0, yes = "+", no = "-")
        }
        
        g <- ggplot(data = d)+
            geom_sina(aes(x = sample, y = mean, col = specificity, 
                          size = size), stroke = 0.5, alpha = 0.4)+
            theme_bw(base_size = 10)+
            ggtitle(label = '', 
                    subtitle = paste0("Ag=", ag_name, ", ", ag_type, 
                                      " [DB=", db, 
                                      ", chains=", chain, "]"))+
            ylab(label = expression(beta))+
            xlab(label = "repertoire")+
            scale_radius(name = "# cells", 
                         breaks = scales::pretty_breaks(n = 4),
                         range = c(0.5, 3))+
            scale_color_manual(name = "specificity", 
                               values = c("steelblue", "orange"))+
            guides(size = guide_legend(order = 2), 
                   colour = guide_legend(order = 1))
        
        return(g)
    }
    
    # input checks, TODO: pack
    if(missing(beta)) {
        stop("beta is missing")
    }
    if(missing(node_summary)) {
        stop("node_summary is missing")
    }
    
    if(missing(ag_species)) {
        warning("ag_species is missing: skip ag_species annotation")
        ag_species <- ''
    }
    if(is.null(ag_species)|any(is.na(ag_species))) {
        warning("ag_species is NULL/NA: skip ag_species annotation")
        ag_species <- ''
    } 
    if(any(is.character(ag_species)==FALSE)) {
        stop("ag_species must be character vector") 
    }
    
    if(missing(ag_genes)) {
        warning("ag_genes is missing: skip ag_genes annotation")
        ag_genes <- ''
    }
    if(is.null(ag_genes)|any(is.na(ag_genes))) {
        warning("ag_genes is NULL/NA: skip ag_genes annotation")
        ag_genes <- ''
    } 
    if(any(is.character(ag_genes)==FALSE)) {
        stop("ag_genes must be character vector") 
    }
    
    if(missing(db)) {
        stop("db is missing")
    }
    if(length(db)!=1) {
        stop("db must be one of 'vdjdb', 'mcpas', or 'tcr3d'")
    }
    if(is.character(db)==FALSE) {
        stop("db must be character")
    }
    dbs <- get_annotation_dbs()
    if(db %in% dbs == FALSE) {
        stop("db must be one of 'vdjdb', 'mcpas', or 'tcr3d'")
    }
    if(missing(chain)) {
        stop("chain is missing")
    }
    if(length(chain)!=1) {
        stop("chain must be one of CDR3a, CDR3b or both")
    }
    if(is.character(chain)==FALSE) {
        stop("chain must be character")
    }
    if(all(chain %in% c("both", "CDR3a", "CDR3b"))==FALSE) {
        stop("chain must be one of: 'CDR3a', 'CDR3b' or 'both' 
        (DBs do not provide annotations for other chain types)")
    }
    check_db_dist(db_dist)
    if(chain == "both") {
        chains <- get_chains(colnames(node_summary))
    } else {
        chains <- chain
    }
    
    # select node_summary columns
    ns <- node_summary[,c(chains, "community", "sample", "clone_size")]
    
    # construct antigen key
    # if not provided -> create only one entry ''
    if(all(ag_species=="")&all(ag_genes=="")) {
        ag_key <- data.frame(ag = unique(ag_species), type = "species")
        ag_key$id <- "I1"
    } else {
        ag_key <- rbind(data.frame(ag = unique(ag_species), type = "species"),
                        data.frame(ag = unique(ag_genes), type = "gene"))
        ag_key <- ag_key[complete.cases(ag_key),]
        ag_key$id <- paste0("I", 1:nrow(ag_key))
    }
    
    ns <- get_db_match(ns = ns, db = db, db_dist = db_dist)
    na <- get_node_ann(ns = ns, ag_key = ag_key, db = db)
    
    d <- na$ann %>% 
        group_by(community, sample) %>% 
        summarise_at(.funs = sum, .vars = c(na$vars, "clone_size"))
    
    d <- merge(x = beta, 
               y = d[, c("community", "sample", na$vars)], 
               by = c("community", "sample"),
               all.x = TRUE)
    
    w <- which(is.na(d), arr.ind = TRUE)
    if(nrow(w)!=0) {
        d[w] <- 0
    }
    
    d <- d[order(d$mean, decreasing = TRUE),]
    
    # violins
    v <- lapply(X = 1:nrow(ag_key), FUN = get_violins, d = d, 
                db = db, ag_key = ag_key, chain = chain)
    
    return(list(node_annotations = na$ann,
                beta_summary = d,
                vars = vars,
                violins = v))
}


# Description:
get_beta_scatterplot <- function(beta,
                                 node_summary,
                                 ag_species = NULL, 
                                 ag_genes = NULL,
                                 db = "vdjdb",
                                 db_dist = 0,
                                 chain = "both") {
    
    get_scatters <- function(x, d, db, ag_key, chain) {
        ag_id <- ag_key$id[x]
        ag_name <- ag_key$ag[x]
        ag_type <- ag_key$type[x]
        
        if(chain=="both") {
            ags <- paste0(ag_id, "_", db, "_summary")
            d$size <- d[, ags]
            d$specificity <- ifelse(d[, ags]>0, yes = "+", no = "-")
        } 
        if(chain == "CDR3a" | chain == "CDR3b") {
            ags <- paste0(ag_id, "_", db, "_", chain)
            d$size <- d[, ags]
            d$specificity <- ifelse(d[, ags]>0, yes = "+", no = "-")
        }
        
        samples <- unique(d$sample)
        gs <- vector(mode = "list", length = length(samples)^2)
        i <- 1
        for(s1 in samples) {
            for(s2 in samples) {
                dx <- d[d$sample == s1, ]
                dy <- d[d$sample == s2, ]
                
                dx <- dx[, c("mean", "specificity", "size", "community")]
                colnames(dx) <- c("mean_x", "specificity_x", 
                                  "size_x", "community")
                
                dy <- dy[, c("mean", "specificity", "size", "community")]
                colnames(dy) <- c("mean_y", "specificity_y", 
                                  "size_y", "community")
                
                dxy <- merge(x = dx, y = dy, by = "community")
                
                dxy$specificity <- ifelse(test = (dxy$specificity_x == "+" | 
                                              dxy$specificity_y == "+"),
                                          yes = "+", no = "-")
                
                gs[[i]] <- ggplot(data = dxy)+
                    geom_point(aes(x = mean_x, y = mean_y, 
                                   col = specificity, 
                                   size = size_x+size_y), alpha = 0.5)+
                    theme_bw(base_size = 10)+
                    ggtitle(label = '', 
                            subtitle = paste0("Ag=", ag_name, ", ", ag_type, 
                                              " ", "[DB=", db, 
                                              ", chains=", chain, "]"))+
                    ylab(label = bquote(beta ~" (" ~ .(s1) ~ ")"))+
                    xlab(label = bquote(beta ~" (" ~ .(s2) ~ ")"))+
                    scale_radius(name = "# cells", 
                                 breaks = scales::pretty_breaks(n = 4),
                                 range = c(0.5, 3))+
                    scale_color_manual(name = "specificity", 
                                       values = c("steelblue", "orange"))+
                    guides(size = guide_legend(order = 2), 
                           colour = guide_legend(order = 1))
                
                i <- i + 1
            }
        }
        return(gs)
    }
    
    # input checks, TODO: pack
    if(missing(beta)) {
        stop("beta is missing")
    }
    if(missing(node_summary)) {
        stop("node_summary is missing")
    }
    
    if(missing(ag_species)) {
        warning("ag_species is missing: skip ag_species annotation")
        ag_species <- ''
    }
    if(is.null(ag_species)|any(is.na(ag_species))) {
        warning("ag_species is NULL/NA: skip ag_species annotation")
        ag_species <- ''
    } 
    if(any(is.character(ag_species)==FALSE)) {
        stop("ag_species must be character vector") 
    }
    
    if(missing(ag_genes)) {
        warning("ag_genes is missing: skip ag_genes annotation")
        ag_genes <- ''
    }
    if(is.null(ag_genes)|any(is.na(ag_genes))) {
        warning("ag_genes is NULL/NA: skip ag_genes annotation")
        ag_genes <- ''
    } 
    if(any(is.character(ag_genes)==FALSE)) {
        stop("ag_genes must be character vector") 
    }
    
    if(missing(db)) {
        stop("db is missing")
    }
    if(length(db)!=1) {
        stop("db must be one of 'vdjdb', 'mcpas', or 'tcr3d'")
    }
    if(is.character(db)==FALSE) {
        stop("db must be character")
    }
    dbs <- get_annotation_dbs()
    if(db %in% dbs == FALSE) {
        stop("db must be one of 'vdjdb', 'mcpas', or 'tcr3d'")
    }
    if(missing(chain)) {
        stop("chain is missing")
    }
    if(length(chain)!=1) {
        stop("chain must be one of CDR3a, CDR3b or both")
    }
    if(is.character(chain)==FALSE) {
        stop("chain must be character")
    }
    if(all(chain %in% c("both", "CDR3a", "CDR3b"))==FALSE) {
        stop("chain must be one of: 'CDR3a', 'CDR3b' or 'both' 
        (DBs do not provide annotations for other chain types)")
    }
    check_db_dist(db_dist)
    if(chain == "both") {
        chains <- get_chains(colnames(node_summary))
    } else {
        chains <- chain
    }
    
    # select node_summary columns
    ns <- node_summary[,c(chains, "community", "sample", "clone_size")]
    
    # if not provided -> create only one entry ''
    if(all(ag_species=="")&all(ag_genes=="")) {
        ag_key <- data.frame(ag = unique(ag_species), type = "species")
        ag_key$id <- "I1"
    } else {
        ag_key <- rbind(data.frame(ag = unique(ag_species), type = "species"),
                        data.frame(ag = unique(ag_genes), type = "gene"))
        ag_key <- ag_key[complete.cases(ag_key),]
        ag_key$id <- paste0("I", 1:nrow(ag_key))
    }
    
    ns <- get_db_match(ns = ns, db = db, db_dist = db_dist)
    na <- get_node_ann(ns = ns, ag_key = ag_key, db = db)
    
    d <- na$ann %>% 
        group_by(community, sample) %>% 
        summarise_at(.funs = sum, .vars = c(na$vars, "clone_size"))
    
    d <- merge(x = beta, 
               y = d[, c("community", "sample", na$vars)], 
               by = c("community", "sample"),
               all.x = TRUE)
    
    w <- which(is.na(d), arr.ind = TRUE)
    if(nrow(w)!=0) {
        d[w] <- 0
    }
    
    d <- d[order(d$mean, decreasing = TRUE),]
    
    # scatterplots
    v <- lapply(X = 1:nrow(ag_key), FUN = get_scatters, d = d, 
                db = db, ag_key = ag_key, chain = chain)
    
    return(list(node_annotations = na$ann,
                beta_summary = d,
                vars = vars,
                scatterplots = v))
}


get_honeycombs <- function(com) {
    n <- ncol(com)
    gs <- vector(mode = "list", length = n*(n-1)/2)
    iter <- 1
    for(i in 1:(ncol(com)-1)) {
        for(j in (i+1):ncol(com)) {
            
            m <- data.frame(x = com[,i], 
                            y = com[,j], 
                            community = 1:nrow(com),
                            sample_1 = colnames(com)[i], 
                            sample_2 = colnames(com)[j],
                            contrast = paste0(colnames(com)[i], "-", colnames(com)[j]),
                            title = paste0(colnames(com)[i], ' (x) vs ', 
                                           colnames(com)[j], " (y)"))
            
            m$x_adj <- ifelse(test = m$x==0, yes = m$x+0.5, no = m$x)
            m$y_adj <- ifelse(test = m$y==0, yes = m$y+0.5, no = m$y)
            
            gs[[iter]] <- ggplot(data = m)+
                facet_wrap(facets = ~title)+
                geom_hex(aes(x = x_adj, y = y_adj, fill=log10(..count..)), 
                         col = "white", bins = 5)+
                geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
                geom_point(aes(x = x_adj, y = y_adj), 
                           alpha = 0.75, size = 0.5, stroke = 0.3)+
                scale_x_continuous(name = "x", trans = "log10")+
                scale_y_continuous(name = "y", trans = "log10")+
                annotation_logticks(base = 10, sides = "lb")+
                scale_fill_gradient(name = "log10(C)",
                                    low = "white", high = "#FFC68A")+
                theme_bw(base_size = 10)+
                theme(legend.position = "right",
                      strip.text.x = element_text(
                          margin = margin(0.04,0,0.04,0, "cm")))
            
            iter <- iter + 1
        }
    }
    
    return(gs)
}




# Description:
get_ag_summary <- function(communities,
                           node_summary,
                           ag_species = NULL, 
                           ag_genes = NULL,
                           db = "vdjdb",
                           db_dist = 0,
                           chain = "both") {
    
    match_db <- function(ns, db, db_dist) {
        
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
                                              "Antigen_species", 
                                              "Antigen_gene", 
                                              "Reference")],
                                     chain = chain,
                                     index = index,
                                     FUN=get_vdjdb_info)))
            }
            if(db_type == "tcr3d") {
                return(unlist(lapply(X = 1:nrow(ns),
                                     ns = ns,
                                     db = db[,c(chain, "Antigen_species", 
                                                "Antigen_gene", 
                                                "Reference")],
                                     chain = chain,
                                     index = index,
                                     FUN = get_tcr3d_info)))
            }
            if(db_type == "mcpas") {
                return(unlist(lapply(X = 1:nrow(ns),
                                     ns = ns,
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
        
        ns$Ag_species <- unlist(lapply(X = a, FUN = function(x){x[["ag_species"]]}))
        ns$Ag_gene <- unlist(lapply(X = a, FUN = function(x) {x[["ag_gene"]]}))
        
        return(ns)
    }
    
    # input checks, TODO: pack
    if(missing(communities)) {
        stop("communities is missing")
    }
    if(missing(node_summary)) {
        stop("node_summary is missing")
    }
    if(missing(ag_species)&missing(ag_genes)) {
        stop("ag_species and ag_genes are missing")
    }
    if(is.character(ag_species)==FALSE&is.character(ag_genes)==FALSE) {
        stop("ag_species and ag_genes must be character vectors")
    }
    if(length(ag_genes)==0) {
        ag_genes <- NA
    }
    if(length(ag_species)==0) {
        ag_species <- NA
    }
    if(missing(db)) {
        stop("db is missing")
    }
    if(length(db)!=1) {
        stop("db must be one of 'vdjdb', 'mcpas', or 'tcr3d'")
    }
    if(is.character(db)==FALSE) {
        stop("db must be character")
    }
    dbs <- get_annotation_dbs()
    if(db %in% dbs == FALSE) {
        stop("db must be one of 'vdjdb', 'mcpas', or 'tcr3d'")
    }
    if(missing(chain)) {
        stop("chain is missing")
    }
    if(length(chain)!=1) {
        stop("chain must be one of CDR3a, CDR3b or both")
    }
    if(is.character(chain)==FALSE) {
        stop("chain must be character")
    }
    if(all(chain %in% c("both", "CDR3a", "CDR3b"))==FALSE) {
        stop("chain must be one of: 'CDR3a', 'CDR3b' or 'both' 
        (DBs do not provide annotations for other chain types)")
    }
    check_db_dist(db_dist)
    if(chain == "both") {
        chains <- get_chains(colnames(node_summary))
    } else {
        chains <- chain
    }
    
    # select node_summary columns
    ns <- node_summary[,c(chains, "community", "sample", "clone_size")]
    
    # keep info about the community size
    n <- aggregate(clone_size~sample, data = ns, FUN = sum)
    n$repertoire_size <- n$clone_size
    n$clone_size <- NULL
    
    ns <- ns[ns$community %in% communities,]
    
    # construct antigen key
    ag_key <- rbind(data.frame(ag = ag_species, type = "species"),
                    data.frame(ag = ag_genes, type = "gene"))
    ag_key <- ag_key[complete.cases(ag_key),]
    
    
    ns <- match_db(ns = ns, db = db, db_dist = db_dist)
    na <- get_node_ann(ns = ns, ag_key = ag_key, db = db)
    
    # Create the list to fill d with 0s
    l <- as.list(rep(0, length(c(na$vars, "clone_size"))))
    names(l) <- c(na$vars, "clone_size")
    d <- na$ann %>% complete(community, sample, fill = l, explicit = TRUE)
    
    d <- d %>% 
        group_by(community, sample) %>% 
        summarise_at(.funs = sum, .vars = c(na$vars, "clone_size"))
    
    w <- which(is.na(d), arr.ind = TRUE)
    if(nrow(w)!=0) {
        d[w] <- 0
    }
    
    d <- merge(x = d, y = n, by = "sample", all.x = TRUE)
    
    return(d)
}




