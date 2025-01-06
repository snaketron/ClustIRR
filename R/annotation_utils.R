

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


get_annotation_dbs <- function() {
    return(c("vdjdb", "mcpas", "tcr3d"))
}


# Description:
# integrate nodes with data from databases: VDJdb, tcr3d, mcpas-tcr
get_node_ann <- function(node_summary, ag_key) {
    
    get_ann <- function(ns, ag, type) {
        dbs <- get_annotation_dbs()
        
        # created column
        agk <- paste0(ag, "_", dbs, "_summary")
        
        # clone size (multiplier)
        clone_size <- ns$clone_size
        
        # columns of interest
        i <- which(regexec(pattern = paste0(
            paste0("info_", dbs, "_"), collapse = "|"), 
            text = colnames(ns))!=-1)
        is <- do.call(rbind, strsplit(x = colnames(ns)[i], 
                                      split = "\\_"))
        is <- data.frame(db = is[,2], chain = is[,3], 
                         info_key = colnames(ns)[i])
        
        # results
        m <- matrix(data = 0, nrow = nrow(ns), ncol = nrow(is)+length(agk))
        colnames(m) <- c(gsub(pattern = "info\\_", 
                              replacement = paste0(ag, "_"), 
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
        for(db in dbs) {
            ags <- paste0(ag, "_", db, "_summary")
            mm <- m[, which(regexpr(pattern = db, text = colnames(m))!=-1 & 
                      regexpr(pattern = ags, text = colnames(m))==-1),
              drop = FALSE]
            m[,ags] <- apply(X = mm, MARGIN = 1, FUN = max)
        }
        m <- data.frame(m)
        return(m)
    }
    
    m <- node_summary
    vars <- c()
    
    for(i in 1:nrow(ag_key)) {
        browser()
        a<- get_ann(ns = node_summary, ag = ag_key$ag[i], type = ag_key$type[i])
        vars <- c(vars, colnames(a))
        m <- cbind(m, a)
    }
   
    return(list(ann = m, vars = vars))
}

# Description:
get_beta_violins <- function(beta,
                             node_summary,
                             ag_species = NULL, 
                             ag_genes = NULL,
                             db = "vdjdb",
                             chain = "both") {
    
    get_violins <- function(x, d, db, chain) {
        
        if(chain=="both") {
            ags <- paste0(x, "_", db, "_summary")
            d$size <- d[, ags]
            d$specificity <- ifelse(d[, ags]>0, yes = "+", no = "-")
        } 
        if(chain == "alpha" | chain == "beta") {
            ags <- paste0(db, "_", chain)
            d$specificity <- ifelse(d[, ags]>0, yes = "+", no = "-")
        }

        g <- ggplot(data = d)+
            geom_sina(aes(x = sample, y = mean, col = specificity, 
                          size = size), alpha = 0.75)+
            theme_bw()+
            ggtitle(label = x, subtitle = paste0(db, ", chains: ", chain))+
            scale_radius(name = "Specific cells", 
                         breaks = scales::pretty_breaks(n = 4),
                         range = c(0.5, 4))
        
        return(g)
    }
    
    if(missing(beta)) {
        stop("beta is missing")
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
        stop("db must be one of vdjdb, mcpas, or tcr3d")
    }
    if(is.character(db)==FALSE) {
        stop("db must be character")
    }
    dbs <- get_annotation_dbs()
    if(db %in% dbs == FALSE) {
        stop("db must be one of vdjdb, mcpas, or tcr3d")
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
    
    # construct antigen key
    ag_key <- rbind(data.frame(ag = ag_species, type = "species"),
                    data.frame(ag = ag_genes, type = "gene"))
    ag_key <- ag_key[complete.cases(ag_key),]
    
    browser()
    na <- get_node_ann(node_summary = node_summary, ag_key = ag_key)
    
    
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
    v <- lapply(X = ag_key$ag, FUN = get_violins, d = d, db = db, chain = chain)
    
    return(list(node_annotations = na$ann,
                beta_summary = d,
                vars = vars,
                violins = v))
}
    