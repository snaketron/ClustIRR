
get_ag_species_hits <- function(node_summary, 
                                db = "vdjdb", 
                                ag_species) {
    db_cols <- which(regexpr(pattern = paste0("info_", db, "_*"), 
                             text = colnames(node_summary))!=-1)
    
    node_summary$cells <- node_summary$clone_size
    node_summary$clones <- 1
    
    table_summary <- c()
    for(i in db_cols) {
        name_col <- gsub(pattern = paste0("info_", db, "_"), 
                         replacement = '', x = colnames(node_summary)[i])
        new_col <- paste0(paste0(ag_species, collapse = '_'), "_", name_col)
        node_summary[,new_col] <- 0
        
        ann <- strsplit(x = node_summary[,i], split = "\\|")
        js <- which(is.na(ann)==FALSE)
        if(length(js)!=0) {
            hits <- vapply(X = ann[js], FUN = function(x, ag_species) {
                x <- gsub(pattern="Antigen_species\\:",replacement='',x=x[[3]])
                x <- unlist(strsplit(x = x, split = "\\;"))
                return(any(x %in% ag_species))
            }, ag_species = ag_species, FUN.VALUE = logical(length = 1))
            js <- js*hits
            if(any(js!=0)) {
                node_summary[js, new_col] <- 1
                
                # table summary
                ts <- node_summary[,c("sample", "community", 
                                      "cells", "clones", new_col)]
                ts$ag <- ts[, new_col]
                
                # compute repertoire properties
                rs <- node_summary %>% 
                    group_by(sample, .drop = TRUE) %>%
                    summarise(rep_cells = sum(cells), 
                              rep_clones = sum(clones),
                              .groups = "drop")
                
                # compute community properties
                cs <- ts %>% 
                    group_by(sample, community, .drop = TRUE) %>%
                    summarise(com_cells = sum(cells), 
                              com_clones = sum(clones),
                              .groups = "drop") %>%
                    complete(sample, community, 
                             fill = list(total_cells = 0, total_clones = 0))
                
                # compute community and ag-specific properties
                as <- ts %>% 
                    group_by(sample, community, ag, .drop = TRUE) %>%
                    summarise(ag_cells = sum(cells), 
                              ag_clones = sum(clones))
                as <- as[as$ag==1,]
                
                # merge stuff
                cs <- merge(x = cs, y = rs, by = c("sample"), all.x = TRUE)
                as <- merge(x = as, y = cs, by = c("sample", "community"), 
                            all.y = TRUE)
                as[which(is.na(as), arr.ind = TRUE)] <- 0
                as$ag <- NULL
                as$ag_key <- new_col
                
                table_summary <- rbind(table_summary, as)
            }
        }
    }
    
    cols <- paste0(paste0(ag_species, collapse = '_'), "_", 
                   gsub(pattern = paste0("info_", db, "_"), 
                        replacement = '', x = colnames(node_summary)[db_cols]))
    
    return(list(node_summary = node_summary,
                new_columns = cols,
                table_summary = table_summary))
}

get_ag_gene_hits <- function(node_summary, 
                             db = "vdjdb", 
                             ag_gene) {
    db_cols <- which(regexpr(pattern = paste0("info_", db, "_*"), 
                             text = colnames(node_summary))!=-1)
    
    node_summary$cells <- node_summary$clone_size
    node_summary$clones <- 1
    
    table_summary <- c()
    for(i in db_cols) {
        name_col <- gsub(pattern = paste0("info_", db, "_"), 
                         replacement = '', x = colnames(node_summary)[i])
        new_col <- paste0(paste0(ag_gene, collapse = '_'), "_", name_col)
        node_summary[,new_col] <- 0
        
        ann <- strsplit(x = node_summary[,i], split = "\\|")
        js <- which(is.na(ann)==FALSE)
        if(length(js)!=0) {
            hits <- vapply(X = ann[js], FUN = function(x, ag_gene) {
                x <- gsub(pattern = "Antigen_gene\\:", replacement='', x=x[[4]])
                x <- unlist(strsplit(x = x, split = "\\;"))
                return(any(x %in% ag_gene))
            }, ag_gene = ag_gene, FUN.VALUE = logical(length = 1))
            js <- js*hits
            if(any(js!=0)) {
                node_summary[js, new_col] <- 1
                
                # table summary
                ts <- node_summary[,c("sample", "community", 
                                      "cells", "clones", new_col)]
                ts$ag <- ts[, new_col]
                
                # compute repertoire properties
                rs <- node_summary %>% 
                    group_by(sample, .drop = TRUE) %>%
                    summarise(rep_cells = sum(cells), 
                              rep_clones = sum(clones),
                              .groups = "drop")
                
                # compute community properties
                cs <- ts %>% 
                    group_by(sample, community, .drop = TRUE) %>%
                    summarise(com_cells = sum(cells), 
                              com_clones = sum(clones),
                              .groups = "drop") %>%
                    complete(sample, community, 
                             fill = list(total_cells = 0, total_clones = 0))
                
                # compute community and ag-specific properties
                as <- ts %>% 
                    group_by(sample, community, ag, .drop = TRUE) %>%
                    summarise(ag_cells = sum(cells), 
                              ag_clones = sum(clones))
                as <- as[as$ag==1,]
                
                # merge stuff
                cs <- merge(x = cs, y = rs, by = c("sample"), all.x = TRUE)
                as <- merge(x = as, y = cs, by = c("sample", "community"), 
                            all.y = TRUE)
                as[which(is.na(as), arr.ind = TRUE)] <- 0
                as$ag <- NULL
                as$ag_key <- new_col
                
                table_summary <- rbind(table_summary, as)
            }
        }
    }
    
    cols <- paste0(paste0(ag_gene, collapse = '_'), "_", 
                   gsub(pattern = paste0("info_", db, "_"), 
                        replacement = '', x = colnames(node_summary)[db_cols]))

    return(list(node_summary = node_summary,
                new_columns = cols,
                table_summary = table_summary))
}
