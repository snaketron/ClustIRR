get_functional_communities <- function(com,
                                       ag_gene,
                                       ag_species){
    
    add_hits <- function(df, col, keys) {

        for (key in keys) {
            col_s <- paste0(col, "_", key)
            df[[col_s]] <- grepl(key, df[[col]], ignore.case = T) 
        }

        new_cols <- paste0(col, "_", keys)
        df$functional_node <- rowSums(df[, new_cols, drop = F]) > 0 |
            df$functional_node
        
        return(df)
    }
    
    node_summary <- com$node_summary
    node_summary$functional_node <- F
    node_summary <- add_hits(node_summary, "Ag_gene", ag_gene)
    node_summary <- add_hits(node_summary, "Ag_species", ag_species)
    
    com$node_summary <- node_summary
    com$functional_communities <- 
        unique(node_summary[node_summary$functional_node, "community"])
    com$community_summary$functional_community <- 
        com$community_summary$community %in% com$functional_communities

    return(com)
    
}








