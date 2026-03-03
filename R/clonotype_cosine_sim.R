get_clonotype_cosine_similarity <- 
    function(clust_irr, 
             chains = c("CDR3b", "CDR3a"),
             title = "Repertoire overlap - Clonotype level"){
        
        l_r <- length(clust_irr$clust_irrs)
        
        if(l_r == 1) {
            stop("clust_irr has to contain at least two repertoires'")
        }
        
        check_clustirr(clust_irr = clust_irr$clust_irrs[[1]])
        
        if (!all(chains %in% c("CDR3a", "CDR3b"))) {
            stop("chains can only contain 'CDR3a' and/or 'CDR3b'")
        }
        
        v_attr <- names(igraph::vertex_attr(clust_irr$graph))
        
        # catch case were only one chain is provided but chains left at default
        if(!all(chains %in% v_attr)){
            if("CDR3a" %in% v_attr) {
                chains <- "CDR3a"
            } else {
                chains <- "CDR3b"
            }
        } 
        
        l <- list()
        
        # extract repertoires
        for(r in 1:l_r){
            s <- clust_irr$clust_irrs[[r]]@inputs[["s"]]
            
            # set dummy for NA chains
            if (nrow(s) == 0 || all(is.na(s[, chains, drop = F]))) {
                print("i am an idiot")
                df <- data.frame(id = "no_id",
                                 val = 1)
            } else {
                s <- s[,c(chains, "clone_size")]
                
                if (length(chains)>1){
                    df <- data.frame(id = paste0(s$CDR3a, "_", s$CDR3b),
                                     val = s$clone_size)
                } else {
                    df <- data.frame(id = s[[chains]],
                                     val = s$clone_size)
                }
            }
            colnames(df)[2] <- names(clust_irr$clust_irrs)[r]
            
            l[[r]] <- df
        }
        
        # deduplicate
        l <- lapply(l, function(df) {
            aggregate(. ~ id, data = df, FUN = sum)
        })
        
        # join as data.tables
        l <- lapply(l, as.data.table)
        l <- Reduce(function(x, y) merge(x, y, by = "id", all = T), l)
        
        l[is.na(l)] <- 0
        l$id <- NULL
        l <- as.matrix(l)
        
        ret <- get_community_cosine_similarity(l)
        ret$g <- ret$g + ggtitle(title)
        
        return(ret)
    }
