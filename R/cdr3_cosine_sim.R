get_cdr3_cosine_similarity <- 
    function(clust_irr, 
             chains){
        
        l_r <- length(clust_irr$clust_irrs)
        
        if(l_r < 2) {
            stop("clust_irr object has to contain at least two repertoires")
        }
        
        check_clustirr(clust_irr = clust_irr$clust_irrs[[1]])

        v_attr <- names(igraph::vertex_attr(clust_irr$graph))
        
        # if no chains argument supplied, set to existing chains
        if(missing(chains) || is.null(chains)) {
            chains <- grep("\\bCDR3[abhlgd]\\b", v_attr, value = T)
        }
        
        # catch case were only one chain is provided but chains left at default
        if(!all(chains %in% v_attr)){
            stop("clust_irr does not contain all CDR3 chains for comparison")
        } 
        
        l <- list()
        i <- 0
        # extract repertoires
        for(r in 1:l_r){
            i <- i+1
            s <- clust_irr$clust_irrs[[r]]@inputs[["s"]]
            
            # set dummy for NA chains
            if (nrow(s) == 0 || all(is.na(s[, chains, drop = F]))) {
                df <- data.frame(id = paste0("no_id_", i),
                                 val = 1)
            } else {
                s <- s[,c(chains, "clone_size")]
                
                if (length(chains)>1){
                    df <- data.frame(id = paste0(s[[chains[1]]], 
                                                 "_", 
                                                 s[[chains[2]]]),
                                     val = s$clone_size)
                } else {
                    df <- data.frame(id = s[[chains]],
                                     val = s$clone_size)
                }
            }
            colnames(df)[2] <- names(clust_irr$clust_irrs)[r]
            
            l[[r]] <- df
        }
        
        # sort alphabetically to match with get_community_cosine_similarity()
        l <- l[order(names(clust_irr$clust_irrs))]
        
        # aggregate duplicated sequences into single count
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
        ret$g <- ret$g
        
        return(ret)
    }
