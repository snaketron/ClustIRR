
merge_clustirrs <- function(x) {
    
    get_ix <- function(xs, ns, chains) {
        ixs <- c()
        for(c in chains) {
            ix <- expand.grid(1:xs, 1:xs, KEEP.OUT.ATTRS = FALSE)
            colnames(ix) <- c("index_i", "index_j")
            ix <- ix[ix$index_i!=ix$index_j,]
            ix$name_i <- ns[ix$index_i]
            ix$name_j <- ns[ix$index_j]
            ix$chain <- c
            ixs <- rbind(ixs, ix)
        }
        return(ixs)
    }
    
    get_graph_cons <- function(x, ix, gs, control) {
        # prepare pair-rep data
        s1_name <- ix$name_i[x]
        s2_name <- ix$name_j[x]
        chain <- ix$chain[x]
        
        s1 <- data.frame(cdr3=vertex_attr(graph = gs[[ix$index_i[x]]], name=chain),
                         chain=chain,
                         id=vertex_attr(graph = gs[[ix$index_i[x]]], name="name"))
        s1 <- s1[is.na(s1$cdr3)==FALSE,]
        if(nrow(s1)==0) {
            return(NULL)
        }
        if(is.null(s1)) {
            return(NULL)
        }
        
        s2 <- data.frame(cdr3 = vertex_attr(graph=gs[[ix$index_j[x]]], name=chain),
                         chain=chain,
                         id=vertex_attr(graph = gs[[ix$index_j[x]]], name="name"))
        s2 <- s2[is.na(s2$cdr3)==FALSE,]
        if(nrow(s2)==0) {
            return(NULL)
        }
        if(is.null(s2)) {
            return(NULL)
        }
        
        b <- get_score_pair(s_from = s1, s_to = s2, control = control)
        
        if(is.null(b)==FALSE && nrow(b)!=0) {
            b$chain <- chain
            b$sample <- paste0(s1_name, "|", s2_name)
            b$type <- "between-repertoire"
            return(b)
        }
        
        return(NULL)
    }
    
    check_controls <- function(x) {
        controls <- lapply(X = x, FUN = function(x) {
            if(length(x$clust_irrs)==1) {
                return(get_clustirr_inputs(x$clust_irrs)$control)
            }
            return(get_clustirr_inputs(x$clust_irrs[[1]])$control)
        })
        check <- logical(length = length(controls))
        check[1] <- TRUE
        for(i in 2:length(controls)) {
            check[i] <- identical(controls[[1]], controls[[i]]);
        }
        if(all(check)) {
            return(list(check = TRUE, control = controls[[1]]))
        } 
        return(list(check = FALSE, control = NULL))
    }
    
    
    
    # compare controls
    compare_controls <- check_controls(x)
    control <- compare_controls$control
    if(compare_controls$check==FALSE) {
        stop("mismatching controls")
    }
    
    # get graphs
    gs <- lapply(X = x, FUN = function(x) {
        return(x$graph)
    })
    g <- disjoint_union(gs)
    
    # get chains
    chains <- get_chains(x = vertex_attr_names(g))
    
    # indices 
    ix <- get_ix(xs = length(gs), ns = seq_len(length(gs)), chains = chains)
    
    ige <- lapply(X = seq_len(nrow(ix)),
                  ix = ix,
                  FUN = get_graph_cons,
                  gs = gs,
                  control = control)
    ige <- do.call(rbind, ige)
    
    
    
    # these are the cols we want to keep in this order
    cols <- c("from", "to", "weight", "cweight", "nweight", "ncweight",
              "type", "chain")
    
    df_v <- as_data_frame(x = g, what = "vertices")
    df_e <- as_data_frame(x = g, what = "edges")
    
    if(nrow(df_e)!=0) {
        if(is.null(ige)==FALSE && nrow(ige)!=0) {
            df_e <- rbind(df_e[, cols], ige[, cols])
        }
    } else {
        if(is.null(ige)==FALSE && nrow(ige)!=0) {
            df_e <- df_e[, cols]
        }
    }
    
    key <- data.frame(pmin(df_e$from, df_e$to), 
                      pmax(df_e$from, df_e$to), df_e$chain)
    df_e <- df_e[!duplicated(key), ]
    df_e <- df_e[is.na(df_e$from)==FALSE, ]
    
    # clean unused vars
    rm(ige)
    
    # build joint graph
    g <- graph_from_data_frame(df_e, directed = FALSE, vertices = df_v)
    
    return(list(graph = g, multigraph = TRUE))
}
