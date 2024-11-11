
get_graph <- function(clust_irr) {
    
    get_edges <- function(clust_irr, cs) {
        eg <- vector(mode="list", length = length(get_clustirr_clust(clust_irr)))
        names(eg) <- names(get_clustirr_clust(clust_irr))
        for(chain in names(get_clustirr_clust(clust_irr))) {
            g <- get_clustirr_clust(clust_irr)[[chain]]
            if(is.null(g) == FALSE && nrow(g) != 0) {
                g <- merge(x = g, y = cs[, c(chain, "clone_id")], 
                           by.x = "from_cdr3", by.y = chain, all.x = TRUE)
                g$from_clone_id <- g$clone_id
                g$clone_id <- NULL
                
                g <- merge(x = g, y = cs[, c(chain, "clone_id")], 
                           by.x = "to_cdr3", by.y = chain, all.x = TRUE)
                g$to_clone_id <- g$clone_id
                g$clone_id <- NULL
                
                # remove duplicated ids
                g <- g[g$from_clone_id != g$to_clone_id, ]
                g$key <- apply(X = g[, c("from_clone_id", "to_clone_id")], 
                               MARGIN = 1, FUN = function(x) {
                                   return(paste0(sort(x),collapse='|'))})
                g <- g[duplicated(g$key)==FALSE,]
                g$key <- NULL
                
                out <- data.frame(from_cdr3 = g[,"from_clone_id"], 
                                  to_cdr3 = g[,"to_clone_id"], 
                                  weight = g[,"weight"],
                                  nweight = g[,"nweight"],
                                  cweight = g[,"cweight"],
                                  ncweight = g[,"ncweight"],
                                  max_len = g[,"max_len"],
                                  max_clen = g[,"max_clen"],
                                  motif = NA,
                                  chain = chain,
                                  type = "global")
                
                eg[[chain]] <- out
            }
        }
        eg <- do.call(rbind, eg)
        return(eg)
    }
    
    get_edge_order <- function(x) {
        return(paste0(sort(c(x[1], x[2])), collapse = '-'))
    }
    
    build_graph <- function(e, cs, sample_id, chains) {
        
        add_edges <- function(e, ig, sample_id, chain) {
            
            get_e <- function(x, sample_id) {
                return(paste0(sample_id, '|', x))
            }
            
            ve <- as.vector(apply(X = e[, c("from_cdr3", "to_cdr3")], 
                                  MARGIN = 1, FUN = get_e, 
                                  sample_id = sample_id))
            ig <- igraph::add_edges(graph = ig, 
                                    edges = ve, 
                                    weight = e$weight,
                                    cweight = e$cweight,
                                    nweight = e$nweight,
                                    ncweight = e$ncweight,
                                    max_len = e$max_len,
                                    max_clen = e$max_clen,
                                    type = "within-repertoire",
                                    chain = chain,
                                    clustering = "global")
            
            return(ig)
        }
        
        ig <- graph_from_data_frame(d = data.frame(from = cs$name[1], 
                                                   to = cs$name[1]),
                                    directed = FALSE,
                                    vertices = cs)
        ig <- delete_edges(ig, edges = 1)
        
        # add global edges
        message("(", sample_id, "): adding edges... \n")
        if(is.null(e)==FALSE && nrow(e)!=0) {
            for(chain in chains) {
                chain_ge <- e[e$chain == chain, ]
                if(nrow(chain_ge)!=0) {
                    ig <- add_edges(e = chain_ge, ig = ig, 
                                    sample_id = sample_id,
                                    chain = chain)
                }
            }
        }
        
        return(ig)
    }
    
    check_clustirr(clust_irr = clust_irr)
    
    # get chains
    chains <- get_chains(x = colnames(get_clustirr_inputs(clust_irr)$s))
    
    # cells
    s <- get_clustirr_inputs(clust_irr)$s
    
    # setting up the sample id
    sample_id <- s$sample[1]
    
    # get clones
    cs <- get_clones(sample_id = sample_id, x = s)
    
    # get edges between clones
    e <- get_edges(clust_irr = clust_irr, cs = cs)
    
    # build graph with only vertices
    if(is.null(e)) {
        ig <- graph_from_data_frame(d = data.frame(from=cs$name[1], 
                                                   to=cs$name[1]), 
                                    directed = FALSE, vertices = cs)
        ig <- delete_edges(ig, edges = 1)
        
        return(list(graph = ig, clones = cs))
    }
    
    # build graph
    message("building graph... \n")
    ig <- build_graph(e = e, cs = cs, sample_id = sample_id, chains = chains)
    
    return(list(graph = ig, clones = cs, joint_graph = FALSE))
}


get_joint_graph <- function(clust_irrs, cores = 1) {
    
    check_input <- function(clust_irrs) {
        if(missing(clust_irrs)==TRUE) {
            stop("clust_irrs input missing")
        }
        if(is.list(clust_irrs)==FALSE) {
            stop("clust_irrs must be a list of clust_irr objects")
        }
        if(length(clust_irrs)<=1) {
            stop("get_joint_graph needs >= 2 clust_irr outputs")
        }
        for(i in 1:length(clust_irrs)) {
            check_clustirr(clust_irr = clust_irrs[[i]])
        }
        
        # check if same chain names
        cs <- colnames(get_clustirr_inputs(clust_irrs[[1]])$s)
        for(i in 2:length(clust_irrs)) {
            if(any(colnames(get_clustirr_inputs(clust_irrs[[i]])$s)!=cs)) {
                stop("different chains in graphs")
            }
        }
    }
    
    get_clust_irrs_names <- function(clust_irrs) {
        clust_irrs_names <- names(clust_irrs)
        if(is.null(clust_irrs_names)) {
            names(clust_irrs) <- paste0("s", 1:length(clust_irrs))
        }
        
        for(i in 1:length(clust_irrs)) {
            clust_irrs[[i]]@inputs[["sample_id"]] <- names(clust_irrs)[i]
        }
        
        return(clust_irrs)
    }
    
    get_v_e <- function(x, what) {
        return(as_data_frame(x$graph, what = what))
    }
    
    # check input
    check_input(clust_irrs = clust_irrs)
    
    # check cores
    check_cores(cores = cores)
    
    # get clust_irrs names
    clust_irrs <- get_clust_irrs_names(clust_irrs = clust_irrs)
    
    # get joint controls
    control <- get_joint_controls(clust_irrs = clust_irrs)
    
    # get igs
    future::plan(future::multisession, workers = I(cores))
    igs <- future_lapply(X = clust_irrs, FUN = get_graph, future.seed = TRUE)
    names(igs) <- names(clust_irrs)
    
    # get chains
    chains <- get_chains(x = colnames(get_clustirr_inputs(clust_irrs[[1]])$s))
    
    ige <- get_intergraph_edges(igs = igs,
                                chains = chains,
                                cores = cores,
                                trim_flank_aa = control$trim_flank_aa,
                                gmi = control$gmi)
    
    # get the vertices/edges of the graph
    df_v <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "vertices"))
    df_e <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "edges"))
    
    # these are the cols we want to keep in this order
    cols <- c("from", "to", "weight", "cweight", "nweight", "ncweight",
              "max_len", "max_clen", "type", "chain", "clustering")
    
    if(nrow(df_e)!=0) {
        if(is.null(ige)==FALSE && nrow(ige)!=0) {
            df_e <- rbind(df_e[, cols], ige[, cols])
        }
    } 
    else {
        if(is.null(ige)==FALSE && nrow(ige)!=0) {
            df_e <- ige[, cols]
        }
    }
    
    # build joint graph
    g <- graph_from_data_frame(df_e, directed = FALSE, vertices = df_v)
    
    return(list(graph = g, clust_irrs = clust_irrs, joint_graph = TRUE))
}


plot_graph <- function(g,
                       select_by = "Ag_species",
                       as_visnet = FALSE, 
                       show_singletons = TRUE,
                       node_opacity = 1) {
    
    check_graph <- function(g) {
        if(missing(g)) {
            stop("missing input g")
        }
        if(is.list(g)==FALSE) {
            stop("missing input g")
        }
        if(is_igraph(g$graph)==FALSE) {
            stop("wrong input g")
        }
        if(is.null(g$graph)) {
            stop("wrong input g")
        }
        if(is.logical(g$joint_graph)==FALSE) {
            stop("wrong input g")
        }
    }
    
    check_graph(g = g)
    check_as_visnet(as_visnet = as_visnet)
    check_show_singletons(show_singletons = show_singletons)
    check_select_by(select_by = select_by)
    check_node_opacity(node_opacity = node_opacity)
    
    # unpack g
    ig <- g$graph
    cs <- as_data_frame(ig, what = "vertices")
    is_jg <- g$joint_graph
    
    if(!show_singletons){
        k <- which(degree(ig) == 0 & V(ig)$clone_size <= 1)
        if(length(k)!=0) {
            ig <- delete_vertices(ig, k)
        }
    }
    
    ig <- config_vertices_plot(g = ig, is_jg = is_jg, 
                               node_opacity = node_opacity)
    if(as_visnet == FALSE) {
        plot(ig, vertex.label = NA)
    }
    if(as_visnet == TRUE) {
        V(ig)$size <- V(ig)$size*5
        if(length(E(ig))==0) {
            message("no edges found in the graph, dummy self-edge added")
            # apparently if no edges, visnetwork can't plot
            ig <- add_edges(graph = ig, edges = c(V(ig)[1],V(ig)[1]))
        } 
        else {
            if(all(E(ig)$width == 1)==FALSE) {
                E(ig)$width <- 2*(1/(1+exp(-(-6-0.1*E(ig)$weight))))
            }
        }
        
        vi <- visIgraph(igraph = ig,
                        idToLabel = TRUE,
                        layout = "layout_components",
                        randomSeed = 1234,
                        physics = FALSE,
                        smooth = FALSE,
                        type = "square")
        
        # set vertex titles -> CDR sequences
        chains <- get_chains(x = colnames(vi$x$nodes))
        vi$x$nodes$title <- apply(X = vi$x$nodes[, chains, drop = FALSE], 
                                  y = chains, 
                                  MARGIN = 1, FUN = function(x, y) {
                                      paste0(paste0("<b>", y, "</b>", ':', x), 
                                             collapse = '<br>')})
        
        
        u <- unique(unlist(strsplit(x=unique(cs[, select_by]),split = ',')))
        vi <- vi %>% visOptions(selectedBy = list(variable = select_by,
                                                  values = u,
                                                  multiple = TRUE))
        
        return(vi)
    }
}

get_intergraph_edges <- function(igs,
                                 chains, 
                                 cores, 
                                 trim_flank_aa,
                                 gmi) {
    
    get_identical_hits <- function(a, b) {
        js <- intersect(a$Seq, b$Seq)
        js <- js[is.na(js)==FALSE]
        if(length(js)==0) {
            return(NULL)
        }
        is <- lapply(X = js, a = a, b = b, FUN = function(x, a, b) {
            return(expand.grid(a$Id[which(a$Seq==x)], 
                               b$Id[which(b$Seq==x)]))
        })
        is <- do.call(rbind, is)
        colnames(is) <- c("QueryId", "TargetId")
        return(is)
    }
    
    get_bscore_trim <- function(x, s1, s2, bm, d, trim_flank_aa) {
        
        a <- s1$Seq[d$QueryId[x]]
        b <- s2$Seq[d$TargetId[x]]
        na <- nchar(a)
        nb <- nchar(b)
        if(is.na(na)|is.na(nb)) {
            return(NA)
        }
        if((na-2*trim_flank_aa)<=0 | (nb-2*trim_flank_aa)<=0) {
            return(NA)
        }
        
        a <- substr(x=a, start = trim_flank_aa+1, stop = nchar(a)-trim_flank_aa)
        b <- substr(x=b, start = trim_flank_aa+1, stop = nchar(b)-trim_flank_aa)
        
        if(is.na(a)|is.na(b)) {
            return(NA)
        }
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
    
    get_bscore <- function(x, s1, s2, bm, d) {
        a <- s1$Seq[d$QueryId[x]]
        b <- s2$Seq[d$TargetId[x]]
        
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
    
    get_blastr <- function(s1, s2, chain, trim_flank_aa, gmi) {
        s1 <- data.frame(Id = 1:nrow(s1), Seq = s1[,chain], name = s1$name,
                         len = nchar(s1[, chain]))
        s2 <- data.frame(Id = 1:nrow(s2), Seq = s2[,chain], name = s2$name,
                         len = nchar(s2[, chain]))
        
        o <- blast(query = s1, 
                   db = s2,
                   maxAccepts = 10^4,
                   minIdentity = gmi,
                   alphabet = "protein",
                   output_to_file = FALSE)
        
        o <- o[is.na(o$QueryMatchSeq)==FALSE&
                   is.na(o$TargetMatchSeq)==FALSE,]
        
        # if empty stop
        if(nrow(o)==0) {
            return(NULL)
        } else {
            # we only need these columns
            o <- o[, c("QueryId", "TargetId")]
            
            # here make sure that short identical sequences are not thrown out
            # if duplicated -> remove
            oi <- get_identical_hits(a = s1, b = s2)
            if(is.null(oi)==FALSE) {
                o <- rbind(o, oi)
                o <- o[duplicated(o)==FALSE, ]
                rm(oi)
            }
        }
        
        # get blosum matrix from pwalign
        data_env <- new.env(parent = emptyenv())
        data("BLOSUM62", envir = data_env, package = "pwalign")
        # ensure: min(BLOSUM62)=0, max=15
        data_env[["BLOSUM62"]] <- data_env[["BLOSUM62"]] + 4
        
        # compute BLSOUM62 score for matches
        o$bs <- sapply(X = 1:nrow(o), FUN = get_bscore, d = o, 
                       s1 = s1, s2 = s2, bm = data_env[["BLOSUM62"]])
        
        # compute BLSOUM62 score for matches
        o$core_bs <- o$bs
        if(trim_flank_aa > 0) {
            o$core_bs <- vapply(X = 1:nrow(o), 
                                FUN = get_bscore_trim, 
                                s1 = s1,
                                s2 = s2,
                                d = o,
                                bm = data_env[["BLOSUM62"]],
                                trim_flank_aa = trim_flank_aa,
                                FUN.VALUE = numeric(1))
        }
        
        out <- data.frame(from = s1$name[o$QueryId],
                          to = s2$name[o$TargetId],
                          weight = -o$bs,
                          cweight = -o$core_bs)
        
        len_s1 <- s1$len[o$QueryId]
        len_s2 <- s2$len[o$TargetId]
        
        max_len <- ifelse(test = len_s1 >= len_s2, yes = len_s1, no = len_s2)
        out$max_len <- max_len
        out$max_clen <- max_len-(2*trim_flank_aa)
        out$max_clen <- ifelse(test=out$max_clen<0, yes = 0, no = out$max_clen)
        
        out$nweight <- out$weight/out$max_len
        out$ncweight <- out$cweight/out$max_clen
        
        return(out)
    }
    
    get_igg <- function(x, ix, igs, trim_flank_aa, gmi) {
        # prepare pair-rep data
        s1_name <- ix$name_i[x]
        s2_name <- ix$name_j[x]
        chain <- ix$chain[x]
        s1 <- igs[[ix$index_i[x]]]$clones
        s2 <- igs[[ix$index_j[x]]]$clones
        rm(igs)
        
        message("joining:", x, ":", s1_name, '|', s2_name, "\n")
        
        # run 
        b <- get_blastr(s1 = s1,
                        s2 = s2,
                        chain = chain,
                        trim_flank_aa = trim_flank_aa,
                        gmi = gmi)
        
        if(is.null(b)==FALSE && nrow(b)!=0) {
            b$chain <- chain
            b$sample <- paste0(s1_name, "|", s2_name)
            b$type <- "between-repertoire"
            b$clustering <- "global"
            return(b)
        }
        
        return(NULL)
    }
    
    get_ix <- function(xs, ns, chains) {
        ix <- c()
        for(i in 1:(xs-1)) {
            for(j in (i+1):xs) {
                for(c in chains) {
                    ix <- rbind(ix, data.frame(index_i = i, 
                                               index_j = j,  
                                               name_i = ns[i],
                                               name_j = ns[j],
                                               chain = c))
                }
            }
        }
        return(ix)
    }
    
    # indices 
    ix <- get_ix(xs = length(igs), ns = names(igs), chains = chains)
    # find global similarities between pairs of clone tables
    message("merging clust_irrs: ", nrow(ix), "\n")
    future::plan(future::multisession, workers = I(cores))
    ige <- future_lapply(X = 1:nrow(ix),
                         ix = ix,
                         FUN = get_igg,
                         igs = igs,
                         trim_flank_aa = trim_flank_aa,
                         gmi = gmi,
                         future.seed = TRUE)
    ige <- do.call(rbind, ige)
    
    return(ige)
}

get_clones <- function(sample_id, x) {
    cs <- x
    cs$id <- NULL
    cs$clone_id <- seq_len(nrow(cs))
    cs$sample <- sample_id
    cs$name <- paste0(sample_id, '|', cs$clone_id)
    cs <- cs[, rev(colnames(cs))]
    return(cs)
}

# Description:
# Compare the control lists of individual clust_irr objects. If they match,
# then extract the control list from one list.
get_joint_controls <- function(clust_irrs) {
    # if these tests pass -> extract controls of any of the elements
    for(i in 1:(length(clust_irrs)-1)) {
        for(j in (i+1):length(clust_irrs)) {
            c_i <- get_clustirr_inputs(clust_irrs[[i]])$control
            c_j <- get_clustirr_inputs(clust_irrs[[j]])$control
            
            if(length(c_i)!=length(c_j)) {
                stop("different controls used in individual clust_irrs")
            }
            
            c_ij <- vapply(X = names(c_i), c_i = c_i, c_j = c_j,
                           FUN = function(x, c_i, c_j) {
                               if(is.null(c_i[[x]]) & is.null(c_j[[x]])) {
                                   return(TRUE)
                               }
                               if(is.na(c_i[[x]]) & is.na(c_j[[x]])) {
                                   return(TRUE)
                               }
                               return(c_i[[x]]==c_j[[x]])}, 
                           FUN.VALUE = logical(1))
            
            if(any(c_ij==FALSE)) {
                stop("different controls used in individual clust_irrs")
            }
        } 
    }
    
    return(get_clustirr_inputs(x = clust_irrs[[1]])$control)
}

