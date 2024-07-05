
get_graph <- function(clust_irr, 
                      sample_id,
                      custom_db = NULL,
                      edit_dist = 0) {
    
    get_local_edges <- function(clust_irr, cs) {
        
        get_motif_to_id <- function(x, cs, lp, chain) {
            y <- lp$cdr3[lp$motif == x]
            if(length(y)>0) {
                return(unique(cs$clone_id[cs[,chain] %in% y]))
            }
            return(NA)
        }
        
        el <- vector(mode="list", length = length(get_clustirr_clust(clust_irr)))
        names(el) <- names(get_clustirr_clust(clust_irr))
        for(chain in names(get_clustirr_clust(clust_irr))) {
            lp <- get_clustirr_clust(clust_irr)[[chain]]$local$lp
            cs <- cs[cs[, chain] %in% lp$cdr3,]
            
            if(is.null(lp) == FALSE && nrow(lp) != 0  && nrow(s) != 0) {
                
                lm <- lapply(X = unique(lp$motif), 
                             FUN = get_motif_to_id, 
                             lp = lp, 
                             cs = cs, 
                             chain = chain)
                names(lm) <-  unique(lp$motif)
                
                if(is.null(lm)==FALSE) {
                    el[[chain]] <- lm
                }
            }
        }
        return(el)
    }
    
    get_global_edges <- function(clust_irr, cs) {
        eg <- vector(mode="list", length = length(get_clustirr_clust(clust_irr)))
        names(eg) <- names(get_clustirr_clust(clust_irr))
        for(chain in names(get_clustirr_clust(clust_irr))) {
            g <- get_clustirr_clust(clust_irr)[[chain]]$global
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
                g$key <- apply(X = g[, c("from_clone_id", "to_clone_id")], MARGIN = 1, 
                               FUN = function(x) {return(paste0(sort(x),collapse='|'))})
                g <- g[duplicated(g$key)==FALSE,]
                g$key <- NULL
                
                out <- data.frame(from_cdr3 = g[,"from_clone_id"], 
                                  to_cdr3 = g[,"to_clone_id"], 
                                  weight = g[,"weight"],
                                  nweight = g[,"nweight"],
                                  cweight = g[,"cweight"],
                                  ncweight = g[,"ncweight"],
                                  max_len = g[, "max_len"],
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
    
    build_graph <- function(le, ge, cs, sample_id, chains) {
        
        add_local_edges <- function(le, ig, sample_id, chain) {
            
            xs <- unlist(lapply(X = le, sample_id = sample_id, chain = chain, 
                                FUN = function(x, sample_id, chain) {
                                    if(length(x)<=1) {
                                        return(NULL)
                                    }
                                    xp <- utils::combn(x = paste0(
                                        sample_id, '|', x), m = 2)
                                    xp <- as.vector(xp)
                                    return(xp)
                                }))
            
            if(length(xs)!=0) {
                ig <- igraph::add_edges(graph = ig,
                                        edges = xs,
                                        weight = 1,
                                        cweight = 1,
                                        nweight = 1,
                                        ncweight = 1,
                                        max_len = NA,
                                        type = "within-repertoire",
                                        chain = chain,
                                        clustering = "local")
            }
            return(ig)
        }
        
        add_global_edges <- function(ge, ig, sample_id, chain) {
            
            get_e <- function(x, sample_id) {
                return(paste0(sample_id, '|', x))
            }
            
            e <- as.vector(apply(X = ge[, c("from_cdr3", "to_cdr3")], 
                                 MARGIN = 1, FUN = get_e, sample_id = sample_id))
            ig <- igraph::add_edges(graph = ig, 
                                    edges = e, 
                                    weight = ge$weight,
                                    cweight = ge$cweight,
                                    nweight = ge$nweight,
                                    ncweight = ge$ncweight,
                                    max_len = ge$max_len,
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
        
        # add local edges
        message("1/2: adding local edges... \n")
        if(length(le)!=0) {
            for(chain in chains) {
                if(length(le[[chain]])!=0) {
                    ig <- add_local_edges(le = le[[chain]], ig = ig, 
                                          sample_id = sample_id, 
                                          chain = chain)
                }
            }
        }
        
        # add global edges
        message("2/2: adding global edges... \n")
        if(is.null(ge)==FALSE && nrow(ge)!=0) {
            for(chain in chains) {
                chain_ge <- ge[ge$chain == chain, ]
                if(nrow(chain_ge)!=0) {
                    ig <- add_global_edges(ge = chain_ge, ig = ig, 
                                           sample_id = sample_id,
                                           chain = chain)
                }
            }
        }
        
        return(ig)
    }
    
    check_clustirr(clust_irr = clust_irr)
    
    custom_db <- check_custom_db(custom_db = custom_db)
    
    check_edit_dist(edit_dist = edit_dist)
    
    # get chains
    chains <- get_chains(x = colnames(get_clustirr_inputs(clust_irr)$s))
    
    # cells
    s <- get_clustirr_inputs(clust_irr)$s
    
    # setting up the sample id
    if(missing(sample_id)==FALSE) {
        sample_id <- sample_id
    } 
    else {
        if("sample_id" %in% names(get_clustirr_inputs(clust_irr))) {
            sample_id <- get_clustirr_inputs(clust_irr)$sample_id
        } else {
            sample_id <- paste0("S", sample(x = 1:10^5, size = 1))
        }
    }
    
    # get clones
    cs <- get_clones(sample_id = sample_id, x = s)
    
    # annotate cs with known CDR3-antigen data 
    cs <- get_db_matches(cs = cs, custom_db = custom_db, 
                         edit_dist = edit_dist)
    
    # get local and global edges between clones
    le <- get_local_edges(clust_irr = clust_irr, cs = cs)
    ge <- get_global_edges(clust_irr = clust_irr, cs = cs)
    
    # build graph with only vertices
    if(length(le)==0 & is.null(ge)) {
        ig <- graph_from_data_frame(d = data.frame(from=cs$name[1], 
                                                   to=cs$name[1]),
                                    directed = FALSE, vertices = cs)
        ig <- delete_edges(ig, edges = 1)
        
        return(list(graph = ig, clones = cs))
    }
    
    # build graph
    message("building graph... \n")
    ig <- build_graph(le = le, ge = ge, cs = cs, 
                      sample_id = sample_id, chains = chains)
    
    
    return(list(graph = ig, clones = cs, joint_graph = FALSE))
}

get_joint_graph <- function(clust_irrs, 
                            cores = 1,
                            custom_db = NULL,
                            edit_dist = 0) {
    
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
        gmd <- numeric(length = length(clust_irrs))
        for(i in 1:length(clust_irrs)) {
            check_clustirr(clust_irr = clust_irrs[[i]])
            gmd[i]<-get_clustirr_inputs(clust_irrs[[i]])$control$global_max_hdist
        }
        if(length(unique(gmd))!=1) {
            stop("all global_max_hdist should be equal")
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
    
    # check edit dist
    check_edit_dist(edit_dist = edit_dist)
    
    # get clust_irrs names
    clust_irrs <- get_clust_irrs_names(clust_irrs = clust_irrs)
    
    # get joint controls
    ctrl <- get_joint_controls(clust_irrs = clust_irrs)
    
    # get igs
    message("start: graph creation... \n")
    igs <- bplapply(X = clust_irrs,
                    FUN = get_graph,
                    custom_db = custom_db,
                    edit_dist = edit_dist,
                    BPPARAM = MulticoreParam(workers = cores))
    names(igs) <- names(clust_irrs)
    
    # get chains
    chains <- get_chains(x = colnames(get_clustirr_inputs(clust_irrs[[1]])$s))
    # get intergraph edges (global)
    if(ctrl$global_smart==FALSE) {
        ige <- get_intergraph_edges_hamming(
            igs = igs, 
            global_max_hdist=ctrl$global_max_hdist, 
            chains = chains, 
            cores = cores)
    }
    else {
        ige <- get_intergraph_edges_blosum(
            igs = igs,
            chains = chains,
            cores = cores,
            trim_flank_aa = ctrl$trim_flank_aa,
            global_min_identity=ctrl$global_min_identity)
    }
    # get the vertices/edges of the graph
    df_v <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "vertices"))
    df_e <- do.call(rbind, lapply(X = igs, FUN = get_v_e, what = "edges"))
    
    # these are the cols we want to keep in this order
    cols <- c("from", "to", "weight", "cweight", "nweight", "ncweight",
              "max_len", "type", "chain", "clustering")
    
    
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
    
    return(list(graph = g, clones = df_v, joint_graph = TRUE))
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
        if(length(g)!=3) {
            stop("missing input g with two objects")
        }
        if(all(names(g)==c("graph", "clones", "joint_graph"))==FALSE) {
            stop("wrong input g")
        }
        if(is.data.frame(g$clones)==FALSE) {
            stop("wrong input g")
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
    cs <- g$clones
    ig <- g$graph
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


get_intergraph_edges_hamming <- function(igs, 
                                         global_max_hdist, 
                                         chains, 
                                         cores) {
    
    get_igg <- function(x, i, igs, global_max_hdist, chain) {
        
        get_hd_row <- function(x, 
                               id_x, 
                               id_y, 
                               seq_x, 
                               seq_y,
                               sample_x,
                               sample_y,
                               global_max_hdist) {
            d <- stringdist(a = seq_x[x], b = seq_y, method = "hamming")
            js <- which(d <= global_max_hdist)
            if(length(js) == 0) {
                return(NULL)
            }
            return(data.frame(from = id_x[x], 
                              to = id_y[js],
                              sample_x = sample_x, 
                              sample_y = sample_y))
        }
        
        get_hd <- function(x, 
                           id_x, 
                           id_y, 
                           seq_x, 
                           seq_y, 
                           len_x, 
                           len_y, 
                           sample_x,
                           sample_y,
                           global_max_hdist) {
            
            is_x <- which(len_x == x)
            is_y <- which(len_y == x)
            
            if(length(is_x)==0|length(is_y)==0) {
                return(NULL)
            }
            
            hd <- lapply(X = seq_along(is_x),
                         FUN = get_hd_row,
                         id_x = id_x[is_x], 
                         id_y = id_y[is_y], 
                         seq_x = seq_x[is_x], 
                         seq_y = seq_y[is_y],
                         sample_x = sample_x,
                         sample_y = sample_y,
                         global_max_hdist = global_max_hdist)
            hd <- do.call(rbind, hd)
            return(hd)
        }
        
        s1_name <- names(igs)[i]
        s2_name <- names(igs)[x]
        
        s1 <- igs[[i]]$clones
        s2 <- igs[[x]]$clones
        
        seq_x <- s1[,chain]
        seq_y <- s2[,chain]
        id_x <- s1[,"name"] 
        id_y <- s2[,"name"]
        len_x <- nchar(seq_x)
        len_y <- nchar(seq_y)
        
        hd <- lapply(X = unique(c(len_x, len_y)),
                     FUN = get_hd,
                     id_x = id_x,
                     id_y = id_y,
                     seq_x = seq_x,
                     seq_y = seq_y,
                     len_x = len_x,
                     len_y = len_y,
                     sample_x = s1_name,
                     sample_y = s2_name,
                     global_max_hdist = global_max_hdist)
        
        hd <- do.call(rbind, hd)
        if(is.null(hd)==FALSE && nrow(hd)!=0) {
            hd$chain <- chain
            hd$sample <- paste0(hd$sample_x, "|", hd$sample_y)
            hd$sample_x <- NULL
            hd$sample_y <- NULL
            hd$weight <- 1
            hd$cweight <- 1
            hd$nweight <- 1
            hd$ncweight <- 1
            hd$max_len <- NA
            hd$type <- "between-repertoire"
            hd$chain <- chain
            hd$clustering <- "global"
            return(hd)
        }
        return(NULL)
    }
    
    # find global similarities between pairs of clone tables
    ige <- vector(mode = "list", length = length(chains)*(length(igs)-1))
    
    count <- 1
    for(i in 1:(length(igs)-1)) {
        message("merging clust_irr index: ", i, "/", (length(igs)-1), "\n")
        for(chain in chains) {
            ige[[count]] <- do.call(rbind,
                                    lapply(X = (i+1):length(igs), 
                                           i = i,
                                           FUN = get_igg,
                                           igs = igs,
                                           chain = chain,
                                           global_max_hdist = global_max_hdist))
            count <- count + 1
        }
    }
    ige <- do.call(rbind, ige)
    
    return(ige)
}

get_intergraph_edges_blosum <- function(igs,
                                        chains, 
                                        cores, 
                                        trim_flank_aa,
                                        global_min_identity) {
    
    get_identical_hits <- function(a, b) {
        js <- intersect(a$Seq, b$Seq)
        if(length(js)==0) {
            return(NULL)
        }
        is <- lapply(X = js, a = a, b = b, FUN = function(x, a, b) {
            return(expand.grid(a$Id[a$Seq==x], b$Id[b$Seq==x]))
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
        if((na-2*trim_flank_aa)<=0 | (nb-2*trim_flank_aa)<=0) {
            return(NA)
        }
        
        a <- substr(x = a, start = trim_flank_aa+1, stop = nchar(a)-trim_flank_aa)
        b <- substr(x = b, start = trim_flank_aa+1, stop = nchar(b)-trim_flank_aa)
        
        return(stringDist(x = c(a, b),
                          method = "substitutionMatrix", 
                          type = "global", 
                          substitutionMatrix = bm, 
                          gapOpening = 10,
                          gapExtension = 4))
    }
    
    get_bscore <- function(x, s1, s2, bm, d) {
        return(stringDist(x = c(s1$Seq[d$QueryId[x]], s2$Seq[d$TargetId[x]]),
                          method = "substitutionMatrix", 
                          type = "global", 
                          substitutionMatrix = bm, 
                          gapOpening = 10,
                          gapExtension = 4))
    }
    
    get_blastr <- function(s1, s2, chain, trim_flank_aa, global_min_identity) {
        s1 <- data.frame(Id = 1:nrow(s1), Seq = s1[,chain], name = s1$name, 
                         len = nchar(s1[, chain]))
        s2 <- data.frame(Id = 1:nrow(s2), Seq = s2[,chain], name = s2$name,
                         len = nchar(s2[, chain]))
        
        o <- blast(query = s1, 
                   db = s2,
                   maxAccepts = 10^4,
                   minIdentity = global_min_identity,
                   alphabet = "protein", 
                   output_to_file = FALSE)
        
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
        
        # data("BLOSUM62", envir = data_env, package = "pwalign")
        data("BLOSUM62", envir = data_env, package = "Biostrings")
        
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
        out$nweight <- out$weight/out$max_len
        out$ncweight <- out$cweight/out$max_len
        
        return(out)
    }
    
    get_igg <- function(x, i, igs, chain, trim_flank_aa, global_min_identity) {
        s1_name <- names(igs)[i]
        s2_name <- names(igs)[x]
        
        b <- get_blastr(s1 = igs[[i]]$clones, 
                        s2 = igs[[x]]$clones, 
                        chain = chain, 
                        trim_flank_aa = trim_flank_aa,
                        global_min_identity = global_min_identity)
        
        if(is.null(b)==FALSE && nrow(b)!=0) {
            b$chain <- chain
            b$sample <- paste0(s1_name, "|", s2_name)
            b$type <- "between-repertoire"
            b$clustering <- "global"
            return(b)
        }
        
        return(NULL)
    }
    
    # find global similarities between pairs of clone tables
    ige <- vector(mode = "list", length = length(chains)*(length(igs)-1))
    
    count <- 1
    for(i in 1:(length(igs)-1)) {
        message("merging clust_irr index: ", i, "/", (length(igs)-1), "\n")
        for(chain in chains) {
            # ige[[count]] <- do.call(rbind, lapply(
            #     X = (i+1):length(igs), 
            #     i = i,
            #     FUN = get_igg,
            #     igs = igs,
            #     trim_flank_aa = trim_flank_aa,
            #     global_min_identity = global_min_identity,
            #     chain = chain))
            ige[[count]] <- do.call(rbind, bplapply(
                X = (i+1):length(igs),
                i = i,
                FUN = get_igg,
                igs = igs,
                trim_flank_aa = trim_flank_aa,
                global_min_identity = global_min_identity,
                chain = chain,
                BPPARAM = MulticoreParam(workers = cores)))
            count <- count + 1
        }
    }
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
# integrate clustirr with data from databases: VDJdb, tcr3d, mcpas-tcr
get_db_matches <- function(cs, custom_db, edit_dist) {
    
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
    
    if(is.null(custom_db)==FALSE) {
        db[["custom"]] <- custom_db
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
                                      d = edit_dist,
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

