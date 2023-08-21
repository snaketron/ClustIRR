get_edges <- function(clust_irr) {
    check_clustirr(clust_irr = clust_irr)

    le <- get_local_edges(clust_irr = clust_irr)
    ge <- get_global_edges(clust_irr = clust_irr)
    e <- rbind(le, ge)

    if (is.null(e) || nrow(e) == 0) {
        warning("No local or global edges found \n")
        return(NULL)
    }

    me <- vector(mode = "list", length = length(clust_irr$clust))
    names(me) <- names(clust_irr$clust)
    for (chain in names(clust_irr$clust)) {
        tmp_d <- clust_irr$inputs$s[, c(chain, "id")]
        # remove clonal expanded cdr3s for version 1
        if (clust_irr$inputs$version != 2) {
            tmp_d <- tmp_d[!duplicated(tmp_d[[chain]]), ]
        }
        tmp_e <- e[e$chain == chain, ]

        if (nrow(tmp_e) != 0 & nrow(tmp_d) != 0) {
            colnames(tmp_d) <- c("CDR3", "id")
            me[[chain]] <- merge(
                x = merge(
                    x = tmp_e, y = tmp_d, by.x = "from_cdr3", by.y = "CDR3"
                ),
                y = tmp_d, by.x = "to_cdr3", by.y = "CDR3"
            )
        }
    }
    me <- do.call(rbind, me)
    if (is.null(me) == FALSE && nrow(me) != 0) {
        me$from <- me$id.x
        me$to <- me$id.y
        me <- me[, c(
            "from", "to", "from_cdr3", "to_cdr3",
            "motif", "type", "chain"
        )]
        # remove self-reference connections to prevent circles (only v2)
        me <- me[!(me$from == me$to), ]
        return(me)
    }
    warning("No local or global edges found \n")
    return(NULL)
}


get_local_edges <- function(clust_irr) {
    get_lp <- function(x, lp, chain) {
        if (sum(lp$motif == x) > 1) {
            p <- t(combn(x = lp$cdr3[lp$motif == x], m = 2))
            return(data.frame(
                from_cdr3 = p[, 1], to_cdr3 = p[, 2],
                motif = x, type = "local",
                chain = chain
            ))
        }
        return(NULL)
    }

    edges_local <- vector(mode = "list", length = length(clust_irr$clust))
    names(edges_local) <- names(clust_irr$clust)
    for (chain in names(clust_irr$clust)) {
        lp <- clust_irr$clust[[chain]]$local$lp
        if (is.null(lp) == FALSE && nrow(lp) != 0) {
            edges_local[[chain]] <- do.call(rbind, lapply(
                X = unique(lp$motif), FUN = get_lp, lp = lp,
                chain = chain
            ))
        }
    }
    edges_local <- do.call(rbind, edges_local)

    return(unique(edges_local))
}


get_global_edges <- function(clust_irr) {
    get_diff_str <- function(m) {
        c1 <- strsplit(m[1], "")[[1]]
        c2 <- strsplit(m[2], "")[[1]]
        c1[c1 != c2] <- "-"
        return(paste(c1, collapse = ""))
    }

    get_gp <- function(gp, chain) {
        return(data.frame(
            from_cdr3 = gp[, 1], to_cdr3 = gp[, 2],
            motif = apply(gp, 1, get_diff_str),
            type = "global",
            chain = chain
        ))
    }

    edges_global <- vector(mode = "list", length = length(clust_irr$clust))
    names(edges_global) <- names(clust_irr$clust)
    for (chain in names(clust_irr$clust)) {
        g <- clust_irr$clust[[chain]]$global
        if (is.null(g) == FALSE && nrow(g) != 0) {
            edges_global[[chain]] <- get_gp(gp = g, chain = chain)
        }
    }
    edges_global <- do.call(rbind, edges_global)
    return(edges_global)
}


get_graph <- function(clust_irr) {
    check_clustirr(clust_irr = clust_irr)

    edges <- get_edges(clust_irr = clust_irr)

    if (is.null(edges)) {
        warning("No local or global edges to build igraph from \n")
        return(NULL)
    }

    ig <- graph_from_data_frame(edges,
        directed = TRUE
    )
    return(ig)
}


plot_graph <- function(clust_irr, expand_clones = FALSE) {
    check_clustirr(clust_irr = clust_irr)
    ig <- get_graph(clust_irr = clust_irr)
    if (is.null(ig)) {
        warning("No graph to plot \n")
        return(NULL)
    }
    edges <- as_data_frame(ig, what = "edges")
    nodes <- as_data_frame(ig, what = "vertices")
    chains <- names(clust_irr$clust)
    types <- unique(edges$type)
    edges <- configure_edges(edges = edges, chains = chains, types = types)
    nodes <- configure_nodes(
        nodes = nodes, edges = edges, chains = chains,
        types = types, s = clust_irr$inputs$s
    )
    if (!expand_clones) {
        nodes <- nodes[!duplicated(nodes$label), ]
        edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, ]
    }
    if(nrow(edges) != 0){
        ledges <- data.frame(
            color = unique(edges$color),
            label = unique(edges$type),
            arrows = "", width = 4
        )
    } else {
        ledges <- data.frame()
    }
    if (length(chains) > 1) {
        chains <- append(chains, "Both chains")
    }
    lnodes <- data.frame(
        label = chains,
        color = "", shape = "dot", size = 8
    )
    lnodes$color <- ifelse(test = lnodes$label == "CDR3b" |
        lnodes$label == "CDR3g" |
        lnodes$label == "CDR3h",
    yes = "blue",
    no = "yellow"
    )
    lnodes$color <- ifelse(test = lnodes$label == "Both chains",
        yes = "green",
        no = lnodes$color
    )
    shapes <- c("dot", "diamond")
    labels <- c("Expanded", "Singleton")
    for (i in seq_along(shapes)) {
        if (shapes[i] %in% unique(nodes$shape)) {
            lnodes <- rbind(
                lnodes,
                data.frame(
                    label = labels[i],
                    color = "black",
                    shape = shapes[i],
                    size = 8
                )
            )
        }
    }
    return(configure_network(
        nodes = nodes, edges = edges,
        lnodes = lnodes, ledges = ledges
    ))
}



configure_edges <- function(edges, chains, types) {
    get_edge_id <- function(x) {
        t <- c(as.numeric(x["from"]), as.numeric(x["to"]))
        return(paste(min(t), max(t)))
    }
    edges$edge_id <- apply(X = edges, MARGIN = 1, FUN = get_edge_id)
    edges <- unique(edges[, c("edge_id", "motif", "type", "chain")])
    m <- data.frame(edge_id = unique(edges$edge_id))
    for (chain in chains) {
        for (type in types) {
            df <- edges[
                edges$chain == chain & edges$type == type,
                c("edge_id", "motif")
            ]
            t <- paste(chain, type, sep = "_")
            c <- paste(chain, type, "count", sep = "_")
            if (nrow(df) > 0) {
                df <- do.call(
                    cbind.data.frame,
                    aggregate(motif ~ edge_id,
                        data = df,
                        FUN = function(x) {
                            count <- length(x)
                            c(
                                motifs = paste(x, collapse = ", "),
                                counts = count
                            )
                        }
                    )
                )
            } else {
                df <- data.frame(
                    a = character(), b = character(),
                    c = character()
                )
            }

            colnames(df) <- c("edge_id", t, c)
            m <- merge(m, df, by = "edge_id", all = TRUE)

            if (t %in% colnames(m)) {
                m[[t]][is.na(m[[t]])] <- "-"
            }
            if (c %in% colnames(m)) {
                m[[c]][is.na(m[[c]])] <- 0
                m[[c]] <- as.numeric(m[[c]])
            }
        }
    }
    edges <- m
    l <- colnames(edges)[grepl("local_count", colnames(edges))]
    edges$local_count <- rowSums(edges[, l, drop = FALSE])
    g <- colnames(edges)[grepl("global_count", colnames(edges))]
    edges$global_count <- rowSums(edges[, g, drop = FALSE])
    edges$total_count <- edges$global_count + edges$local_count
    edges$from <-
        as.numeric(apply(X = edges, MARGIN = 1, function(x) {
            strsplit(x["edge_id"], " ")[[1]][1]
        }))
    edges$to <-
        as.numeric(apply(X = edges, MARGIN = 1, function(x) {
            strsplit(x["edge_id"], " ")[[1]][2]
        }))
    edges$length <- 15
    edges$width <- edges$global_count * 5 + edges$local_count
    edges$type <- ifelse(test = (edges$global_count > 0),
        yes = "global",
        no = "local"
    )
    edges$type <- ifelse(test = (edges$global_count > 0 &
        edges$local_count > 0),
    yes = "local & global",
    no = edges$type
    )

    edges$color <- ifelse(test = (edges$type == "local"),
        yes = "gray",
        no = "orange"
    )
    edges$color <- ifelse(test = (edges$type == "local & global"),
        yes = "#9A0000",
        no = edges$color
    )

    edges$title <- apply(
        X = edges, MARGIN = 1, FUN = get_motif_list,
        chains = chains
    )
    edges$arrows <- ""
    edges$dashes <- FALSE
    edges$smooth <- FALSE
    edges$shadow <- FALSE
    return(edges)
}


get_motif_list <- function(x, chains) {
    tbg <- tbg_g <- tbg_l <- tad <- tad_g <- tad_l <- t_gl <- t_lo <- ret <- ""
    bg <- bg_g <- bg_l <- ad <- ad_g <- ad_l <- gl <- lo <- FALSE
    m_y <- "<mark style=\"background-color: yellow; color: black;\">"
    m_b <- "<mark style=\"background-color: blue; color: white;\">"
    for (c in chains) {
        gl <- lo <- FALSE
        m <- m_b
        if (c %in% c("CDR3a", "CDR3d", "CDR3l")) {
            m <- m_y
        }
        g <- paste(c, "global", sep = "_")
        l <- paste(c, "local", sep = "_")
        if (g %in% names(x)) {
            if (x[g] != "-") {
                t_gl <-
                    paste0("<b>Global:</b><br>", m, x[g], "</mark><br>")
                gl <- TRUE
            }
        }
        if (l %in% names(x)) {
            if (x[l] != "-") {
                t_lo <-
                    paste0("<b>Local:</b><br>", m, x[l], "</mark><br>")
                lo <- TRUE
            }
        }
        if (gl | lo) {
            if (c %in% c("CDR3b", "CDR3g", "CDR3h")) {
                tbg_g <- t_gl
                tbg_l <- t_lo
                if (c == "CDR3b") {
                    tbg <- "<b>CDR3-&beta;</b><br>"
                }
                if (c == "CDR3g") {
                    tbg <- "<b>CDR3-&gamma;</b><br>"
                }
                if (c == "CDR3h") {
                    tbg <- "<b>CDR3-H</b><br>"
                }
                bg <- TRUE
            }
            if (c %in% c("CDR3a", "CDR3d", "CDR3l")) {
                tad_g <- t_gl
                tad_l <- t_lo
                if (c == "CDR3a") {
                    tad <- "<b>CDR3-&alpha;</b><br>"
                }
                if (c == "CDR3d") {
                    tad <- "<b>CDR3-&delta;</b><br>"
                }
                if (c == "CDR3l") {
                    tad <- "<b>CDR3-L</b><br>"
                }
                ad <- TRUE
            }
        }
    }
    if (bg & ad) {
        return(paste0(tbg, tbg_g, tbg_l, "<br>", tad, tad_g, tad_l))
    }
    if (bg & !ad) {
        return(paste0(tbg, tbg_g, tbg_l))
    }
    if (!bg & ad) {
        return(paste0(tad, tad_g, tad_l))
    }
    return("")
}


configure_nodes <- function(nodes, edges, chains, types, s) {
    names(nodes) <- "id"
    nodes <- merge(nodes, s, by = "id")
    nodes[is.na(nodes)] <- "<NA>"
    s[is.na(s)] <- "<NA>"
    min_size <- 20
    clone_boost <- 10
    if (length(chains) > 1) {
        dec <- TRUE
        if (any(chains %in% c("CDR3h", "CDR3l"))) {
            dec <- FALSE
        }
        chains <- sort(chains, decreasing = dec)
        bg <- substr(chains[1], 5, 5)
        ad <- substr(chains[2], 5, 5)
        nodes$label <- paste(nodes[[chains[1]]], " (", bg, ") - ",
            nodes[[chains[2]]], " (", ad, ")",
            sep = ""
        )
        nodes$clone_count <- apply(X = nodes, MARGIN = 1, function(x) {
            sum(s[chains[1]] == x[chains[1]] &
                          s[chains[2]] == x[chains[2]])
        })
    } else {
        nodes$label <- nodes[[chains]]
        nodes$clone_count <- apply(X = nodes, MARGIN = 1, function(x) {
            sum(s[chains] == x[chains])
        })
    }
    nodes$size <- log2(nodes$clone_count) * clone_boost + min_size
    nodes$color.border <- "black"
    nodes$color.highlight <- "red"
    l <- list()
    for (c in chains) {
        l <- c(l, paste0(c, "_", types))
    }
    nodes$group <- apply(
        X = nodes, MARGIN = 1, FUN = get_u_motifs,
        edges = edges, l = l
    )
    for (i in seq_along(l)) {
        nodes[l[[i]]] <- unlist(lapply(nodes$group, `[`, i))
        nodes[nodes == ""] <- "-"
    }
    nodes$group <- apply(X = nodes, MARGIN = 1, function(x) {
        u <- unlist(x["group"])
        paste(u[nzchar(u)], collapse = ", ")
    })
    nodes$shape <- ifelse(test = nodes$clone_count > 1,
        yes = "dot", no = "diamond"
    )
    nodes$color.background <- apply(
        X = nodes, MARGIN = 1, FUN = get_color,
        l = l
    )
    nodes$title <- apply(
        X = nodes, MARGIN = 1, FUN = set_node_title,
        chains = chains
    )

    nodes$shadow <- FALSE
    nodes <- nodes[order(nodes$label), ]
    return(nodes)
}


get_u_motifs <- function(x, edges = edges, l = l) {
    get_unique_str <- function(x) {
        res <- unique(e[x][e[x] != "-"])
        return(
            paste(
                unique(
                    unlist(strsplit(res, ", "))), collapse = ", "))
    }

    id <- as.numeric(x["id"])
    e <- edges[edges$from == id | edges$to == id, ]
    return(lapply(X = l, FUN = get_unique_str))
}


get_color <- function(x, l) {
    t <- vector(mode = "character")

    for (c in l) {
        if (x[[c]] != "-") {
            t <- c(t, c)
        }
    }

    bg <- "CDR3b_global" %in% t | "CDR3b_local" %in% t |
        "CDR3g_global" %in% t | "CDR3g_local" %in% t |
        "CDR3h_global" %in% t | "CDR3h_local" %in% t

    ad <- "CDR3a_global" %in% t | "CDR3a_local" %in% t |
        "CDR3d_global" %in% t | "CDR3d_local" %in% t |
        "CDR3l_global" %in% t | "CDR3l_local" %in% t

    if (ad & bg) {
        return("green")
    }
    if (bg) {
        return("blue")
    }
    return("yellow")
}


set_node_title <- function(x, chains) {
    m <- "<mark style=\"background-color: white; color: black;\">"
    mb <- "</mark><mark style=\"background-color: blue; color: white;\">"
    my <- "</mark><mark style=\"background-color: yellow; color: black;\">"
    mb_b <- paste0(mb, "<i>(&beta;)</i></mark><br>", m)
    my_a <- paste0(my, "<i>(&alpha;)</i>")
    mb_g <- paste0(mb, "<i>(&gamma;)</i></mark><br>", m)
    my_d <- paste0(my, "<i>(&delta;)</i>")
    mb_h <- paste0(mb, "<i>(H)</i></mark><br>", m)
    my_l <- paste0(my, "<i>(L)</i>")

    l <- x[["label"]]
    l <- gsub(pattern = " \\(b\\) - ", replacement = mb_b, x = l)
    l <- gsub(pattern = " \\(a\\)", replacement = my_a, x = l)
    l <- gsub(pattern = " \\(g\\) - ", replacement = mb_g, x = l)
    l <- gsub(pattern = " \\(d\\)", replacement = my_d, x = l)
    l <- gsub(pattern = " \\(h\\) - ", replacement = mb_h, x = l)
    l <- gsub(pattern = " \\(l\\)", replacement = my_l, x = l)
    l <- paste0(m, "<b>", l, "</b></mark><br><br>")

    cc <- x[["clone_count"]]
    cc <- paste0("<b>Clone count:", cc, "</b><br>")

    m <- get_motif_list(x, chains)
    return(paste0(l, cc, "<br>", m))
}


configure_network <- function(nodes, edges, ledges, lnodes) {
    id_style <- paste(
        "'width: ",
        max(nchar(nodes$label)) * 7.5,
        "px;'"
    )
    return(
        visNetwork(nodes = nodes, edges = edges) %>%
            visIgraphLayout(
                layout = "layout_components",
                randomSeed = 1234
            ) %>%
            visOptions(
                highlightNearest =
                    list(
                        enabled = TRUE,
                        degree = 1,
                        algorithm = "hierarchical"
                    ),
                selectedBy = list(
                    variable = "group",
                    multiple = TRUE,
                    sort = TRUE
                ),
                manipulation = FALSE,
                nodesIdSelection = list(
                    enabled = TRUE,
                    style = id_style
                )
            ) %>%
            visLegend(
                addEdges = ledges, addNodes = lnodes,
                useGroups = FALSE, position = "right",
                width = 0.15, zoom = FALSE
            )
    )
}

