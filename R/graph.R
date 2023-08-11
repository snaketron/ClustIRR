get_edges <- function(clust_irr) {
    check_clustirr(clust_irr = clust_irr)

    le <- get_local_edges(clust_irr = clust_irr)
    ge <- get_global_edges(clust_irr = clust_irr)
    e <- base::rbind(le, ge)

    if (base::is.null(e) || base::nrow(e) == 0) {
        base::warning("No local or global edges found \n")
        return(NULL)
    }

    me <- base::vector(mode = "list", length = base::length(clust_irr$clust))
    base::names(me) <- base::names(clust_irr$clust)
    for (chain in base::names(clust_irr$clust)) {
        tmp_d <- clust_irr$inputs$s[, c(chain, "id")]
        # remove clonal expanded cdr3s for version 1 and 2
        if (clust_irr$inputs$version != 3) {
            tmp_d <- tmp_d[!base::duplicated(tmp_d[[chain]]), ]
        }
        tmp_e <- e[e$chain == chain, ]

        if (base::nrow(tmp_e) != 0 & base::nrow(tmp_d) != 0) {
            base::colnames(tmp_d) <- c("CDR3", "id")
            me[[chain]] <- base::merge(
                x = base::merge(
                    x = tmp_e, y = tmp_d, by.x = "from_cdr3", by.y = "CDR3"
                ),
                y = tmp_d, by.x = "to_cdr3", by.y = "CDR3"
            )
        }
    }
    me <- base::do.call(base::rbind, me)
    if (base::is.null(me) == FALSE && base::nrow(me) != 0) {
        me$from <- me$id.x
        me$to <- me$id.y
        me <- me[, c(
            "from", "to", "from_cdr3", "to_cdr3",
            "motif", "type", "chain"
        )]
        # remove self-reference connections to prevent circles (only v3)
        me <- me[!(me$from == me$to), ]
        return(me)
    }
    base::warning("No local or global edges found \n")
    return(NULL)
}


get_local_edges <- function(clust_irr) {
    get_lp <- function(x, lp, chain) {
        if (base::sum(lp$motif == x) > 1) {
            p <- base::t(utils::combn(x = lp$cdr3[lp$motif == x], m = 2))
            return(base::data.frame(
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
        if (base::is.null(lp) == FALSE && base::nrow(lp) != 0) {
            edges_local[[chain]] <- base::do.call(base::rbind, base::lapply(
                X = base::unique(lp$motif), FUN = get_lp, lp = lp,
                chain = chain
            ))
        }
    }
    edges_local <- base::do.call(base::rbind, edges_local)

    return(base::unique(edges_local))
}


get_global_edges <- function(clust_irr) {
    get_diff_str <- function(m) {
        c1 <- base::strsplit(m[1], "")[[1]]
        c2 <- base::strsplit(m[2], "")[[1]]
        c1[c1 != c2] <- "-"
        return(base::paste(c1, collapse = ""))
    }

    get_gp <- function(gp, chain) {
        return(base::data.frame(
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
        if (base::is.null(g) == FALSE && base::nrow(g) != 0) {
            edges_global[[chain]] <- get_gp(gp = g, chain = chain)
        }
    }
    edges_global <- base::do.call(base::rbind, edges_global)
    return(edges_global)
}


get_graph <- function(clust_irr) {
    check_clustirr(clust_irr = clust_irr)

    edges <- get_edges(clust_irr = clust_irr)

    if (base::is.null(edges)) {
        base::warning("No local or global edges to build igraph from \n")
        return(NULL)
    }

    ig <- igraph::graph_from_data_frame(edges,
        directed = TRUE
    )
    return(ig)
}


plot_graph <- function(clust_irr, expand_clones = FALSE) {
    check_clustirr(clust_irr = clust_irr)
    ig <- get_graph(clust_irr = clust_irr)
    if (base::is.null(ig)) {
        base::warning("No graph to plot \n")
        return(NULL)
    }
    edges <- igraph::as_data_frame(ig, what = "edges")
    nodes <- igraph::as_data_frame(ig, what = "vertices")
    chains <- base::names(clust_irr$clust)
    types <- base::unique(edges$type)
    edges <- configure_edges(edges = edges, chains = chains, types = types)
    nodes <- configure_nodes(
        nodes = nodes, edges = edges, chains = chains,
        types = types, s = clust_irr$inputs$s
    )
    if (!expand_clones) {
        nodes <- nodes[!base::duplicated(nodes$label), ]
        edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, ]
    }
    if(base::nrow(edges) != 0){
        ledges <- base::data.frame(
            color = base::unique(edges$color),
            label = base::unique(edges$type),
            arrows = "", width = 4
        )
    } else {
        ledges <- base::data.frame()
    }
    if (base::length(chains) > 1) {
        chains <- base::append(chains, "Both chains")
    }
    lnodes <- base::data.frame(
        label = chains,
        color = "", shape = "dot", size = 8
    )
    lnodes$color <- base::ifelse(test = lnodes$label == "CDR3b" |
        lnodes$label == "CDR3g" |
        lnodes$label == "CDR3h",
    yes = "blue",
    no = "yellow"
    )
    lnodes$color <- base::ifelse(test = lnodes$label == "Both chains",
        yes = "green",
        no = lnodes$color
    )
    shapes <- base::c("dot", "diamond")
    labels <- base::c("Expanded", "Singleton")
    for (i in base::seq_along(shapes)) {
        if (shapes[i] %in% base::unique(nodes$shape)) {
            lnodes <- base::rbind(
                lnodes,
                base::data.frame(
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
        t <- base::c(base::as.numeric(x["from"]), base::as.numeric(x["to"]))
        return(base::paste(base::min(t), base::max(t)))
    }
    edges$edge_id <- base::apply(X = edges, MARGIN = 1, FUN = get_edge_id)
    edges <- base::unique(edges[, c("edge_id", "motif", "type", "chain")])
    m <- base::data.frame(edge_id = unique(edges$edge_id))
    for (chain in chains) {
        for (type in types) {
            df <- edges[
                edges$chain == chain & edges$type == type,
                c("edge_id", "motif")
            ]
            t <- base::paste(chain, type, sep = "_")
            c <- base::paste(chain, type, "count", sep = "_")
            if (base::nrow(df) > 0) {
                df <- base::do.call(
                    base::cbind.data.frame,
                    stats::aggregate(motif ~ edge_id,
                        data = df,
                        FUN = function(x) {
                            count <- base::length(x)
                            c(
                                motifs = base::paste(x, collapse = ", "),
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

            base::colnames(df) <- c("edge_id", t, c)
            m <- base::merge(m, df, by = "edge_id", all = TRUE)

            if (t %in% base::colnames(m)) {
                m[[t]][base::is.na(m[[t]])] <- "-"
            }
            if (c %in% base::colnames(m)) {
                m[[c]][base::is.na(m[[c]])] <- 0
                m[[c]] <- base::as.numeric(m[[c]])
            }
        }
    }
    edges <- m
    l <- base::colnames(edges)[grepl("local_count", base::colnames(edges))]
    edges$local_count <- base::rowSums(edges[, l, drop = FALSE])
    g <- base::colnames(edges)[grepl("global_count", base::colnames(edges))]
    edges$global_count <- base::rowSums(edges[, g, drop = FALSE])
    edges$total_count <- edges$global_count + edges$local_count
    edges$from <-
        base::as.numeric(base::apply(X = edges, MARGIN = 1, function(x) {
            base::strsplit(x["edge_id"], " ")[[1]][1]
        }))
    edges$to <-
        base::as.numeric(base::apply(X = edges, MARGIN = 1, function(x) {
            base::strsplit(x["edge_id"], " ")[[1]][2]
        }))
    edges$length <- 15
    edges$width <- edges$global_count * 5 + edges$local_count
    edges$type <- base::ifelse(test = (edges$global_count > 0),
        yes = "global",
        no = "local"
    )
    edges$type <- base::ifelse(test = (edges$global_count > 0 &
        edges$local_count > 0),
    yes = "local & global",
    no = edges$type
    )

    edges$color <- base::ifelse(test = (edges$type == "local"),
        yes = "gray",
        no = "orange"
    )
    edges$color <- base::ifelse(test = (edges$type == "local & global"),
        yes = "#9A0000",
        no = edges$color
    )

    edges$title <- base::apply(
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
        if (c %in% base::c("CDR3a", "CDR3d", "CDR3l")) {
            m <- m_y
        }
        g <- base::paste(c, "global", sep = "_")
        l <- base::paste(c, "local", sep = "_")
        if (g %in% base::names(x)) {
            if (x[g] != "-") {
                t_gl <-
                    base::paste0("<b>Global:</b><br>", m, x[g], "</mark><br>")
                gl <- TRUE
            }
        }
        if (l %in% base::names(x)) {
            if (x[l] != "-") {
                t_lo <-
                    base::paste0("<b>Local:</b><br>", m, x[l], "</mark><br>")
                lo <- TRUE
            }
        }
        if (gl | lo) {
            if (c %in% base::c("CDR3b", "CDR3g", "CDR3h")) {
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
            if (c %in% base::c("CDR3a", "CDR3d", "CDR3l")) {
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
        return(base::paste0(tbg, tbg_g, tbg_l, "<br>", tad, tad_g, tad_l))
    }
    if (bg & !ad) {
        return(base::paste0(tbg, tbg_g, tbg_l))
    }
    if (!bg & ad) {
        return(base::paste0(tad, tad_g, tad_l))
    }
    return("")
}


configure_nodes <- function(nodes, edges, chains, types, s) {
    base::names(nodes) <- "id"
    nodes <- base::merge(nodes, s, by = "id")
    nodes[base::is.na(nodes)] <- "<NA>"
    s[base::is.na(s)] <- "<NA>"
    min_size <- 20
    clone_boost <- 10
    if (base::length(chains) > 1) {
        dec <- TRUE
        if (base::any(chains %in% base::c("CDR3h", "CDR3l"))) {
            dec <- FALSE
        }
        chains <- base::sort(chains, decreasing = dec)
        bg <- substr(chains[1], 5, 5)
        ad <- substr(chains[2], 5, 5)
        nodes$label <- base::paste(nodes[[chains[1]]], " (", bg, ") - ",
            nodes[[chains[2]]], " (", ad, ")",
            sep = ""
        )
        nodes$clone_count <- base::apply(X = nodes, MARGIN = 1, function(x) {
            base::sum(s[chains[1]] == x[chains[1]] & 
                          s[chains[2]] == x[chains[2]])
        })
    } else {
        nodes$label <- nodes[[chains]]
        nodes$clone_count <- apply(X = nodes, MARGIN = 1, function(x) {
            base::sum(s[chains] == x[chains])
        })
    }
    nodes$size <- base::log2(nodes$clone_count) * clone_boost + min_size
    nodes$color.border <- "black"
    nodes$color.highlight <- "red"
    l <- list()
    for (c in chains) {
        l <- base::c(l, base::paste0(c, "_", types))
    }
    nodes$group <- base::apply(
        X = nodes, MARGIN = 1, FUN = get_u_motifs,
        edges = edges, l = l
    )
    for (i in base::seq_along(l)) {
        nodes[l[[i]]] <- base::unlist(lapply(nodes$group, `[`, i))
        nodes[nodes == ""] <- "-"
    }
    nodes$group <- base::apply(X = nodes, MARGIN = 1, function(x) {
        u <- base::unlist(x["group"])
        base::paste(u[base::nzchar(u)], collapse = ", ")
    })
    nodes$shape <- base::ifelse(test = nodes$clone_count > 1,
        yes = "dot", no = "diamond"
    )
    nodes$color.background <- base::apply(
        X = nodes, MARGIN = 1, FUN = get_color,
        l = l
    )
    nodes$title <- base::apply(
        X = nodes, MARGIN = 1, FUN = set_node_title,
        chains = chains
    )

    nodes$shadow <- FALSE
    nodes <- nodes[order(nodes$label), ]
    return(nodes)
}


get_u_motifs <- function(x, edges = edges, l = l) {
    get_unique_str <- function(x) {
        res <- base::unique(e[x][e[x] != "-"])
        return(
            base::paste(
                base::unique(
                    base::unlist(base::strsplit(res, ", "))), collapse = ", "))
    }

    id <- base::as.numeric(x["id"])
    e <- edges[edges$from == id | edges$to == id, ]
    return(base::lapply(X = l, FUN = get_unique_str))
}


get_color <- function(x, l) {
    t <- base::vector(mode = "character")

    for (c in l) {
        if (x[[c]] != "-") {
            t <- base::c(t, c)
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
    mb_b <- base::paste0(mb, "<i>(&beta;)</i></mark><br>", m)
    my_a <- base::paste0(my, "<i>(&alpha;)</i>")
    mb_g <- base::paste0(mb, "<i>(&gamma;)</i></mark><br>", m)
    my_d <- base::paste0(my, "<i>(&delta;)</i>")
    mb_h <- base::paste0(mb, "<i>(H)</i></mark><br>", m)
    my_l <- base::paste0(my, "<i>(L)</i>")

    l <- x[["label"]]
    l <- base::gsub(pattern = " \\(b\\) - ", replacement = mb_b, x = l)
    l <- base::gsub(pattern = " \\(a\\)", replacement = my_a, x = l)
    l <- base::gsub(pattern = " \\(g\\) - ", replacement = mb_g, x = l)
    l <- base::gsub(pattern = " \\(d\\)", replacement = my_d, x = l)
    l <- base::gsub(pattern = " \\(h\\) - ", replacement = mb_h, x = l)
    l <- base::gsub(pattern = " \\(l\\)", replacement = my_l, x = l)
    l <- base::paste0(m, "<b>", l, "</b></mark><br><br>")

    cc <- x[["clone_count"]]
    cc <- base::paste0("<b>Clone count:", cc, "</b><br>")

    m <- get_motif_list(x, chains)
    return(base::paste0(l, cc, "<br>", m))
}


configure_network <- function(nodes, edges, ledges, lnodes) {
    id_style <- base::paste(
        "'width: ",
        base::max(base::nchar(nodes$label)) * 7.5,
        "px;'"
    )
    return(
        visNetwork::visNetwork(nodes = nodes, edges = edges) %>%
            visNetwork::visIgraphLayout(
                layout = "layout_components",
                randomSeed = 1234
            ) %>%
            visNetwork::visOptions(
                highlightNearest =
                    base::list(
                        enabled = TRUE,
                        degree = 1,
                        algorithm = "hierarchical"
                    ),
                selectedBy = base::list(
                    variable = "group",
                    multiple = TRUE,
                    sort = TRUE
                ),
                manipulation = FALSE,
                nodesIdSelection = base::list(
                    enabled = TRUE,
                    style = id_style
                )
            ) %>%
            visNetwork::visLegend(
                addEdges = ledges, addNodes = lnodes,
                useGroups = FALSE, position = "right",
                width = 0.15, zoom = FALSE
            )
    )
}


check_clustirr <- function(clust_irr) {
    if (!methods::is(clust_irr, "clust_irr")) {
        base::stop("Input has to be object of class clust_irr")
    }
}
