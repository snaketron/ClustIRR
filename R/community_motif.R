
get_cdr3_motifs <- function(node_summary, community_id, gap_remove_prob = 0.8) {
    
    
    chains <- get_chains(colnames(node_summary))
    
    i <- which(node_summary$community == community_id)
    if(length(i)==0) {
        stop("community not found")
    }
    else if(length(i)==1) {
        stop(paste0("community has only one clonotype: ", 
                    node_summary[, c(chains, clone_size)]))
    } else {
        gs <- lapply(X = chains, FUN = get_motif, 
                     node_summary = node_summary[i, ],
                     gap_remove_prob = gap_remove_prob)
        return(patchwork::wrap_plots(gs, nrow = 1))
    }
}

get_motif <- function(chain, node_summary, gap_remove_prob) {
    cdrs <- node_summary[,chain]
    cells <- node_summary[,"clone_size"]
    
    i <- which(is.na(cdrs))
    if(length(i)!=0) {
        cdrs <- cdrs[-i]
        cells <- cells[-i]
        if(length(cdrs)==0) {
            warning("no CDRs after removing NAs")
            return(ggplot() + theme_bw())
        }
    }
    
    aln <- msa::msaClustalOmega(inputSeqs = AAStringSet(cdrs), type = "protein")
    q <- data.frame(s = as.character(aln))
    q <- gglogo::ggfortify(q, s, method = "frequency")
    q <- q[!q$position %in% q$position[which(q$info >= gap_remove_prob 
                                             & q$element == "-")], ]
    q$pos <- as.numeric(as.factor(as.numeric(as.character(q$position))))
    
    i <- which.max(cells)
    top <- data.frame(chain = chain, cdr = cdrs[i], cells = cells[i])
    
    g <- ggplot(data = q)+
        gglogo::geom_logo(aes(x = position, y = info, group = element, 
                              label = element), fill = "black", col = "black",
                          alpha = 0.1, position = "classic")+
        theme_bw(base_size = 10)+
        theme(legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.x=element_blank())+
        ylab(label = "prob.")+
        scale_fill_manual(values = c("white", "red"))+
        scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))+
        ggtitle(label = paste0(chain, ": ", top$cdr, 
                               " (cells=", top$cells, ")"))
    
    return(g)
}
