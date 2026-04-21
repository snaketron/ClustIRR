
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
    
    # remove gaps if they exceed specific threshold on position
    aln <- msa::msaClustalOmega(inputSeqs = AAStringSet(cdrs), type = "protein")
    aln <- as.character(aln)
    alnm <- do.call(rbind, strsplit(aln, split = ""))
    gap_rm <- apply(X = alnm, MARGIN = 2, gap = gap_remove_prob, 
                    FUN = function(x, gap) {
        return(sum(x=="-")/length(x) >= gap)
    })
    if(any(gap_rm)) {
        alnm <- alnm[, which(gap_rm), drop = FALSE]
    }
    aln <- apply(alnm, 1, paste0, collapse = "")
    
    g <- ggseqlogo(aln, seq_type = 'aa', method = "probability")+
        theme_logo(base_size = 10)+
        theme(legend.position = "none")+
        scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1))
    
    return(g)
}
