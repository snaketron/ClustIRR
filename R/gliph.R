# Description:
# This is rewritten code for GLIPH v1 and v2 
# Details are missing, TODO
# data_sample: must be a data.frame that has the following columns
# (CDR3b, TRBV, TRBJ, CDR3a, TRAV, TRAJ, sample_id), bare minimum
# (CDR3b, sample_id)
gliph <- function(data_sample,
                  data_ref,
                  version = 2,
                  ks = c(2, 3, 4),
                  cores = 1,
                  control = NULL) {
    
    # 1. parameter check
    
    
    # a. control check
    control <- get_control(control_in = control)
    
    
    
    # gliph v1 and v2 use nonredundant CDR3s for downstream analysis
    
    
    
    # b. filter sample data based on inputs (e.g. remove CDR3 with small size, 
    # remove CDR3s without C and F, ...,)
    
    # check pars
    # check data_sample & data_ref
    # check ks
    # check B
    # check cores
    
    # filter data based on input
    # 
    

    # get chains to be analyzed
    chains <- get_chains(colnames(data_sample))
    
    # run analysis for each chain (if available)
    clust <- vector(mode = "list", length = length(chains))
    names(clust) <- chains
    
    edges <- vector(mode = "list", length = length(chains))
    names(edges) <- chains

    for(chain in chains) {
        
        if(control$trim_flanks) {
            # NAs ignored by qqgram, how about global dist? TODO
            data_sample[, chain] <- get_trimmed_flanks(
                x = data_sample[, chain],
                flank_size = control$flank_size)
            data_ref[, chain] <- get_trimmed_flanks(
                x = data_ref[, chain],
                flank_size = control$flank_size)
        }
        
        # run local + global clustering
        if(version==1) {
            clust[[chain]] <- get_chain_run_v1(
                cdr3 = unique(data_sample[, chain]), # unique -> critical
                cdr3_ref = unique(data_ref[, chain]), # unique -> critical
                ks = ks, 
                cores = cores, 
                B = B,
                control = control)
        } 
        if(version==2) {
            clust[[chain]] <- get_chain_run_v2(
                cdr3 = unique(data_sample[, chain]), # unique -> critical
                cdr3_ref = unique(data_ref[, chain]), # unique -> critical
                ks = ks, 
                cores = cores, 
                control = control)
        }
        if(version==3) {
            # v3 = v2 but we do not use non-redundant set of CDR3s
            clust[[chain]] <- get_chain_run_v2(
                cdr3 = data_sample[, chain], # no-unique
                cdr3_ref = data_ref[, chain], # no-unique
                ks = ks, 
                cores = cores, 
                control = control)
        }
        
        # create edges of a graph
        edges[[chain]] <- get_edges(
            local_pairs = clust[[chain]]$local_pairs,
            global_pairs = clust[[chain]]$global_pairs, 
            cdr3 = ifelse(test = version==3, 
                          yes = data_sample[, chain],
                          no = unique(data_sample[, chain])),
            chain = chain)
    }
    
    return(list(clust = clust, 
                edges = edges,
                data_sample = data_sample))
}



# Description:
# This function will perform the scoring of the clusters based on
# e.g. HLA, V, J, clone size TODO
score <- function(gliph,
                  data_sample,
                  data_meta) {
    
}



# Description:
# This function will create an igraph or visGraph object. TODO
plot_graph <- function(gliph_score) {
    
}