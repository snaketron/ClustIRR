# KK: commented out for test coverage checks (currently not used)

# Description:
# Create an igraph object from the gliph_output and export strongly connected
# components and/or modules. Return a map between cell IDs and component and/or
# module IDs.
# get_graph <- function(gliph_output,
#                       edge_type = "local+global",
#                       chain_type = "CDR3a+CDR3b") {
#
#     # Description
#     # extracts required edges from gliph output object
#     get_graph_edges <- function(gliph_output,
#                                 edge_type,
#                                 chain_type) {
#
#         chains <- unlist(strsplit(x = chain_type, split = "\\+"))
#         edges <- unlist(strsplit(x = edge_type, split = "\\+"))
#
#         # gliph_output$edges = list over chains, rbind chain-spec.data.frames
#         edge_data <- do.call(rbind, gliph_output$edges)
#
#         # make selection
#         j <- which(edge_data$chain %in% chains
#                    & edge_data$edge_type %in% edges)
#         if (length(j) == 0) {
#             return(NULL)
#         }
#         return(edge_data[j,])
#     }
#
#     # Description:
#     # Get the strongly connected components in a graph
#     get_components <- function(graph) {
#
#         c <- igraph::components(graph = graph,
#                                 mode = "strong")
#
#         # cell ID to cluster ID map
#         c <- data.frame(ID = names(c$membership),
#                         component_ID = as.numeric(c$membership))
#
#         return(c)
#     }
#
#
#
#
#     # Description:
#     # Use graph clustering to find densely connected subgraphs
#     get_modules <- function(graph) {
#
#     }
#
#
#     # edge_type = "local+global" == "global+local"
#     # edge_type = "local"
#     # edge_type = "global"
#
#     # chain_type = "CDR3a+CDR3b"
#     # chain_type = "CDR3a"
#     # chain_type = "CDR3b"
#
#     # check if chains are part of gliph output
#     # check if edge_type is part of gliph ouput -> warn
#     # check if elements of edge_type and chain_type are allowed values
#
#     edges <- get_graph_edges(
#         gliph_output = gliph_output,
#         edge_type = edge_type,
#         chain_type = chain_type)
#
#     # no edges -> no graph
#     if (is.null(edges) == TRUE) {
#         return(list(graph = NA,
#                     components = NA,
#                     modules = NA))
#     }
#
#
#     # if edges are available: create igraph from data.frame
#     n <- data.frame(name = gliph_output$data_sample$ID)
#
#     # remove self-edges (loops)
#     edges <- edges[edges$from != edges$to, ]
#
#     # create graph
#     g <- igraph::graph_from_data_frame(d = edges,
#                                        directed = FALSE,
#                                        vertices = n)
#
#     # extract components
#     c <- get_components(graph = g)
#
#     # extract modules
#     # m <- get_modules(graph = g)
#     m <- NA
#
#     return(list(graph = g,
#                 components = c,
#                 modules = m))
# }
#
#
# # Description:
# # This function will perform the scoring of the components (graph
# # clusters/modules) based on e.g. HLA, V, J
#
# score_by_cellvar <- function(graph_data,
#                              data_sample,
#                              labels) {
#
#
#
#
# }
#
#
# # Description:
# # This function will perform the scoring of the clusters based on
# # e.g. HLA, V, J
# # TODO clone size
# score_by_subjectvar <- function(graph_data,
#                                 data_sample,
#                                 data_meta) {
#
#     # enrichment of cluster sizes (definition of cluster: local vs. global)
#
#
#     # enrichment of categorical features (such as HLAs, samples, biological
#     # condition, TRV, TRJ, cdr3 length, SHM, etc.)
#
#     # probabilistic evaluation vs.
#
# }


# #Jan's code
#
# cluster_scoring <- function(cluster_list,
#                             cdr3_sequences,
#                             refdb_beta = "gliph_reference",
#                             v_usage_freq = NULL,
#                             cdr3_length_freq = NULL,
#                             ref_cluster_size = "original",
#                             gliph_version = 1,
#                             sim_depth = 1000,
#                             hla_cutoff = 0.1,
#                             n_cores = 1){
#
#
#     doParallel::registerDoParallel(n_cores)
#
#     actCluster <- NULL
#     res <- foreach::foreach(
#         actCluster = base::seq_along(cluster_list)) %dopar% {
#
#         ### Get sequences and information of current cluster
#         act_seq_infos <- cluster_list[[actCluster]]
#         num_members <- base::nrow(act_seq_infos) # number of ALL members
#         ori_num_members <- base::length(base::unique(
#             act_seq_infos$CDR3b)) # number of all unique CDR3b sequences
#         all_scores <- base::c()
#
#         ### Get network size score from lookup table
#         score_network_size <- 1
#         nearest_sample_size <- base::order(base::abs(1-base::as.numeric(
#             base::colnames(ref_cluster_sizes)[-1])/
#                 base::nrow(cdr3_sequences)))[1]
#         if(ori_num_members > 100){
#             score_network_size <-
#                 ref_cluster_sizes[100,nearest_sample_size+1]
#         } else {
#             score_network_size <- ref_cluster_sizes[
#                 ori_num_members,nearest_sample_size+1]
#         }
#         all_scores <- base::c(all_scores, score_network_size)
#
#         ### Enrichment of CDR3 length (spectratype) within cluster
#         score_cdr3_length <- base::c()
#         # calculate score of sample (product of all frequencies)
#         pick_freqs <- base::data.frame(base::table(base::nchar(base::unique(
#             act_seq_infos$CDR3b))))
#         base::colnames(pick_freqs) <- base::c("object", "probs")
#         pick_freqs$probs <- pick_freqs$probs/ori_num_members
#         sample_score <- base::round(base::prod(pick_freqs$probs), digits = 3)
#
#         # generate random subsamples
#         random_subsample <- base::list()
#         for(i in base::seq_len(sim_depth)){
#             random_subsample[[i]] <- base::sample.int(
#                 n = base::length(cdr3_length_ref_frequencies),
#                 size = ori_num_members,
#                 prob = cdr3_length_ref_frequencies, replace = TRUE)
#         }
#
#         # calculate score of subsamples (product of all frequencies)
#         pick_freqs <- stringdist::seq_qgrams(
#             .list = random_subsample)[,-1]/ori_num_members
#         pick_freqs[pick_freqs == 0] <- 1
#         pick_freqs <- base::round(base::exp(base::colSums(base::log(
#             pick_freqs))), digits = 3)
#         # vectorized way to calculate the product of each column
#         if(gliph_version == 1){
#             score_cdr3_length <- base::sum(pick_freqs
#                                            >= sample_score)/sim_depth
#         } else {
#             score_cdr3_length <- base::sum(pick_freqs
#                                            > sample_score)/sim_depth
#         }
#         if(score_cdr3_length == 0) score_cdr3_length <- 1/sim_depth
#         # minimum score of 1/sim_depth
#
#         all_scores <- base::c(all_scores, score_cdr3_length)
#
#         ### Enrichment of v genes within cluster
#         score_vgene <- base::c()
#         if(vgene_info == TRUE){
#
#             # calculate score of sample (product of all frequencies)
#             pick_freqs <- base::data.frame(base::table(act_seq_infos$TRBV))
#             base::colnames(pick_freqs) <- base::c("object", "probs")
#             pick_freqs$probs <- pick_freqs$probs/num_members
#             sample_score <- base::round(
#                 base::prod(pick_freqs$probs), digits = 3)
#
#             # generate random subsamples
#             random_subsample <- base::list()
#             for(i in base::seq_len(sim_depth)){
#                 random_subsample[[i]] <- base::sample.int(
#                     n = base::length(vgene_ref_frequencies),
#                     size = num_members,
#                     prob = vgene_ref_frequencies, replace = TRUE)
#             }
#
#             # calculate score of subsamples (product of all frequencies)
#             pick_freqs <- stringdist::seq_qgrams(
#                 .list = random_subsample)[,-1]/num_members
#             pick_freqs[pick_freqs == 0] <- 1
#             pick_freqs <- base::round(base::exp(base::colSums(base::log(
#                 pick_freqs))), digits = 3)
#             # vectorized way to calculate the product of each column
#             if(gliph_version == 1){
#                 score_vgene <- base::sum(pick_freqs >= sample_score)/
#                     sim_depth
#             } else {
#                 score_vgene <- base::sum(pick_freqs > sample_score)/
#                     sim_depth
#             }
#
#             if(score_vgene == 0) score_vgene <- 1/sim_depth
#             # minimum score of 1/sim_depth
#             all_scores <- base::c(all_scores, score_vgene)
#         }
#
#         ### Enrichment of clonal expansion within cluster
#         score_clonal_expansion <- base::c()
#         if(counts_info == TRUE){
#             sample_score <- base::sum(base::as.numeric(
#                 act_seq_infos$counts))/
#                 num_members
#             counter <- 0
#             for(i in base::seq_len(sim_depth)){
#                 random_subsample <- base::sample(
#                     x = cdr3_sequences$counts,
#                     size = num_members, replace = FALSE)
#                 test_score <- base::sum(base::as.numeric(
#                     random_subsample))/num_members
#                 if(test_score>=sample_score) counter <- counter+1
#             }
#             if(counter == 0)
#                 score_clonal_expansion <- 1/sim_depth
#             else
#                 score_clonal_expansion <- counter/sim_depth
#             score_clonal_expansion <- base::round(
#                 score_clonal_expansion, digits = 3)
#
#             all_scores <- base::c(all_scores, score_clonal_expansion)
#         }
#
#         ### Enrichment of common HLA among donor TCR contributors in cluster
#         score_hla <- base::c()
#         lowest_hla <- ""
#         if(hla_info == TRUE && patient_info == TRUE){
#             act_seq_infos <- act_seq_infos[act_seq_infos$HLA != "" &
#                                           !base::is.na(act_seq_infos$HLA),]
#
#             if(base::nrow(act_seq_infos) > 0){
#                 act_seq_infos$patient <-
#                          gsub(":.*", "",act_seq_infos$patient)
#
#                 score_hla <- 1
#                 for(act_hla in base::seq_len(num_HLAs)){
#                     crg_patient_count <- base::length(base::unique(
#                         act_seq_infos$patient))
#                     crg_patient_hla_count <- base::sum(base::unlist(
#                         base::lapply(
#                         all_patient_hlas[base::unique(
#                                          act_seq_infos$patient)],
#                         function(x){
#                             if(all_hlas$HLA[act_hla] %in% x)
#                                 1
#                             else
#                                 0
#                         })))
#                     if(crg_patient_hla_count > 1){
#                         act_Prob <- base::sum(base::choose(
#                             all_hlas$counts[act_hla],
#                             crg_patient_hla_count:crg_patient_count)*
#                                 base::choose(
#                                     num_patients-all_hlas$counts[act_hla],
#                                     crg_patient_count-
#                                     crg_patient_hla_count:crg_patient_count)/
#                                base::choose(num_patients,
#                                             crg_patient_count))
#                         if(act_Prob<score_hla) score_hla <- act_Prob
#                         if(act_Prob < hla_cutoff){
#                             if(lowest_hla == ""){
#                                 lowest_hla <-
#                                     base::paste(all_hlas$HLA[act_hla],
#                                                 " [(", crg_patient_hla_count,
#                                                 "/",
#                                                 crg_patient_count,
#                                                 ") vs (",
#                                                 all_hlas$counts[act_hla],
#                                                 "/",
#                                                 num_patients,
#                                                 ") = ",
#                                                 base::formatC(act_Prob,
#                                                               digits = 1,
#                                                               format = "e"),
#                                                 "]",
#                                                 sep = "" )
#                             } else {
#                                 lowest_hla <-
#                                     base::paste(lowest_hla, ", ",
#                                                 all_hlas$HLA[act_hla],
#                                                 " [(", crg_patient_hla_count,
#                                                 "/",
#                                                 crg_patient_count, ") vs (",
#                                                 all_hlas$counts[act_hla],
#                                                 "/",
#                                                 num_patients,
#                                                 ") = ",
#                                                 base::formatC(act_Prob,
#                                                               digits = 1,
#                                                               format = "e"),
#                                                 "]",
#                                                 sep = "" )
#                             }
#                         }
#                     }
#
#                 }
#             } else {
#                 score_hla <- 1
#             }
#
#             all_scores <- base::c(all_scores, score_hla)
#         }
#
#         ### Total score
#         if(gliph_version == 1) score_final <- base::prod(all_scores)*0.001*64
#         else if(gliph_version == 2) score_final <- base::prod(all_scores)
#
#         ### Output
#         all_scores <- c(score_final, all_scores)
#         all_scores <- base::formatC(all_scores, digits = 1, format = "e")
#         output <- base::c(base::names(cluster_list)[actCluster], all_scores)
#         if(hla_info == TRUE && patient_info == TRUE){
#             output <- base::c(output, lowest_hla)
#         }
#         output
#     }
#
#     doParallel::stopImplicitCluster()
#
#     res <- base::data.frame(base::matrix(base::unlist(res),
#                                          ncol = 2+base::length(
#         score_names), byrow = TRUE))
#     base::colnames(res) <- base::c("leader.tag", "total.score", score_names)
#
#     for(i in base::c("total.score", score_names))
#         if(i != "lowest.hlas") res[,i] <- base::as.numeric(res[,i])
#
#     # set all theoretical numeric values (act. characters) to numeric values
#     if(base::is.data.frame(res)){
#         for(i in base::seq_len(base::ncol(res))){
#             if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(
#                 res[,i])))) == FALSE) res[,i] <- base::as.numeric(res[,i])
#         }
#     }
#
#     # Closing time!
#     base::return(res[,-1])
# }
