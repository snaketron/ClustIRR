
# Description:
# Now exported to user, but used within function gliph.R
parameter_check <- function(data_sample,
                            data_ref,
                            version,
                            ks,
                            cores,
                            control) {

    # add individual tests from below here

}


check_data_sample <- function(data_sample) {

}

check_data_ref <- function(data_ref) {

}

check_version <- function(version) {

}

check_ks_local_min_ove <- function(ks, local_min_ove) {

}

check_cores <- function(cores) {

}

check_B <- function(B) {

}

check_global_max_dist <- function(global_max_dist) {

}

check_local_min_p <- function(local_min_p) {

}

check_local_min_o <- function(local_min_o) {

}

check_trim_flanks <- function(trim_flanks) {

}

check_flank_size <- function(flank_size) {

}

check_trim_flanks_flank_size <- function(trim_flanks,
                                         flank_size) {

}

check_global_pairs <- function(global_pairs, data_sample) {
    # global_pairs must be integer matrix with two columns and u rows.
    # In each entry global_pairs will store an index (integer) i = 1, ..., n,
    # pointing to a specific CDR3 sequence from data_sample which is a
    # data.frame with nrow = n.
    # row in global pairs gives us the indices of two CDR3s that are globally
    # similar (e.g. Hamming dist < 1 or distance computed using external tool)
}




# functions from Kai/Jan's script -> we should use them to populate code above
#
# ch_result_folder <- function(result_folder, sim_depth, lcminp, lcminove){
#     if(!base::is.character(result_folder))
#         base::stop("result_folder has to be a character object")
#     if(base::length(result_folder) > 1)
#         base::stop("result_folder has to be a single path")
#     save_results <- FALSE
#     warn <- "WARNING: Saving the results is cancelled. The following file
#   already exists in the result folder:\n"
#     if(result_folder != ""){
#         if(base::substr(result_folder,base::nchar(result_folder),
#                         base::nchar(result_folder)) != "/")
#             result_folder <- base::paste0(result_folder,"/")
#         if(!base::dir.exists(result_folder)) base::dir.create(result_folder)
#         save_results <- TRUE
#         if(base::file.exists(base::paste0(result_folder,
#                                           "kmer_resample_",sim_depth,"_log.txt"))){
#             save_results <- FALSE
#             base::cat(base::paste(warn, base::paste0(result_folder,"kmer_resample_",
#                                                      sim_depth,"_log.txt"),"\n"))}
#         if(base::file.exists(base::paste0(result_folder,"kmer_resample_",
#                                           sim_depth,"_minp",lcminp,"_ove",
#                                           base::paste(lcminove,
#                                                       collapse = "_"),".txt"))){
#             save_results <- FALSE
#             base::cat(base::paste(warn, base::paste0(result_folder,"kmer_resample_",
#                                                      sim_depth,"_minp",lcminp,"_ove",
#                                                      base::paste(lcminove, collapse = "_"),
#                                                      ".txt"),"\n"))}
#         if(base::file.exists(base::paste0(result_folder,"kmer_resample_",
#                                           sim_depth,"_all_motifs.txt"))){
#             save_results <- FALSE
#             base::cat(base::paste(warn, base::paste0(result_folder,"kmer_resample_",
#                                                      sim_depth,"_all_motifs.txt"),"\n"))}
#         if(base::file.exists(base::paste0(result_folder, "clone_network.txt"))){
#             save_results <- FALSE
#             base::cat(base::paste(warn, base::paste0(result_folder,
#                                                      "clone_network.txt"),"\n"))}
#         if(base::file.exists(base::paste0(result_folder,"convergence_groups.txt"))){
#             save_results <- FALSE
#             base::cat(base::paste(warn, base::paste0(result_folder,
#                                                      "convergence_groups.txt"),"\n"))}
#         if(base::file.exists(base::paste0(result_folder,"parameter.txt"))){
#             save_results <- FALSE
#             base::cat(base::paste(warn, base::paste0(result_folder,"parameter.txt"),
#                                   "\n"))}
#     }
#     return(save_results)
# }
#
#
# ch_refdb_beta <- function(refdb_beta) {
#     # if(!(refdb_beta %in% base::c("gliph_reference", "human_v1.0_CD4",
#     # "human_v1.0_CD8", "human_v1.0_CD48", "human_v2.0_CD4",
#     #                        "human_v2.0_CD8", "human_v2.0_CD48",
#     #"mouse_v1.0_CD4", "mouse_v1.0_CD8", "mouse_v1.0_CD48")) &&
#     #    !base::is.data.frame(refdb_beta)){
#     if(!base::is.data.frame(refdb_beta)){
#         if(base::length(refdb_beta) != 1 || !is.character(refdb_beta)){
#             base::stop("refdb_beta has to be a data frame (containing CDR3b sequences
#       in the first column and optional V-gene information in the
#                  second column) or the value 'gliph_reference'")
#         } else if(!(refdb_beta %in% base::c("gliph_reference"))){
#             base::stop("refdb_beta has to be a data frame (containing CDR3b sequences
#                  in the first column and optional V-gene information in the
#                  second column) or the value 'gliph_reference'")
#         }
#     }
# }
#
# ch_v_usage_freq <- function(v_usage_freq){
#     if(!base::is.null(v_usage_freq)){
#         if(base::is.data.frame(v_usage_freq)){
#             if(base::ncol(v_usage_freq) < 2)
#                 base::stop("v_usage_freq has to be a data frame containing
#                    V-gene information in the first column and the corresponding
#                    frequency in a naive  T-cell repertoire in the second
#                    column.")
#             if(base::nrow(v_usage_freq) < 1) base::stop("v_usage_freq has to
#                                                     contain at least one row.")
#             if(base::suppressWarnings(
#                 base::any(base::is.na(base::as.numeric(v_usage_freq[,2])))) == TRUE){
#                 base::stop("The second column of v_usage_freq must contain the frequency
#                    of the corresponding V-gene in the first column in a naive
#                    T-cell repertoire.")
#             } else v_usage_freq[,2] <- as.numeric(v_usage_freq[,2])
#
#         } else {base::stop("v_usage_freq has to be a data frame containing V-gene
#                        information in the first column and the corresponding
#                        frequency in a naive T-cell repertoire in the second
#                        column.")}
#     }
# }
#
# ch_cdr3_length_freq <- function(cdr3_length_freq){
#     if(!base::is.null(cdr3_length_freq)){
#         if(base::is.data.frame(cdr3_length_freq)){
#             if(base::ncol(cdr3_length_freq) < 2)
#                 base::stop("cdr3_length_freq has to be a data frame containing
#                    CDR3 lengths in the first column and the corresponding
#                    frequency in a naive  T-cell repertoire in the
#                    second column.")
#             if(base::nrow(cdr3_length_freq) < 1)
#                 base::stop("cdr3_length_freq has to contain at least one row.")
#             if(base::suppressWarnings(base::any(
#                 base::is.na(base::as.numeric(cdr3_length_freq[,2])))) == TRUE){
#                 base::stop("The second column of cdr3_length_freq must contain the
#                    frequency of the corresponding CDR3 length in the first
#                    column in a naive T-cell repertoire.")
#             } else cdr3_length_freq[,2] <- as.numeric(cdr3_length_freq[,2])
#
#         } else {base::stop("cdr3_length_freq has to be a data frame containing CDR3
#                        lengths in the first column and the corresponding
#                        frequency in a naive T-cell repertoire in the
#                        second column.")}
#     }
# }
#
# ch_ref_cluster_size <- function(ref_cluster_size){
#     if(!(ref_cluster_size %in% base::c("original", "simulated") ||
#          !base::is.character(ref_cluster_size) ||
#          base::length(ref_cluster_size) > 1)){
#         base::stop("ref_cluster_size has to be either 'original' or 'simulated'.")
#     }
# }
#
# ch_sim_depth <- function(sim_depth){
#     if(!base::is.numeric(sim_depth)) base::stop("sim_depth has to be numeric")
#     if(base::length(sim_depth) > 1)
#         base::stop("sim_depth has to be a single number")
#     if(sim_depth < 1) base::stop("sim_depth must be at least 1")
#     sim_depth <- base::round(sim_depth)
#     return(sim_depth)
# }
#
# ch_lcminp <- function(lcminp){
#     if(!base::is.numeric(lcminp)) base::stop("lcminp has to be numeric")
#     if(base::length(lcminp) > 1) base::stop("lcminp has to be a single number")
#     if(lcminp <= 0) base::stop("lcminp must be greater than 0")
# }
#
# ch_lcminove <- function(lcminove, motif_length){
#     if(!base::is.numeric(lcminove)) base::stop("lcminove has to be numeric")
#     if(base::length(lcminove) > 1 &&
#        base::length(lcminove) != base::length(motif_length))
#         base::stop("lcminove has to be a single number or of same length as
#                motif_length")
#     if(base::any(lcminove < 1)) base::stop("lcminove must be at least 1")
# }
#
#
# ch_kmer_mindepth <- function(kmer_mindepth){
#     if(!base::is.numeric(kmer_mindepth))
#         base::stop("kmer_mindepth has to be numeric")
#     if(base::length(kmer_mindepth) > 1)
#         base::stop("kmer_mindepth has to be a single number")
#     if(kmer_mindepth < 1) base::stop("kmer_mindepth must be at least 1")
#     kmer_mindepth <- base::round(kmer_mindepth)
# }
#
#
# ch_accept_sequences_with_C_F_start_end <-
#     function(accept_sequences_with_C_F_start_end){
#         if(!base::is.logical(accept_sequences_with_C_F_start_end))
#             base::stop("accept_sequences_with_C_F_start_end has to be logical")
#     }
#
# ch_min_seq_length_a <- function(min_seq_length){
#     if(!base::is.numeric(min_seq_length))
#         base::stop("min_seq_length has to be numeric")
#     if(base::length(min_seq_length) > 1)
#         base::stop("min_seq_length has to be a single number")
#     if(min_seq_length < 0) base::stop("min_seq_length must be at least 0")
#     min_seq_length <- base::round(min_seq_length)
#     return(min_seq_length)
# }
#
# ch_gccutoff <- function(gccutoff){
#     if(!base::is.null(gccutoff) && !base::is.numeric(gccutoff))
#         base::stop("gccutoff has to be NULL or numeric")
#     if(!base::is.null(gccutoff) && base::length(gccutoff)>1)
#         base::stop("gccutoff has to be NULL or a single number")
#     if(!base::is.null(gccutoff) && gccutoff < 0)
#         base::stop("gccutoff must be at least 0")
# }
#
# ch_structboundaries <- function(structboundaries){
#     if(!base::is.logical(structboundaries))
#         base::stop("structboundaries has to be logical")
# }
#
# ch_boundary_size <- function(boundary_size, structboundaries, min_seq_length){
#     if(!base::is.numeric(boundary_size))
#         base::stop("boundary_size has to be numeric")
#     if(base::length(boundary_size) > 1)
#         base::stop("boundary_size has to be a single number")
#     if(boundary_size < 0) base::stop("boundary_size must be at least 0")
#     boundary_size <- base::round(boundary_size)
#     if(structboundaries == TRUE)
#         min_seq_length <- base::max(min_seq_length, ((boundary_size * 2)+1))
#     return(list(boundary_size, min_seq_length))
# }
#
# ch_motif_length <- function(motif_length){
#     if(!base::is.numeric(motif_length))
#         base::stop("motif_length has to be numeric")
#     if(base::any(motif_length < 1))
#         base::stop("values of motif_length must be at least 1")
#     motif_length <- base::round(motif_length)
# }
#
# ch_discontinuous <- function(discontinuous){
#     if(!base::is.logical(discontinuous))
#         base::stop("discontinuous has to be logical")
# }
#
# ch_make_depth_fig <- function(make_depth_fig){
#     if(!base::is.logical(make_depth_fig))
#         base::stop("make_depth_fig has to be logical")
# }
#
# ch_local_similarities <- function(local_similarities){
#     if(!base::is.logical(local_similarities))
#         base::stop("local_similarities has to be logical")
# }
#
# ch_global_similarities <- function(global_similarities, local_similarities){
#     if(!base::is.logical(global_similarities))
#         base::stop("global_similarities has to be logical")
#     if(local_similarities == FALSE && global_similarities == FALSE)
#         base::stop("Either local_similarities or global_similarities
#                have to be TRUE")
# }
#
# ch_global_vgene <- function(global_vgene){
#     if(!base::is.logical(global_vgene))
#         base::stop("global_vgene has to be logical")
# }
#
# ch_positional_motifs <- function(positional_motifs){
#     if(!base::is.logical(positional_motifs))
#         base::stop("positional_motifs has to be logical")
# }
#
# ch_cdr3_len_stratify <- function(cdr3_len_stratify){
#     if(!base::is.logical(cdr3_len_stratify))
#         base::stop("cdr3_len_stratify has to be logical")
# }
#
# ch_vgene_stratify <- function(vgene_stratify){
#     if(!base::is.logical(vgene_stratify))
#         base::stop("vgene_stratify has to be logical")
# }
#
# ch_public_tcrs <- function(public_tcrs){
#     if(!base::is.logical(public_tcrs))
#         base::stop("public_tcrs has to be logical")
# }
#
# ch_cluster_min_size <- function(cluster_min_size){
#     if(!base::is.numeric(cluster_min_size))
#         base::stop("cluster_min_size has to be numeric")
#     if(base::length(cluster_min_size) > 1)
#         base::stop("cluster_min_size has to be a single number")
#     if(cluster_min_size < 1)
#         base::stop("cluster_min_size must be at least 1")
#     cluster_min_size <- base::round(cluster_min_size)
#     return(cluster_min_size)
# }
#
# ch_hla_cutoff <- function(hla_cutoff){
#     if(!base::is.numeric(hla_cutoff))
#         base::stop("hla_cutoff has to be numeric")
#     if(base::length(hla_cutoff) > 1)
#         base::stop("hla_cutoff has to be a single number")
#     if(hla_cutoff > 1 || hla_cutoff < 0)
#         base::stop("hla_cutoff must be between 0 and 1")
# }
#
# ch_n_cores <- function(n_cores){
#     if(!base::is.null(n_cores))
#     {
#         if(!base::is.numeric(n_cores)) base::stop("n_cores has to be numeric")
#         if(base::length(n_cores) > 1)
#             base::stop("n_cores has to be a single number")
#         if(n_cores < 1) base::stop("n_cores must be at least 1")
#         n_cores <- base::round(n_cores)
#     }
#     return(n_cores)
# }
#
# ch_min_seq_length_b <- function(structboundaries, min_seq_length,
#                                 boundary_size){
#     if(structboundaries == TRUE)
#         min_seq_length <- base::max(min_seq_length, 2*boundary_size)
#     return(min_seq_length)
# }

