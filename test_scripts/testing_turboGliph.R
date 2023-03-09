# This script contains first benchmarking results. Not reproducible, yet. 
# More tests needed with seed.


#### Inputs ####
require(turboGliph)
data("hs_CD8_ref")

data_sample <- hs_CD8_ref[sample(x = 1:nrow(hs_CD8_ref), size = 2000, replace = F), 1:3]
data_ref <- hs_CD8_ref[, 1:3]
ks <- c(2, 3, 4)
B <- 1000
cores <- 1

control_input <- list(
    B = 1000,
    global_max_dist = 1,
    local_min_p = 0.05, 
    local_min_ove = c(10^3, 10^2, 10^1), 
    local_min_o = 3,
    trim_flanks = FALSE,
    flank_size = 3)

source("R/util_v1.R")
source("R/util_v2.R")
source("R/util_v1_v2.R")
source("R/gliph.R")


#### Run ####
out_v1 <- gliph(data_sample = data_sample,
                data_ref = data_ref,
                ks = ks,
                cores = cores,
                version = 1,
                control = control_input)

out_v2 <- gliph(data_sample = data_sample,
                data_ref = data_ref,
                ks = ks,
                cores = cores,
                version = 2,
                control = control_input)

out_v3 <- gliph(data_sample = data_sample,
                data_ref = data_ref,
                ks = ks,
                cores = cores,
                version = 3,
                control = control_input)


jan_v1 <- turboGliph::turbo_gliph(
    cdr3_sequences = data_sample,
    result_folder = "/home/sktron/Desktop/tmp/",
    refdb_beta = data_ref,
    lcminp = 0.05,
    gccutoff = 1,
    structboundaries = F,
    boundary_size = 0,
    cluster_min_size = 1,
    accept_sequences_with_C_F_start_end = FALSE)


jan_v2 <- turboGliph::gliph2(
    cdr3_sequences = data_sample,
    result_folder = "/home/sktron/Desktop/tmp/",
    refdb_beta = data_ref,
    lcminp = 0.05,
    accept_sequences_with_C_F_start_end = FALSE)



#### Compare results #####
e_v1 <- out_v1$clust$CDR3b$motif_enrichment
e_v2 <- out_v2$clust$CDR3b$motif_enrichment
e_v3 <- out_v3$clust$CDR3b$motif_enrichment

e_jan_v1 <- jan_v1$motif_enrichment$all_motifs
e_jan_v2 <- jan_v2$motif_enrichment$all_motifs


e_v1 <- e_v1[e_v1$filter == TRUE,]
e_v2 <- e_v2[e_v2$filter == TRUE,]
e_v3 <- e_v3[e_v3$filter == TRUE,]

e_jan_v1 <- jan_v1$motif_enrichment$selected_motifs
e_jan_v2 <- jan_v2$motif_enrichment$selected_motifs

