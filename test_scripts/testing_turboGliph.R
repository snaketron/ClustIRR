# This script contains first benchmarking results. Not reproducible, yet.
# More tests needed with seed.


#### Inputs ####
require(turboGliph)

# initial input data
# data("hs_CD8_ref")
# ref <- base::grep(pattern = "^C.*F$",x = hs_CD8_ref$CDR3b,
#                   perl = TRUE,value = TRUE)
# ref <- data.frame(CDR3b = ref, TRBV = NA, TRBJ = NA)
#
# data_sample <- ref[sample(x = 1:nrow(ref), size = 2000, replace = F), 1:3]
# data_ref <- ref[, 1:3]

## idea: use the same input as jan did
## this way maybe comparison with gliph2 version

# input sequences and related information
ex_path = "test_scripts/jans_tests/Input_demo_gliph2.csv"
# reference database
re_path = "test_scripts/jans_tests/ref_CD4_v1.0.txt"

data_sample <- read.csv2(file = ex_path, sep = "\t")
data_ref <- read.csv2(file = re_path, sep = "\t", header = FALSE)
colnames(data_ref) <- c("CDR3b", "TRBV", "TRBJ")

# prepare for gliphR input
ref <- base::grep(pattern = "^C.*F$",x = data_ref$CDR3b,
                  perl = TRUE,value = TRUE)
ref <- data.frame(CDR3b = ref, TRBV = NA, TRBJ = NA)
data_ref <- ref[, 1:3]
data_sample <- data_sample[, 1:3]

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
source("R/score.R")



#### Run ####
out_v1 <- gliph(data_sample = data_sample,
                data_ref = data_ref,
                ks = ks,
                cores = cores,
                version = 1,
                control = control_input)

# only working with initial input
# out_v2 <- gliph(data_sample = data_sample,
#                 data_ref = data_ref,
#                 ks = ks,
#                 cores = cores,
#                 version = 2,
#                 control = control_input)


# out_v3 <- gliph(data_sample = data_sample,
#                 data_ref = data_ref,
#                 ks = ks,
#                 cores = cores,
#                 version = 3,
#                 control = control_input)


###
gliph_output <- out_v1
ge <- get_graph(gliph_output = gliph_output,
                chain = "CDR3a+CDR3b",
                edge_type = "local+global")

ge$graph
ge$components


# plot(ge)
# E(ge)$edge_type
###




# this results in an error at scoring stage when using initial input data:
# "Fehler in { : task 1 failed - "ungültiges erstes Argument""
jan_v1 <- turboGliph::turbo_gliph(
    cdr3_sequences = data_sample,
    #result_folder = "test_scripts/test_results/",
    refdb_beta = data_ref,
    lcminp = 0.05,
    gccutoff = 1,
    #structboundaries = F, # this breaks if left in with new example data
    boundary_size = 0,
    cluster_min_size = 1,
    accept_sequences_with_C_F_start_end = T)

## this also does not work with initial input data, this time stopping at part 2:
# Fehler in base::seq_len(base::min(base::max(reference_seqs$nchar),
# base::max(sample_seqs$nchar))) :
#    Argument muss sich in eine nicht-negative ganze Zahl umwandeln lassen
# Zusätzlich: Warnmeldung:
#    In base::max(sample_seqs$nchar) :
#    kein nicht-fehlendes Argument für max; gebe -Inf zurück

jan_v2 <- turboGliph::gliph2(
    cdr3_sequences = data_sample,
    #result_folder = "test_scripts/test_results/",
    refdb_beta = data_ref,
    lcminp = 0.05,
    v_usage_freq = NULL,
    cdr3_length_freq = NULL,
    #ref_cluster_size = NULL,
    #sim_depth = 1000,
    #motif_distance_cutoff = 0,
    accept_sequences_with_C_F_start_end = T,
    #structboundaries = F,  # this breaks if left in with new example data
    min_seq_length = 1)


#e_v2
jan_v1$motif_enrichment$selected_motifs
jan_v2$motif_enrichment$selected_motifs

#### Compare results #####
e_v1 <- out_v1$clust$CDR3b$motif_enrichment
# e_v2 <- out_v2$clust$CDR3b$motif_enrichment
#. _v3 <- out_v3$clust$CDR3b$motif_enrichment

# e_jan_v1 <- jan_v1$motif_enrichment$all_motifs
# e_jan_v2 <- jan_v2$motif_enrichment$all_motifs

# x <- merge(x = e_v2, y = e_jan_v2, by = "motif")
# x$d <- x$num_in_sample-x$f_sample

e_v1 <- e_v1[e_v1$filter == TRUE,]
# e_v2 <- e_v2[e_v2$filter == TRUE,]
# e_v3 <- e_v3[e_v3$filter == TRUE,]

e_jan_v1 <- jan_v1$motif_enrichment$selected_motifs
e_jan_v2 <- jan_v2$motif_enrichment$selected_motifs

