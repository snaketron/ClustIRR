### Recreation of Jan's comparison of gliph2
### turboGliph vs server version
### also tested with turboGliph package version of gliph2
### nearly same differences, count difference exactly the same
### tested with local executable of gliph2:
### - global: no difference in count
### - local:  46 / 49 instead of 31 / 49, also 46 included in package
###           but p-values seem to differ even more


## data import & gliph2 online results (precomputed)

# input sequences and related information
ex_path = "test_scripts/jans_tests/Input_demo_gliph2.csv"
# # results of the GLIPH2 algorithm for the input
# se_path =
#     "test_scripts/jans_tests/Input_demo_gliph2_original_webserve_results.csv"
# results of the GLIPH2 algorithm, computed with local exec
se_path =
    "test_scripts/gliph2_exec_comparison/Comparison_cluster.csv"
gl_perl_path =
    "test_scripts/perl_gliph_output/Input_demo_gliph2.csv-kmer_resample_1000_minp0.001_ove10.txt"

# reference database
re_path = "test_scripts/jans_tests/ref_CD4_v1.0.txt"

example_input <- read.csv2(file = ex_path, sep = "\t")
server_output <- read.csv2(file = se_path, sep = ",")
gliph_perl_output <- read.csv2(file = gl_perl_path, sep = "\t")
example_reference <- read.csv2(file = re_path, sep = "\t", header = FALSE)
colnames(example_reference) <- c("CDR3b", "TRBV", "TRBJ")


# reduce to cluster information
server_output_cluster <- unique(server_output[, 1:12])
for(i in 3:11) {
    # convert columns to numeric values
    server_output_cluster[,i] <- as.numeric(server_output_cluster[,i])
}




## turboGliph
library(turboGliph)
library(doParallel)
source("test_scripts/jans_tests/oldGliph2.R")
# jans modified version of gliph2
package_output <- oldGliph2(cdr3_sequences = example_input,
                            refdb_beta = example_reference,
                            accept_sequences_with_C_F_start_end = FALSE,
                            lcminove = 10,
                            min_seq_length = 0)
# package version of turboGliph
package_output_tg <- turboGliph::turbo_gliph(cdr3_sequences = example_input,
                            refdb_beta = example_reference,
                            accept_sequences_with_C_F_start_end = FALSE,
                            lcminove = 10)

tg_output = package_output_tg[["motif_enrichment"]][["selected_motifs"]]

# gliph2 package output
package_output_cluster <-
    package_output$cluster_properties
package_output_cluster <-
    package_output_cluster[, c(1,2,7,3:5,8,11,12,10,9,6)]
package_output_cluster <-
    package_output_cluster[order(package_output_cluster$fisher.score),]

cat("gliph2 local executable output")
head(server_output_cluster)

cat("turboGliph gliph2 output")
head(package_output_cluster, n=10)

cat("gliph perl version output")
gliph_perl_output[order(-gliph_perl_output$Counts),]

cat("turboGliph gliph version output")
head(tg_output[order(-tg_output$counts),], n=10)

## comparison - global similarity
global_cluster_server <-
    server_output_cluster[grep("%", server_output_cluster$pattern),]
global_cluster_package <-
    package_output_cluster[package_output_cluster$type == "global",]

cat("Number of global clusters\n")
cat("Server: ", nrow(global_cluster_server), "\n")   # server = local: 97
cat("Package: ", nrow(global_cluster_package), "\n") # packge: 115


# Simplify tag for comparison of shared cluster count
global_cluster_package$tag.simple <-
    gsub("_.*", "", global_cluster_package$tag)

# Get shared cluster results
global_cluster_server_shared <-
    global_cluster_server[global_cluster_server$pattern %in%
                              global_cluster_package$tag.simple,]
global_cluster_package_shared <-
    global_cluster_package[global_cluster_package$tag.simple %in%
                               global_cluster_server$pattern,]

cat("Number of shared global clusters\n")
cat("Server: ", nrow(global_cluster_server_shared), "\n")   # server=local: 97
cat("Package: ", nrow(global_cluster_package_shared), "\n") # package: 97


# Unite fisher score information in one date frame
global_fisher_dataframe <-
    data.frame(pattern = global_cluster_server_shared$pattern,
               fisherServer = global_cluster_server_shared$Fisher_score)
global_fisher_dataframe$fisherPackage <-
    sapply(1:nrow(global_fisher_dataframe),
           function(x)
               return(global_cluster_package_shared$fisher.score[
                   global_cluster_package_shared$tag.simple ==
                       global_fisher_dataframe$pattern[x]]))

# Print cluster with different fisher score
global_fisher_dataframe[global_fisher_dataframe$fisherServer !=
                            global_fisher_dataframe$fisherPackage,]


global_cluster_package_unique <-
    global_cluster_package[!(global_cluster_package$tag.simple %in%
                                 global_cluster_server$pattern),]
global_cluster_package_unique



## local similarities
local_cluster_server <-
    server_output_cluster[!grepl("%", server_output_cluster$pattern) &
                              server_output_cluster$pattern != "single",]
local_cluster_package <-
    package_output_cluster[package_output_cluster$type == "local",]

cat("Number of local clusters\n")
cat("Server: ", nrow(local_cluster_server), "\n")   # server: 31 # local: 46
cat("Package: ", nrow(local_cluster_package), "\n") # package: 49


# Simplify tag
local_cluster_package$tag.simple <- gsub("_.*", "", local_cluster_package$tag)

# Get shared cluster results
local_cluster_server_shared <-
    local_cluster_server[local_cluster_server$pattern %in%
                             local_cluster_package$tag.simple,]
local_cluster_package_shared <-
    local_cluster_package[local_cluster_package$tag.simple %in%
                              local_cluster_server$pattern,]

cat("Number of shared local clusters\n")
cat("Server: ", nrow(local_cluster_server_shared), "\n")   # s: 31 # l: 46
cat("Package: ", nrow(local_cluster_package_shared), "\n") # p: 31 # p: 46


# Unite fisher score information in one date frame
local_fisher_dataframe <-
    data.frame(pattern = local_cluster_server_shared$pattern,
               fisherServer = local_cluster_server_shared$Fisher_score)
local_fisher_dataframe$fisherPackage <-
    sapply(1:nrow(local_fisher_dataframe),
           function(x)
               return(local_cluster_package_shared$fisher.score[
                   local_cluster_package_shared$tag.simple ==
                       local_fisher_dataframe$pattern[x]][1]))

# "GLVL" occurs two times in data frame, add second fisher score manually
local_fisher_dataframe$fisherPackage[
    local_fisher_dataframe$pattern == "GLVL"][2] <-
    local_cluster_package_shared$fisher.score[
        local_cluster_package_shared$tag.simple == "GLVL"][2]

# Print cluster with different fisher score
local_fisher_dataframe[local_fisher_dataframe$fisherServer !=
                           local_fisher_dataframe$fisherPackage,]

# clusters only detected by the package

local_cluster_package_unique <-
    local_cluster_package[!(local_cluster_package$tag.simple %in%
                                local_cluster_server$pattern),]
local_cluster_package_unique
