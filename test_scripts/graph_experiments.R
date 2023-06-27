# tests with graphs (just a collection of code. not necessarily working)

library(ClustIRR)
library(knitr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(visNetwork)
library(turboGliph)
library(igraph)

# load package input data
data("CD8")

# set the random seed for sampling reproducibility
set.seed(987)

# sample 500 sequences from the reference dataset as sample dataset
data_sample <- data.frame(CDR3b = CD8[sample(x = 1:nrow(CD8), 
                                             size = 100, 
                                             replace = FALSE),])

# artificially enrich motif 'RQWW' inside sample dataset
substr(x = data_sample$CDR3b[1:5], start = 6, stop = 9) <- "RQWW"

# add the artificial clonal expansion of two sequences to the sample dataset
clones <- data.frame(CDR3b = rep(x = c("CATSRAAKPDGLAALETQYF",
                                       "CATSRAAKPDGLAALSTQYF"),
                                 times = 15))
data_sample <- rbind(data_sample, clones)


## clustirr

# ClustIRR version 2
clustirr_output_v2 <- cluster_irr(data_sample = data_sample,  
                                  data_ref = CD8,
                                  version = 2,
                                  ks = 4,
                                  cores = 1,
                                  control = list(
                                   B = 1000,
                                   global_max_dist = 1,
                                   local_max_fdr = 0.05,
                                   local_min_ove = 2,
                                   local_min_o = 3,
                                   trim_flank_aa = 3,
                                   low_mem = FALSE,
                                   global_pairs = NULL))

# ClustIRR version 3
clustirr_output_v3 <- cluster_irr(data_sample = data_sample,  
                                  data_ref = CD8,
                                  version = 3,
                                  ks = 4,
                                  cores = 1,
                                  control = list(
                                   B = 1000,
                                   global_max_dist = 1,
                                   local_max_fdr = 0.05,
                                   local_min_ove = 2,
                                   local_min_o = 3,
                                   trim_flank_aa = 3,
                                   low_mem = FALSE,
                                   global_pairs = NULL))


## v2 visualisation

# quick and dirty
gl_v2 <- data.frame(cdr3=c("CASSYVLRAADGYTF", "CASSYVLRAADGDTF",
                           "CATSRAAKPDGLAALSTQYF", "CATSRAAKPDGLAALETQYF"),
                    type="global", 
                    group=c(7,7,8,8))

lc_v2 <- clustirr_output_v2$clust$CDR3b$local$lp
lc_v2$group <- as.numeric(factor(lc_v2$motif))
lc_v2$type <- "local"
nodes <- bind_rows(lc_v2, gl_v2)
nodes$id <- seq.int(nrow(nodes))
nodes$title <- nodes$cdr3

edges <- subset(nodes, select = c(id, type, group))
edges <- transform(edges, color = ifelse(type=="local", "#E5E5E5", "#4D4D4D"))
edges$from <- edges$id
edges$to <- 0
edges$width <- 10

j = 0
for(i in 1:nrow(edges)){
 if (j == 0){
  first <- i
 }
 j = i+1
 if(isTRUE(edges$group[i] == edges$group[j])){
  edges$to[i] <- j
  
 } else {
  edges$to[i] <- first
  j = 0
 }
}

# fix circle / self ref
for(i in 1:nrow(edges)){
 dst <- edges$to[i]
 if(!dst == 0){
  if(edges$to[dst] == i){
   edges$to[dst] <- 0
  }
 }
 if (edges$from[i] == edges$to[i])
 {
  edges$to[i] = 0
 }
}

visNetwork(nodes, edges) 


## v3

# quick and dirty
gl_v3 <- data.frame(cdr3=c("CASSYVLRAADGYTF", "CASSYVLRAADGDTF",
                           "CATSRAAKPDGLAALSTQYF", "CATSRAAKPDGLAALETQYF"),
                    type="global", 
                    group=c(21,21,22,22))

lc_v3 <- clustirr_output_v3$clust$CDR3b$local$lp
lc_v3$group <- as.numeric(factor(lc_v3$motif))
lc_v3$type <- "local"
nodes <- bind_rows(lc_v3, gl_v3)
nodes$id <- seq.int(nrow(nodes))
nodes$title <- nodes$cdr3

edges <- subset(nodes, select = c(id, type, group))
edges <- transform(edges, color = ifelse(type=="local", "#E5E5E5", "#4D4D4D"))
edges$from <- edges$id
edges$to <- 0
edges$width <- 10

j = 0
for(i in 1:nrow(edges)){
 if (j == 0){
  first <- i
 }
 j = i+1
 if(isTRUE(edges$group[i] == edges$group[j])){
  edges$to[i] <- j
  
 } else {
  edges$to[i] <- first
  j = 0
 }
}

# fix circle / self ref
for(i in 1:nrow(edges)){
 dst <- edges$to[i]
 if(!dst == 0){
  if(edges$to[dst] == i){
   edges$to[dst] <- 0
  }
 }
 if (edges$from[i] == edges$to[i])
 {
  edges$to[i] = 0
 }
}

visNetwork(nodes, edges) 

## test get edges

gv2 <- get_edges(clustirr_output_v2)
nodes <- data.frame(id = 1:5)
nodes$title <- nodes$id
edges <- data.frame(from = gv2$from_ID, to = gv2$to_ID)

visNetwork(nodes, edges)

## try with igraph
# out_df$color <- ifelse(test = out_df$type == "global", 
#                        yes = "grey30", no = "grey90")
# 
# g <- graph_from_data_frame(out_df)
# 
# plot(g,
#      edge.width = 2,
#      edge.arrow.size = 0.1,
#      vertex.label = V(g)$counts,
#      vertex.label.cex = 1,
#      edge.color = get.edge.attribute(graph = g, name = "color"),
#      main = paste0("test"))




## turboGliph tests
tg <- turbo_gliph(cdr3_sequences = data_sample, 
                  refdb_beta = CD8, 
                  lcminove = 10, 
                  motif_length = 4,
                  gccutoff = 1)

g2 <- gliph2(cdr3_sequences = data_sample, 
             refdb_beta = CD8, 
             lcminove = 1000, 
             motif_length = 4)

plot_network(tg)

plot_network(g2)



## igraph
data_sample <- data.frame(CDR3b = CD8[1:5000, "CDR3b"])
data_ref <- CD8

# run analysis
out <- cluster_irr(data_sample = data_sample,
                   data_ref = data_ref,
                   version = 3,
                   ks = 4,
                   cores = 1,
                   control = list(
                     B = 1000,
                     global_max_dist = 1,
                     local_max_fdr = 0.05,
                     local_min_ove = 2,
                     local_min_o = 1,
                     trim_flank_aa = 3,
                     global_pairs = NULL,
                     low_mem = FALSE))

# create data.frame
out_df <- get_edges(out)

out_df$color <- ifelse(test = out_df$type == "global", 
                     yes = "grey30", no = "grey90")

# visualize data with igraph
g <- graph_from_data_frame(out_df)



plot(g,
     edge.width = 2,
     edge.arrow.size = 0.1,
     vertex.label = V(g)$counts,
     vertex.label.cex = 1,
     edge.color = get.edge.attribute(graph = g, name = "color"),
     main = paste0("test"))

plot(g)
