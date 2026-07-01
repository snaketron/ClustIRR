# First one has to download the raw data from ParseBiosciences:

# Performance of Evercode™ TCR v4 in Antigen-Stimulated Human T Cells, 
# https://www.parsebiosciences.com/datasets/performance-of-evercode-tcr-v4-in-antigen-stimulated-human-t-cells/; 
# Parse Biosciences; Accessed 2026 June 30

# Then extract two folders:
# CMV HLA-A0201 
# Influenza M1 v4 

# Then process the data and extract:
# require(Seurat)


t <- read.csv("data/CMV HLA-A0201 v4/clonotype_frequency.tsv", sep = "\t")
t$group <- "CMV"

b <- read.csv("data/CMV HLA-A0201 v4/barcode_report.tsv", sep = "\t")
b <- b[is.na(b$clonotype_id)==F, ]
b <- aggregate(cell_id~clonotype_id, data = b, FUN = paste0, collapse = ',')

t <- merge(t, b, by = "clonotype_id")

# Read in cell meta data
cell_meta <- read.csv("data/CMV HLA-A0201 v4/cell_metadata.csv", row.names = 1)
d <- ReadParseBio(data.dir = "data/CMV HLA-A0201 v4/")
d <- CreateSeuratObject(d, min.features = 100, min.cells = 100, names.field = 0, meta.data = cell_meta)
d@meta.data$cell_id <- rownames(d@meta.data)
d <- subset(d, subset = nFeature_RNA < 5000 & nCount_RNA < 20000)
d <- subset(d, subset = cell_id %in% unlist(strsplit(t$cell_id, split = "\\,")))
d <- NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)

genes_activated_t_cells <- c("CD69", "CD44", "SELL", "CD38","HLA-DRA", "CTLA4",
                             "GZMB", "PRF1", "IFNG", "IL2RA", "IL2RB", "MKI67")

e <- t(as.matrix(d@assays$RNA$data[
    rownames(d@assays$RNA$data) %in% c(genes_activated_t_cells), ]))
gex <- data.frame(e)
gex$cell_id <- rownames(e)
gex$sample <- "CMV"

s <- t[,c("TRA", "TRB", "group", "count")]
colnames(s) <- c("CDR3a", "CDR3b", "sample", "clone_size")

cmv_A02 <- list(s = s, meta = t, gex = gex)




t <- read.csv("data/CMV HLA-B35 v4/clonotype_frequency.tsv", sep = "\t")
t$group <- "CMV_B35"

b <- read.csv("data/CMV HLA-B35 v4/barcode_report.tsv", sep = "\t")
b <- b[is.na(b$clonotype_id)==F, ]
b <- aggregate(cell_id~clonotype_id, data = b, FUN = paste0, collapse = ',')

t <- merge(t, b, by = "clonotype_id")

# Read in cell meta data
cell_meta <- read.csv("data/CMV HLA-B35 v4/cell_metadata.csv", row.names = 1)
d <- ReadParseBio(data.dir = "data/CMV HLA-B35 v4/")
d <- CreateSeuratObject(d, min.features = 100, min.cells = 100, names.field = 0, meta.data = cell_meta)
d@meta.data$cell_id <- rownames(d@meta.data)
d <- subset(d, subset = nFeature_RNA < 5000 & nCount_RNA < 20000)
d <- subset(d, subset = cell_id %in% unlist(strsplit(t$cell_id, split = "\\,")))
d <- NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)

genes_activated_t_cells <- c("CD69", "CD44", "SELL", "CD38","HLA-DRA", "CTLA4",
                             "GZMB", "PRF1", "IFNG", "IL2RA", "IL2RB", "MKI67")

e <- t(as.matrix(d@assays$RNA$data[
    rownames(d@assays$RNA$data) %in% c(genes_activated_t_cells), ]))
gex <- data.frame(e)
gex$cell_id <- rownames(e)
gex$sample <- "CMV_B35"

s <- t[,c("TRA", "TRB", "group", "count")]
colnames(s) <- c("CDR3a", "CDR3b", "sample", "clone_size")

cmv_B35 <- list(s = s, meta = t, gex = gex)





t <- read.csv("data/Influenza M1 v4/clonotype_frequency.tsv", sep = "\t")
t$group <- "Flu"

b <- read.csv("data/Influenza M1 v4/barcode_report.tsv", sep = "\t")
b <- b[is.na(b$clonotype_id)==F, ]
b <- aggregate(cell_id~clonotype_id, data = b, FUN = paste0, collapse = ',')

t <- merge(t, b, by = "clonotype_id")

# Read in cell meta data
cell_meta <- read.csv("data/Influenza M1 v4/cell_metadata.csv", row.names = 1)
d <- ReadParseBio(data.dir = "data/Influenza M1 v4/")
d <- CreateSeuratObject(d, min.features = 100, min.cells = 100, names.field = 0, meta.data = cell_meta)
d@meta.data$cell_id <- rownames(d@meta.data)
d <- subset(d, subset = nFeature_RNA < 5000 & nCount_RNA < 20000)
d <- subset(d, subset = cell_id %in% unlist(strsplit(t$cell_id, split = "\\,")))
d <- NormalizeData(d, normalization.method = "LogNormalize", scale.factor = 10000)

genes_activated_t_cells <- c("CD69", "CD44", "SELL", "CD38","HLA-DRA", "CTLA4", 
                             "GZMB", "PRF1", "IFNG", "IL2RA", "IL2RB", "MKI67")


e <- t(as.matrix(d@assays$RNA$data[
    rownames(d@assays$RNA$data) %in% c(genes_activated_t_cells), ]))
gex <- data.frame(e)
gex$cell_id <- rownames(e)
gex$sample <- "Flu"

s <- t[,c("TRA", "TRB", "group", "count")]
colnames(s) <- c("CDR3a", "CDR3b", "sample", "clone_size")

flu <- list(s = s, meta = t, gex = gex)


d <- list(s = rbind(cmv_A02$s, cmv_B35$s, flu$s), 
          meta = rbind(cmv_A02$meta, cmv_B35$meta, flu$meta), 
          gex = rbind(cmv_A02$gex, cmv_B35$gex, flu$gex))
D_CMV_FLU <- d
save(D_CMV_FLU, file = "D_CMV_FLU.RData", compress = TRUE)


