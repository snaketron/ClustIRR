# This script describes how to download and prepare the small example
# dataset of human CD8+ T cells. The original data is publicly available
# at http://50.255.35.37:8080/

# create temporary directory and download the scRNA-seq data:
if (dir.exists("temp_folder") == FALSE) {
    dir.create(path = "temp_folder")
}

# download data (last tested on 02. March 2023)
utils::download.file(
    url = "http://50.255.35.37:8080/downloads/human_v2.0.zip",
    destfile = "temp_folder/human_v2.0.zip")


# unzip
utils::unzip(zipfile = "temp_folder/human_v2.0.zip", exdir = "temp_folder/")

# human (hs) CD8
CD8 <- read.csv(file = "temp_folder/human_v2.0/ref_CD8_v2.0.txt", sep = "\t",
                header = FALSE)
set.seed(seed = 1234321)
CD8 <- CD8[sample(x = 1:nrow(CD8), size = 10^4, replace = FALSE), ]
colnames(CD8) <- c("CDR3b", "TRBV", "TRBJ")
CD8$TRBV <- NULL
CD8$TRBJ <- NULL
save(CD8, file = "data/CD8.RData", compress = "xz")

