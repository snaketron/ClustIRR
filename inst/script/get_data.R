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
CDR3b <- CD8[sample(x = 1:nrow(CD8), size = 10^4, replace = FALSE), ]
colnames(CDR3b) <- c("CDR3b", "TRBV", "TRBJ")
CDR3b$TRBV <- NULL
CDR3b$TRBJ <- NULL

set.seed(seed = 1245421)
CDR3a <- CD8[sample(x = 1:nrow(CD8), size = 10^4, replace = FALSE), ]
colnames(CDR3a) <- c("CDR3a", "TRBV", "TRBJ")
CDR3a$TRBV <- NULL
CDR3a$TRBJ <- NULL

# create dummy CDR3ab
CDR3ab <- data.frame(CDR3a = CDR3a$CDR3a, CDR3b = CDR3b$CDR3b)
save(CDR3ab, file = "data/CDR3ab.RData", compress = "xz")

