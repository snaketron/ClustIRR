# This script describes how to download and prepare the small example
# dataset of human CD8+ T cells. The original data is publicly available
# at http://50.255.35.37:8080/

# download data (last tested on 09 January 2025)
utils::download.file(url = "http://50.255.35.37:8080/downloads/human_v2.0.zip",
                     destfile = "human_v2.0.zip")

# unzip
utils::unzip(zipfile = "human_v2.0.zip", exdir = ".")

# human (hs) CD8
CD8 <- read.csv(file = "human_v2.0/ref_CD8_v2.0.txt", sep = "\t",
                header = FALSE)

set.seed(seed = 1234321)
CDR3b <- CD8[sample(x = 1:nrow(CD8), size = 10^4, replace = FALSE), ]
colnames(CDR3b) <- c("CDR3b", "TRBV", "TRBJ")

set.seed(seed = 1245421)
CDR3a <- CD8[sample(x = 1:nrow(CD8), size = 10^4, replace = FALSE), ]
colnames(CDR3a) <- c("CDR3a", "TRAV", "TRAJ")

# create dummy CDR3ab
CDR3ab <- cbind(CDR3a, CDR3b)
CDR3ab$TRAV <- gsub(pattern = "TRB", replacement = "TRA", x = CDR3ab$TRAV)
CDR3ab$TRAJ <- gsub(pattern = "TRB", replacement = "TRA", x = CDR3ab$TRAJ)
save(CDR3ab, file = "data/CDR3ab.RData", compress = "xz")

