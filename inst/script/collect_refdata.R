# This script describes how to download and prepare the reference (naive) 
# CD8 and CD4 TCR repertoires of humans and mice. Both datasets are provided
# by the maintainers of GLIPH2 (http://50.255.35.37:8080/).

# create temporary directory and download the scRNA-seq data:
if(dir.exists("temp_folder")==FALSE) {
    dir.create(path = "temp_folder")
}

# download data (last tested on 02. March 2023) 
utils::download.file(
    url = "http://50.255.35.37:8080/downloads/human_v2.0.zip",
    destfile = "temp_folder/human_v2.0.zip")

utils::download.file(
    url = "http://50.255.35.37:8080/downloads/mouse_v1.0.zip",
    destfile = "temp_folder/mouse_v1.0.zip")

# unzip
utils::unzip(zipfile = "temp_folder/human_v2.0.zip", exdir = "temp_folder/")
utils::unzip(zipfile = "temp_folder/mouse_v1.0.zip", exdir = "temp_folder/")


# human (hs) CD8
CD8 <- read.csv(file = "temp_folder/human_v2.0/ref_CD8_v2.0.txt", sep = "\t", 
                header = FALSE)
colnames(CD8) <- c("CDR3b", "TRBV", "TRBJ")
CD8$CDR3a <- NA
CD8$TRAV <- NA
CD8$TRAJ <- NA
hs_CD8_ref <- CD8
rm(CD8)
save(hs_CD8_ref, file = "data/hs_CD8_ref.RData")


# human (hs) CD4
CD4 <- read.csv(file = "temp_folder/human_v2.0/ref_CD4_v2.0.txt", sep = "\t", 
                header = FALSE)
colnames(CD4) <- c("CDR3b", "TRBV", "TRBJ")
CD4$CDR3a <- NA
CD4$TRAV <- NA
CD4$TRAJ <- NA
hs_CD4_ref <- CD4
rm(CD4)
save(hs_CD4_ref, file = "data/hs_CD4_ref.RData")



# mouse (mm) CD8
CD8 <- read.csv(file = "temp_folder/mouse_v1.0/ref_CD8_ms.txt", sep = "\t", 
                header = FALSE)
colnames(CD8) <- c("CDR3b", "TRBV", "TRBJ")
CD8$CDR3a <- NA
CD8$TRAV <- NA
CD8$TRAJ <- NA
mm_CD8_ref <- CD8
rm(CD8)
save(mm_CD8_ref, file = "data/mm_CD8_ref.RData")


# mouse (mm) CD4
CD4 <- read.csv(file = "temp_folder/mouse_v1.0/ref_CD4_ms.txt", sep = "\t", 
                header = FALSE)
colnames(CD4) <- c("CDR3b", "TRBV", "TRBJ")
CD4$CDR3a <- NA
CD4$TRAV <- NA
CD4$TRAJ <- NA
mm_CD4_ref <- CD4
rm(CD4)
save(mm_CD4_ref, file = "data/mm_CD4_ref.RData")