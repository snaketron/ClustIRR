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


##### Dataset D1: TCRab repertoires 'a', 'b', and 'c' with n=500 #####
# repertoire size
n <- 500

# a
a <- data.frame(CDR3a = CDR3ab$CDR3a[1:n], 
                CDR3b = CDR3ab$CDR3b[1:n], 
                sample = "a")

# metadata a
ma <- data.frame(clone_id = 1:n, 
                 cell = sample(x = c("CD8", "CD4"), replace = TRUE, size = n),
                 HLA_A = "HLA-A∗24",
                 age = 24,
                 TRAV = CDR3ab$TRAV[1:n],
                 TRAJ = CDR3ab$TRAJ[1:n],
                 TRBV = CDR3ab$TRBV[1:n],
                 TRBJ = CDR3ab$TRBJ[1:n])
                 
# b
b <- data.frame(CDR3a = CDR3ab$CDR3a[1:n], 
                CDR3b = CDR3ab$CDR3b[1:n], 
                sample = "b")

# metadata b
mb <- data.frame(clone_id = 1:n, 
                 cell = sample(x = c("CD8", "CD4"), replace = TRUE, size = n),
                 HLA_A = "HLA-A∗02",
                 age = 30,
                 TRAV = CDR3ab$TRAV[1:n],
                 TRAJ = CDR3ab$TRAJ[1:n],
                 TRBV = CDR3ab$TRBV[1:n],
                 TRBJ = CDR3ab$TRBJ[1:n])

# c
c <- data.frame(CDR3a = CDR3ab$CDR3a[1:n], 
                CDR3b = CDR3ab$CDR3b[1:n], 
                sample = "c")

# metadata c
mc <- data.frame(clone_id = 1:n, 
                 cell = sample(x = c("CD8", "CD4"), replace = TRUE, size = n),
                 HLA_A = "HLA-A∗11",
                 age = 40,
                 TRAV = CDR3ab$TRAV[1:n],
                 TRAJ = CDR3ab$TRAJ[1:n],
                 TRBV = CDR3ab$TRBV[1:n],
                 TRBJ = CDR3ab$TRBJ[1:n])


get_clonal_expansion <- function(n, p_expanded) {
    s <- sample(x = c(0, 1), size = n, prob = c(1-p_expanded, 
                                                p_expanded), replace = T)
    y <- vapply(X = s, FUN.VALUE = numeric(1), FUN = function(x) {
        if(x == 0) {
            return(rpois(n = 1, lambda = 0.5))
        }
        return(rpois(n = 1, lambda = 50))
    })
    return(y)
}

# simulate expansion of specific communities
set.seed(1243)
clone_size <- rpois(n = n, lambda = 3)+1
expansion_factor <- get_clonal_expansion(n = n, p_expanded = 0.02)

a$clone_size <- clone_size
b$clone_size <- clone_size+expansion_factor*1
c$clone_size <- clone_size+expansion_factor*2

D1 <- list(a = a, b = b, c = c, ma = ma, mb = mb, mc = mc)
save(D1, file = "data/D1.RData", compress = "xz")
rm(a, b, c, ma, mb, mc)
