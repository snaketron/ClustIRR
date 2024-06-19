# Download TCR3d data
# https://tcr3d.ibbr.umd.edu/cancerseqs
# https://tcr3d.ibbr.umd.edu/virusseqs
t <- read.csv("tcr3d_cancer.csv", sep = ",")
t <- t[, c("CDR3.alpha.", "CDR3.beta.", "Cancer.BR.Type", "Antigen", "Reference")]
colnames(t) <- c("CDR3a", "CDR3b", "Antigen_species", "Antigen_gene", "Reference")
t_c <- t

t <- read.csv("tcr3d_virus.csv", sep = ",")
t <- t[, c("CDR3.alpha.", "CDR3.beta.", "Virus.BR.Type", "Antigen", "Reference")]
colnames(t) <- c("CDR3a", "CDR3b", "Antigen_species", "Antigen_gene", "Reference")
t_v <- t

t <- rbind(t_c, t_v)
t$CDR3_species <- "Human"
t$CDR3g <- NA
t$CDR3d <- NA
t$CDR3h <- NA
t$CDR3l <- NA

t <- t[, c("CDR3a", "CDR3b",
           "CDR3g", "CDR3d", 
           "CDR3h", "CDR3l", 
           "CDR3_species", "Antigen_species", 
           "Antigen_gene", "Reference")]

tcr3d <- t
save(tcr3d, file = "data/tcr3d.RData", compress = T)
