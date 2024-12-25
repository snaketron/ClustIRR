# Download McPAS-TCR data
# McPAS-TCR - http://friedmanlab.weizmann.ac.il/McPAS-TCR/
t <- read.csv("McPAS-TCR.csv", sep = ",")
t <- t[, c("CDR3.alpha.aa", "CDR3.beta.aa", "Species", "Pathology",
           "Antigen.protein", "PubMed.ID")]
t <- t[duplicated(t)==FALSE, ]

colnames(t) <- c("CDR3a", "CDR3b", "CDR3_species", "Antigen_species",
                 "Antigen_gene", "Reference")

t$CDR3g <- NA
t$CDR3d <- NA
t$CDR3h <- NA
t$CDR3l <- NA

t <- t[, c("CDR3a", "CDR3b",
           "CDR3g", "CDR3d", 
           "CDR3h", "CDR3l", 
           "CDR3_species", "Antigen_species", 
           "Antigen_gene", "Reference")]

mcpas <- t
mcpas$Antigen_gene <- gsub(pattern = "\\|", replacement = ',', 
                           x = mcpas$Antigen_gene)
save(mcpas, file = "data/mcpas.RData", compress = T)
