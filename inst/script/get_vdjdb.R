# Download VDJdb .tsv file from: https://vdjdb.cdr3.net/
t <- read.csv(file = "SearchTable-2025-01-09 15_15_58.769.tsv", sep = "\t")

tb <- t[t$Gene=="TRB", c("CDR3", "Species", "Epitope.species", 
                         "Epitope.gene", "Reference")]
tb$CDR3a <- NA
tb$CDR3b <- tb$CDR3
tb$CDR3 <- NULL
tb <- tb[, c("CDR3a", "CDR3b", "Species", "Epitope.species", 
             "Epitope.gene", "Reference")]

ta <- t[t$Gene=="TRA", c("CDR3", "Species", "Epitope.species", 
                         "Epitope.gene", "Reference")]
ta$CDR3a <- ta$CDR3
ta$CDR3b <- NA
ta$CDR3 <- NULL
ta <- ta[, c("CDR3a", "CDR3b", "Species", "Epitope.species", 
             "Epitope.gene", "Reference")]

t <- rbind(ta, tb)
t <- t[duplicated(t)==F,]

colnames(t) <- c("CDR3a", "CDR3b", "CDR3_species",
                 "Antigen_species", "Antigen_gene",
                 "Reference")

t$CDR3g <- NA
t$CDR3d <- NA
t$CDR3h <- NA
t$CDR3l <- NA

t <- t[, c("CDR3a", "CDR3b",
           "CDR3g", "CDR3d", 
           "CDR3h", "CDR3l", 
           "CDR3_species", "Antigen_species", 
           "Antigen_gene", "Reference")]

vdjdb <- t
save(vdjdb, file = "data/vdjdb.RData", compress = T)
