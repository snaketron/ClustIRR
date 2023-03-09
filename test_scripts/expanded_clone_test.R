require(turboGliph)
require(ggplot2)
require(ggforce)
require(patchwork)

data("hs_CD8_ref")

set.seed(seed = 4123)
data("hs_CD8_ref")
data_sample <- hs_CD8_ref[sample(x = 1:nrow(hs_CD8_ref), 
                                 size = 2000, replace = F), 1:3]
data_ref <- hs_CD8_ref[, 1:3]
ks <- c(2, 3, 4)
B <- 1000
cores <- 1

# enrich the first CDR3 by +9
for(i in 1:9) {
  data_sample <- rbind(data_sample, data_sample[1, ])
}


# control <- NULL
# control <- get_control(control_in = control)



source("R/gliph_v1.R")
source("R/gliph_v2.R")
source("R/util_v1.R")
source("R/util_v2.R")
source("R/util_v1_v2.R")



o1 <- gliph_v1(
  data_sample = data_sample,
  data_ref = data_ref,
  ks = ks,
  B = B,
  cores = cores,
  control = NULL)

o2 <- gliph_v2(
  data_sample = data_sample,
  data_ref = data_ref,
  ks = ks,
  cores = cores,
  control = NULL)

o4 <- turboGliph::gliph2(
  cdr3_sequences = data_sample$CDR3b,
  result_folder = "/home/sktron/Desktop/tmp/",
  refdb_beta = data_ref$CDR3b,
  lcminp = 0.05)



# we would expect that motifs found in the expanded clone 
# are found as enriched
# CASSLEVGKSRVTGELFF
# "VGKS" "SRVT" "LEVG" "RVTG" "EVGK" "GKSR" "KSRV"
x <- o1$motif_enrichment
y <- o2$motif_enrichment
z <- o4$motif_enrichment$all_motifs

motifs_in_clone <- unique(
    c(names(stringdist::qgrams("CASSLEVGKSRVTGELFF", q = 2)[1,]),
      names(stringdist::qgrams("CASSLEVGKSRVTGELFF", q = 3)[1,]),
      names(stringdist::qgrams("CASSLEVGKSRVTGELFF", q = 4)[1,])))

x$in_clone <- ifelse(test = x$motif %in% motifs_in_clone, yes = T, no = F)
y$in_clone <- ifelse(test = y$motif %in% motifs_in_clone, yes = T, no = F)
z$in_clone <- ifelse(test = z$motif %in% motifs_in_clone, yes = T, no = F)
z$k <- nchar(z$motif)
z$filter <- ifelse(test = ((z$fisher.score<=0.05&z$num_in_sample>=3)
                           &((z$num_fold>=1000&z$k==2)
                             |(z$num_fold>=100&z$k==3)
                             |(z$num_fold>=10&z$k==4))), 
                   yes = T, no = F)

x$fdr <-p.adjust(p = x$p, method = "fdr")
y$fdr <-p.adjust(p = y$p_value, method = "fdr")
z$fdr <-p.adjust(p = z$fisher.score, method = "fdr")

# Motifs of large clone are naturally not identified by 
# the TurboGliph implementation

# Top row: Gliph 1 returns p=0 for outliers (result of limited B), hence
# -log10(0) -> not shown as point in top-row, left panel
(ggplot()+
        geom_sina(data = x, 
                  aes(x = "Gliph_1", y = -log10(p), col = in_clone, 
                      shape = filter, group = in_clone))+
        geom_sina(data = y, 
                  aes(x = "Gliph_2", y = -log10(p_value), col = in_clone, 
                      shape = filter, group = in_clone))+
        geom_sina(data = z, 
                  aes(x = "Gliph_2 (Jan)", y = -log10(fisher.score), col = in_clone, 
                      shape = filter, group = in_clone))+
        theme_bw()+
        scale_shape_manual(values = c(21, 19))+
        theme(legend.position = "top"))/
    (ggplot()+
         geom_sina(data = x, 
                   aes(x = "Gliph_1", y = ove, col = in_clone, 
                       shape = filter, group = in_clone))+
         geom_sina(data = y, 
                   aes(x = "Gliph_2", y = ove, col = in_clone, 
                       shape = filter, group = in_clone))+
         geom_sina(data = z, 
                   aes(x = "Gliph_2 (Jan)", y = num_fold, col = in_clone,
                       shape = filter, group = in_clone))+
         theme(legend.position = "top")+
         theme_bw()+
         scale_shape_manual(values = c(21, 19))+
         theme(legend.position = "top")+
         scale_y_log10())


x_f <- x[x$filter==TRUE,]
y_f <- y[y$filter==TRUE,]
z_f <- z[z$filter==TRUE,]

sort(intersect(x_f$motif, z_f$motif))
sort(intersect(y_f$motif, z_f$motif))
sort(intersect(x_f$motif, y_f$motif))


u <- merge(x = x, y = y, by = "motif", all = T)
u <- merge(x = o1$motif_enrichment, 
           y = o2$motif_enrichment, 
           by = "motif", 
           all = T)

ggplot(data = u)+
    geom_point(aes(y = -log10(p_value), x = -log10(p)))+
    geom_abline(slope = 1, intercept = 0, col = "red")+
    theme_bw()

