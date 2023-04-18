require(turboGliph)
require(ggplot2)
require(ggforce)
require(patchwork)
require(gliphR)


data("hs_CD8_ref")
ref <- base::grep(pattern = "^C.*F$",x = hs_CD8_ref$CDR3b,
                  perl = TRUE,value = TRUE)
ref <- hs_CD8_ref[hs_CD8_ref$CDR3b %in% ref, ]
rm(hs_CD8_ref)

set.seed(seed = 1298)
data_sample <- ref[sample(x = 1:nrow(ref), size = 5000, replace = F), 1:3]
data_ref <- ref[, 1:3]
ks <- c(2, 3, 4)
cores <- 1
rm(ref)

# enrich the first CDR3 by +19 = 20fold
for(i in 1:19) {
  data_sample <- rbind(data_sample, data_sample[1, ])
}

control_input <- list(
  B = 1000,
  global_max_dist = 1,
  local_max_fdr = 0.05,
  local_min_ove = 2,
  local_min_o = 1,
  trim_flanks = FALSE,
  flank_size = 3)



#### Run ####
out_v1 <- gliph(data_sample = data_sample,
                data_ref = data_ref,
                ks = ks,
                cores = cores,
                version = 1,
                control = control_input)


out_v2 <- gliph(data_sample = data_sample,
                data_ref = data_ref,
                ks = ks,
                cores = cores,
                version = 2,
                control = control_input)


out_v3 <- gliph(data_sample = data_sample,
                data_ref = data_ref,
                ks = ks,
                cores = cores,
                version = 3,
                control = control_input)


# save data for benchmarking
dir.create(path = "test_scripts/benchmark_expansion")
save(out_v1, file = "test_scripts/benchmark_expansion/out_v1.RData")
save(out_v2, file = "test_scripts/benchmark_expansion/out_v2.RData")
save(out_v3, file = "test_scripts/benchmark_expansion/out_v3.RData")

# save data for benchmarking
save(data_ref, file = "test_scripts/benchmark_expansion/data_ref.RData")
save(data_sample, file = "test_scripts/benchmark_expansion/data_sample.RData")

# create ting-compatible data
data_sample$TRBV <- gsub(pattern = "TRB", replacement = '', x = data_sample$TRBV)
data_sample$TRBJ <- gsub(pattern = "TRB", replacement = '', x = data_sample$TRBJ)
write.table(x = data_sample, file = "test_scripts/benchmark_expansion/data_sample.tsv",
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)

data_ref$TRBV <- gsub(pattern = "TRB", replacement = '', x = data_ref$TRBV)
data_ref$TRBJ <- gsub(pattern = "TRB", replacement = '', x = data_ref$TRBJ)
write.table(x = data_ref, file = "test_scripts/benchmark_expansion/data_ref.tsv",
            append = F, quote = F, sep = "\t", row.names = F, col.names = T)


# ting -t data_sample.tsv -r data_ref.tsv -k "ting_bench_kmers" -o "ting_output"
# ting -t data_sample.tsv -r data_ref.tsv -k "ting_bench_kmers" -o "ting_output" -p 1000000000

#
# CDR3b		TRBV	TRBJ	CDR3a		TRAV		TRAJ	Sample-ID
# CAADTSSGANVLTF	TRBV30	TRBJ2-6	CALSDEDTGRRALTF	TRAV19		TRAJ5	09/02171
# CAATGGDRAYEQYF	TRBV2	TRBJ2-7	CAASSGANSKLTF	TRAV13-1	TRAJ56	03/04922
# CAATQQGETQYF	TRBV2	TRBJ2-5	CAASYGGSARQLTF	TRAV13-1	TRAJ22	02/02591
# CACVSNTEAFF	TRBV28	TRBJ1-1	CAGDLNGAGSYQLTF	TRAV25		TRAJ28	PBMC8631
# CAGGKGNSPLHF	TRBV2	TRBJ1-6	CVVLRGGSQGNLIF	TRAV12-1	TRAJ42	02/02071
# CAGQILAGSDTQYF	TRBV6-4	TRBJ2-3	CATASGNTPLVF	TRAV17		TRAJ29	09/00181
# CAGRTGVSTDTQYF	TRBV5-1	TRBJ2-3	CAVTPGGGADGLTF	TRAV41		TRAJ45	02/02591
# CAGYTGRANYGYTF	TRBV2	TRBJ1-2	CVVNGGFGNVLHC	TRAV12-1	TRAJ35	01/08733

v1 <- out_v1$clust$CDR3b$local$m
v2 <- out_v2$clust$CDR3b$local$m
v3 <- out_v3$clust$CDR3b$local$m

table(v1$pass)
table(v2$pass)
table(v3$pass)

w <- merge(x = v2, y = v3, by = "motif")
w$ove.change <- w$ove.x/w$ove.y
ggplot(data = w)+
  geom_point(aes(x = ove.x, y = ove.y))+
  geom_abline(slope = 1, intercept = 0)

ggplot(data = w)+
  geom_point(aes(x = p_value.x, y = p_value.y))+
  geom_abline(slope = 1, intercept = 0)+
  scale_x_log10()+
  scale_y_log10()


u <- out_v2$clust$CDR3b$motif_enrichment
u <- out_v3$clust$CDR3b$motif_enrichment
u[u$motif %in% c("REVQ", "YIN", "EGLP", "IGRA", "LREA", "QVAR"),]
out_v2$clust$CDR3b$motif_enrichment[which(out_v2$clust$CDR3b$motif_enrichment$filter),]
out_v3$clust$CDR3b$motif_enrichment[which(out_v3$clust$CDR3b$motif_enrichment$filter),]

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

