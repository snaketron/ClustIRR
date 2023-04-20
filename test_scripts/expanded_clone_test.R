# require(turboGliph)
require(ggplot2)
require(ggforce)
require(patchwork)
require(gliphR)


data("hs_CD8_ref")
ref <- base::grep(pattern = "^C.*F$",
                  x = hs_CD8_ref$CDR3b,
                  perl = TRUE,
                  value = TRUE)
ref <- hs_CD8_ref[hs_CD8_ref$CDR3b %in% ref, ]
rm(hs_CD8_ref)

set.seed(seed = 1298)
data_sample <- ref[sample(x = 1:nrow(ref), size = 2000, replace = F), 1:3]
data_ref <- ref[, 1:3]
ks <- c(2, 3, 4)
cores <- 1
rm(ref)

# enrich the first CDR3 by +19 = 20fold
for(i in 1:49) {
  data_sample <- rbind(data_sample, data_sample[1, ])
}

control_input <- list(
  B = 10000,
  global_max_dist = 1,
  local_max_fdr = 0.05,
  local_min_ove = 2,
  local_min_o = 1,
  trim_flanks = TRUE,
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


v1 <- out_v1$clust$CDR3b$local$m
v2 <- out_v2$clust$CDR3b$local$m
v3 <- out_v3$clust$CDR3b$local$m

table(v1$pass)
table(v2$pass)
table(v3$pass)



w <- merge(x = v1, y = v2, by = "motif")
w$ove.change <- w$ove.x/w$ove.y
g1 <- ((ggplot(data = w)+
          geom_point(aes(x = ove.x, y = ove.y, col = as.character(k.x)), size = 1)+
          geom_abline(slope = 1, intercept = 0)+
          scale_color_discrete(name = "k")+
          theme_bw())|
         (ggplot(data = w)+
            geom_point(aes(x = fdr.x, y = fdr.y, col = as.character(k.x)), size = 1)+
            geom_abline(slope = 1, intercept = 0)+
            scale_x_log10()+
            scale_y_log10()+
            scale_color_discrete(name = "k")+
            theme_bw()))



w <- merge(x = v2, y = v3, by = "motif")
w$ove.change <- w$ove.x/w$ove.y
g2 <- ((ggplot(data = w)+
          geom_point(aes(x = ove.x, y = ove.y, col = as.character(k.x)), size = 1)+
          geom_abline(slope = 1, intercept = 0)+
          scale_color_discrete(name = "k")+
          theme_bw())|
      (ggplot(data = w)+
          geom_point(aes(x = fdr.x, y = fdr.y, col = as.character(k.x)), size = 1)+
          geom_abline(slope = 1, intercept = 0)+
          scale_x_log10()+
          scale_y_log10()+
         scale_color_discrete(name = "k")+
          theme_bw()))

g1/g2



v3$motif[which(v3$fdr<=10^(-10))]
# "EY"   "QP"   "YQ"
# "EYQ"  "GEY"  "PRG"  "QPR"  "SQP"  "YQP"
# "EYQP" "GEYQ" "GGEY" "PRGG" "QPRG" "RGGE" "SQPR"
# CASSQPRGGEYQPQHF
#    SQPR
#       RGGE
#     QPRG
#    SQPR
#      PRGG
#         GEYQ
#        GGEY
#          EYQP










#### sanity check about OvE: OK ####

# V3
sample_trim <- get_trimmed_flanks(x = data_sample$CDR3b, flank_size = 3)
ref_trim <- get_trimmed_flanks(x = data_ref$CDR3b, flank_size = 3)

k <- 4
ms <- v3$motif[v2$k == k]
n_s <- sum(stringdist::qgrams(sample_trim, q = k)[1,])
n_f <- sum(stringdist::qgrams(ref_trim, q = k)[1,])

y_s <- stringdist::qgrams(sample_trim, q = k)[1,]
y_f <- stringdist::qgrams(ref_trim, q = k)[1,]

ove <- c()
for(m in ms) {
  cat(m, "\n")
  a <- as.numeric(y_s[m])
  if(length(a)==0) {
    a <- 0
  }
  b <- as.numeric(y_f[m])
  if(length(b)==0) {
    b <- 0
  }
  a <- a/n_s
  b <- b/n_f
  ove <- c(ove, a/b)
}
plot(v3$ove[v3$k == k], ove);abline(0,1)
rm(a, b, m)




# V2
sample_trim <- get_trimmed_flanks(x = unique(data_sample$CDR3b), flank_size = 3)
ref_trim <- get_trimmed_flanks(x = unique(data_ref$CDR3b), flank_size = 3)

k <- 4
ms <- v2$motif[v2$k == k]
n_s <- sum(stringdist::qgrams(sample_trim, q = k)[1,])
n_f <- sum(stringdist::qgrams(ref_trim, q = k)[1,])

y_s <- stringdist::qgrams(sample_trim, q = k)[1,]
y_f <- stringdist::qgrams(ref_trim, q = k)[1,]

ove <- c()
for(m in ms) {
  cat(m, "\n")
  a <- as.numeric(y_s[m])
  if(length(a)==0) {
    a <- 0
  }
  b <- as.numeric(y_f[m])
  if(length(b)==0) {
    b <- 0
  }
  a <- a/n_s
  b <- b/n_f
  ove <- c(ove, a/b)
}
plot(v2$ove[v2$k == k], ove);abline(0,1)
rm(a, b, m)





# V1
sample_trim <- get_trimmed_flanks(x = unique(data_sample$CDR3b), flank_size = 3)
ref_trim <- get_trimmed_flanks(x = unique(data_ref$CDR3b), flank_size = 3)

k <- 4
ms <- v1$motif[v1$k == k]
n_s <- sum(stringdist::qgrams(sample_trim, q = k)[1,])
n_f <- sum(stringdist::qgrams(ref_trim, q = k)[1,])

y_s <- stringdist::qgrams(sample_trim, q = k)[1,]
y_f <- stringdist::qgrams(ref_trim, q = k)[1,]

ove <- c()
for(m in ms) {
  cat(m, "\n")
  a <- as.numeric(y_s[m])
  if(length(a)==0) {
    a <- 0
  }
  b <- as.numeric(y_f[m])
  if(length(b)==0) {
    b <- 0
  }
  a <- a/n_s
  b <- b/n_f
  ove <- c(ove, a/b)
}
plot(v1$ove[v1$k == k], ove);abline(0,1)
rm(a, b, m)


