colors <- c("black", "white", "black")

kmer_1 <- matrix(rep(c(-1,1)), nrow=5, ncol=4)
kmer_2 <- matrix(rep(c(-1,1)), nrow=25, ncol=16)
kmer_3 <- matrix(rep(c(-1,1)), nrow=125, ncol=64)
kmer_4 <- matrix(rep(c(-1,1)), nrow=625, ncol=256)
kmer_5 <- matrix(rep(c(-1,1)), nrow=3125, ncol=1024)

kmer_1[2,3] <- 0
kmer_2[10,6] <- 0
kmer_3[60,18] <- 0
kmer_4[100,60] <- 0
kmer_5[1220,600] <- 0

image(kmer_1, col = colors, main = "1mer")
image(kmer_2, col = colors, main = "2mer")
image(kmer_3, col = colors, main = "3mer")
image(kmer_4, col = colors, main = "4mer")
image(kmer_5, col = colors, main = "5mer")


