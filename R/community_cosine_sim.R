
get_cosine_similarity <- function(com) {
    cos_sim <- matrix(data = 0, nrow = ncol(com), ncol = ncol(com))
    for(i in 1:ncol(com)) {
        for(j in 1:ncol(com)) {
            v1 <- com[,i]
            v2 <- com[,j]
            cos_sim[i,j] <- sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
        }
    }
    colnames(cos_sim) <- colnames(com)
    rownames(cos_sim) <- colnames(com)
    cos_sim <- reshape2::melt(cos_sim)
    colnames(cos_sim) <- c("i", "j", "CS")
    
    g <- ggplot(data = cos_sim)+
        geom_tile(aes(x = i, y = j, fill = CS), col = "white")+
        geom_text(aes(x = i, y = j, label = round(CS, digits = 2)), 
                  size = 2.75, col = "black")+
        scale_fill_distiller(palette = "Spectral", 
                             limits = c(-0.01, 1.01),
                             breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                             labels = c(0, 0.25, 0.5, 0.75, 1.0))+
        theme(legend.position = "right")+
        guides(fill = guide_colourbar(barheight = 5, barwidth = 0.5))+
        xlab(label = '')+
        ylab(label = '')
    
    return(list(g = g, cs = cos_sim))
}
