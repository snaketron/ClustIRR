
get_community_cosine_similarity <- function(com) {
    n <- sqrt(colSums(com^2))
    cos_sim <- crossprod(com)/(n %o% n)
    colnames(cos_sim) <- colnames(com)
    rownames(cos_sim) <- colnames(com)
    cos_sim <- reshape2::melt(cos_sim, varnames = c("i", "j"), value.name="CS")
    
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
