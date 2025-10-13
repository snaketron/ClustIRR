
get_honeycombs <- function(com) {
    n <- ncol(com)
    gs <- vector(mode = "list", length = n*(n-1)/2)
    iter <- 1
    for(i in 1:(ncol(com)-1)) {
        for(j in (i+1):ncol(com)) {
            
            m <- data.frame(x = com[,i], 
                            y = com[,j], 
                            community = 1:nrow(com),
                            sample_1 = colnames(com)[i], 
                            sample_2 = colnames(com)[j],
                            contrast = paste0(colnames(com)[i], "-", colnames(com)[j]),
                            title = paste0(colnames(com)[i], ' (x) vs ', 
                                           colnames(com)[j], " (y)"))
            
            m$x_adj <- ifelse(test = m$x==0, yes = m$x+0.5, no = m$x)
            m$y_adj <- ifelse(test = m$y==0, yes = m$y+0.5, no = m$y)
            
            gs[[iter]] <- ggplot(data = m)+
                facet_wrap(facets = ~title)+
                geom_hex(aes(x = x_adj, y = y_adj, fill=log10(..count..)), 
                         col = "white", bins = 5)+
                geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
                geom_point(aes(x = x_adj, y = y_adj), 
                           alpha = 0.75, size = 0.5, stroke = 0.3)+
                scale_x_continuous(name = "x", trans = "log10")+
                scale_y_continuous(name = "y", trans = "log10")+
                annotation_logticks(base = 10, sides = "lb")+
                scale_fill_gradient(name = "log10(C)",
                                    low = "white", high = "#FFC68A")+
                theme_bw(base_size = 10)+
                theme(legend.position = "right",
                      strip.text.x = element_text(
                          margin = margin(0.04,0,0.04,0, "cm")))
            
            iter <- iter + 1
        }
    }
    
    return(gs)
}
