
get_beta_violin <- function(beta,
                            node_summary,
                            ag,
                            ag_species = TRUE,
                            db = "vdjdb") {
    if(ag_species) {
        h <- get_ag_species_hits(node_summary = node_summary, 
                                 db = db, 
                                 ag_species = ag)
    } else {
        h <- get_ag_gene_hits(node_summary = node_summary, 
                              db = db, 
                              ag_gene = ag)
    }
    
    
    h$node_summary$spec <- ifelse(
        test = apply(X=h$node_summary[,h$new_columns], MARGIN=1, FUN=sum)==0,
        yes = "-", no = "+")
    h <- h$node_summary
    h <- aggregate(clone_size~sample+spec+community, data = h, FUN = sum)
    dup <- which(h$community %in% h$community[h$spec=="+"] & h$spec=="-")
    if(length(dup)!=0) {
        h <- h[-dup,]
    }
    
    h$size <- ifelse(test = h$spec == "+", yes = h$clone_size, no = 1)
    o <- merge(x = h, y = beta, by = c("sample", "community"), all.x = TRUE)
    
    g <- ggplot(data = o)+
        geom_sina(aes(x = sample, y = mean, col = spec, size = size), 
                  stroke = 0.02, alpha = 0.4)+
        theme_bw(base_size = 10)+
        ggtitle(label = '', 
                subtitle = paste0("Ag=", ag, ", DB=", db))+
        ylab(label = expression(beta))+
        xlab(label = "repertoire")+
        scale_radius(name = "# cells", 
                     breaks = scales::pretty_breaks(n = 4),
                     range = c(1, 4))+
        scale_color_manual(name = "specificity", 
                           values = c("steelblue", "orange"))+
        guides(size = guide_legend(order = 2), 
               colour = guide_legend(order = 1))
    
    return(g)
    
}
