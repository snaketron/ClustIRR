# Description:
#'
#' @return data.frame
#' @export
get_edges <- function(clust_irr) {

    if(!base::class(clust_irr) == "clust_irr"){
        base::stop("Input has to be object of class clust_irr")
    }
    
    le <- get_local_edges(clust_irr = clust_irr)
    ge <- get_global_edges(clust_irr = clust_irr)
    e <- base::rbind(le, ge)

    if(base::is.null(e) || base::nrow(e)==0) {
        return(NULL)
    }

    me <- base::vector(mode = "list", length = base::length(clust_irr$clust))
    base::names(me) <- base::names(clust_irr$clust)
    for(chain in base::names(clust_irr$clust)) {
        tmp_d <- clust_irr$inputs$data_sample[, c(chain, "ID")]
        tmp_e <- e[e$chain == chain, ]

        if(base::nrow(tmp_e)!=0 & base::nrow(tmp_d)!=0) {
            base::colnames(tmp_d) <- c("CDR3", "ID")
            me[[chain]] <- base::merge(x = base::merge(
                x = tmp_e, y = tmp_d, by.x = "from", by.y = "CDR3"),
                y = tmp_d, by.x = "to", by.y = "CDR3")
        }
    }
    me <- base::do.call(base::rbind, me)
    if(base::is.null(me)==FALSE && base::nrow(me)!=0) {
        me$from_ID <- me$ID.x
        me$to_ID <- me$ID.y
        me <- me[,c("from_ID", "to_ID", "from", "to", "motif", "type", "chain")]
        return(me)
    }
    return(NULL)
}

get_local_edges <- function(clust_irr) {

    get_lp <- function(x, lp, chain) {
        if(base::sum(lp$motif == x)>1) {
            p <- base::t(utils::combn(x = lp$cdr3[lp$motif == x], m = 2))
            return(base::data.frame(from = p[,1], to = p[,2],
                                    motif = x, type = "local",
                                    chain = chain))
        }
        return(NULL)
    }

    edges_local <- vector(mode = "list", length = length(clust_irr$clust))
    names(edges_local) <- names(clust_irr$clust)
    for(chain in names(clust_irr$clust)) {
        lp <- clust_irr$clust[[chain]]$local$lp
        if(base::is.null(lp)==FALSE && base::nrow(lp)!=0) {
            edges_local[[chain]] <- base::do.call(base::rbind, base::lapply(
                X = base::unique(lp$motif), FUN = get_lp, lp = lp,
                chain = chain))
        }
    }
    edges_local <- base::do.call(base::rbind, edges_local)
    return(edges_local)
}

get_global_edges <- function(clust_irr) {

    get_gp <- function(gp, chain) {
        return(base::data.frame(from = gp[,1], to = gp[,2],
                                motif = NA, type = "global",
                                chain = chain))
    }

    edges_global <- vector(mode = "list", length = length(clust_irr$clust))
    names(edges_global) <- names(clust_irr$clust)
    for(chain in names(clust_irr$clust)) {
        g <- clust_irr$clust[[chain]]$global
        if(base::is.null(g)==FALSE && base::nrow(g)!=0) {
            edges_global[[chain]] <- get_gp(gp = g, chain = chain)
        }
    }
    edges_global <- base::do.call(base::rbind, edges_global)
    return(edges_global)
}
