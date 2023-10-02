# Description:
# Given input data, get the number of chains that ought to be analyzed
# one by one.
# x = columns(data_sample)
get_chains <- function(x) {
    js <- c(which(regexpr(pattern = "CDR3b", text = x) != -1),
            which(regexpr(pattern = "CDR3a", text = x) != -1),
            which(regexpr(pattern = "CDR3g", text = x) != -1),
            which(regexpr(pattern = "CDR3d", text = x) != -1),
            which(regexpr(pattern = "CDR3h", text = x) != -1),
            which(regexpr(pattern = "CDR3l", text = x) != -1))
    return(x[js])
}



# Description:
# cut the left/right flanks of each CDR3 sequence by flank_size amino
# acids
get_trimmed_flanks <- function(x, flank_size) {
    n <- ifelse(test = deparse(substitute(x)) == "cdr3",
                yes = "sample",
                no = "reference")
    t <- flank_size*2
    l <- nchar(x)
    if(max(l, na.rm = TRUE) <= t) {
        stop("all input CDR3s are shorter than 2 x trim_flank_aa")
    }
    
    x <- substr(x = x, start = flank_size + 1, stop = nchar(x) - flank_size)
    
    x[x == ""] <- NA
    
    if(any(is.na(x))) {
        warning("some input CDR3s are shorter than 2 x trim_flank_aa")
    }
    
    return(x)
}


# This function packs the results into clust_irr object
get_clustirr_output_obj <- function(clust, 
                                    s,
                                    r,
                                    version, 
                                    ks, 
                                    cores, 
                                    control) {
    
    
    return(new("clust_irr", clust = clust,
               inputs = list(s = s, 
                             r = r, 
                             version = version, 
                             ks = ks, 
                             cores = cores,
                             control = control)))
}



# This is the class of the outputs produced by function clust_irr. Object
# from this class are used as input of plot_graph
setClass("clust_irr", representation(clust = "list", inputs = "list"))

# Accessors
setGeneric(name = "get_clustirr_clust", 
           def = function(x) {standardGeneric("get_clustirr_clust")})
setGeneric(name = "get_clustirr_inputs", 
           def = function(x) {standardGeneric("get_clustirr_inputs")})

setMethod(f = "get_clustirr_clust", 
          signature = signature("clust_irr"), 
          definition = function(x) {x@clust})
setMethod(f = "get_clustirr_inputs", 
          signature = signature("clust_irr"), 
          definition = function(x) {x@inputs})

