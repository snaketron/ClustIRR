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
