# Description:
# Given input data, get the number of chains that ought to be analyzed
# one by one.
# x = columns(data_sample)
get_chains <- function(x) {
    js <- c(
        which(regexpr(pattern = "CDR3b", text = x) != -1),
        which(regexpr(pattern = "CDR3a", text = x) != -1)
    )
    return(x[js])
}



# Description:
# cut the left/right flanks of each CDR3 sequence by flank_size amino
# acids
get_trimmed_flanks <- function(x,
                               flank_size) {
    x <- base::substr(
        x = x,
        start = flank_size + 1,
        stop = base::nchar(x) - flank_size
    )
    x[x == ""] <- NA
    return(x)
}
