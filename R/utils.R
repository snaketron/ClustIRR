# Description:
# Given input data, get the number of chains that ought to be analyzed
# one by one.
# x = columns(data_sample)
get_chains <- function(x) {
    js <- c(
        which(regexpr(pattern = "CDR3b", text = x) != -1),
        which(regexpr(pattern = "CDR3a", text = x) != -1),
        which(regexpr(pattern = "CDR3g", text = x) != -1),
        which(regexpr(pattern = "CDR3d", text = x) != -1),
        which(regexpr(pattern = "CDR3h", text = x) != -1),
        which(regexpr(pattern = "CDR3l", text = x) != -1)
    )
    return(x[js])
}



# Description:
# cut the left/right flanks of each CDR3 sequence by flank_size amino
# acids
get_trimmed_flanks <- function(x,flank_size) {
    n <- ifelse(test = deparse(substitute(x)) == "cdr3",
        yes = "sample",
        no = "reference"
    )
    t <- flank_size * 2
    l <- nchar(x)
    s <- "trim_flank_aa too high, no sequences left to cluster after trimming"
    if (max(l, na.rm = TRUE) <= t) {
        stop(s)
    }

    x <- substr(
        x = x,
        start = flank_size + 1,
        stop = nchar(x) - flank_size
    )

    x[x == ""] <- NA

    if (any(is.na(x))) {
        w <- paste0(
            "cdr3 sequences shorter than ",
            t,
            " (trim_flank_aa*2) of the ",
            n,
            " dataset \n were trimmed completely before local clustering \n \n"
        )
        warning(w)
    }

    return(x)
}
