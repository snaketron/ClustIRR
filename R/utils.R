# Description:
# Given input data, get the number of chains that ought to be analyzed
# one by one.
# x = columns(data_sample)
get_chains <- function(x) {
    js <- c(
        base::which(base::regexpr(pattern = "CDR3b", text = x) != -1),
        base::which(base::regexpr(pattern = "CDR3a", text = x) != -1),
        base::which(base::regexpr(pattern = "CDR3g", text = x) != -1),
        base::which(base::regexpr(pattern = "CDR3d", text = x) != -1),
        base::which(base::regexpr(pattern = "CDR3h", text = x) != -1),
        base::which(base::regexpr(pattern = "CDR3l", text = x) != -1)
    )
    return(x[js])
}



# Description:
# cut the left/right flanks of each CDR3 sequence by flank_size amino
# acids
get_trimmed_flanks <- function(x,
                               flank_size) {
  
  n <- base::ifelse(test = base::deparse(base::substitute(x)) == "cdr3",
                    yes = "sample",
                    no = "reference")
  t <- flank_size*2
  l <- base::nchar(x)

    if(base::max(l) <= t){
    s <- "trim_flank_aa too high, no sequences left to cluster after trimming"
    base::stop(s)
  }
  
    x <- base::substr(
        x = x,
        start = flank_size + 1,
        stop = base::nchar(x) - flank_size
    )
    
    x[x == ""] <- NA
    
    if(base::any(base::is.na(x))) {
      w <- base::paste0(
        "cdr3 sequences shorter than ",
        t, 
        " (trim_flank_aa*2) of the ",
        n,
        " dataset \n were trimmed completely before local clustering \n \n"
      )
      base::warning(w)
    }
    
    return(x)
}
