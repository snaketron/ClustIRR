# Description:
# Now exported to user, but used within function gliph.R
#' Input checks for gliphR parameters
#'
#' @param data_sample, data.frame: TCR sample
#' @param data_ref, data.frame: reference database
#' @param version, integer: version = 1, 2 or 3, gliph version to use
#' @param ks, vector of integers: motif lengths to use (default ks=(2,3,4)),
#' @param cores, integer: number of CPU cores to use
#' @param control, list: auxiliary input parameters
#' @noRd
parameter_check <- function(data_sample,
                            data_ref,
                            version,
                            ks,
                            cores,
                            control) {
    check_data_sample(data_sample)
    check_data_ref(dat_ref)
    check_data_sample_and_ref(data_sample, data_ref)
    check_version(version)
    check_ks(ks)
    check_local_min_ove(control$local_min_ove)
    # we do not need this anymore
    # check_ks_and_local_min_ove(ks, control$local_min_ove)
    check_cores(cores)
    check_B(control$B)
    check_global_max_dist(control$global_max_dist)
    check_local_min_p(control$local_min_p)
    #check_local_min_o(local_min_o)
    #check_trim_flanks(trim_flanks)
    #check_flank_size(flank_size)
    #check_trim_flanks_flank_size(trim_flanks,flank_size)
    #check_global_pairs(global_pairs, data_sample)

}

check_data_sample <- function(data_sample) {
    if(base::missing(data_sample)){
        base::stop("missing input parameter data_sample")
        }
    if(!base::is.data.frame(data_sample)){
        base::stop("data_sample has to be of type data frame")
        }
    if(base::nrow(data_sample) == 0){
        base::stop("data_sample contains zero rows")
        }
    if(!base::any(base::colnames(data_sample) %in% base::c("CDR3a", "CDR3b"))){
        base::stop("data_sample is missing a 'CDR3a' or 'CDR3b' column")
        }
    if(base::any(base::colnames(data_sample) %in% "CDR3a") &&
       !base::is.character(data_sample$CDR3a)){
        base::stop("CDR3a column has to of type character")
        }
    if(base::any(base::colnames(data_sample) %in% "CDR3b") &&
       !base::is.character(data_sample$CDR3b)){
        base::stop("CDR3b column has to of type character")
        }
    if(base::any(base::is.na(data_sample))){
        base::warning("data_sample contains NA values")
        }
    if(base::any(data_sample=="")){
        base::warning("data_sample contains empty values")
        }
    # checks for optional v-gene / j-gene columns tba if necessary
}

check_data_ref <- function(data_ref) {
    if(base::missing(data_ref)){
        base::stop("missing input parameter data_ref")
        }
    if(!base::is.data.frame(data_ref) &&
       !(data_ref %in% base::c("gliph_reference"))){
        base::stop("data_ref has to be of type data frame or 'gliph_reference'")
        }
    if(base::nrow(data_ref) == 0){
        base::stop("data_ref contains zero rows")
        }
    if(!base::any(base::colnames(data_ref) %in% base::c("CDR3a", "CDR3b"))){
        base::stop("data_ref is missing a 'CDR3a' or 'CDR3b' column")
        }
    if(base::any(base::colnames(data_ref) %in% "CDR3a") &&
       !base::is.character(data_ref$CDR3a)){
        base::stop("CDR3a column has to of type character")
        }
    if(base::any(base::colnames(data_ref) %in% "CDR3b") &&
       !base::is.character(data_ref$CDR3b)){
        base::stop("CDR3b column has to of type character")
        }
    if(base::any(base::is.na(data_ref))){
        base::warning("data_ref contains NA values")
        }
    if(base::any(data_ref=="")){
        base::warning("data_ref contains empty values")
        }
    # checks for optional v-gene / j-gene columns tba if necessary
}

check_data_sample_and_ref <- function(data_sample, data_ref) {
    if(base::any(base::colnames(data_sample) %in% base::c("CDR3a")) &&
       !base::any(base::colnames(data_ref) %in% base::c("CDR3a"))){
        base::stop("data_sample contains CDR3a column, but data_ref does not")
        }
    if(base::any(base::colnames(data_ref) %in% base::c("CDR3a")) &&
       !base::any(base::colnames(data_sample) %in% base::c("CDR3a"))){
        base::stop("data_ref contains CDR3a column, but data_sample does not")
        }
    if(base::any(base::colnames(data_sample) %in% base::c("CDR3b")) &&
       !base::any(base::colnames(data_ref) %in% base::c("CDR3b"))){
        base::stop("data_sample contains CDR3b column, but data_ref does not")
        }
    if(base::any(base::colnames(data_ref) %in% base::c("CDR3b")) &&
       !base::any(base::colnames(data_sample) %in% base::c("CDR3b"))){
        base::stop("data_ref contains CDR3b column, but data_sample does not")
        }
}

check_version <- function(version) {
    # not necessary because of default value
    # if(base::missing(version)){
    #     base::stop("missing input parameter version")
    #     }
    if(!base::is.numeric(version)){
        base::stop("version has to be numeric")
        }
    if(base::length(version)!=1){
        base::stop("version has to be a single number")
        }
    if(!(version %in% base::c(1, 2, 3))){
        base::stop("version has to be 1, 2 or 3")
        }
}

check_ks <- function(ks) {
    if(!base::is.numeric(ks)){
        base::stop("ks has to be numeric")
    }
    if(!is_wholenumber(ks)){
        base::stop("ks has to be a whole number")
    }
    if(base::any(ks < 1)){
        base::stop("ks must be at least 1")
    }
    if(base::any(is.infinite(ks))) {
        base::stop("ks must be an integer > 1")
    }
}

check_local_min_ove <- function(local_min_ove) {
    if(!base::is.numeric(local_min_ove)){
        base::stop("local_min_ove has to be numeric")
    }
    if(!is_wholenumber(local_min_ove)){
        base::stop("local_min_ove has to be a whole number")
    }
    # while I agree with the next line, I would still give the user full
    # autonomy on filtering thresholds.
    # if(local_min_ove < 1){
    #     base::stop("local_min_ove must be at least 1")
    # }
    if(base::is.infinite(local_min_ove)){
        base::stop("local_min_ove must be an integer")
    }
}

# I have changed local_min_ove to be a single number, and not a vector
# associated with the different ks -> the check is thus not necessary
# check_ks_and_local_min_ove <- function(ks, local_min_ove) {
#     if(base::length(local_min_ove) > 1 &&
#        base::length(local_min_ove) != base::length(ks)){
#         base::stop("local_min_ove has to be a single number
#                     or the same length as ks")
#     }
#
# }

check_cores <- function(cores) {
    if(!base::is.numeric(cores)){
        base::stop("cores has to be numeric")
    }
    if(!is_wholenumber(cores)){
        base::stop("cores has to be a whole number")
    }
    if(base::any(cores < 1)){
        base::stop("cores must be at least 1")
    }
}

check_B <- function(B) { # is B == simulation_depth? yes
    if(!base::is.numeric(B)){
        base::stop("B has to be numeric")
    }
    if(!is_wholenumber(B)){
        base::stop("B has to be a whole number")
    }
    if(base::any(B < 1)){
        base::stop("B must be at least 1")
    }
    if(base::any(B > 10000)){ # to be tested. SK: I think this was an issue with Jan's code, this is fixed now.
        base::warning("B > 10000 can lead to unstable behaviour")
    }
}

check_global_max_dist <- function(global_max_dist) {
    if(!base::is.numeric(global_max_dist)){
        base::stop("global_max_dist has to be numeric")
    }
    if(!is_wholenumber(global_max_dist)){
        base::stop("global_max_dist has to be a whole number")
    }
    if(base::any(global_max_dist < 1)){
        base::stop("global_max_dist must be at least 1")
    }
}

check_local_min_p <- function(local_min_p) {
    if(!base::is.numeric(local_min_p)){
        base::stop("local_min_p has to be numeric")
    }
    if(base::any(local_min_p <= 0)){
        base::stop("local_min_p must be > 0")
    }
    if(base::any(local_min_p > 1)){
        base::stop("local_min_p must be <= 1")
    }
}

check_local_min_o <- function(local_min_o) {
    # which parameter of turboGliph/gliph2 is this?
    # kmer_mindepth
}

check_trim_flanks <- function(trim_flanks) {
    # which parameter of turboGliph/gliph2 is this?
    # structboundaries
}

check_flank_size <- function(flank_size) {
    # which parameter of turboGliph/gliph2 is this?
    # boundary_size
}

check_trim_flanks_flank_size <- function(trim_flanks,
                                         flank_size) {

}

check_global_pairs <- function(global_pairs, data_sample) {
    # global_pairs must be integer matrix with two columns and u rows.
    # In each entry global_pairs will store an index (integer) i = 1, ..., n,
    # pointing to a specific CDR3 sequence from data_sample which is a
    # data.frame with nrow = n.
    # row in global pairs gives us the indices of two CDR3s that are globally
    # similar (e.g. Hamming dist < 1 or distance computed using external tool)
}

is_wholenumber <- function(x) {
    return(abs(x - round(x)) < .Machine$double.eps^0.5)
}

