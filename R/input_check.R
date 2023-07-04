# Description:
# Check user provided input and generate errors and warnings, if necessary
input_check <- function(s,
                        r,
                        version,
                        ks,
                        cores,
                        control) {
    check_s(s)
    check_r(r)
    check_s_and_r(s, r)
    check_version(version)
    check_ks(ks)
    check_local_min_ove(control$local_min_ove)
    check_cores(cores)
    check_B(control$B)
    check_global_max_dist(control$global_max_dist)
    check_local_max_fdr(control$local_max_fdr)
    check_local_min_o(control$local_min_o)
    check_trim_flank_aa(control$trim_flank_aa)
    check_low_mem(control$low_mem)
    check_global_pairs(control$global_pairs, s)
}

check_s <- function(s) {
    check_missing(s)
    check_dataframe(s)
    check_rowcount(s)
    check_dataframe_colnames(s, base::c("CDR3a", "CDR3b"))
    if(base::any(base::colnames(s) %in% "CDR3a") &&
        !base::is.character(s$CDR3a)) {
        base::stop("CDR3a column has to of type character")
    }
    if(base::any(base::colnames(s) %in% "CDR3b") &&
        !base::is.character(s$CDR3b)) {
        base::stop("CDR3b column has to of type character")
    }
    check_dataframe_na(s)
    check_dataframe_empty(s)
}

check_r <- function(r) {
    check_missing(r)
    check_dataframe(r)
    check_rowcount(r)
    check_dataframe_colnames(r, base::c("CDR3a", "CDR3b"))
    if(base::any(base::colnames(r) %in% "CDR3a") &&
        !base::is.character(r$CDR3a)) {
        base::stop("CDR3a column has to of type character")
    }
    if(base::any(base::colnames(r) %in% "CDR3b") &&
        !base::is.character(r$CDR3b)) {
        base::stop("CDR3b column has to of type character")
    }
    check_dataframe_na(r)
    check_dataframe_empty(r)
}

check_s_and_r <- function(s, r) {
    if(base::any(base::colnames(s) %in% base::c("CDR3a")) &&
        !base::any(base::colnames(r) %in% base::c("CDR3a"))) {
        base::stop("s contains CDR3a column, but r does not")
    }
    if(base::any(base::colnames(r) %in% base::c("CDR3a")) &&
        !base::any(base::colnames(s) %in% base::c("CDR3a"))) {
        base::stop("r contains CDR3a column, but s does not")
    }
    if(base::any(base::colnames(s) %in% base::c("CDR3b")) &&
        !base::any(base::colnames(r) %in% base::c("CDR3b"))) {
        base::stop("s contains CDR3b column, but r does not")
    }
    if(base::any(base::colnames(r) %in% base::c("CDR3b")) &&
        !base::any(base::colnames(s) %in% base::c("CDR3b"))) {
        base::stop("r contains CDR3b column, but s does not")
    }
}

check_version <- function(version) {
    check_numeric(version)
    check_singlevalue(version)
    if(!(version %in% base::c(1, 2, 3))) {
        base::stop("version has to be 1, 2 or 3")
    }
}

check_ks <- function(ks) {
    check_infinity(ks)
    check_numeric(ks)
    check_wholenumber(ks)
    check_lessthan(ks, 1)
}

check_local_min_ove <- function(local_min_ove) {
    check_infinity(local_min_ove)
    check_numeric(local_min_ove)
    check_singlevalue(local_min_ove)
}

check_cores <- function(cores) {
    check_infinity(cores)
    check_numeric(cores)
    check_wholenumber(cores)
    check_singlevalue(cores)
    check_lessthan(cores, 1)
}

check_B <- function(B) {
    check_infinity(B)
    check_numeric(B)
    check_wholenumber(B)
    check_singlevalue(B)
    check_lessthan(B, 1)
}

check_global_max_dist <- function(global_max_dist) {
    check_infinity(global_max_dist)
    check_numeric(global_max_dist)
    check_wholenumber(global_max_dist)
    check_singlevalue(global_max_dist)
    check_lessthan(global_max_dist, 1)
}

check_local_max_fdr <- function(local_max_fdr) {
    check_infinity(local_max_fdr)
    check_numeric(local_max_fdr)
    check_singlevalue(local_max_fdr)
    check_lessthan(local_max_fdr, 0)
    check_greaterthan(local_max_fdr, 1)
}

check_local_min_o <- function(local_min_o) { # kmer_mindepth
    check_infinity(local_min_o)
    check_numeric(local_min_o)
    check_wholenumber(local_min_o)
    check_singlevalue(local_min_o)
}

check_trim_flank_aa <- function(trim_flank_aa) { # boundary_size
    check_singlevalue(trim_flank_aa)
    check_infinity(trim_flank_aa)
    check_numeric(trim_flank_aa)
    check_wholenumber(trim_flank_aa)
    check_positive(trim_flank_aa)
}

check_global_pairs <- function(global_pairs, s) {
    if(!base::is.null(global_pairs)) {
        check_rowcount(global_pairs)
        check_matrix(global_pairs)
        check_matrix_type(global_pairs, type = "character")
        check_matrix_column_count(global_pairs, 3)
        if(all(global_pairs[, c(1,2)] %in% s)==FALSE) {
            stop("not all CDR3s from global_pair are found in s")
        }
    }
}

check_low_mem <- function(low_mem) {
    check_singlevalue(x = low_mem)
    check_logical(x = low_mem)
}

# Description:
# Setup control list.
# control_in: user generated list (if missing -> use default)
get_control <- function(control_in) {
    control <- base::list(
        B = 1000,
        global_max_dist = 1,
        local_max_fdr = 0.05,
        local_min_ove = 2,
        local_min_o = 1,
        trim_flank_aa = 0,
        global_pairs = NULL,
        low_mem = FALSE
    )

    # if missing control_in -> use default values
    if(base::missing(control_in) || base::is.null(control_in)) {
        return(control)
    }
    if(base::is.list(control_in) == FALSE) {
        base::stop("control must be a list")
    }
    if(base::all(base::names(control_in) %in% base::names(control)) == FALSE){
        base::stop("unrecognized elements found in control")
    }

    ns <- names(control_in)
    for(i in seq_len(length(control_in))) {
        control[[ns[i]]] <- control_in[[ns[i]]]
    }
    return(control)
}


#### Helper functions ####

check_dataframe <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be of type data frame"
    )
    if(!base::is.data.frame(x)) {
        base::stop(w)
    }
}

check_dataframe_colnames <- function(x, c) {
    w <- base::paste0(
        base::paste0(
            base::deparse(base::substitute(x)),
            " has to contain the following columns: "
        ),
        base::paste0(c, collapse = " or ")
    )
    if(!base::any(base::colnames(x) %in% c)) {
        base::stop(w)
    }
}

check_dataframe_empty <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " contains empty values"
    )
    if(base::any(x == "", na.rm = TRUE)) {
        base::warning(w)
    }
}

check_dataframe_na <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " contains NA value"
    )
    if(base::any(base::is.na(x))) {
        base::warning(w)
    }
}

check_greaterthan <- function(x, v) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be <= ",
        v
    )
    if(base::any(x > v)) {
        base::stop(w)
    }
}

check_infinity <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be a finite number"
    )
    if(base::any(base::is.infinite(x))) {
        base::stop(w)
    }
}

check_lessthan <- function(x, v) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be >= ",
        v
    )
    if(base::any(x < v)) {
        base::stop(w)
    }
}

check_logical <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be logical"
    )
    if(base::any(base::is.na(x))) {
        base::stop(w)
    }
    if(!base::is.logical(x)) {
        base::stop(w)
    }
}

check_matrix <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be of type matrix"
    )
    if(!base::is.matrix(x)) {
        base::stop(w)
    }
}

check_matrix_column_count <- function(x, c) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to have ", c, " columns"
    )
    if(base::ncol(x) != c) {
        base::stop(w)
    }
}

check_matrix_type <- function(x, type) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be a numeric matrix"
    )

    if(type=="numeric") {
        if(is.numeric(x)==FALSE) {
            base::stop(w)
        }
    }
    if(type=="character") {
        if(is.character(x)==FALSE) {
            base::stop(w)
        }
    }
}

check_missing <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " parameter is missing"
    )
    if(base::missing(x) || base::is.null(x)) {
        base::stop(w)
    }
}

check_numeric <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be numeric"
    )
    if(base::any(!base::is.numeric(x))) {
        base::stop(w)
    }
}

check_rowcount <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " contains zero rows"
    )
    if(base::nrow(x) == 0) {
        base::stop(w)
    }
}

check_singlevalue <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be a single value"
    )
    if(base::any(base::is.na(x))) {
        base::stop(w)
    }
    if(base::length(x) != 1) {
        base::stop(w)
    }
}

check_wholenumber <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be a whole number"
    )
    if(base::any(!(abs(x - round(x)) < .Machine$double.eps^0.5))) {
        base::stop(w)
    }
}

check_positive <- function(x) {
    w <- base::paste0(
        base::deparse(base::substitute(x)),
        " has to be positive number"
    )
    if(x < 0) {
        base::stop(w)
    }
}
