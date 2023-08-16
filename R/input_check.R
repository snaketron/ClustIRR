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
    check_ks_and_trim_flank_aa(ks, control$trim_flank_aa, s, r)
}

check_s <- function(s) {
    check_missing(s)
    check_dataframe(s)
    check_rowcount(s)
    check_dataframe_colnames(s)
    check_dataframe_na(s)
    check_dataframe_empty(s)
    check_aa(s)
}

check_r <- function(r) {
    check_missing(r)
    check_dataframe(r)
    check_rowcount(r)
    check_dataframe_colnames(r)
    check_dataframe_na(r)
    check_dataframe_empty(r)
    check_aa(r)
}

check_s_and_r <- function(s, r) {
    if (!all(sort(colnames(s)) ==
        sort(colnames(r)))) {
        stop("s has to contain the same columns as r")
    }
}

check_version <- function(version) {
    check_numeric(version)
    check_singlevalue(version)
    if (!(version %in% c(1, 2, 3))) {
        stop("version has to be 1, 2 or 3")
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
    if (!is.null(global_pairs)) {
        check_rowcount(global_pairs)
        check_matrix(global_pairs)
        check_matrix_type(global_pairs, type = "character")
        check_matrix_column_count(global_pairs, 3)
        if (all(global_pairs[, c(1, 2)] %in% s) == FALSE) {
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
    control <- list(
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
    if (missing(control_in) || is.null(control_in)) {
        return(control)
    }
    if (is.list(control_in) == FALSE) {
        stop("control must be a list")
    }
    if (all(names(control_in) %in% names(control)) == FALSE) {
        stop("unrecognized elements found in control")
    }

    ns <- names(control_in)
    for (i in seq_len(length(control_in))) {
        control[[ns[i]]] <- control_in[[ns[i]]]
    }
    return(control)
}

check_ks_and_trim_flank_aa <- function(ks, trim, s, r){
    k_max <- max(ks)
    trim <- trim*2
    s_max <- max(
        apply(s, 2, function(x) max(nchar(x))))
    r_max <- max(
        apply(r, 2, function(x) max(nchar(x))))
    rs_minmax <- min(s_max, r_max, na.rm = TRUE)

    if((rs_minmax - trim) < k_max){
        stop("ks has to be smaller than the biggest trimmed sequence")
    }

}

#### Helper functions ####

check_aa <- function(x) {
    aa <- "[^ACDEFGHIKLMNPQRSTVWY]"
    w <- paste0(deparse(substitute(x)),
                      " contains non-standard or lowercase amino acid codes")
    res <- lapply(x, function(y) length(grep(aa, y)) != 0)
    if (any(unlist(res))) {
        stop(w)
    }

}

check_dataframe <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be of type data frame"
    )
    if (!is.data.frame(x)) {
        stop(w)
    }
}

check_dataframe_colnames <- function(x) {
    c <- c("CDR3a", "CDR3b", "CDR3d", "CDR3g", "CDR3l", "CDR3h")
    w <- paste0(
        paste0(
            deparse(substitute(x)),
            " has to contain at least one of the following columns: "
        ),
        "CDR3a and/or CDR3b or CDR3d and/or CDR3g or CDR3h and/or CDR3l"
    )
    if (!any(colnames(x) %in% c)) {
        stop(w)
    }
    for (n in colnames(x)) {
        if (!is.character(x[[n]])) {
            s <- paste0(n, " column has to of type character")
            stop(s)
        }
    }
    if (any(colnames(x) %in% c("CDR3a", "CDR3b")) &
        any(colnames(x) %in% c("CDR3d", "CDR3g"))) {
        stop("CDR3a/b can't be mixed with CDR3d/g columns")
    }
    if (any(colnames(x) %in% c("CDR3a", "CDR3b")) &
        any(colnames(x) %in% c("CDR3h", "CDR3l"))) {
        stop("CDR3a/b can't be mixed with CDR3l/h columns")
    }
    if (any(colnames(x) %in% c("CDR3d", "CDR3g")) &
        any(colnames(x) %in% c("CDR3h", "CDR3l"))) {
        stop("CDR3d/g can't be mixed with CDR3l/h columns")
    }
}

check_dataframe_empty <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " contains empty values"
    )
    if (any(x == "", na.rm = TRUE)) {
        warning(w)
    }
}

check_dataframe_na <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " contains NA value"
    )
    if (any(is.na(x))) {
        warning(w)
    }
}

check_greaterthan <- function(x, v) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be <= ",
        v
    )
    if (any(x > v)) {
        stop(w)
    }
}

check_infinity <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be a finite number"
    )
    if (any(is.infinite(x))) {
        stop(w)
    }
}

check_lessthan <- function(x, v) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be >= ",
        v
    )
    if (any(x < v)) {
        stop(w)
    }
}

check_logical <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be logical"
    )
    if (any(is.na(x))) {
        stop(w)
    }
    if (!is.logical(x)) {
        stop(w)
    }
}

check_matrix <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be of type matrix"
    )
    if (!is.matrix(x)) {
        stop(w)
    }
}

check_matrix_column_count <- function(x, c) {
    w <- paste0(
        deparse(substitute(x)),
        " has to have ", c, " columns"
    )
    if (ncol(x) != c) {
        stop(w)
    }
}

check_matrix_type <- function(x, type) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be a numeric matrix"
    )

    if (type == "numeric") {
        if (is.numeric(x) == FALSE) {
            stop(w)
        }
    }
    if (type == "character") {
        if (is.character(x) == FALSE) {
            stop(w)
        }
    }
}

check_missing <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " parameter is missing"
    )
    if (missing(x) || is.null(x)) {
        stop(w)
    }
}

check_numeric <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be numeric"
    )
    if (any(!is.numeric(x))) {
        stop(w)
    }
}

check_rowcount <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " contains zero rows"
    )
    if (nrow(x) == 0) {
        stop(w)
    }
}

check_singlevalue <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be a single value"
    )
    if (any(is.na(x))) {
        stop(w)
    }
    if (length(x) != 1) {
        stop(w)
    }
}

check_wholenumber <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be a whole number"
    )
    if (any(!(abs(x - round(x)) < .Machine$double.eps^0.5))) {
        stop(w)
    }
}

check_positive <- function(x) {
    w <- paste0(
        deparse(substitute(x)),
        " has to be positive number"
    )
    if (x < 0) {
        stop(w)
    }
}

