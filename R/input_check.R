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
input_check <- function(data_sample,
                        data_ref,
                        version,
                        ks,
                        cores,
                        control) {
    check_data_sample(data_sample)
    check_data_ref(data_ref)
    check_data_sample_and_ref(data_sample, data_ref)
    check_version(version)
    check_ks(ks)
    check_local_min_ove(control$local_min_ove)
    check_cores(cores)
    check_B(control$B)
    check_global_max_dist(control$global_max_dist)
    check_local_max_fdr(control$local_max_fdr)
    check_local_min_o(control$local_min_o)
    check_trim_flanks(control$trim_flanks)
    check_flank_size(control$flank_size)
    check_low_mem(control$low_mem)
    check_trim_flanks_flank_size(control$trim_flanks,control$flank_size)
    check_global_pairs(control$global_pairs, data_sample)
}

check_data_sample <- function(data_sample) {
    check_missing(data_sample)
    check_dataframe(data_sample)
    check_rowcount(data_sample)
    check_dataframe_colnames(data_sample, base::c("CDR3a", "CDR3b"))
    if (base::any(base::colnames(data_sample) %in% "CDR3a") &&
        !base::is.character(data_sample$CDR3a)) {
        base::stop("CDR3a column has to of type character")
    }
    if (base::any(base::colnames(data_sample) %in% "CDR3b") &&
        !base::is.character(data_sample$CDR3b)) {
        base::stop("CDR3b column has to of type character")
    }
    check_dataframe_na(data_sample)
    check_dataframe_empty(data_sample)
    # checks for optional v-gene / j-gene columns tba if necessary
}

check_data_ref <- function(data_ref) {
    if (!base::deparse(base::substitute(data_ref)) == "data_ref") {
        check_missing(data_ref)
        check_dataframe(data_ref)
        check_rowcount(data_ref)
        check_dataframe_colnames(data_ref, base::c("CDR3a", "CDR3b"))
        if (base::any(base::colnames(data_ref) %in% "CDR3a") &&
            !base::is.character(data_ref$CDR3a)) {
            base::stop("CDR3a column has to of type character")
        }
        if (base::any(base::colnames(data_ref) %in% "CDR3b") &&
            !base::is.character(data_ref$CDR3b)) {
            base::stop("CDR3b column has to of type character")
        }
        check_dataframe_na(data_ref)
        check_dataframe_empty(data_ref)
    }
    # checks for optional v-gene / j-gene columns tba if necessary
}

check_data_sample_and_ref <- function(data_sample, data_ref) {
    if (base::any(base::colnames(data_sample) %in% base::c("CDR3a")) &&
        !base::any(base::colnames(data_ref) %in% base::c("CDR3a"))) {
        base::stop("data_sample contains CDR3a column, but data_ref does not")
    }
    if (base::any(base::colnames(data_ref) %in% base::c("CDR3a")) &&
        !base::any(base::colnames(data_sample) %in% base::c("CDR3a"))) {
        base::stop("data_ref contains CDR3a column, but data_sample does not")
    }
    if (base::any(base::colnames(data_sample) %in% base::c("CDR3b")) &&
        !base::any(base::colnames(data_ref) %in% base::c("CDR3b"))) {
        base::stop("data_sample contains CDR3b column, but data_ref does not")
    }
    if (base::any(base::colnames(data_ref) %in% base::c("CDR3b")) &&
        !base::any(base::colnames(data_sample) %in% base::c("CDR3b"))) {
        base::stop("data_ref contains CDR3b column, but data_sample does not")
    }
}

check_version <- function(version) {
    check_numeric(version)
    check_singlevalue(version)
    if (!(version %in% base::c(1, 2, 3))) {
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
    #check_wholenumber(local_min_ove) # SK: min_ove can be a real number
    check_singlevalue(cores) # SK: contrary to Jan's implementation I
    # programmed min_ove as a single number (threshold)
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
    check_lessthan(global_max_dist,1)
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

check_trim_flanks <- function(trim_flanks) { # structboundaries
    check_logical(trim_flanks)
    check_singlevalue(trim_flanks)
}

check_flank_size <- function(flank_size) { # boundary_size
    check_infinity(flank_size)
    check_numeric(flank_size)
    check_wholenumber(flank_size)
    check_singlevalue(flank_size)
}

check_trim_flanks_flank_size <- function(trim_flanks,
                                         flank_size) {
    # When testing the inputs you could check what happens downstream if the
    # user sets trim_flanks = TRUE & flank_size = 3 (that is 6 amino acids
    # are supposed to be trimmed from both flanks of the CDR3), but a given
    # CDR3 sequence is 4 amino acids long.


    # SK: this is the behavior we want. We should not produce such a warning
    # if(trim_flanks){
    #     if(flank_size > 2){
    #         base::warning(base::paste0("flank_size is set to ",
    #                              flank_size,
    #                              ". Bigger CRD3 sequences could potentially ",
    #                              "get cut off at flanks"))
    #     }
    # }
}

check_global_pairs <- function(global_pairs, data_sample) {
    # global_pairs has to be integer matrix with two columns and u rows.
    # In each entry global_pairs will store an index (integer) i = 1, ..., n,
    # pointing to a specific CDR3 sequence from data_sample which is a
    # data.frame with nrow = n.
    # row in global pairs gives us the indices of two CDR3s that are globally
    # similar (e.g. Hamming dist < 1 or distance computed using external tool)

    # global_pairs -> optional user provided input
    # (e.g. by smarter global clustering)
    # there is no equivalent parameter in gliph/gliph2
    if (!base::is.null(global_pairs)) {
        check_rowcount(global_pairs)
        check_matrix(global_pairs)
        check_matrix_type(global_pairs)
        check_matrix_column_count(global_pairs, 2)
        # TODO max index should not exceed data sample
    }
}

check_low_mem <- function(low_mem) {
    check_missing(x = low_mem)
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
        local_min_o = 3,
        trim_flanks = FALSE,
        flank_size = 3,
        global_pairs = NULL,
        low_mem = FALSE)

    # if missing control_in -> use default values
    if(base::missing(control_in)|base::is.null(control_in)) {
        return(control)
    }
    if(base::is.list(control_in)==FALSE) {
        base::stop("control must be a list")
    }
    if(base::is.data.frame(control_in)==TRUE) {
        base::stop("control must be a list")
    }
    if(base::all(base::names(control_in) %in% base::names(control))==FALSE) {
        base::stop("unrecognized elements found in control")
    }

    # edit control by user-defined control_in
    # There are packages to update list given another list, but here we can
    # live with the following "inefficiency" as the list is generally small
    # (~5 elements)
    ns <- names(control_in)
    for(i in 1:length(control_in)) {
        control[[ns[i]]] <- control_in[[ns[i]]]
    }
    return(control)
}


#Helper functions---------------------------------------------------------------

check_dataframe <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be of type data frame")
    if (!base::is.data.frame(x)) {
        base::stop(w)
    }
}

check_dataframe_colnames <- function(x, c){
    w <- base::paste0(base::paste0(base::deparse(base::substitute(x)),
                          " has to contain the following column(s): "),
                      base::paste0(c, collapse = " or "))
    if (!base::any(base::colnames(x) %in% c)) {
        base::stop(w)
    }
}

check_dataframe_empty <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " contains empty values")
    if (base::any(x == "", na.rm = TRUE)) {
        base::warning(w)
    }
}

check_dataframe_na <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " contains NA value")
    if (base::any(base::is.na(x))) {
        base::warning(w)
    }
}

check_greaterthan <- function(x, v){
    w <- base::paste0(base::deparse(base::substitute(x)),
                 " has to be <= ",
                 v)
    if (base::any(x > v)) {
        base::stop(w)
    }
}

check_infinity <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be a finite number")
    if (base::any(base::is.infinite(x))) {
        base::stop(w)
    }
}

check_lessthan <- function(x, v){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be >= ",
                      v)
    if (base::any(x < v)) {
        base::stop(w)
    }
}

check_logical <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be logical")
    if (base::any(base::is.na(x))) {
        base::stop(w)
    }
    if (!base::is.logical(x)) {
        base::stop(w)
    }
}

check_matrix <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be of type matrix")
    if (!base::is.matrix(x)) {
        base::stop(w)
    }
}

check_matrix_column_count <- function(x, c){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be have ", c, " columns")
    if (base::ncol(x) != c) {
        base::stop(w)
    }
}

check_matrix_type <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be an integer type matrix")
    if (!base::typeof(x) == "integer") {
        base::stop(w)
    }
}

check_missing <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " parameter is missing")
    if (base::missing(x)) {
        base::stop(w)
    }
}

check_numeric <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be numeric")
    if (base::any(!base::is.numeric(x))) {
        base::stop(w)
    }
}

check_rowcount <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " contains zero rows")
    if (base::nrow(x) == 0) {
        base::stop(w)
    }
}

check_singlevalue <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be a single value")
    if (base::is.na(x)) {
        base::stop(w)
    }
    if (base::length(x) != 1) {
        base::stop(w)
    }
}

check_wholenumber <- function(x){
    w <- base::paste0(base::deparse(base::substitute(x)),
                      " has to be a whole number")
    if (base::any(!(abs(x - round(x)) < .Machine$double.eps^0.5))) {
        base::stop(w)
    }
}

