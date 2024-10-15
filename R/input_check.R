# Description:
# Check user provided input and generate errors and warnings, if necessary
input_check <- function(s, control) {
  
  check_s(s = s)
  check_global_max_hdist(control$global_max_hdist)
  check_trim_flank_aa(control$trim_flank_aa)
  check_low_mem(control$low_mem)
  check_global_hamming(control$global_hamming)
  check_gmi(control$gmi)
}

check_s <- function(s) {
  check_missing(s)
  check_dataframe(s)
  check_s_cols(s)
  check_rowcount(s)
  check_dataframe_na(s)
  check_dataframe_empty(s)
  check_aa(s)
}

check_cores <- function(cores) {
  check_infinity(cores)
  check_numeric(cores)
  check_wholenumber(cores)
  check_singlevalue(cores)
  check_lessthan(cores, 1)
}

check_global_max_hdist <- function(global_max_hdist) {
  check_infinity(global_max_hdist)
  check_numeric(global_max_hdist)
  check_wholenumber(global_max_hdist)
  check_singlevalue(global_max_hdist)
  check_lessthan(global_max_hdist, 1)
}

check_trim_flank_aa <- function(trim_flank_aa) { # boundary_size
  check_singlevalue(trim_flank_aa)
  check_infinity(trim_flank_aa)
  check_numeric(trim_flank_aa)
  check_wholenumber(trim_flank_aa)
  check_positive(trim_flank_aa)
}

check_low_mem <- function(low_mem) {
  check_singlevalue(x = low_mem)
  check_logical(x = low_mem)
}

check_global_hamming <- function(global_hamming) {
  check_singlevalue(x = global_hamming)
  check_logical(x = global_hamming)
}

check_as_visnet <- function(as_visnet) {
  if(missing(as_visnet)) {
    stop("as_visnet missing")
  }
  check_singlevalue(x = as_visnet)
  check_logical(x = as_visnet)
}

check_show_singletons <- function(show_singletons) {
  if(missing(show_singletons)) {
    stop("show_singletons missing")
  }
  check_singlevalue(x = show_singletons)
  check_logical(x = show_singletons)
}

check_select_by <- function(select_by) {
  if(missing(select_by)) {
    stop("select_by is missing")
  }
  check_singlevalue(x = select_by)
  if(is.null(select_by) || is.na(select_by)) {
    stop("select_by is NULL or NA")
  }
  if(is.character(select_by)==FALSE) {
    stop("type of select_by must be character")
  }
  if(select_by %in% c("Ag_species", "Ag_gene") == FALSE) {
    stop("select_by can be either Ag_species or Ag_gene")
  }
}

check_node_opacity <- function(node_opacity) {
    check_infinity(node_opacity)
    check_numeric(node_opacity)
    check_singlevalue(node_opacity)
    check_probability(node_opacity)
}

check_gmi <- function(gmi) {
  check_infinity(gmi)
  check_numeric(gmi)
  check_singlevalue(gmi)
  check_probability(gmi)
}

check_edit_dist <- function(edit_dist) {
    check_infinity(edit_dist)
    check_numeric(edit_dist)
    check_wholenumber(edit_dist)
    check_singlevalue(edit_dist)
    check_lessthan(edit_dist, 0)
}

# Description:
# Setup control list.
# control_in: user generated list (if missing -> use default)
get_control <- function(control_in) {
  control = list(gmi = 0.7,
                 trim_flank_aa = 3,
                 global_hamming = FALSE,
                 global_max_hdist = 1,
                 low_mem = FALSE)
  
  # if missing control_in -> use default values
  if(missing(control_in) || is.null(control_in)) {
    return(control)
  }
  if(is.list(control_in) == FALSE) {
    stop("control must be a list")
  }
  if(all(names(control_in) %in% names(control)) == FALSE) {
    stop("unrecognized elements found in control")
  }
  
  ns <- names(control_in)
  for (i in seq_len(length(control_in))) {
    control[[ns[i]]] <- control_in[[ns[i]]]
  }
  return(control)
}


check_clustirr <- function(clust_irr) {
  # if missing control_in -> use default values
  if(missing(clust_irr) || is.null(clust_irr)) {
    stop("input clust_irr is empty")
  }
  
  if(is(clust_irr, class2 = "clust_irr")==FALSE) {
    stop("input clust_irr is not class clust_irr")
  }
}


check_custom_db <- function(custom_db) {
  if(missing(custom_db)) {
    return(NULL)
  }
  if(is.null(custom_db)) {
    return(NULL)
  }
  
  if(is.data.frame(custom_db)==FALSE) {
    stop("custom_db must be a data.frame") 
  }
  if(nrow(custom_db)<=0) {
    stop("custom_db must be a data.frame") 
  }
  
  e <- new.env()
  name <- data("tcr3d", package = "ClustIRR", envir = e)[1]
  e <- e[[name]]
  if(all(colnames(custom_db) == colnames(e))==FALSE) {
    stop("wrong columns in custom_db, see internal data tcr3d.")
  }
  return(custom_db)
}


#### Helper functions ####

check_aa <- function(x) {
  aa <- "[^ACDEFGHIKLMNPQRSTVWY]"
  for(col in get_chains(colnames(x))) {
    res <- lapply(x[,col], function(y) length(grep(aa, y)) != 0)
    if(any(unlist(res))) {
      stop("non-standard amino acid symbols in input CDR")
    }
  }
}

check_dataframe <- function(x) {
  w <- paste0(deparse(substitute(x)), " has to be a dataframe")
  if(!is.data.frame(x)) {
    stop(w)
  }
}

check_s_cols <- function(x) {
  chains <- paste0("CDR3", c("a","b","g","d","l","h"))
  if(!any(colnames(x) %in% c("clone_size", "sample", chains))) {
    stop(paste0("allowed columns in s: ", chains, " clone_size and sample"))
  }
  if((any(colnames(x) %in% chains[1:2])==TRUE &
      all(colnames(x) %in% c(chains[1:2], "clone_size", "sample"))==FALSE)|
     (any(colnames(x) %in% chains[3:4])==TRUE &
      all(colnames(x) %in% c(chains[3:4], "clone_size", "sample"))==FALSE)|
     (any(colnames(x) %in% chains[5:6])==TRUE &
      all(colnames(x) %in% c(chains[5:6], "clone_size", "sample"))==FALSE)) {
    stop(paste0("allowed chain pairs: CDR3a-CDR3b CDR3d-CDR3g CDR3l-CDR3h"))
  }
  
  if(any(colnames(x) == "clone_size")) {
    if(any(is.numeric(x$clone_size)==FALSE)) {
      stop("clone_size must be numeric")
    }
  }
  if(any(colnames(x) == "sample")) {
    if(any(is.character(x$sample)==FALSE)) {
      stop("sample must be character")
    }
    if(length(unique(x$sample))!=1) {
      stop("multiple sample ids found")
    }
  }
}

check_dataframe_empty <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " contains empty values"
  )
  if(any(x == "", na.rm = TRUE)) {
    warning(w)
  }
}

check_dataframe_na <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " contains NA value"
  )
  if(any(is.na(x))) {
    warning(w)
  }
}

check_greaterthan <- function(x, v) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be <= ",
    v
  )
  if(any(x > v)) {
    stop(w)
  }
}

check_infinity <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be a finite number"
  )
  if(any(is.infinite(x))) {
    stop(w)
  }
}

check_lessthan <- function(x, v) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be >= ",
    v
  )
  if(any(x < v)) {
    stop(w)
  }
}

check_logical <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be logical"
  )
  if(any(is.na(x))) {
    stop(w)
  }
  if(!is.logical(x)) {
    stop(w)
  }
}

check_matrix <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be of type matrix"
  )
  if(!is.matrix(x)) {
    stop(w)
  }
}

check_matrix_column_count <- function(x, c) {
  w <- paste0(
    deparse(substitute(x)),
    " has to have ", c, " columns"
  )
  if(ncol(x) != c) {
    stop(w)
  }
}

check_matrix_type <- function(x, type) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be a numeric matrix"
  )
  
  if(type == "numeric") {
    if(is.numeric(x) == FALSE) {
      stop(w)
    }
  }
  if(type == "character") {
    if(is.character(x) == FALSE) {
      stop(w)
    }
  }
}

check_missing <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " is missing"
  )
  if(missing(x) || is.null(x)) {
    stop(w)
  }
}

check_numeric <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be numeric"
  )
  if(any(!is.numeric(x))) {
    stop(w)
  }
}

check_rowcount <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " contains zero rows"
  )
  if(nrow(x) == 0) {
    stop(w)
  }
}

check_singlevalue <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be a single value"
  )
  if(any(is.na(x))) {
    stop(w)
  }
  if(length(x) != 1) {
    stop(w)
  }
}

check_wholenumber <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be a whole number"
  )
  if(any(!(abs(x - round(x)) < .Machine$double.eps^0.5))) {
    stop(w)
  }
}

check_positive <- function(x) {
  w <- paste0(
    deparse(substitute(x)),
    " has to be positive number"
  )
  if(x < 0) {
    stop(w)
  }
}

check_probability <- function(x) {
    w <- paste0(deparse(substitute(x)), " must be probability")
    if(x < 0 | x > 1) {
        stop(w)
    }
}
