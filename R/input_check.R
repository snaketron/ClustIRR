# Description:
# Check user provided input and generate errors and warnings, if necessary
input_check <- function(s,
                        r,
                        ks,
                        cores,
                        control) {
  
  check_s_r(s = s, r = r)
  check_ks(ks)
  check_cores(cores)
  check_global_max_hdist(control$global_max_hdist)
  check_local_max_fdr(control$local_max_fdr)
  check_local_min_o(control$local_min_o)
  check_trim_flank_aa(control$trim_flank_aa)
  check_low_mem(control$low_mem)
  check_global_smart(control$global_smart)
  check_global_pairs(control$global_pairs, s)
}

check_s_r <- function(s, r) {
  check_missing(s)
  check_dataframe(s)
  check_r_s_cols(s)
  check_rowcount(s)
  check_dataframe_na(s)
  check_dataframe_empty(s)
  check_aa(s)
  
  if(missing(r)||is.null(r)) {
    message("missing input r, global clustering mode only")
  } 
  else {
    check_dataframe(r)
    check_r_s_cols(r)
    check_rowcount(r)
    check_dataframe_na(r)
    check_dataframe_empty(r)
    check_aa(r)
    
    if(!all(sort(colnames(s)) == sort(colnames(r)))) {
      stop("s has to contain the same columns as r")
    }
  }
}

check_ks <- function(ks) {
  check_infinity(ks)
  check_numeric(ks)
  check_wholenumber(ks)
  check_lessthan(ks, 1)
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
  if(missing(global_pairs)) {
    stop("global_pairs is missing")
  }
  if(is.null(global_pairs)==FALSE) {
    if(is.data.frame(global_pairs)==FALSE) {
      stop("global_pairs must be a data.frame")
    }
    if(nrow(global_pairs)==0) {
      stop("global_pairs must be a data.frame with at least one row")
    }
    if(ncol(global_pairs)!=4) {
      stop("global_pairs must have 4 columns: from_cdr3, to_cdr3 weight, chain")
    }
    if(any(colnames(global_pairs) %in% c("from_cdr3", "to_cdr3", 
                                         "weight", "chain"))==FALSE) {
      stop("global_pairs must have 4 columns: from_cdr3, to_cdr3 weight, chain")
    }
    
    if(is.character(global_pairs$from_cdr3)==FALSE) {
      stop("from_cdr3 must be a character vector")
    }
    if(is.character(global_pairs$to_cdr3)==FALSE) {
      stop("to_cdr3 must be a character vector")
    }
    if(any(is.na(global_pairs$to_cdr3) | is.null(global_pairs$to_cdr3))) {
      stop("some chains in global_pairs$to_cdr3 are NA or NULL")
    }
    if(any(is.na(global_pairs$from_cdr3) | is.null(global_pairs$from_cdr3))) {
      stop("some chains in global_pairs$from_cdr3 are NA or NULL")
    }
    
    if(is.numeric(global_pairs$weight)==FALSE) {
      stop("weight must be a numeric vector")
    }
    if(global_pairs$weight < 0 | global_pairs$weight > 1) {
      stop("weight must be a number between 0 and 1")
    }
    if(is.na(global_pairs$weight) | is.null(global_pairs$weight)) {
      stop("some weights are NA or NULL")
    }
    
    if(is.character(global_pairs$chain)==FALSE) {
      stop("chain must be a character vector")
    }
    if(any(global_pairs$chain %in% get_chains(colnames(s)))==FALSE) {
      stop("some chains in global_pairs$chain are not part of s")
    }
    if(any(is.na(global_pairs$chain) | is.null(global_pairs$chain))) {
      stop("some chains in global_pairs$chain are NA or NULL")
    }
    for(chain in unique(global_pairs$chain)) {
      if(any(s[,chain] %in% global_pairs$from_cdr3[
        global_pairs$chain==chain])==FALSE) {
        stop("some CDR3 sequences in global_pair are not found in s")
      }
      if(any(s[,chain] %in% global_pairs$to_cdr3[
        global_pairs$chain==chain])==FALSE) {
        stop("some CDR3 sequences in global_pair are not found in s")
      }
    }
  }
}

check_low_mem <- function(low_mem) {
  check_singlevalue(x = low_mem)
  check_logical(x = low_mem)
}

check_global_smart <- function(global_smart) {
  check_singlevalue(x = global_smart)
  check_logical(x = global_smart)
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

# Description:
# Setup control list.
# control_in: user generated list (if missing -> use default)
get_control <- function(control_in) {
  control <- list(global_smart = TRUE,
                  global_max_hdist = 1,
                  local_max_fdr = 0.05,
                  local_min_o = 1,
                  trim_flank_aa = 0,
                  global_pairs = NULL,
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
  w <- paste0(deparse(substitute(x)), " has to be of type data frame")
  if(!is.data.frame(x)) {
    stop(w)
  }
}

check_r_s_cols <- function(x) {
  if(!any(colnames(x) %in% c("clone_size", 
                             paste0("CDR3", c("a","b","g","d","l","h"))))) {
    s <- paste0("unallowed columns in s/r, allowed are ",
                "CDR3a, CDR3b, CDR3d, CDR3g, CDR3l, CDR3h and clone_size")  
    stop(s)
  }
  if((any(colnames(x) %in% c("CDR3a", "CDR3b"))==TRUE &
      all(colnames(x) %in% c("CDR3a", "CDR3b", "clone_size"))==FALSE)|
     (any(colnames(x) %in% c("CDR3g", "CDR3d"))==TRUE &
      all(colnames(x) %in% c("CDR3g", "CDR3d", "clone_size"))==FALSE)|
     (any(colnames(x) %in% c("CDR3l", "CDR3h"))==TRUE &
      all(colnames(x) %in% c("CDR3l", "CDR3h", "clone_size"))==FALSE)) {
    s <- paste0("mixed chains, allowed chain combinations are ",
                "CDR3a x CDR3b, CDR3d x CDR3g, CDR3l x CDR3h")
    stop(s)
  }
  
  if(any(colnames(x) == "clone_size")) {
    if(any(is.numeric(x$clone_size)==FALSE)) {
      stop("clone_size must be numeric")
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
    " parameter is missing"
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
