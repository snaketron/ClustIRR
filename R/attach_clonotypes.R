attach_clonotypes <- function(s, node_summary, control, keep_best = TRUE) {

    validate_inputs <- function(s, node_summary, control, keep_best) {
        check_missing(s)
        check_dataframe(s)
        check_rowcount(s)
        check_dataframe_na(s)
        check_dataframe_empty(s)
        check_s_cols(s)
        check_aa(s)

        check_missing(node_summary)
        check_dataframe(node_summary)
        check_rowcount(node_summary)
        check_dataframe_na(node_summary)
        check_dataframe_empty(node_summary)

        if(!any(colnames(node_summary) %in% "community")) {
            stop("node_summary must include community")
        }

        if(!all(get_chains(colnames(s)) %in% colnames(node_summary))) {
            stop("node_summary must include the chain columns used in s")
        }

        check_aa(node_summary)

        if(missing(control) || is.null(control)) {
            stop("control is missing")
        }
        if(is.list(control) == FALSE) {
            stop("control must be a list")
        }
        if(all(names(control) %in% names(get_control())) == FALSE) {
            stop("unrecognized elements found in control")
        }

        check_trim_flank_aa(control$trim_flank_aa)
        check_blast_gmi(control$blast_gmi)
        check_blast_cores(control$blast_cores)

        check_singlevalue(keep_best)
        check_logical(keep_best)
    }

    validate_inputs(s = s, node_summary = node_summary, control = control, keep_best = keep_best)

    attach <- function(s, node_summary, chain, control, keep_best) {
        # rblast query
        db_from <- AAStringSet(s[,chain])
        names(db_from) <- 1:nrow(s)

        # rblast db
        db_to <- AAStringSet(node_summary[,chain])
        names(db_to) <- 1:nrow(node_summary)

        fasta_path <- tempfile(pattern = "fs_", fileext = ".fasta", tmpdir = ".")
        db_path <- tempfile(tmpdir = ".")

        writeXStringSet(db_to, filepath = fasta_path)
        makeblastdb(db_name = db_path, dbtype = "prot",
                    file = fasta_path, verbose = FALSE)
        bdb <- blast(db_path, type = "blastp")
        blast_fmt <- paste("qseqid sseqid pident length qstart qend",
                           "sstart send evalue score qseq sseq qlen slen")
        blast_args <- paste0("-num_threads ", control$blast_cores,
                             " -mt_mode 1",
                             " -evalue 10 -gapopen 11 -gapextend 2",
                             " -word_size 2 -matrix BLOSUM62 -comp_based_stats 0")
        o <- predict(bdb, db_from, BLAST_args = blast_args,
                     custom_format = blast_fmt,
                     keep_tmp = FALSE, verbose = FALSE)

        # clean up
        unlink(fasta_path, force = TRUE)

        # Remove BLAST database files sharing db_path as prefix
        unlink(paste0(db_path, ".*"), force = TRUE)

        # remove low gmi matches
        o <- o[o$pident >= control$blast_gmi*100,]
        if(nrow(o)==0) {
            return(NULL)
        }

        # remove partial hits
        j <- which(o$qstart!=1 | o$sstart!=1 | o$qend!=o$qlen | o$send!=o$slen)
        if(length(j)!=0) {
            o <- o[-j,]
        }
        # if empty and no duplicates found stop
        if(nrow(o)==0) {
            return(NULL)
        }

        out <- c()
        if(nrow(o)>0) {
            out <- data.frame(from = o$qseqid,
                              to = o$sseqid,
                              len = o$length,
                              weight = o$score,
                              nweight = o$score/o$length)

            data_env <- get_blosum62()
            bs <- t(vapply(X = seq_len(nrow(o)),
                           FUN.VALUE = numeric(3),
                           FUN = get_bscore,
                           o = o,
                           gap_o = -11,
                           gap_e = -2,
                           trim = control$trim_flank_aa,
                           b = data_env[["BLOSUM62"]]))
            out$cweight <- bs[,2]
            out$clen <- bs[,3]

            out$nweight <- out$weight/out$len
            out$ncweight <- out$cweight/out$clen
            out$chain <- chain
            out$community <- node_summary$community[out$to]
            out <- out[order(out$nweight, decreasing = TRUE),]
            if(keep_best) {
                out <- out[duplicated(out$from)==FALSE,]
            }
        }
        return(out)
    }

    # get chains to be analyzed
    chains <- get_chains(colnames(s))

    return(lapply(X = chains,
                  FUN = attach,
                  s = s,
                  node_summary = node_summary,
                  control = control,
                  keep_best = keep_best))
}
