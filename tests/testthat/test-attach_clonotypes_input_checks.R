skip_if_not(rBLAST::has_blast(), message = "BLAST+ is not installed.")

test_that("attach_clonotypes validates input arguments", {
    data("CDR3ab", package = "ClustIRR")

    s <- data.frame(
        CDR3b = CDR3ab[1:20, "CDR3b"],
        CDR3a = CDR3ab[1:20, "CDR3a"],
        clone_size = 1,
        sample = "a"
    )

    node_summary <- data.frame(
        community = rep(1L, nrow(s)),
        CDR3b = s$CDR3b,
        CDR3a = s$CDR3a,
        clone_size = 1L
    )

    control <- list(
        blast_gmi = 0.8,
        blast_cores = 1,
        trim_flank_aa = 3,
        db_dist = 0,
        db_custom = NULL,
        knn = FALSE,
        k = 30
    )

    expect_error(attach_clonotypes(),
                 regexp = "s is missing")

    bad_s <- s
    bad_s$CDR3b[1] <- "X"
    expect_error(
        attach_clonotypes(s = bad_s, node_summary = node_summary, control = control),
        regexp = "non-standard amino acid symbols in input CDR"
    )

    bad_node_summary <- node_summary
    bad_node_summary$community <- NULL
    expect_error(
        attach_clonotypes(s = s, node_summary = bad_node_summary, control = control),
        regexp = "node_summary must include community"
    )

    bad_control <- control
    bad_control$blast_gmi <- 1.5
    expect_error(
        attach_clonotypes(s = s, node_summary = node_summary, control = bad_control),
        regexp = "blast_gmi must be probability"
    )

    bad_control <- control
    bad_control$blast_cores <- 0
    expect_error(
        attach_clonotypes(s = s, node_summary = node_summary, control = bad_control),
        regexp = "blast_cores has to be >= "
    )

    expect_error(
        attach_clonotypes(s = s, node_summary = node_summary,
                          control = control, keep_best = c(TRUE, FALSE)),
        regexp = "keep_best has to be a single value"
    )
})
