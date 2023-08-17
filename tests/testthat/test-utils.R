test_that("get_trimmed_flanks() leaves sequences to analyse after trimming", {
    # input sequences. min(nchar(s)) = 3, max(nchar(s)) = 15
    s <- base::c("CDRCDRCDRCDR", "CDRCDRCDRCDRCDR", "CDRCDR", "CDRC")

    # trim by 8 aa from each side to trim all sequences completely
    t <- 8
    e <- "all input CDR3s are shorter than 2 x trim_flank_aa"
    expect_error(get_trimmed_flanks(s, t), regexp = e)

    # trim by 2 aa from each side to trim all sequences <= 4 completely
    t <- 2

    # check warning for trimmed sample dataset
    cdr3 <- s
    n <- "sample"
    ws <- "some input CDR3s are shorter than 2 x trim_flank_aa"
    expect_warning(get_trimmed_flanks(cdr3, t), regexp = ws, fixed = TRUE)

    # check warning for trimmed reference dataset
    cdr3_ref <- s
    wr <- "some input CDR3s are shorter than 2 x trim_flank_aa"
    expect_warning(get_trimmed_flanks(cdr3_ref, t), regexp = wr, fixed = TRUE)
})
