test_that("check_singlevalue() runs as expected", {
    expect_error(check_singlevalue(NA),
        regexp = "NA has to be a single value"
    )
})

test_that("check_missing() runs as expected", {
    expect_error(check_missing(),
        regexp = " parameter is missing"
    )
    expect_error(check_missing(NULL),
        regexp = "NULL parameter is missing"
    )
})

test_that("check_logical() runs as expected", {
    expect_error(check_logical(NA),
        regexp = "NA has to be logical"
    )
})

test_that("check_dataframe_na() runs as expected", {
    df <- base::data.frame(CDR3b = NA)
    expect_warning(check_dataframe_na(df),
        regexp = "df contains NA value"
    )
})

test_that("check_dataframe_empty() runs as expected", {
    df <- base::data.frame(CDR3b = "")
    expect_warning(check_dataframe_empty(df),
        regexp = "df contains empty values"
    )
})

test_that("check_dataframe_colnames() runs as expected", {
    df <- base::data.frame(CDR3 = base::character())
    e <- base::paste0(
        "df has to contain at least one of the following columns: ",
        "CDR3a and/or CDR3b or CDR3d and/or CDR3g or CDR3h and/or CDR3l"
    )
    expect_error(check_dataframe_colnames(df), regexp = e, fixed = TRUE)
})

test_that("check_dataframe_colnames() takes only character inputs", {
    df <- base::data.frame("CDR3a" = base::c(1, 2))
    expect_error(check_dataframe_colnames(df),
        regexp = "CDR3a column has to of type character"
    )
    df <- base::data.frame("CDR3b" = NA)
    expect_error(check_dataframe_colnames(df),
        regexp = "CDR3b column has to of type character"
    )
    df <- base::data.frame("CDR3d" = base::c(0.005, 0.99))
    expect_error(check_dataframe_colnames(df),
        regexp = "CDR3d column has to of type character"
    )
    df <- base::data.frame("CDR3g" = FALSE)
    expect_error(check_dataframe_colnames(df),
        regexp = "CDR3g column has to of type character"
    )
    df <- base::data.frame("CDR3l" = base::c(0.005, 0.99))
    expect_error(check_dataframe_colnames(df),
        regexp = "CDR3l column has to of type character"
    )
    df <- base::data.frame("CDR3h" = FALSE)
    expect_error(check_dataframe_colnames(df),
        regexp = "CDR3h column has to of type character"
    )
})

test_that("check_dataframe_colnames() takes only valid chain combinations", {
    m <- "CDR3a/b can't be mixed with CDR3d/g columns"
    s <- base::data.frame(
        "CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF"),
        "CDR3d" = c("CATSRLEFAQYF", "CATSRLEARATYF")
    )
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3a", "CDR3g")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3b", "CDR3d")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3b", "CDR3g")
    expect_error(check_dataframe_colnames(s), regexp = m)

    m <- "CDR3a/b can't be mixed with CDR3l/h columns"
    base::names(s) <- base::c("CDR3a", "CDR3l")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3a", "CDR3h")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3b", "CDR3l")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3b", "CDR3h")
    expect_error(check_dataframe_colnames(s), regexp = m)

    m <- "CDR3d/g can't be mixed with CDR3l/h columns"
    base::names(s) <- base::c("CDR3d", "CDR3l")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3d", "CDR3h")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3g", "CDR3l")
    expect_error(check_dataframe_colnames(s), regexp = m)
    base::names(s) <- base::c("CDR3g", "CDR3h")
    expect_error(check_dataframe_colnames(s), regexp = m)
})

test_that("check_s_and_r() runs as expected", {
    s <- base::data.frame("CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF"))
    r <- base::data.frame("CDR3b" = c("CATSRLEFAQYF", "CATSRLEARATYF"))
    expect_error(check_s_and_r(s, r),
        regexp = "s has to contain the same columns as r"
    )
    s <- base::data.frame(
        "CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF"),
        "CDR3b" = c("CATSRLEFAQYF", "CATSRLEARATYF")
    )
    r <- base::data.frame(
        "CDR3b" = c("CATSRLEFAQYF", "CATSRLEARATYF"),
        "CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF")
    )
    expect_no_error(check_s_and_r(s, r))
})


test_that("get_control() runs as expected", {
    control <- base::list(
        global_max_dist = 1,
        local_max_fdr = 0.05,
        local_min_ove = 2,
        local_min_o = 1,
        trim_flank_aa = 0,
        global_pairs = NULL,
        low_mem = FALSE
    )

    expect_equal(get_control(NULL), control) # return control if NULL
    expect_error(get_control("test"),
        regexp = "control must be a list"
    ) # no list as input
    control$add <- "one too much"
    expect_error(get_control(control),
        regexp = "unrecognized elements found in control"
    )
})

test_that("check_aa() runs as expected", {
    s <- base::data.frame(CDR3b = "CSSEDNDSSS12")
    expect_error(
        check_aa(s), 
        regexp = "s contains non-standard or lowercase amino acid codes"
    )
    s <- base::data.frame(CDR3b = "CSSEDNDSSSAS", CDR3a = "CSSEDNDSBSZT")
    expect_error(
        check_aa(s), 
        regexp = "s contains non-standard or lowercase amino acid codes"
    )
    s <- base::data.frame(CDR3b = base::c("ACDEFGHIKLMN", "ACDEFGHdKLMN"))
    expect_error(
        check_aa(s), 
        regexp = "s contains non-standard or lowercase amino acid codes"
    )
    s <- base::data.frame(CDR3b = "ACDEFGHIKLMN", CDR3a = "PQRSTVWY")
    expect_no_error(check_aa(s))
})


test_that("check_ks_and_trim_flank_aa() works as expected", {
    data("CDR3ab")
    s <- base::data.frame(CDR3b = CDR3ab[1:10, "CDR3b"]) # max(nchar(...)) 15
    r <- base::data.frame(CDR3b = CDR3ab[1:100, "CDR3b"]) # max(nchar(...)) 24
    
    # 16 > 15
    expect_error(
        check_ks_and_trim_flank_aa(16, 0, s, r),
        regexp = "ks has to be smaller than the biggest trimmed sequence"
    ) 
    expect_no_error(check_ks_and_trim_flank_aa(15, 0, s, r)) 
    
    # 14 > 13 (15-2*1)
    expect_error(
        check_ks_and_trim_flank_aa(14, 1, s, r),
        regexp = "ks has to be smaller than the biggest trimmed sequence"
    )
    expect_no_error(check_ks_and_trim_flank_aa(13, 1, s, r)) 
    
    # 10 > 9 (15-2*3)
    expect_error(
        check_ks_and_trim_flank_aa(10, 3, s, r),
        regexp = "ks has to be smaller than the biggest trimmed sequence"
    ) 
    expect_no_error(check_ks_and_trim_flank_aa(9, 3, s, r)) 
    
})

