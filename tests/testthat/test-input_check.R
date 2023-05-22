test_that("check_singlevalue() runs as expected", {
    expect_error(check_singlevalue(NA))
})

test_that("check_missing() runs as expected", {
    expect_error(check_missing())
    expect_error(check_missing(NULL))
})

test_that("check_logical() runs as expected", {
    expect_error(check_logical(NA))
})

test_that("check_dataframe_na() runs as expected", {
    df <- data.frame(CDR3b = NA)
    expect_warning(check_dataframe_na(df))
})

test_that("check_dataframe_empty() runs as expected", {
    df <- data.frame(CDR3b = "")
    expect_warning(check_dataframe_empty(df))
})


test_that("check_dataframe_colnames() runs as expected", {
    df <- data.frame(CDR3 = character())
    expect_error(check_dataframe_colnames(df, base::c("CDR3a", "CDR3b")))
})

test_that("get_control() runs as expected", {
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

    expect_equal(get_control(NULL), control) # return control if NULL
    expect_error(get_control("test")) # no list as input
    control$add <- "one too much"
    expect_error(get_control(control))
})
