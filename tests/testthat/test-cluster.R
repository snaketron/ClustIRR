# load minimal input data
data("CDR3ab")
s <- base::data.frame(CDR3b = CDR3ab[1:500, "CDR3b"])
r <- base::data.frame(CDR3b = CDR3ab[1:1000, "CDR3b"])

# detect cores
cores <- future::availableCores()

# set ks and control input parameters
ks <- c(2, 3, 4)
control_input <- list(
    B = 100,
    global_max_dist = 1,
    local_max_fdr = 0.05,
    local_min_ove = 2,
    local_min_o = 1,
    trim_flank_aa = 3,
    low_mem = FALSE,
    global_pairs = NULL
)


# check everything for all three versions of the algorithm
for (version in c(1, 2, 3)) {
    test_that("s parameter takes only valid input", {
        expect_error(cluster_irr( # missing input
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "s parameter is missing")
        expect_error(cluster_irr(
            s = as.matrix(s), # matrix
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "s has to be of type data frame")
        expect_error(cluster_irr(
            s = s[0:0, ], # 0 row df
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "s has to be of type data frame")
        df <- s
        df["CDR3a"] <- 42
        expect_error(cluster_irr(
            s = df, # numerical CDR3a column
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "CDR3a column has to of type character")
        df <- s
        df["CDR3b"] <- 42
        expect_error(cluster_irr(
            s = df, # numerical CDR3b column
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "CDR3b column has to of type character")
    })


    test_that("r parameter takes only valid input", {
        expect_error(cluster_irr(
            s = s,
            # missing input
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "r parameter is missing")
        expect_error(cluster_irr(
            s = s,
            r = as.matrix(s), # matrix
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "r has to be of type data frame")
        expect_error(cluster_irr(
            s = s,
            r = r[0:0, ], # 0 row df
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "r has to be of type data frame")
        df <- r
        df["CDR3a"] <- 42
        expect_error(cluster_irr(
            s = s,
            r = df, # numerical CDR3a column
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "CDR3a column has to of type character")
        df <- r
        df["CDR3b"] <- 42
        expect_error(cluster_irr(
            s = s,
            r = df, # numerical CDR3b column
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "CDR3b column has to of type character")
    })

    test_that("different combinations of s and r run as expected", {
        s_cdr3a <- s
        r_cdr3a <- r
        base::colnames(s_cdr3a)[1] <- "CDR3a"
        base::colnames(r_cdr3a)[1] <- "CDR3a"
        # cdr3a column only in s
        expect_error(check_s_and_r(
            s_cdr3a,
            r
        ), regexp = "s has to contain the same columns as r")
        # cdr3a column only in r
        expect_error(check_s_and_r(
            s,
            r_cdr3a
        ), regexp = "s has to contain the same columns as r")
        base::colnames(r_cdr3a)[1] <- "CDR3c"
        # cdr3b column only in s
        expect_error(check_s_and_r(
            s,
            r_cdr3a
        ), regexp = "s has to contain the same columns as r")
        base::colnames(s_cdr3a)[1] <- "CDR3c"
        # cdr3b column only in s
        expect_error(check_s_and_r(
            s_cdr3a,
            r
        ), regexp = "s has to contain the same columns as r")
    })

    test_that("version parameter takes only valid input", {
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = "three", # non-numeric
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "version has to be numeric")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = c(1, 2, 3), # multiple values
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "version has to be a single value")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = 1.6, # float
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "version has to be 1, 2 or 3")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = 4, # non-existing version
            ks = ks,
            cores = cores,
            control = control_input
        ), regexp = "version has to be 1, 2 or 3")
    })

    test_that("ks parameter takes only valid input", {
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = Inf, # infinite
            cores = cores,
            control = control_input
        ), regexp = "ks has to be a finite number")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = "hello", # non-numeric
            cores = cores,
            control = control_input
        ), regexp = "ks has to be numeric")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = 1.7, # float
            cores = cores,
            control = control_input
        ), regexp = "ks has to be a whole number")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = 0, # < 1
            cores = cores,
            control = control_input
        ), regexp = "ks has to be >= 1")
    })

    test_that("cores parameter takes only valid input", {
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = Inf, # infinite
            control = control_input
        ), regexp = "cores has to be a finite number")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = "all of them", # non-numeric
            control = control_input
        ), regexp = "cores has to be numeric")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = 1.7, # float
            control = control_input
        ), regexp = "cores has to be a whole number")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = c(1, 2, 3), # multiple values
            control = control_input
        ), regexp = "cores has to be a single value")
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = 0, # < 1
            control = control_input
        ), regexp = "cores has to be >= 1")
    })

    test_that("control_input$B parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$B <- Inf # infinity
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "B has to be a finite number")
        control_input_tmp$B <- "Everyone" # non-numeric
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "B has to be numeric")
        control_input_tmp$B <- 1.7 # float
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "B has to be a whole number")
        control_input_tmp$B <- c(1, 2, 3) # multiple values
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "B has to be a single value")
        control_input_tmp$B <- 0 # < 1
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "B has to be >= 1")
    })

    test_that("control_input$global_max_dist param takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$global_max_dist <- Inf # infinity
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_max_dist has to be a finite number")
        control_input_tmp$global_max_dist <- "Everyone" # non-numeric
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_max_dist has to be numeric")
        control_input_tmp$global_max_dist <- 1.7 # float
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_max_dist has to be a whole number")
        control_input_tmp$global_max_dist <- c(1, 2) # multiple values
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_max_dist has to be a single value")
        control_input_tmp$global_max_dist <- 0 # < 1
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_max_dist has to be >= 1")
    })

    test_that("control_input$local_max_fdr param takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$local_max_fdr <- Inf # infinity
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_max_fdr has to be a finite number")
        control_input_tmp$local_max_fdr <- "Everyone" # non-numeric
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_max_fdr has to be numeric")
        control_input_tmp$local_max_fdr <- -1 # < 0
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_max_fdr has to be >= 0")
        control_input_tmp$local_max_fdr <- 2 # > 1
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_max_fdr has to be <= 1")
    })

    test_that("control_input$local_min_ove parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$local_min_ove <- Inf # infinity
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_min_ove has to be a finite number")
        control_input_tmp$local_min_ove <- "Everyone" # non-numeric
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_min_ove has to be numeric")
        control_input_tmp$local_min_ove <- c(1, 2, 3) # multiple values
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_min_ove has to be a single value")
    })

    test_that("control_input$local_min_o parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$local_min_o <- Inf # infinity
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_min_o has to be a finite number")
        control_input_tmp$local_min_o <- "Everyone" # non-numeric
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_min_o has to be numeric")
        control_input_tmp$local_min_o <- 1.7 # float
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_min_o has to be a whole number")
        control_input_tmp$local_min_o <- c(1, 2, 3) # multiple values
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "local_min_o has to be a single value")
    })

    test_that("control_input$trim_flank_aa parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$trim_flank_aa <- Inf # infinity
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "trim_flank_aa has to be a finite number")
        control_input_tmp$trim_flank_aa <- "Everyone" # non-numeric
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "trim_flank_aa has to be numeric")
        control_input_tmp$trim_flank_aa <- 1.7 # float
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "trim_flank_aa has to be a whole number")
        control_input_tmp$trim_flank_aa <- -2 # positive
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "trim_flank_aa has to be positive number")
        control_input_tmp$trim_flank_aa <- c(1, 2, 3) # single value
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "trim_flank_aa has to be a single value")
    })

    test_that("control_input$low_mem parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$low_mem <- c(1, 2, 3) # multiple values
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "low_mem has to be a single value")
        control_input_tmp$low_mem <- "TRUE" # not logical
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "low_mem has to be logical")
    })

    test_that("control_input$global_pairs parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$global_pairs <- matrix(
            data = NA,
            nrow = 0,
            ncol = 1
        ) # rowcount zero
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_pairs contains zero rows")
        control_input_tmp$global_pairs <- base::data.frame(c(1, 2)) # rowcount 0
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_pairs has to be of type matrix")
        control_input_tmp$global_pairs <- matrix(
            data = "NA",
            nrow = 1,
            ncol = 1
        ) # non-integer matrix
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_pairs has to have 3 columns")
        control_input_tmp$global_pairs <- matrix(
            data = 17L,
            nrow = 1,
            ncol = 3
        ) # ≠ 2 columns
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_pairs has to be a numeric matrix")
        n <- nrow(s) + 1
        control_input_tmp$global_pairs <- matrix(
            data = as.integer(n),
            nrow = 30,
            ncol = 2
        ) # wrong index
        expect_error(cluster_irr(
            s = s,
            r = r,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ), regexp = "global_pairs has to be a numeric matrix")
    })


    # test all versions with correct input
    test_that("cluster_irr works with correct input", {
        control_input_tmp <- control_input
        control_input_tmp$trim_flank_aa <- 0
        expect_no_error(cluster_irr(
            s = s,
            r = r,
            ks = ks,
            cores = cores,
            version = version,
            control = control_input_tmp
        ))
    })

    # test all versions with correct input and low_mem = true
    test_that("cluster_irr works with correct input in low_mem mode", {
        control_input_tmp <- control_input
        control_input_tmp$low_mem <- TRUE
        control_input_tmp$trim_flank_aa <- 0
        expect_no_error(cluster_irr(
            s = s,
            r = r,
            ks = ks,
            cores = cores,
            version = version,
            control = control_input_tmp
        ))
    })
    
    # na checks
    test_that("cluster_irr works with NA included in s and/or r", {
        nas <- base::data.frame(CDR3b = c(NA, NA))
        s_na <- base::rbind(s, nas)
        r_na <- base::rbind(r, nas)
        base::suppressWarnings({
            expect_no_error(cluster_irr(
                s = s_na,
                r = r,
                ks = ks,
                cores = cores,
                version = version,
                control = control_input
            ))
            expect_no_error(cluster_irr(
                s = s,
                r = r_na,
                ks = ks,
                cores = cores,
                version = version,
                control = control_input
            ))
            expect_no_error(cluster_irr(
                s = s_na,
                r = r_na,
                ks = ks,
                cores = cores,
                version = version,
                control = control_input
            ))
        })
    })
    
}

# test get_clust() functions with global_pairs input
test_that("get_clust functions works with global_pairs input", {
    control_input_tmp <- control_input
    control_input_tmp$global_pairs <- 
        base::matrix(data = 17L, nrow = 10, ncol = 2)
    control_input_tmp$trim_flank_aa <- 0
    expect_no_error(get_clust_v1(
        cdr3 = s,
        cdr3_r = r,
        ks = ks,
        cores = cores,
        control = control_input_tmp
    ))
    expect_no_error(get_clust_v23(
        cdr3 = s,
        cdr3_r = r,
        ks = ks,
        cores = cores,
        control = control_input_tmp
    ))
})



