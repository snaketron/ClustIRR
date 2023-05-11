# create valid input parameters
cores <- parallel::detectCores()
data("hs_CD8_ref")
# minimal input, only 1000 rows + 50 samples
data_ref <- hs_CD8_ref[1:1000, 1:3]
data_sample <- hs_CD8_ref[sample(
    x = 1:nrow(data_ref),
    size = 50, replace = FALSE
), 1:3]
ks <- c(2, 3, 4)
control_input <- list(
    B = 100,
    global_max_dist = 1,
    local_max_fdr = 0.05,
    local_min_ove = 2,
    local_min_o = 3,
    trim_flank_aa = 3,
    low_mem = FALSE,
    global_pairs = NULL
)

# check everything for all three gliphR versions
for (version in c(1, 2, 3)) {
    test_that("data_sample parameter takes only valid input", {
        expect_error(gliph( # missing input
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = as.matrix(data_sample), # matrix
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample[0:0, ], # 0 row df
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        df <- data_sample
        df["CDR3a"] <- 42
        expect_error(gliph(
            data_sample = df, # numerical CDR3a column
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        df <- data_sample
        df["CDR3b"] <- 42
        expect_error(gliph(
            data_sample = df, # numerical CDR3b column
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
    })


    test_that("data_ref parameter takes only valid input", {
        expect_error(gliph(
            data_sample = data_sample,
            # missing input
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = as.matrix(data_sample), # matrix
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref[0:0, ], # 0 row df
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        df <- data_ref
        df["CDR3a"] <- 42
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = df, # numerical CDR3a column
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
        df <- data_ref
        df["CDR3b"] <- 42
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = df, # numerical CDR3b column
            version = version,
            ks = ks,
            cores = cores,
            control = control_input
        ))
    })

    test_that("different combinations of sample and ref run as expected", {
        data_sample_cdr3a <- data_sample
        data_ref_cdr3a <- data_ref
        colnames(data_sample_cdr3a)[1] <- "CDR3a"
        colnames(data_ref_cdr3a)[1] <- "CDR3a"
        # cdr3a column only in data_sample
        expect_error(check_data_sample_and_ref(
            data_sample_cdr3a,
            data_ref
        ))
        # cdr3a column only in data_ref
        expect_error(check_data_sample_and_ref(
            data_sample,
            data_ref_cdr3a
        ))
        colnames(data_ref_cdr3a)[1] <- "CDR3c"
        # cdr3b column only in data_sample
        expect_error(check_data_sample_and_ref(
            data_sample,
            data_ref_cdr3a
        ))
        colnames(data_sample_cdr3a)[1] <- "CDR3c"
        # cdr3b column only in data_sample
        expect_error(check_data_sample_and_ref(
            data_sample_cdr3a,
            data_ref
        ))
    })

    test_that("version parameter takes only valid input", {
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = "three", # non-numeric
            ks = ks,
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = c(1, 2, 3), # multiple values
            ks = ks,
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = 1.6, # float
            ks = ks,
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = 4, # non-existing version
            ks = ks,
            cores = cores,
            control = control_input
        ))
    })

    test_that("ks parameter takes only valid input", {
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = Inf, # infinite
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = "hello", # non-numeric
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = 1.7, # float
            cores = cores,
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = 0, # < 1
            cores = cores,
            control = control_input
        ))
    })

    test_that("cores parameter takes only valid input", {
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = Inf, # infinite
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = "all of them", # non-numeric
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = 1.7, # float
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = c(1, 2, 3), # multiple values
            control = control_input
        ))
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = 0, # < 1
            control = control_input
        ))
    })

    test_that("control_input$B parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$B <- Inf # infinity
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$B <- "Everyone" # non-numeric
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$B <- 1.7 # float
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$B <- c(1, 2, 3) # multiple values
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$B <- 0 # < 1
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })

    test_that("control_input$global_max_dist param takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$global_max_dist <- Inf # infinity
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$global_max_dist <- "Everyone" # non-numeric
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$global_max_dist <- 1.7 # float
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$global_max_dist <- c(1, 2) # multiple values
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$global_max_dist <- 0 # < 1
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })

    test_that("control_input$local_max_fdr param takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$local_max_fdr <- Inf # infinity
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_max_fdr <- "Everyone" # non-numeric
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_max_fdr <- 1.7 # float
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_max_fdr <- -1 # < 0
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_max_fdr <- 2 # > 1
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })

    test_that("control_input$local_min_ove parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$local_min_ove <- Inf # infinity
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_min_ove <- "Everyone" # non-numeric
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_min_ove <- c(1, 2, 3) # multiple values
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })

    test_that("control_input$local_min_o parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$local_min_o <- Inf # infinity
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_min_o <- "Everyone" # non-numeric
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_min_o <- 1.7 # float
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$local_min_o <- c(1, 2, 3) # multiple values
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })

    test_that("control_input$trim_flank_aa parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$trim_flank_aa <- Inf # infinity
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$trim_flank_aa <- "Everyone" # non-numeric
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$trim_flank_aa <- 1.7 # float
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$trim_flank_aa <- -2 # positive
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$trim_flank_aa <- c(1, 2, 3) # single value
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })

    test_that("control_input$low_mem parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$low_mem <- c(1, 2, 3) # multiple values
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$low_mem <- "TRUE" # not logical
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })

    test_that("control_input$global_pairs parameter takes only valid input", {
        control_input_tmp <- control_input
        control_input_tmp$global_pairs <- matrix(
            data = NA,
            nrow = 0,
            ncol = 1
        ) # rowcount zero
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$global_pairs <- data.frame(c(1, 2)) # rowcount zero
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$global_pairs <- matrix(
            data = "NA",
            nrow = 1,
            ncol = 1
        ) # non-integer matrix
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        control_input_tmp$global_pairs <- matrix(
            data = 17L,
            nrow = 1,
            ncol = 3
        ) # ≠ 2 columns
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
        n <- nrow(data_sample) + 1
        control_input_tmp$global_pairs <- matrix(
            data = as.integer(n),
            nrow = 30,
            ncol = 2
        ) # wrong index
        expect_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            version = version,
            ks = ks,
            cores = cores,
            control = control_input_tmp
        ))
    })


    # test all versions with correct input
    test_that("gliph works with correct input", {
        expect_no_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            ks = ks,
            cores = cores,
            version = version,
            control = control_input
        ))
    })

    # test all versions with correct input and low_mem = true
    test_that("gliph works with correct input in low_mem mode", {
        control_input_tmp <- control_input
        control_input_tmp$low_mem <- TRUE
        expect_no_error(gliph(
            data_sample = data_sample,
            data_ref = data_ref,
            ks = ks,
            cores = cores,
            version = version,
            control = control_input_tmp
        ))
    })
}

# test get_clust() functions with global_pairs input
test_that("get_clust functions works with global_pairs input", {
    control_input_tmp <- control_input
    control_input_tmp$global_pairs <- matrix(data = 17L, nrow = 10, ncol = 2)
    expect_no_error(get_clust_v1(
        cdr3 = data_sample,
        cdr3_ref = data_ref,
        ks = ks,
        cores = cores,
        control = control_input_tmp
    ))
    expect_no_error(get_clust_v23(
        cdr3 = data_sample,
        cdr3_ref = data_ref,
        ks = ks,
        cores = cores,
        control = control_input_tmp
    ))
})
