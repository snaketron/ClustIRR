# create valid input parameters
cores <- parallel::detectCores()
data("hs_CD8_ref")
data_ref <- hs_CD8_ref[, 1:3]
data_sample <- hs_CD8_ref[sample(x = 1:nrow(data_ref),
                                 size = 500, replace = FALSE), 1:3]
ks <- c(2, 3, 4)
B <- 100
control_input <- list(
    B = B,
    global_max_dist = 1,
    local_max_fdr = 0.05,
    local_min_ove = 2,
    local_min_o = 3,
    trim_flank_aa = 3,
    low_mem = FALSE,
    global_pairs = NULL)

# check everything for all three gliphR versions
for (version in 1:3){

    test_that("data_sample parameter takes only valid input", {
        expect_error(gliph(data_sample = NULL, # missing input
                           data_ref = data_ref,
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        expect_error(gliph(data_sample = as.matrix(data_sample), # matrix
                           data_ref = data_ref,
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        expect_error(gliph(data_sample = data_sample[0:0,], # 0 row df
                           data_ref = data_ref,
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        df <- data_sample
        df['CDR3b']=42
        expect_error(gliph(data_sample = df, # numerical CDR3b column
                           data_ref = data_ref,
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        df <- data_sample
        df['TRBV']=NA
        expect_warning(gliph(data_sample = df, # contains NA
                             data_ref = data_ref,
                             ks = ks,
                             cores = cores,
                             version = version,
                             control = control_input))
        df <- data_sample
        df['TRBV']=""
        expect_warning(gliph(data_sample = df, # contains ""
                             data_ref = data_ref,
                             ks = ks,
                             cores = cores,
                             version = version,
                             control = control_input))

    })


    test_that("data_ref parameter takes only valid input", {
        expect_error(gliph(data_sample = data_sample,
                           data_ref = NULL, # missing input
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        expect_error(gliph(data_sample = data_sample,
                           data_ref = as.matrix(data_sample), # matrix
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        expect_error(gliph(data_sample = data_sample,
                           data_ref = data_ref[0:0,], # 0 row df
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        df <- data_ref
        df['CDR3b']=42
        expect_error(gliph(data_sample = data_sample,
                           data_ref = df, # numerical CDR3b column
                           ks = ks,
                           cores = cores,
                           version = version,
                           control = control_input))
        df <- data_ref
        df['TRBV']=NA
        expect_warning(gliph(data_sample = data_sample,
                             data_ref = df, # contains NA
                             ks = ks,
                             cores = cores,
                             version = version,
                             control = control_input))
        df <- data_ref
        df['TRBV']=""
        expect_warning(gliph(data_sample = data_sample,
                             data_ref = df, # contains ""
                             ks = ks,
                             cores = cores,
                             version = version,
                             control = control_input))

    })

    # # test with correct input (maybe not necessary)
    # test_that("gliph works with correct input", {
    #     expect_no_error(gliph(data_sample = data_sample,
    #                        data_ref = data_ref,
    #                        ks = ks,
    #                        cores = cores,
    #                        version = version,
    #                        control = control_input))
    # })


}
