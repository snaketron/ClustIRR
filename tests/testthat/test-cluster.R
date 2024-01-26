# load minimal input data
data("CDR3ab")
s <- base::data.frame(CDR3b = CDR3ab[1:500, "CDR3b"])
r <- base::data.frame(CDR3b = CDR3ab[1:1000, "CDR3b"])

# detect cores
cores <- 1

# set ks and control input parameters
ks <- c(2, 3, 4)
control_input <- list(
  global_max_dist = 1,
  local_max_fdr = 0.05,
  local_min_ove = 2,
  local_min_o = 1,
  trim_flank_aa = 3,
  low_mem = FALSE,
  global_pairs = NULL,
  global_smart = FALSE)


# check everything for all three versions of the algorithm
test_that("s parameter takes only valid input", {
  expect_error(cluster_irr( # missing input
    r = r,
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "s parameter is missing")
  expect_error(cluster_irr(
    s = as.matrix(s), # matrix
    r = r,
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "s has to be of type data frame")
  expect_error(cluster_irr(
    s = s[0:0, ], # 0 row df
    r = r,
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "s has to be of type data frame")
  df <- s
  df["CDR3a"] <- 42
  expect_error(cluster_irr(
    s = df, # numerical CDR3a column
    r = r,
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "non\\-character columns in s\\/r")
  df <- s
  df["CDR3b"] <- 42
  expect_error(cluster_irr(
    s = df, # numerical CDR3b column
    r = r,
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "non\\-character columns in s\\/r")
})


test_that("r parameter takes only valid input", {
  expect_message(cluster_irr(
    s = s,
    # missing r
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "missing input r, global clustering mode only")
  expect_error(cluster_irr(
    s = s,
    r = as.matrix(s), # matrix
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "r has to be of type data frame")
  expect_error(cluster_irr(
    s = s,
    r = r[0:0, ], # 0 row df
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "r has to be of type data frame")
  df <- r
  df["CDR3a"] <- 42
  expect_error(cluster_irr(
    s = s,
    r = df, # numerical CDR3a column
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "non\\-character columns in s\\/r")
  df <- r
  df["CDR3b"] <- 42
  expect_error(cluster_irr(
    s = s,
    r = df, # numerical CDR3b column
    ks = ks,
    cores = cores,
    control = control_input
  ), regexp = "non\\-character columns in s\\/r")
})

test_that("different combinations of s and r run as expected", {
  s_cdr3a <- s
  r_cdr3a <- r
  base::colnames(s_cdr3a)[1] <- "CDR3a"
  base::colnames(r_cdr3a)[1] <- "CDR3a"
  # cdr3a column only in s
  expect_error(check_s_r(
    s_cdr3a,
    r
  ), regexp = "s has to contain the same columns as r")
  # cdr3a column only in r
  expect_error(check_s_r(s, r_cdr3a), 
               regexp = "s has to contain the same columns as r")
  base::colnames(r_cdr3a)[1] <- "CDR3c"
  # cdr3b column only in s
  expect_error(check_s_r(s, r_cdr3a), 
               regexp = paste0("unallowed columns in s\\/r, allowed are ",
                               "CDR3a CDR3b CDR3d CDR3g CDR3l CDR3h"))
  base::colnames(s_cdr3a)[1] <- "CDR3c"
  # cdr3b column only in s
  expect_error(check_s_r(s_cdr3a, r), 
               regexp = paste0("unallowed columns in s/r, allowed are ",
                               "CDR3a CDR3b CDR3d CDR3g CDR3l CDR3h"))
})

test_that("ks parameter takes only valid input", {
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = Inf, # infinite
    cores = cores,
    control = control_input
  ), regexp = "ks has to be a finite number")
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = "hello", # non-numeric
    cores = cores,
    control = control_input
  ), regexp = "ks has to be numeric")
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = 1.7, # float
    cores = cores,
    control = control_input
  ), regexp = "ks has to be a whole number")
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = 0, # < 1
    cores = cores,
    control = control_input
  ), regexp = "ks has to be >= 1")
})

test_that("cores parameter takes only valid input", {
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = Inf, # infinite
    control = control_input
  ), regexp = "cores has to be a finite number")
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = "all of them", # non-numeric
    control = control_input
  ), regexp = "cores has to be numeric")
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = 1.7, # float
    control = control_input
  ), regexp = "cores has to be a whole number")
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = c(1, 2, 3), # multiple values
    control = control_input
  ), regexp = "cores has to be a single value")
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = 0, # < 1
    control = control_input
  ), regexp = "cores has to be >= 1")
})

test_that("control_input$global_max_dist param takes only valid input", {
  control_input_tmp <- control_input
  control_input_tmp$global_max_dist <- Inf # infinity
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "global_max_dist has to be a finite number")
  control_input_tmp$global_max_dist <- "Everyone" # non-numeric
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "global_max_dist has to be numeric")
  control_input_tmp$global_max_dist <- 1.7 # float
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "global_max_dist has to be a whole number")
  control_input_tmp$global_max_dist <- c(1, 2) # multiple values
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "global_max_dist has to be a single value")
  control_input_tmp$global_max_dist <- 0 # < 1
  expect_error(cluster_irr(
    s = s,
    r = r,
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
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_max_fdr has to be a finite number")
  control_input_tmp$local_max_fdr <- "Everyone" # non-numeric
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_max_fdr has to be numeric")
  control_input_tmp$local_max_fdr <- -1 # < 0
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_max_fdr has to be >= 0")
  control_input_tmp$local_max_fdr <- 2 # > 1
  expect_error(cluster_irr(
    s = s,
    r = r,
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
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_min_ove has to be a finite number")
  control_input_tmp$local_min_ove <- "Everyone" # non-numeric
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_min_ove has to be numeric")
  control_input_tmp$local_min_ove <- c(1, 2, 3) # multiple values
  expect_error(cluster_irr(
    s = s,
    r = r,
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
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_min_o has to be a finite number")
  control_input_tmp$local_min_o <- "Everyone" # non-numeric
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_min_o has to be numeric")
  control_input_tmp$local_min_o <- 1.7 # float
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "local_min_o has to be a whole number")
  control_input_tmp$local_min_o <- c(1, 2, 3) # multiple values
  expect_error(cluster_irr(
    s = s,
    r = r,
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
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "trim_flank_aa has to be a finite number")
  control_input_tmp$trim_flank_aa <- "Everyone" # non-numeric
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "trim_flank_aa has to be numeric")
  control_input_tmp$trim_flank_aa <- 1.7 # float
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "trim_flank_aa has to be a whole number")
  control_input_tmp$trim_flank_aa <- -2 # positive
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "trim_flank_aa has to be positive number")
  control_input_tmp$trim_flank_aa <- c(1, 2, 3) # single value
  expect_error(cluster_irr(
    s = s,
    r = r,
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
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "low_mem has to be a single value")
  control_input_tmp$low_mem <- "TRUE" # not logical
  expect_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ), regexp = "low_mem has to be logical")
})


test_that("cluster works with correct input", {
  control_input_tmp <- control_input
  control_input_tmp$trim_flank_aa <- 0
  expect_no_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ))
})

# test all versions with correct input and low_mem = true
test_that("cluster works with correct input in low_mem mode", {
  control_input_tmp <- control_input
  control_input_tmp$low_mem <- TRUE
  control_input_tmp$trim_flank_aa <- 0
  expect_no_error(cluster_irr(
    s = s,
    r = r,
    ks = ks,
    cores = cores,
    control = control_input_tmp
  ))
})

# na checks
test_that("cluster works with NA included in s and/or r", {
  nas <- base::data.frame(CDR3b = c(NA, NA))
  s_na <- base::rbind(s, nas)
  r_na <- base::rbind(r, nas)
  base::suppressWarnings({
    expect_no_error(cluster_irr(
      s = s_na,
      r = r,
      ks = ks,
      cores = cores,
      control = control_input
    ))
    expect_no_error(cluster_irr(
      s = s,
      r = r_na,
      ks = ks,
      cores = cores,
      control = control_input
    ))
    expect_no_error(cluster_irr(
      s = s_na,
      r = r_na,
      ks = ks,
      cores = cores,
      control = control_input
    ))
  })
})
