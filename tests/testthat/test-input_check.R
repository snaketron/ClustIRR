test_that("check_singlevalue() runs as expected", {
  expect_error(check_singlevalue(NA), regexp = "NA has to be a single value")
})

test_that("check_missing() runs as expected", {
  expect_error(check_missing(), regexp = " parameter is missing")
  expect_error(check_missing(NULL), regexp = "NULL parameter is missing")
})

test_that("check_logical() runs as expected", {
  expect_error(check_logical(NA), regexp = "NA has to be logical")
})

test_that("check_dataframe_na() runs as expected", {
  df <- data.frame(CDR3b = NA)
  expect_warning(check_dataframe_na(df), regexp = "df contains NA value")
})

test_that("check_dataframe_empty() runs as expected", {
  df <- data.frame(CDR3b = "")
  expect_warning(check_dataframe_empty(df), regexp = "df contains empty values")
})

test_that("check_r_s_cols() runs as expected", {
  df <- data.frame(CDR3 = character())
  e <- paste0("unallowed columns in s/r, allowed are ",
              "CDR3a, CDR3b, CDR3d, CDR3g, CDR3l, CDR3h ",
              "and clone_size")
  expect_error(check_r_s_cols(df), regexp = e, fixed = TRUE)
})

test_that("check_s_r() takes only character inputs", {
  df <- data.frame("CDR3a" = c(1, 2))
  expect_error(check_s_r(df),
               regexp = "non-standard amino acid symbols in input CDR")
  df <- data.frame("CDR3b" = NA)
  expect_warning(check_s_r(df),
               regexp = "s contains NA value")
  df <- data.frame("CDR3d" = c(0.005, 0.99))
  expect_error(check_s_r(df),
               regexp = "non-standard amino acid symbols in input CDR")
  df <- data.frame("CDR3l" = c(0.005, 0.99))
  expect_error(check_s_r(df),
               regexp = "non-standard amino acid symbols in input CDR")
})

test_that("check_r_s_cols() takes only valid chain combinations", {
  m <- paste0("mixed chains, allowed chain combinations are ",
              "CDR3a x CDR3b, CDR3d x CDR3g, CDR3l x CDR3h")
  
  s <- data.frame("CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF"),
                  "CDR3d" = c("CATSRLEFAQYF", "CATSRLEARATYF"))
  
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3a", "CDR3g")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3b", "CDR3d")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3b", "CDR3g")
  expect_error(check_r_s_cols(s), regexp = m)
  
  names(s) <- c("CDR3a", "CDR3l")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3a", "CDR3h")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3b", "CDR3l")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3b", "CDR3h")
  expect_error(check_r_s_cols(s), regexp = m)
  
  names(s) <- c("CDR3d", "CDR3l")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3d", "CDR3h")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3g", "CDR3l")
  expect_error(check_r_s_cols(s), regexp = m)
  names(s) <- c("CDR3g", "CDR3h")
  expect_error(check_r_s_cols(s), regexp = m)
})

test_that("check_s_r() runs as expected", {
  s <- data.frame("CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF"))
  r <- data.frame("CDR3b" = c("CATSRLEFAQYF", "CATSRLEARATYF"))
  expect_error(check_s_r(s,r), regexp= "s has to contain the same columns as r")
  s <- data.frame(
    "CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF"),
    "CDR3b" = c("CATSRLEFAQYF", "CATSRLEARATYF"))
  r <- data.frame(
    "CDR3b" = c("CATSRLEFAQYF", "CATSRLEARATYF"),
    "CDR3a" = c("CATSRLEFAQYF", "CATSRLEARATYF"))
  expect_no_error(check_s_r(s, r))
})

test_that("get_control() runs as expected", {
  control <- list(global_hamming = FALSE,
                  global_max_hdist = 1,
                  global_min_identity = 0.7,
                  local_max_fdr = 0.05,
                  local_min_o = 1,
                  trim_flank_aa = 0,
                  low_mem = FALSE)
  
  expect_equal(get_control(NULL), control)
  expect_error(get_control("test"), regexp = "control must be a list")
  
  control$add <- "one too much"
  expect_error(get_control(control), 
               regexp = "unrecognized elements found in control")
})

test_that("check_aa() runs as expected", {
  s <- data.frame(CDR3b = "CSSEDNDSSS12")
  expect_error(check_aa(s), 
               regexp = "non-standard amino acid symbols in input CDR")
  s <- data.frame(CDR3b = "CSSEDNDSSSAS", CDR3a = "CSSEDNDSBSZT")
  expect_error(check_aa(s), 
               regexp = "non-standard amino acid symbols in input CDR")
  s <- data.frame(CDR3b = c("ACDEFGHIKLMN", "ACDEFGHdKLMN"))
  expect_error(check_aa(s), 
               regexp = "non-standard amino acid symbols in input CDR")
  s <- data.frame(CDR3b = "ACDEFGHIKLMN", CDR3a = "PQRSTVWY")
  expect_no_error(check_aa(s))
})

