test_that("get_trimmed_flanks() leaves sequences to analyse after trimming", {
  # input sequences. min(nchar(s)) = 3, max(nchar(s)) = 15
  s <- base::c("CDRCDRCDRCDR", "CDRCDRCDRCDRCDR", "CDRCDR", "CDRC")

  # trim by 8 aa from each side to trim all sequences completely
  t <- 8
  e <- "trim_flank_aa too high, no sequences left to cluster after trimming"
  expect_error(get_trimmed_flanks(s, t), regexp = e)
  
  # trim by 2 aa from each side to trim all sequences <= 4 completely
  t <- 2
  
  # check warning for trimmed sample dataset
  cdr3 <- s
  n <- "sample"
  ws <- base::paste0(
    "cdr3 sequences shorter than 4 (trim_flank_aa*2) of the ",
    "sample dataset \n were trimmed completely before local clustering")
  expect_warning(get_trimmed_flanks(cdr3, t), regexp = ws, fixed = TRUE)
  
  # check warning for trimmed reference dataset
  cdr3_ref <- s
  wr <- base::paste0(
    "cdr3 sequences shorter than 4 (trim_flank_aa*2) of the ",
    "reference dataset \n were trimmed completely before local clustering")
  expect_warning(get_trimmed_flanks(cdr3_ref, t), regexp = wr, fixed = TRUE)
  
})
