test_that("get_trimmed_flanks() leaves sequences to analyse after trimming", {
  # input sequences. min(nchar(s)) = 3, max(nchar(s)) = 15
  s <- c("CDRCDRCDRCDR", "CDRCDRCDRCDRCDR", "CDRCDR", "CDRC")

  # trim by 8 aa from each side to trim all sequences completely
  t <- 8
  e <- "trim_flank_aa too high, no sequences left to cluster after trimming"
  expect_error(get_trimmed_flanks(s, t), regexp = e)
  
  # trim by 2 aa from each side to trim all sequences <= 4 completely
  t <- 2
  w <- base::paste0(
    "cdr3 sequences shorter than ",
    t*2, 
    " (trim_flank_aa*2) of the ",
    base::deparse(base::substitute(n)),
    " dataset \n were trimmed completely before local clustering \n \n"
  )
  expect_warning(get_trimmed_flanks(s, t), regexp = "w")
  
})
