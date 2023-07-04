test_that("get_edges() takes only clust_irr object as input", {
    expect_error(get_edges(NA),
                 regexp = "Input has to be object of class clust_irr")
})

test_that("get_edges() takes clust_irr object as input", {
  data("CDR3ab")
  s <- data.frame(CDR3b = CDR3ab[1:50, "CDR3b"])
  r <- data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"])

  # run clustirr
  x <- cluster_irr(s, r)

  expect_no_error(get_edges(x))
})
