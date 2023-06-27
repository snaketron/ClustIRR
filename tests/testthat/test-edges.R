test_that("get_edges() takes only clust_irr object as input", {
    expect_error(get_edges(NA),
                 regexp = "Input has to be object of class clust_irr")
})

test_that("get_edges() takes clust_irr object as input", {
  # load ref dataset
  data("CD8")
  # sample 500 sequences from the reference dataset as sample dataset
  data_sample <- data.frame(CDR3b = CD8[sample(x = 1:nrow(CD8), 
                                               size = 50, 
                                               replace = FALSE),])
  # run clustirr
  x <- cluster_irr(data_sample, CD8)
  
  expect_no_error(get_edges(x))
})
