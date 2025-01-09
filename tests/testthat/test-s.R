###### testing s with weird inputs #######

test_that("s parameter takes only valid input", {
  expect_error(cluster_irr(), 
               regexp = "s is missing")
  
  expect_error(cluster_irr(s = NULL), 
               regexp = "s is missing")
  
  expect_error(cluster_irr(s = "C"), 
               regexp = "s has to be a dataframe")
  
  expect_error(cluster_irr(s = TRUE), 
               regexp = "s has to be a dataframe")
  
  expect_error(cluster_irr(s = 1), 
               regexp = "s has to be a dataframe")
  
  expect_error(cluster_irr(s = NA), 
               regexp = "s has to be a dataframe")
})



###### testing s with complex inputs #######

test_that("s parameter takes only valid input", {
  data("CDR3ab", package = "ClustIRR")
  s <- CDR3ab[1:100,c("CDR3a", "CDR3b")]
  s$clone_size <- 1
  s$sample <- "a"
  expect_no_error(cluster_irr(s))
  
  
  # test s
  s <- CDR3ab[1:100,c("CDR3a", "CDR3b")]
  s$clone_size <- 1
  s$sample <- "a"
  s$CDR3b[1] <- NA
  expect_warning(cluster_irr(s = s))
  
  
  # test sample
  data("CDR3ab", package = "ClustIRR")
  s <- CDR3ab[1:100,c("CDR3a", "CDR3b")]
  s$clone_size <- 1
  s$sample <- NA
  expect_error(cluster_irr(s = s),
               regexp = "sample must be character")
  
  data("CDR3ab", package = "ClustIRR")
  s <- CDR3ab[1:100,c("CDR3a", "CDR3b")]
  s$clone_size <- 1
  s$sample <- 1
  expect_error(cluster_irr(s = s),
               regexp = "sample must be character")
  
  data("CDR3ab", package = "ClustIRR")
  s <- CDR3ab[1:100,c("CDR3a", "CDR3b")]
  s$clone_size <- 1
  s$sample <- "a"
  s$sample <- as.factor(s$sample)
  expect_error(cluster_irr(s = s),
               regexp = "sample must be character")
  
  data("CDR3ab", package = "ClustIRR")
  s <- CDR3ab[1:100,c("CDR3a", "CDR3b")]
  s$clone_size <- 1
  s$sample <- T
  expect_error(cluster_irr(s = s),
               regexp = "sample must be character")
  
  data("CDR3ab", package = "ClustIRR")
  s <- CDR3ab[1:100,c("CDR3a", "CDR3b")]
  s$clone_size <- 1
  s$sample <- "a"
  s$sample[1] <- "b"
  expect_error(cluster_irr(s = s),
               regexp = "multiple sample ids found")
})



# more coming soon ...
