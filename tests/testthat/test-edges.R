test_that("get_edges() takes only clust_irr object as input", {
    expect_error(get_edges(NA),
                 regexp = "Input has to be object of class clust_irr")
})
