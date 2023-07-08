data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:50, "CDR3b"])
r <- data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"])
x <- cluster_irr(s, r)

test_that("get_edges() takes enriched clust_irr object as input", {
  s2 <- s
  substr(s2 = s2$CDR3b[1:20], start = 6, stop = 9) <- "LEAR"
  xi <- cluster_irr(s2, r)
  expect_no_error(get_edges(xi))
})

test_that("get_edges() takes only clust_irr object as input", {
    expect_error(get_edges(NA),
                 regexp = "Input has to be object of class clust_irr")
})

test_that("get_edges() warns if no local or global edges are found", {
  expect_warning(get_edges(x),
                 regexp = "No local or global edges found \n")
})

test_that("get_graph() takes only clust_irr object as input", {
  expect_error(get_graph(NA),
               regexp = "Input has to be object of class clust_irr")
})

test_that("get_graph() warns if no local or global edges are found", {
  expect_warning(get_graph(x),
                 regexp = "No local or global edges to build igraph from \n")
})

test_that("plot_graph() takes only clust_irr object as input", {
  expect_error(plot_graph(NA),
               regexp = "Input has to be object of class clust_irr")
})

test_that("plot_graph() warns if no local or global edges are found", {
  expect_warning(plot_graph(x),
                 regexp = "No graph to plot \n")
})

