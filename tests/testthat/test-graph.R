data("CDR3ab")
s <- data.frame(CDR3b = CDR3ab[1:50, "CDR3b"])
r <- data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"])
x <- suppressWarnings(cluster_irr(s, r))
y <- suppressWarnings(cluster_irr(r, s))

clust_irrs <- vector(mode = "list", length = 2)
clust_irrs[[1]] <- x
clust_irrs[[2]] <- y

test_that("plot_joint_graph() takes only clust_irr object as input", {
  expect_error(plot_joint_graph(NA),
               regexp = "clust_irrs must be a list of clust_irr objects")
})

test_that("plot_joint_graph() works as expected", {
  expect_no_error(plot_joint_graph(clust_irrs))
  expect_no_error(plot_joint_graph(clust_irrs, as_visnet = TRUE))
  expect_no_error(plot_joint_graph(clust_irrs, show_singletons = FALSE))
})

test_that("plot_graph() takes only clust_irr object as input", {
    expect_error(plot_graph(NA),
        regexp = "input clust_irr is not class clust_irr")
})

test_that("plot_graph() warns if no local or global edges are found", {
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
    expect_no_error(plot_graph(x, as_visnet = TRUE, show_singletons = FALSE))
})

test_that("plot_graph() can handle edges of only one type", {
    # only alpha edges
    s <- data.frame(CDR3b = c("CASYSGANVLTF", "CASSSTGLRDSRDTQYF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
    # only beta edges 
    s <- data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFSDYF"))
    x <- cluster_irr(s, r)
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
})

test_that("plot_graph() can handle a single clone without external edges", {
    s <- data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
})

test_that("plot_graph() can handle NA values in s and/or r", {
    s_e <- rbind(s, 
                     data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG")))
    s_na <- data.frame(CDR3b = c("CASSSEFEG", "CASSSEFEG", "CASSYSDPNEQFF",
                                 NA, NA))
    r_na <- rbind(r, data.frame(CDR3b = c(NA, NA)))
    suppressWarnings({x <- cluster_irr(s_na, r)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    
    suppressWarnings({x <- cluster_irr(s_e, r_na)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    
    suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
})

test_that("plot_graph() can handle NA values in s and/or r with 2 chains", {
    s_ba <- data.frame(CDR3b = CDR3ab[1:50, "CDR3b"],
                             CDR3a = CDR3ab[1:50, "CDR3a"])
    s_ba <- rbind(s_ba, 
                        data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG"),
                                         CDR3a = c("CISSEDFS", "CISSEDFS")))
    r_ba <- data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"],
                             CDR3a = CDR3ab[1:10000, "CDR3a"])
    s_na <- data.frame(CDR3b = c("CASSEDFDG", "CASSEDFDG", NA, NA, 
                                 "CASSYSDPNEQFF"),
                             CDR3a = c("CSSSEDFDG", NA, "CSIEDFDG", NA,
                                             "CASSYSDPNEQFF"))
    r_na <- rbind(r_ba, data.frame(CDR3b = NA, CDR3a = NA))
    suppressWarnings({x <- cluster_irr(s_na, r_ba)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    
    suppressWarnings({x <- cluster_irr(s_ba, r_na)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    
    suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
})



test_that("get_graph() takes only clust_irr object as input", {
    expect_error(get_graph(NA),
                 regexp = "input clust_irr is not class clust_irr")
})

test_that("get_graph() no problem", {
    expect_no_error(get_graph(x))
})

test_that("get_graph() can handle edges of only one type", {
    # only alpha edges
    s <- data.frame(CDR3b = c("CASYSGANVLTF", "CASSSTGLRDSRDTQYF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(get_graph(x))
    # only beta edges 
    s <- data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFSDYF"))
    x <- cluster_irr(s, r)
    expect_no_error(get_graph(x))
})

test_that("get_graph() can handle a single clone without external edges", {
    s <- data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(get_graph(x))
})

test_that("get_graph() can handle NA values in s and/or r", {
    s_e <- rbind(s, 
                       data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG")))
    s_na <- data.frame(CDR3b = c("CASSSEFEG", "CASSSEFEG", "CASSYSDPNEQFF",
                                 NA, NA))
    r_na <- rbind(r, data.frame(CDR3b = c(NA, NA)))
    suppressWarnings({x <- cluster_irr(s_na, r)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    suppressWarnings({x <- cluster_irr(s_e, r_na)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
})

test_that("get_graph() can handle NA values in s and/or r with 2 chains", {
    s_ba <- data.frame(CDR3b = CDR3ab[1:50, "CDR3b"],
                             CDR3a = CDR3ab[1:50, "CDR3a"])
    s_ba <- rbind(s_ba, 
                        data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG"),
                                         CDR3a = c("CISSEDFS", "CISSEDFS")))
    r_ba <- data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"],
                             CDR3a = CDR3ab[1:10000, "CDR3a"])
    s_na <- data.frame(CDR3b = c("CASSEDFDG", "CASSEDFDG", NA, NA,
                                             "CASSYSDPNEQFF"),
                             CDR3a = c("CSSSEDFDG", NA, "CSIEDFDG", NA,
                                             "CASSYSDPNEQFF"))
    r_na <- rbind(r_ba, data.frame(CDR3b = NA, CDR3a = NA))
    suppressWarnings({x <- cluster_irr(s_na, r_ba)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    suppressWarnings({x <- cluster_irr(s_ba, r_na)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
})
