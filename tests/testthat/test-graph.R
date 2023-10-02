data("CDR3ab")
s <- base::data.frame(CDR3b = CDR3ab[1:50, "CDR3b"])
r <- base::data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"])
x <- cluster_irr(s, r)


test_that("plot_graph() takes only clust_irr object as input", {
    expect_error(plot_graph(NA),
        regexp = "input clust_irr is not class clust_irr")
})

test_that("plot_graph() warns if no local or global edges are found", {
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
})

test_that("plot_graph() can handle edges of only one type", {
    # only alpha edges
    s <- base::data.frame(CDR3b = c("CASYSGANVLTF", "CASSSTGLRDSRDTQYF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- base::data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
    # only beta edges 
    s <- base::data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFSDYF"))
    x <- cluster_irr(s, r)
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
})

test_that("plot_graph() can handle a single clone without external edges", {
    s <- base::data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- base::data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(plot_graph(x))
    expect_no_error(plot_graph(x, as_visnet = TRUE))
})

test_that("plot_graph() can handle NA values in s and/or r", {
    s_e <- base::rbind(s, 
                     base::data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG")))
    s_na <- base::data.frame(CDR3b = base::c("CASSSEFEG", "CASSSEFEG", NA, NA))
    r_na <- base::rbind(r, base::data.frame(CDR3b = c(NA, NA)))
    base::suppressWarnings({x <- cluster_irr(s_na, r)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    
    base::suppressWarnings({x <- cluster_irr(s_e, r_na)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    
    base::suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- plot_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
})

test_that("plot_graph() can handle NA values in s and/or r with 2 chains", {
    s_ba <- base::data.frame(CDR3b = CDR3ab[1:50, "CDR3b"],
                             CDR3a = CDR3ab[1:50, "CDR3a"])
    s_ba <- base::rbind(s_ba, 
                        base::data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG"),
                                         CDR3a = c("CISSEDFS", "CISSEDFS")))
    r_ba <- base::data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"],
                             CDR3a = CDR3ab[1:10000, "CDR3a"])
    s_na <- base::data.frame(CDR3b = base::c("CASSEDFDG", "CASSEDFDG", NA, NA,
                                             "CASSEDFDG"),
                             CDR3a = base::c("CSSSEDFDG", NA, "CSIEDFDG", NA,
                                             "CSIEDFDG"))
    r_na <- base::rbind(r_ba, base::data.frame(CDR3b = NA, CDR3a = NA))
    base::suppressWarnings({x <- cluster_irr(s_na, r_ba)})
    p <- plot_graph(x)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
    
    base::suppressWarnings({x <- cluster_irr(s_ba, r_na)})
    p <- plot_graph(x)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
    
    base::suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- plot_graph(x)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
    p <- plot_graph(x, as_visnet = TRUE)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
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
    s <- base::data.frame(CDR3b = c("CASYSGANVLTF", "CASSSTGLRDSRDTQYF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- base::data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(get_graph(x))
    # only beta edges 
    s <- base::data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFSDYF"))
    x <- cluster_irr(s, r)
    expect_no_error(get_graph(x))
})

test_that("get_graph() can handle a single clone without external edges", {
    s <- base::data.frame(CDR3b = c("CASYSGANVLTF", "CASYSGANVLTF"),
                          CDR3a = c("CAGNTGNQFYF", "CAGNTGNQFYF"))
    r <- base::data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                          CDR3a = CDR3ab[1:100, "CDR3a"])
    x <- cluster_irr(s, r)
    expect_no_error(get_graph(x))
})

test_that("get_graph() can handle NA values in s and/or r", {
    s_e <- base::rbind(s, 
                       base::data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG")))
    s_na <- base::data.frame(CDR3b = base::c("CASSSEFEG", "CASSSEFEG", NA, NA))
    r_na <- base::rbind(r, base::data.frame(CDR3b = c(NA, NA)))
    base::suppressWarnings({x <- cluster_irr(s_na, r)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    base::suppressWarnings({x <- cluster_irr(s_e, r_na)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
    base::suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- get_graph(x)
    expect_true(sum(is.na(p$x$nodes$clone_count)) == 0)
})

test_that("get_graph() can handle NA values in s and/or r with 2 chains", {
    s_ba <- base::data.frame(CDR3b = CDR3ab[1:50, "CDR3b"],
                             CDR3a = CDR3ab[1:50, "CDR3a"])
    s_ba <- base::rbind(s_ba, 
                        base::data.frame(CDR3b = c("CASSSEDFDG", "CASSSEDFDG"),
                                         CDR3a = c("CISSEDFS", "CISSEDFS")))
    r_ba <- base::data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"],
                             CDR3a = CDR3ab[1:10000, "CDR3a"])
    s_na <- base::data.frame(CDR3b = base::c("CASSEDFDG", "CASSEDFDG", NA, NA,
                                             "CASSEDFDG"),
                             CDR3a = base::c("CSSSEDFDG", NA, "CSIEDFDG", NA,
                                             "CSIEDFDG"))
    r_na <- base::rbind(r_ba, base::data.frame(CDR3b = NA, CDR3a = NA))
    base::suppressWarnings({x <- cluster_irr(s_na, r_ba)})
    p <- get_graph(x)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
    base::suppressWarnings({x <- cluster_irr(s_ba, r_na)})
    p <- get_graph(x)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
    base::suppressWarnings({x <- cluster_irr(s_na, r_na)})
    p <- get_graph(x)
    expect_true(base::sum(base::is.na(p$x$nodes$clone_count)) == 0)
})
