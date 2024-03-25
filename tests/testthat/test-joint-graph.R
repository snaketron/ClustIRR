data("CDR3ab")
s_1 <- data.frame(CDR3b = CDR3ab[1:50, "CDR3b"])
s_2 <- data.frame(CDR3b = CDR3ab[51:100, "CDR3b"])

r <- data.frame(CDR3b = CDR3ab[1:10000, "CDR3b"])
x_1 <- suppressWarnings(cluster_irr(s_1, r))
x_2 <- suppressWarnings(cluster_irr(s_2, r))

clust_irrs <- vector(mode = "list", length = 2)
names(clust_irrs) <- c("S1", "S0")
clust_irrs[[1]] <- x_1
clust_irrs[[2]] <- x_2


test_that("get_joint_graph() takes only clust_irr object as input", {
    expect_error(get_joint_graph(NA, NA),
                 regexp = "clust_irrs must be a list of clust_irr objects")
    expect_error(plot_joint_graph(NA, NA),
                 regexp = "clust_irrs must be a list of clust_irr objects")
    expect_error(plot_joint_graph(NA, NA, as_visnet = TRUE),
                 regexp = "clust_irrs must be a list of clust_irr objects")
    
})

test_that("get_joint_graph() no problem", {
    expect_no_error(get_joint_graph(clust_irrs = clust_irrs))
    expect_no_error(plot_joint_graph(clust_irrs = clust_irrs))
    expect_no_error(plot_joint_graph(clust_irrs = clust_irrs, 
                                     as_visnet = TRUE))
    
})
 
