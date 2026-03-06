test_that("get_functional_communities works as expected", {

# test data
data("CDR3ab", package = "ClustIRR")
s <- data.frame(CDR3b = CDR3ab[, "CDR3b"],
                CDR3a = CDR3ab[, "CDR3a"],
                clone_size = 1,
                sample = rep(c("a", "b", "c", "d"), 125))[1:500,]

# run ClustIRR analysis
c <- clustirr(s = s)

com <- detect_communities(graph = c$graph,
                          algorithm = "leiden",
                          resolution = 1, 
                          weight = "ncweight",
                          metric = "average",
                          chains = c("CDR3b", "CDR3a"))

ag_gene <- c("PB1", "pp65")
ag_species <- c("Influenza", "YFV")

expect_no_error(fc <- get_functional_communities(com, ag_gene, ag_species))

# test that 6 functional communities are detected
expect_true(length(fc$functional_communities) == 6)

# test if works with no matches
expect_no_error(get_functional_communities(com, "ag_gene", "ag_species"))


})
