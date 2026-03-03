test_that("get_clonotype_cosine_similarity works as expected", {
    
    data("CDR3ab", package = "ClustIRR")
    s <- data.frame(CDR3b = CDR3ab[1:5, "CDR3b"],
                    CDR3a = CDR3ab[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = "a")
    
    # run ClustIRR analysis
    c <- clustirr(s = s)
    
    # only one clust_irr should throw an error
    expect_error(get_clonotype_cosine_similarity(c))
    
    s <- data.frame(CDR3b = CDR3ab[1:5, "CDR3b"],
                    CDR3a = CDR3ab[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = c("a", "a", "b", "c", "d"))
    
    c <- clustirr(s = s)
    
    expect_no_error(get_clonotype_cosine_similarity(c))
    expect_no_error(get_clonotype_cosine_similarity(c, "CDR3a"))
    expect_no_error(get_clonotype_cosine_similarity(c, "CDR3b"))
    expect_no_error(get_clonotype_cosine_similarity(c, c("CDR3a", "CDR3b")))
    expect_no_error(get_clonotype_cosine_similarity(c, c("CDR3b", "CDR3a")))
    
    # wrong chains input
    expect_error(get_clonotype_cosine_similarity(c, "CDR2a"))
    
    # non-clustirr input
    expect_error(get_clonotype_cosine_similarity(2))
    
    s$CDR3b[3] <- NA
    s$CDR3a[4] <- NA

    c <- expect_warning(clustirr(s = s))
    
    # repertoires with missing CDR3a / CDR3b
    expect_no_error(get_clonotype_cosine_similarity(c))
    expect_no_error(get_clonotype_cosine_similarity(c, "CDR3a"))
    expect_no_error(get_clonotype_cosine_similarity(c, "CDR3b"))
    
    # only CDR3a chain as input and no chains argument
    s <- data.frame(CDR3a = CDR3ab[1:100, "CDR3a"], 
                    sample = rep(c("A", "B"), 50), 
                    clone_size = 1)
    
    # run analysis
    c <- clustirr(s = s)
    
    # Compute similarities and plot
    expect_no_error(get_clonotype_cosine_similarity(c))
    
    # only CDR3b chain as input and no chains argument
    s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                    sample = rep(c("A", "B"), 50), 
                    clone_size = 1)
    
    # run analysis
    c <- clustirr(s = s)
    
    # Compute similarities and plot
    expect_no_error(get_clonotype_cosine_similarity(c))
    
})