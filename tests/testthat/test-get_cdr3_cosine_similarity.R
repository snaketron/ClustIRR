test_that("get_cdr3_cosine_similarity works as expected", {
    
    data("CDR3ab", package = "ClustIRR")
    s <- data.frame(CDR3b = CDR3ab[1:5, "CDR3b"],
                    CDR3a = CDR3ab[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = "a")
    
    # run ClustIRR analysis
    c <- clustirr(s = s)
    
    # only one clust_irr should throw an error
    expect_error(get_cdr3_cosine_similarity(c))
    
    # test that it works as expected with correct input
    s <- data.frame(CDR3b = CDR3ab[1:5, "CDR3b"],
                    CDR3a = CDR3ab[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = c("a", "a", "b", "c", "d"))
    
    c <- clustirr(s = s)
    
    expect_no_error(get_cdr3_cosine_similarity(c))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3a"))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3b"))
    expect_no_error(get_cdr3_cosine_similarity(c, c("CDR3a", "CDR3b")))
    expect_no_error(get_cdr3_cosine_similarity(c, c("CDR3b", "CDR3a")))
    
    # wrong chains input
    expect_error(get_cdr3_cosine_similarity(c, "CDR2a"))
    
    # non-clustirr input
    expect_error(get_cdr3_cosine_similarity(2))
    
    s$CDR3b[3] <- NA
    s$CDR3a[4] <- NA
    
    c <- expect_warning(clustirr(s = s))
    
    # repertoires with missing CDR3a / CDR3b
    expect_no_error(get_cdr3_cosine_similarity(c))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3a"))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3b"))
    
    # only CDR3a chain as input and no chains argument
    s <- data.frame(CDR3a = CDR3ab[1:100, "CDR3a"], 
                    sample = rep(c("A", "B"), 50), 
                    clone_size = 1)
    
    c <- clustirr(s = s)
    
    # chains provided not matching chains to plot
    expect_error(get_cdr3_cosine_similarity(c, c("CDR3a", "CDR3b")))
    expect_error(get_cdr3_cosine_similarity(c, "CDR3b"))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3a"))
    
    
    # only CDR3b chain as input and no chains argument
    s <- data.frame(CDR3b = CDR3ab[1:100, "CDR3b"], 
                    sample = rep(c("A", "B"), 50), 
                    clone_size = 1)
    
    c <- clustirr(s = s)
    
    expect_error(get_cdr3_cosine_similarity(c, c("CDR3a", "CDR3b")))
    expect_error(get_cdr3_cosine_similarity(c, "CDR3a"))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3b"))
    
    
    # test duplicate removal
    s <- data.frame(CDR3a = CDR3ab[1:100, "CDR3a"],
                    CDR3b = CDR3ab[1:100, "CDR3b"], 
                    sample = c(rep("a", 50), rep("b", 50)),
                    clone_size = 1)
    
    s <- rbind(s, s[1:20,])
    
    c <- clustirr(s = s)
    
    expect_no_error(get_cdr3_cosine_similarity(c,  c("CDR3a", "CDR3b")))
    
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3b"))
    
    # test with CDR3h-l chains
    s <- data.frame(CDR3h = CDR3ab[1:5, "CDR3b"],
                    CDR3l = CDR3ab[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = c("a", "a", "b", "c", "d"))
    
    c <- clustirr(s = s)
    
    expect_no_error(get_cdr3_cosine_similarity(c))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3h"))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3l"))
    expect_no_error(get_cdr3_cosine_similarity(c, c("CDR3h", "CDR3l")))
    expect_no_error(get_cdr3_cosine_similarity(c, c("CDR3l", "CDR3h")))
    
    # test with CDR3g-d chains
    s <- data.frame(CDR3g = CDR3ab[1:5, "CDR3b"],
                    CDR3d = CDR3ab[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = c("a", "a", "b", "c", "d"))
    
    c <- clustirr(s = s)
    
    expect_no_error(get_cdr3_cosine_similarity(c))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3d"))
    expect_no_error(get_cdr3_cosine_similarity(c, "CDR3g"))
    expect_no_error(get_cdr3_cosine_similarity(c, c("CDR3d", "CDR3g")))
    expect_no_error(get_cdr3_cosine_similarity(c, c("CDR3g", "CDR3d")))
    
    
    
    
})
