test_that("get_blosum does not delete duplicate edges", {
    
    data("CDR3ab", package = "ClustIRR")
    s <- data.frame(CDR3b = CDR3ab[1:5, "CDR3b"],
                    CDR3a = CDR3ab[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = "a")
    
    
    # insert beta chain
    cdr3b <- s[3,]$CDR3b # CSALTPGLIYNEQFF
    cdr3a <- CDR3ab[6, "CDR3a"] # CSLGPSKSQYF
    
    # add one identical clone
    s <- rbind (s, data.frame(CDR3b = cdr3b,
                              CDR3a = cdr3a,
                              clone_size = 1,
                              sample = "a"))
    
    
    # run ClustIRR analysis
    c <- clustirr(s = s)
    
    
    df_g_e <- igraph::as_data_frame(c$graph, what = "edges") 
    
    expect_true(nrow(df_g_e)>0)
    
    
})