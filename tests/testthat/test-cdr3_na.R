skip_if_not(rBLAST::has_blast(), message = "BLAST+ is not installed.")

test_that("samples with one chain only containing NA do not stop clustirr", {
    
    data("D1", package = "ClustIRR")
    s <- data.frame(CDR3b = D1$s[1:5, "CDR3b"],
                    CDR3a = D1$s[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = c(rep("a", 3), rep("b", 2)))
    
    # set 1 CDR3a chain of sample a and all CDR3a chains of sample b to NA 
    s$CDR3a[c(1, 4, 5)] <- NA
    
    # check clustering run - should run and give out warning
    expect_warning(c <- clustirr(s = s),
                   "s contains NA value")
    
    # check edges
    df_e <- igraph::as_data_frame(c$graph, what = "edges")
    
    # there should be no edges between the NA values
    expect_true(nrow(df_e) == 0)
    
})


test_that("test that rBlast tmp files are cleaned up", {
    
    data("D1", package = "ClustIRR")
    s <- data.frame(CDR3b = D1$s[1:5, "CDR3b"],
                    CDR3a = D1$s[1:5, "CDR3a"],
                    clone_size = 1,
                    sample = c(rep("a", 3), rep("b", 2)))
    
    c <- clustirr(s = s)
    
    # there should not be any temp file of rBlast left in the wd  
    expect_false(length(list.files(
        pattern = "\\.(pdb|phd|phi|phr|pin|pog|pot|psq|ptf|pto)$")
    ) > 0 )
    
})
