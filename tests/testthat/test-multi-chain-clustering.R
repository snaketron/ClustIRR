test_that("clustering two chains - within-repertoire", {
    
    # test data
    data("CDR3ab", package = "ClustIRR")
    a1 <- data.frame(CDR3a = CDR3ab[1:100, "CDR3a"],
                     CDR3b = CDR3ab[1:100, "CDR3b"],
                     clone_size = 1,
                     sample = "a")
    a2 <- data.frame(CDR3a = CDR3ab[401:500, "CDR3a"],
                     CDR3b = CDR3ab[401:500, "CDR3b"],
                     clone_size = 1,
                     sample = "a")
    
    # insert clones with two identical chains
    cdr3a <- "CAVEYSGAGSYQLTF"
    cdr3b <- "CASSQEFGSSYEQYF"
    a1_cdr3ab <- data.frame(CDR3a = cdr3a,
                            CDR3b = cdr3b,
                            clone_size = 135,
                            sample = "a")
    a2_cdr3ab <- data.frame(CDR3a = cdr3a,
                            CDR3b = cdr3b,
                            clone_size = 230,
                            sample = "a")
    
    # insert big clones that perfectly overlap on cdr3a chain 
    # set cdr3b of the clones to non-overlapping sequences
    cdr3a_only <- "CASSSSSSSSSSSSF"
    a1_cdr3a <- data.frame(CDR3a = cdr3a_only,
                           CDR3b = "CAGGGGGGGGGGGGF",
                           clone_size = 500,
                           sample = "a")
    a2_cdr3a <- data.frame(CDR3a = cdr3a_only,
                           CDR3b = "CAQQQQQQQQQQQQF",
                           clone_size = 300,
                           sample = "a")
    
    # insert big clones that perfectly overlap on cdr3b chain 
    # set cdr3a of the clones to non-overlapping sequences
    cdr3b_only <- "CAYYYYYYYYYYYYF"
    a1_cdr3b <- data.frame(CDR3a = "CAAAAAAAAAAAAAF",
                           CDR3b = cdr3b_only,
                           clone_size = 150,
                           sample = "a")
    a2_cdr3b <- data.frame(CDR3a = "CAEEEEEEEEEEEEF",
                           CDR3b = cdr3b_only,
                           clone_size = 200,
                           sample = "a")
    
    s <- rbind(a1, a1_cdr3ab, a1_cdr3b, a1_cdr3a, 
               a2, a2_cdr3ab, a2_cdr3b, a2_cdr3a)
    
    # run ClustIRR analysis
    c <- clustirr(s = s)
    
    # check nodes
    df_v <- igraph::as_data_frame(c$graph, what = "vertices")
    
    # check that the 2 nodes for each test pair exist in the graph
    cdr3ab_n <- which(df_v$CDR3b == cdr3b & df_v$CDR3a == cdr3a,)
    cdr3a_n <- which(df_v$CDR3a == cdr3a_only,)
    cdr3b_n <- which(df_v$CDR3b == cdr3b_only,)
    
    # there should be two clones for each test
    expect_true(length(cdr3ab_n) == 2)
    expect_true(length(cdr3a_n) == 2)
    expect_true(length(cdr3b_n) == 2)
    
    # check edges
    df_e <- igraph::as_data_frame(c$graph, what = "edges")
    
    cdr3ab_e_n1 <- df_v$name[cdr3ab_n[1]]
    cdr3ab_e_n2 <- df_v$name[cdr3ab_n[2]]
    cdr3ab_e <- df_e[df_e$from == cdr3ab_e_n1 & df_e$to == cdr3ab_e_n2 |
                         df_e$from == cdr3ab_e_n2 & df_e$to == cdr3ab_e_n1,]
    
    cdr3a_e_n1 <- df_v$name[cdr3a_n[1]]
    cdr3a_e_n2 <- df_v$name[cdr3a_n[2]]
    cdr3a_e <- df_e[df_e$from == cdr3a_e_n1 & df_e$to == cdr3a_e_n2 |
                         df_e$from == cdr3a_e_n2 & df_e$to == cdr3a_e_n1,]
    
    cdr3b_e_n1 <- df_v$name[cdr3b_n[1]]
    cdr3b_e_n2 <- df_v$name[cdr3b_n[2]]
    cdr3b_e <- df_e[df_e$from == cdr3b_e_n1 & df_e$to == cdr3b_e_n2 |
                        df_e$from == cdr3b_e_n2 & df_e$to == cdr3b_e_n1,]
    
    # there should be two edges between the clones with two identical chains
    # and one edge between each clone pair with one identical chain
    expect_true(nrow(cdr3ab_e) == 2)
    expect_true(nrow(cdr3a_e) == 1)
    expect_true(nrow(cdr3b_e) == 1)
    
    
})


test_that("clustering two chains - between-repertoire", {
    
    # test data
    data("CDR3ab", package = "ClustIRR")
    a <- data.frame(CDR3a = CDR3ab[1:100, "CDR3a"],
                     CDR3b = CDR3ab[1:100, "CDR3b"],
                     clone_size = 1,
                     sample = "a")
    b <- data.frame(CDR3a = CDR3ab[401:500, "CDR3a"],
                     CDR3b = CDR3ab[401:500, "CDR3b"],
                     clone_size = 1,
                     sample = "b")
    
    # insert clones with two identical chains
    cdr3a <- "CAVEYSGAGSYQLTF"
    cdr3b <- "CASSQEFGSSYEQYF"
    a_cdr3ab <- data.frame(CDR3a = cdr3a,
                            CDR3b = cdr3b,
                            clone_size = 135,
                            sample = "a")
    b_cdr3ab <- data.frame(CDR3a = cdr3a,
                            CDR3b = cdr3b,
                            clone_size = 230,
                            sample = "b")
    
    # insert big clones that perfectly overlap on cdr3a chain 
    # set cdr3b of the clones to non-overlapping sequences
    cdr3a_only <- "CASSSSSSSSSSSSF"
    a_cdr3a <- data.frame(CDR3a = cdr3a_only,
                           CDR3b = "CAGGGGGGGGGGGGF",
                           clone_size = 500,
                           sample = "a")
    b_cdr3a <- data.frame(CDR3a = cdr3a_only,
                           CDR3b = "CAQQQQQQQQQQQQF",
                           clone_size = 300,
                           sample = "b")
    
    # insert big clones that perfectly overlap on cdr3b chain 
    # set cdr3a of the clones to non-overlapping sequences
    cdr3b_only <- "CAYYYYYYYYYYYYF"
    a_cdr3b <- data.frame(CDR3a = "CAAAAAAAAAAAAAF",
                           CDR3b = cdr3b_only,
                           clone_size = 150,
                           sample = "a")
    b_cdr3b <- data.frame(CDR3a = "CAEEEEEEEEEEEEF",
                           CDR3b = cdr3b_only,
                           clone_size = 200,
                           sample = "b")
    
    s <- rbind(a, a_cdr3ab, a_cdr3b, a_cdr3a, b, b_cdr3ab, b_cdr3b, b_cdr3a)
    
    # run ClustIRR analysis
    c <- clustirr(s = s)
    
    # check nodes
    df_v <- igraph::as_data_frame(c$graph, what = "vertices")
    
    # check that the 2 nodes for each test pair exist in the graph
    cdr3ab_n <- which(df_v$CDR3b == cdr3b & df_v$CDR3a == cdr3a,)
    cdr3a_n <- which(df_v$CDR3a == cdr3a_only,)
    cdr3b_n <- which(df_v$CDR3b == cdr3b_only,)
    
    # there should be two clones for each test
    expect_true(length(cdr3ab_n) == 2)
    expect_true(length(cdr3a_n) == 2)
    expect_true(length(cdr3b_n) == 2)
    
    # check edges
    df_e <- igraph::as_data_frame(c$graph, what = "edges")
    
    cdr3ab_e_n1 <- df_v$name[cdr3ab_n[1]]
    cdr3ab_e_n2 <- df_v$name[cdr3ab_n[2]]
    cdr3ab_e <- df_e[df_e$from == cdr3ab_e_n1 & df_e$to == cdr3ab_e_n2 |
                         df_e$from == cdr3ab_e_n2 & df_e$to == cdr3ab_e_n1,]
    
    cdr3a_e_n1 <- df_v$name[cdr3a_n[1]]
    cdr3a_e_n2 <- df_v$name[cdr3a_n[2]]
    cdr3a_e <- df_e[df_e$from == cdr3a_e_n1 & df_e$to == cdr3a_e_n2 |
                        df_e$from == cdr3a_e_n2 & df_e$to == cdr3a_e_n1,]
    
    cdr3b_e_n1 <- df_v$name[cdr3b_n[1]]
    cdr3b_e_n2 <- df_v$name[cdr3b_n[2]]
    cdr3b_e <- df_e[df_e$from == cdr3b_e_n1 & df_e$to == cdr3b_e_n2 |
                        df_e$from == cdr3b_e_n2 & df_e$to == cdr3b_e_n1,]
    
    # there should be two edges between the clones with two identical chains
    # and one edge between each clone pair with one identical chain
    expect_true(nrow(cdr3ab_e) == 2)
    expect_true(nrow(cdr3a_e) == 1)
    expect_true(nrow(cdr3b_e) == 1)
    
    
})