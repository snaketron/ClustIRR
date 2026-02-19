test_that("save_interactive_graph works and cleans up locally", {
    
    # be sure not to break anything because of missing package
    skip_if_not_installed("htmlwidgets")
    
    # Dummy visNetwork
    n <- data.frame(id = 1:3)
    e <- data.frame(from = c(1,2), to = c(2,3))
    g <- visNetwork::visNetwork(n, e)
    
    file_name <- "test"
    output_folder <- tempdir()
    
    
    # Run function
    expect_message(save_interactive_graph(graph = g, 
                                          file_name = file_name, 
                                          output_folder = output_folder,
                                          overwrite = T), 
                   "exported to")
    
    expect_true(file.exists(file.path(output_folder, 
                                      paste0(file_name, ".html"))))
    

    expect_false(file.exists(paste0(file_name, ".html")))
    
    expect_error(save_interactive_graph(graph = g, 
                                        file_name = file_name, 
                                        output_folder = output_folder,
                                        overwrite = F))
})