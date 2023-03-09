parameter_check_test <- function(...){
  
  l <- list(...)
  
  l$cdr3_sequences = 4
  
  l$save_results = TRUE
  
  return(l)
  
}

cdr3_sequences <- 2
result_folder <- 'profiles'

p <- parameter_check_test(cdr3_sequences=cdr3_sequences,
                    result_folder=result_folder)
print(p)

p$cdr3_sequences


rm(cdr3_sequences, p, result_folder)
   