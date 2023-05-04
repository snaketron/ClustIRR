test_that("check_version works", {
 expect_error(check_version("non numeric input"),
               "version has to be numeric")
 expect_error(check_version(c(2,3,4)),
              "version has to be a single value")
 expect_error(check_version(4),
              "version has to be 1, 2 or 3")
})
