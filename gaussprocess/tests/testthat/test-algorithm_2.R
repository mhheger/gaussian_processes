test_that("non_numerical input", {
  expect_error(predict_gauss(c("a","b"), c(1,2), function(x,y) 1, 1, 1),
               regex= "non_numeric")
  expect_error(predict_gauss(c(1,2),c("a","b"), function(x,y) 1, 1, 1),
               regex= "non_numeric")
  expect_error(predict_gauss(c(1,2),c(1,2), function(x,y) 1, "a", 1),
               regex= "non_numeric")
  expect_error(predict_gauss(c(1,2),c(1,2), function(x,y) 1, 1, "a"),
               regex= "non_numeric")
})


test_that("input data doesn't fit", {
  expect_error(predict_gauss(c(1), c(1,2), function(x,y) 1, 1, 1), regex = "dimension")
})
