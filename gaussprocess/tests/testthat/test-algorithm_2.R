test_that("non-numerical input", {
  expect_error(predict_gauss(c("a","b"), c(1,2), function(x,y) 1, 1, 1),
               regex= "non-numeric")
  expect_error(predict_gauss(c(1,2),c("a","b"), function(x,y) 1, 1, 1),
               regex= "non-numeric")
  expect_error(predict_gauss(c(1,2),c(1,2), function(x,y) 1, "a", 1),
               regex= "non-numeric")
  expect_error(predict_gauss(c(1,2),c(1,2), function(x,y) 1, 1, "a"),
               regex= "non-numeric")
  expect_error(predict_gauss(data.frame(x="a"), 1, function(x,y) 1, 1, 1),
               regexp = "non-numeric")
})


test_that("size of data doesn't fit", {
  expect_error(predict_gauss(c(1), c(1,2), function(x,y) 1, 1, 1),
               regex = "dimension")
  expect_error(predict_gauss(c(1), c(1), function(x,y) 1, 1, c(1,1)),
               regex = "dimension")
  expect_error(predict_gauss(c(1,2,3), c(1,2), function(x,y) 1, 1, 1),
               regex = "dimension")
  expect_error(predict_gauss(c(1), c(1), function(x,y) 1, c(1,1), c(1)),
               regex = "length")
})

test_that("wrong input type", {
  expect_error(predict_gauss(1, 1,1, 1, 1),
               regexp = "function")
})

test_that("different kind of input data", {
  expect_equal(predict_gauss(data.frame(x = 1:10), 1:10, function(x,y) 1, 1, 2.5),
               predict_gauss(1:10, 1:10, function(x,y) 1, 1, 2.5),
               tolerance = testthat_tolerance()
               )
  expect_equal(predict_gauss(data.frame(x = 1:10, y=1:10), 1:10, function(x,y) 1, 1, c(2.5,1)),
               predict_gauss(matrix(c(1:10,1:10), nrow=2, byrow= T ), 1:10, function(x,y) 1, 1, c(2.5,1)),
               tolerance = testthat_tolerance()
  )
})

test_that("input data for list converting", {
  expect_equal(convert_to_list(data.frame(x=1:10,y=1:10),n = 10),
               convert_to_list(matrix(c(1:10,1:10),nrow=2, byrow=T), n=10),
               tolerance = testthat_tolerance()
  )
  expect_equal(convert_to_list(c(1,1,2,2,3,3,4,4),n = 2),
               convert_to_list(matrix(c(1:4,1:4),nrow=2, byrow=T), n=2),
               tolerance = testthat_tolerance()
  )
  expect_error(convert_to_list(array(1:100, dim= c(10,10,10)),n = 2))
})

test_that("input for cov_cross", {
  expect_error(cov_cross(list(list(1,2)), list(c(1,2)), function(x,y) sum(x)+sum(y)))
  expect_error(cov_cross(c(1,2), c(1,2), function(x,y) sum(x)+sum(y)))
  expect_error(cov_cross(c(1,2), list(c(1,2)), function(x,y) sum(x)+sum(y)))
})

