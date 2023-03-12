#----tests for find_mode_laplace---------------

test_that("unallowed matrix inputs", {
  expect_error(find_mode_laplace(matrix(1:6, nrow= 2), c(1,1,1), "probit"),
               regexp = "square")
  expect_error(find_mode_laplace(c("a","b","c"), c(1,1,1), "probit"),
               regexp = "numeric")
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), c(1,1,1,1), "probit"),
               regexp = "dimension")
  expect_error(find_mode_laplace(tibble::tibble(x=1:3, y=1:3), c(1,1,1,1), "probit"),
               regexp = "matrix")
})

test_that("unallowed y inputs", {
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), list("a",2,"b"), "probit"),
               regexp = "numeric")
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), c(NA, NA, NA), "probit"),
               regexp = "numeric")
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), matrix(rep(1,6), nrow=2), "probit"),
               regexp = "dimension")
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), c(1,2,3), "probit"),
               regexp = "1")
})

test_that("unallowed function name", {
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), c(1,-1,1), "wischmopp"),
               regexp = "logit")
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), c(1,-1,1), function(y,f) 1),
               regexp = "character")
  expect_error(find_mode_laplace(matrix(1:9, nrow=3), c(1,-1,1), c("logit", "logit")),
               regexp = "character")
})

#----tests for pred_laplace----------
test_that("unallowed mode inputs", {
  expect_error(pred_laplace(c(1,2,3), c(1,2), c(1,1), function(x,y) exp(-abs(x*y)), "probit", 1))
  expect_error(pred_laplace(c("a","b"), c(1,2), c(1,1), function(x,y) exp(-abs(x*y)), "probit", 1),
               regex = "numeric")
})

test_that("unallowed X_learn inputs", {
  expect_error(pred_laplace(c(1,2), c("a","b"), c(1,1), function(x,y) exp(-abs(x*y)), "probit", 1),
               regex = "numeric")
  expect_error(pred_laplace(c(1,2), c(1,2,3), c(1,1), function(x,y) exp(-abs(x*y)), "probit", 1),
               regex = "dimension")
})

test_that("unallowed y_learn inputs", {
  expect_error(pred_laplace(c(1,2), c(1,1), c("a","b"), function(x,y) exp(-abs(x*y)), "probit", 1),
               regex = "1")
  expect_error(pred_laplace(c(1,2), c(1,1), c(1,2,3), function(x,y) exp(-abs(x*y)), "probit", 1),
               regex = "dimension|1")
})
test_that("unallowed covariance-function", {
  expect_error(pred_laplace(c(1,2), c(1,2), c(1,1), "function(x,y) exp(-abs(x*y))", "probit", 1),
               regex = "function")
  expect_error(pred_laplace(c(1,2), c(1,2), c(1,1), 1, "probit", 1),
               regex = "function")
})
test_that("unallowed function name", {
  expect_error(pred_laplace(c(1,2), c(1,1), c(1,2), function(x,y) exp(-abs(x*y)), "prowin", 1),
               regex = "probit")
  expect_error(pred_laplace(c(1,2), c(1,1), c(1,2), function(x,y) exp(-abs(x*y)), c("probit","logit"), 1),
               regex = "character")
})

#----tests for predict_laplace----------
test_that("unallowed X_learn inputs", {
  expect_error(predict_laplace(c("a","b"), c(1,1), function(x,y) exp(-abs(x*y)), "probit", 1))
  expect_error(predict_laplace(c(1,2,3), c(1,1), function(x,y) exp(-abs(x*y)), "probit", 1),
               regex = "dimension")
})


# tests for find_mode_mc_laplace -------------------------------------------
test_that("unallowed labels", {
  K_list = list(
    matrix(1:4, nrow = 2),
    diag(c(1,2))
  )
  expect_error(find_mode_mc_laplace(K_list, c(2,1)), regex = "1 or 0")
})

# tests to fullfill
test_that("unallowed matrix list", {
  K_list = list(
    matrix(1:4, nrow = 2),
    diag(c(1,2,3))
  )
  expect_error(find_mode_mc_laplace(K_list, c(0,1)), regex = "dimension")
})




