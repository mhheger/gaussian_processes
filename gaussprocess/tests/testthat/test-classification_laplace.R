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


# test for diag_block_matrix ----------------------------------------------

test_that("Check for unallowed list entries", {
  m_list <- list(
    matrix(1:5, nrow = 5 ),
    diag(10),
    diag(c(12,23,2))
    )
  x <-
  expect_error(diag_block_matrix(m_list), regex = "quadratic")
  m_list[4] <- "a"
  m_list[[1]] <- diag(2)
  expect_error(diag_block_matrix(m_list), regex = "matrix")
  m_list[[4]] <- 1:3
  expect_error(diag_block_matrix(m_list), regex = "non-matrix")
}
          )


# test for pred_mc_laplace -----------------------------------------------

test_that("Checks for unallowed inputs",{
  K_list <- list(
    matrix(1:5, nrow = 5 ),
    diag(3),
    diag(c(12,23,2))
  )
  X_learn <- list(1,2,3)
  cov_list <- list(
    function(x) 1,
    function(x) 2,
    function(x) 3
  )
  y <- numeric(9)
  f_mode <- numeric(9)
  x_input <- 1

  expect_error(pred_mc_laplace(X_learn, y, K_list, f_mode, cov_list, x_input ),
               regex = "matrices")
  K_list[[1]] <- diag(3)
  expect_error(pred_mc_laplace(unlist(X_learn), y, K_list, f_mode, cov_list, x_input ),
               regex = "list")
  expect_error(pred_mc_laplace(X_learn, y, K_list, f_mode[1:6], cov_list, x_input ),
               regex = "Length")
  expect_error(pred_mc_laplace(X_learn, y, K_list, f_mode, cov_list, c(x_input,1) ),
               regex = "length")
  cov_list[[2]] <- 2
  expect_error(pred_mc_laplace(X_learn, y, K_list, f_mode, cov_list, x_input ),
               regex = "closure")
})

test_that("Checking if code runs without error", {
  x <- 1:50
  y <-as.integer(x<5)
  y <- c ( y, as.integer(x>5))
  X_learn <- as.list(x)
  covariance_list <- list(
    function(x,y) exp(- 0.5*sqrt(sum((x-y)^2))),
    function(x,y) exp(- 0.5*sqrt(sum((x-y)^2)))
  )

  K <- list(
    cov_cross(X_learn, X_learn,covariance_list[[1]]),
    cov_cross(X_learn, X_learn,covariance_list[[2]])
  )

  f_mode <- find_mode_mc_laplace(K,y)$mode

  expect_no_condition(pred_mc_laplace(X_learn,y, K, f_mode, covariance_list, 20))

})


# testing class methods gp_classification ---------------------------------
test_that("initialization has no problems", {
  expect_no_condition(gp_classification$new(n=2))
  expect_no_condition(gp_classification$new(list("gamma_exp", "linear"), n= 2))
  expect_error(gp_classification$new(list("linear"), n=3), regexp = "Length")
})

test_that("Adding data", {
  p <- gp_classification$new(n= 2)
  expect_no_condition(p$add_data(list(1,2,3), c(1,1,1,0,0,0)))
  expect_error(p$add_data(list(1,2,3), c(1,2,3)), regex = "Length")
  expect_error(p$add_data(list(1,2,3), c(1,2,2,2,2,1)), regex = "label")
  expect_error(p$add_data(list(1,2,3,4), c(1,1,1,1,1,1), regex = "Length"))
  expect_error(p$add_data(list("a",2,3), c(1,1,1,1,1,1), regex = "non-numeric"))
})

test_that("Getting predictions", {
  mcgp <- gp_classification$new(n = 2)
  y <- c(
    c(0,1),       #first point belongs to the first class, second point doesn't
    c(1,0)        #second point belongs to the first class, first point doesn't
  )

  X_learn <- list(
    c(1,2,3),
    c(-1,-2,-3)
  )

  mcgp$add_data(X_learn, y)

  expect_no_condition(mcgp$get_prediction(c(1,2,4)))
  expect_error(mcgp$get_prediction(c(1,2,4,5)), regex = "Length")
  expect_error(mcgp$get_prediction(c(1,2,"a")), regex = "numeric")

})

test_that("Setting parameters", {
  mcgp <- gp_classification$new(n = 2)
  y <- c(
    c(0,1),       #first point belongs to the first class, second point doesn't
    c(1,0)        #second point belongs to the first class, first point doesn't
  )

  X_learn <- list(
    c(1,2,3),
    c(-1,-2,-3)
  )

  mcgp$add_data(X_learn, y)
  expect_no_condition(mcgp$set_parameters(1, "linear", list(sigma = c(1,1,1))))
  expect_error(mcgp$set_parameters(3, "linear", list(sigma = c(1,1,1))), regex = "index")
  expect_error(mcgp$set_parameters(2, "lir", list(sigma = c(1,1,1))), regex = "cov")
  expect_error(mcgp$set_parameters(c(1,2), "linear", list(sigma = c(1,1,1))), regex = "length")
})

