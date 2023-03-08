
# Testing init_cov --------------------------------------------------------
test_that("Handle non numeric inputs", {
  expect_error(init_cov("linear", sigma = list(1,2,3,3)), regex = "numeric")
  expect_error(init_cov("linear", sigma = 12, alpha = 23, gamma = "abba"), regex = "numeric")
})

test_that("Handle false covariance functions", {
  expect_error(
    init_cov("staubsauger", sigma = 0, l = 1, sigma0 = 12),
    regex = "known covariance-function"
  )
  expect_error(
    init_cov(linear, sigma = 0, l = 1, sigma0 = 12),
    regex = "character"
  )
  expect_error(
    init_cov(function(x,y) exp(-x*y), sigma = 0),
    regex = "character"
  )
})


# Testing new.gp ----------------------------------------------------------
test_that("Handling false input for initializing gp object",{
  expect_error(new.gp("staubsauger", cov_fun = "squared_exp", response_fun = "logit"),
               regex = "not one of")
  expect_error(new.gp("regression", cov_fun = "blub", response_fun = "logit"),
               regex = "not one of")
  expect_error(new.gp("regression", cov_fun = "linear", response_fun = "logisch"),
               regex = "not one of")
  expect_error(suppressWarnings(new.gp("regression", cov_fun = function(x){x}, response_fun = "logit"),
               regex = "character"))
  expect_error(new.gp(c("regression", "classification"), cov_fun = "linear", response_fun = "logit"),
               regex = "character")
  expect_error(new.gp(NULL, cov_fun = "linear", response_fun = "logit"),
               regex = "NULL")
})


# Testing add_data --------------------------------------------------------
##Testing error-handling of unallowed inputs

test_that("Handling unallowed inputs", {
  expect_error(add_data(list(), 1:10, 1:10, 0.1),
               regex = "gp")
  expect_error(add_data(NULL, 1:10, 1:10, 0.1),
               regex = "gp")
  expect_error(add_data(new.gp(), "a", 1:10, 0.1),
               regex = "numeric")
  expect_error(add_data(new.gp(), NULL, 1:10, 0.1),
               regex = "numeric")
  expect_error(add_data(new.gp(), 1:10, NULL, 0.1),
               regex = "numeric")
  expect_error(add_data(new.gp(), 1:12, 1:5, 0.1),
               regex = "dimension")
  expect_error(add_data(new.gp(), data.frame(x=1:10, y = 1:10), "f"),
               regex = "name")
  expect_error(add_data(new.gp(), data.frame(x=letters[1:10], y = 1:10), "y"),
               regex = "numeric")
})

##Testing if different ways of input lead to similiar results
test_that("Different ways of input", {
  expect_identical(add_data(new.gp(), data.frame(x=1:10, y = 1:10), "y"),
                   add_data(new.gp(), 1:10, 1:10)
                                )
  expect_identical(add_data(new.gp(), data.frame(x=1:10, y = 1:10), 1:10),
                   add_data(new.gp(), matrix(c(1:10, 1:10), nrow=2), 1:10))
  xx <- lapply(1:10, function(x) c(x,x))
  expect_identical(add_data(new.gp(), data.frame(x=1:10, y = 1:10), 1:10),
                   add_data(new.gp(), xx, 1:10))
  expect_identical(add_data(new.gp(), data.frame(x=1:10, y = 1:10), 1:10),
                   add_data(new.gp(), matrix(unlist(xx), ncol=2), 1:10))
})


# Testing get_cov_name ----------------------------------------------------

test_that("Testing error handling", {
  expect_error(get_cov_name(12), regex= "gp")
  expect_error(get_cov_name(list("Bla")), regex= "gp")
  expect_error(get_cov_name(NULL), regex= "gp")
})

test_that("Checking for correct results", {
  p <- new.gp()
  expect_identical(get_cov_name(p),"squared_exp")
  p <- new.gp(cov_fun = "linear")
  expect_identical(get_cov_name(p),"linear")
})


# Testing get_K -----------------------------------------------------------
test_that("Testing error handling", {
  expect_error(get_K(12), regex= "gp")
  expect_error(get_K(list("Bla")), regex= "gp")
  expect_error(get_K(NULL), regex= "gp")
})


# Testing get_cov ---------------------------------------------------------
test_that("Testing error handling", {
  expect_error(get_cov(12), regex= "gp")
  expect_error(get_cov(list("Bla")), regex= "gp")
  expect_error(get_cov(NULL), regex= "gp")
})

test_that("Checking for correct results", {
  p <- new.gp() %>%
    add_data(1:10, 1:10)
  expect_identical(get_cov(p)$name,"squared_exp")
  expect_identical(get_cov(p)$cov_fun(1,2), exp(-0.5*(1-2)^2))
})


# Testing get_data --------------------------------------------------------
test_that("Testing error handling", {
  expect_error(get_data(12), regex= "gp")
  expect_error(get_data(list("Bla")), regex= "gp")
  expect_error(get_data(NULL), regex= "gp")
})

test_that("Testing if code runs without errors", {
  X_learn <- matrix(1:100, nrow = 25)
  y <- 1:4
  expect_no_condition(p <- new.gp() %>% add_data(X_learn, y) %>% get_data())
})



# Testing get_parameter  -------------------------------------------------
test_that("Testing if output is correct", {
  X_learn <- matrix(1:100, nrow = 25)
  y <- 1:4
  p <- new.gp() %>% add_data(X_learn, y)
  l1 <- list(sigma = rep(0,25), l = 1, alpha = 1, sigma0 = 1, gamma = 1)
  expect_identical(get_parameter(p),l1)
})


# Testing get_mean_fun ----------------------------------------------------
test_that("Testing if output is correct", {
  p <- new.gp()
  add_data(p,12,2)
  p$set_mean_fun(function(x) return(x^2))
  expect_identical(get_mean_fun(p),function(x) return(x^2))
  p$set_mean_fun(10)
  expect_identical(get_mean_fun(p)(12) ,(function(x) return(10))(12))
})


# Testing get_noise -------------------------------------------------------
test_that("Testing error handling", {
  expect_error(get_noise(12), regex= "gp")
  expect_error(get_noise(list("Bla")), regex= "gp")
  expect_error(get_noise(NULL), regex= "gp")
})

test_that("Testing if output is correct", {
  p <- new.gp()
  p$set_noise(12)
  expect_identical(get_noise(p),12)
})


# Testing get_log_marginal_likelihood -------------------------------------
test_that("Testing error handling", {
  expect_error(get_log_marginal_likelihood(12), regex= "gp")
  expect_error(get_log_marginal_likelihood(list("Bla")), regex= "gp")
  expect_error(get_log_marginal_likelihood(NULL), regex= "gp")
})

test_that("Testing if code runs without errors", {
  X_learn <- matrix(1:100, nrow = 25)
  y <- 1:4
  expect_no_condition(p <- new.gp() %>% add_data(X_learn, y) %>% get_log_marginal_likelihood())
})


# Testing get_prediction --------------------------------------------------

test_that("Testing false input data", {
  X_learn <- matrix(1:100, nrow = 25)
  y <- 1:4
  p <- new.gp() %>% add_data(X_learn, y)
  expect_error(get_prediction(p, c(1,23)), regexp = "dimension")
  expect_error(get_prediction(p, LETTERS[1:25]), regexp = "numeric")
  expect_error(get_prediction(p, NULL), regexp = "numeric")
})
test_that("Testing if code runs without errors", {
  X_learn <- matrix(1:100, nrow = 25)
  y <- 1:4
  p <- new.gp() %>% add_data(X_learn, y)
  expect_no_condition(get_prediction(p, rep(1,25)))
})


# Testing set_cov ---------------------------------------------------------

test_that("Testing not allowed input values", {
  p1 <- new.gp()
  expect_error(set_cov(p1, "linear_cov"), regex= "isn't a known covariance")
  expect_error(set_cov(p1, function(x)x^2), regex = "character")
  expect_error(set_cov(p1, NULL), regex= "character")
})

test_that("Reference semantic is working correct", {
  p1 <- new.gp()
  expect_identical(set_cov(p1, "gamma_exp") %>%  get_cov_name(), "gamma_exp")
})


# Testing set_noise -------------------------------------------------------
test_that("Testing not allowed input values", {
  p1 <- new.gp()
  expect_error(set_noise(p1, "a"), regex= "double")
  expect_error(set_noise(p1, NULL), regex = "positive")
  expect_error(set_noise(p1, NA), regex= "positive")
  expect_error(set_noise(p1, c(1,2,3)), regex = "dimension")
})


# Testing set_parameters --------------------------------------------------
test_that("Testing not allowed input values", {
  p1 <- new.gp()
  add_data(p1, data.frame(x=1:10, y=21:30), y = 2:11)
  expect_error(set_parameter(p1, sigma = 1), regex = "length")
  expect_error(set_parameter(p1, l = "a"), regex ="numeric")
  expect_error(set_parameter(p1, list(l = c(1,2))), regex = "length")
})

test_that("Testing if code runs without errors", {
  p <- new.gp()
  add_data(p, 1:10,1:10)
  expect_no_condition(set_parameter(p, sigma = 12, gamma = 1.5))
  par_list <- list( sigma0 = 3, alpha = 2, l = 13)
  expect_no_condition(set_parameter(p, par_list))
  expect_identical(set_parameter(p, l= 1, sigma = 1), set_parameter(p, list(sigma = 1, l = 1)))
})


# Testing set_mean_fun ----------------------------------------------------
test_that("Testing not allowed input values", {
  p1 <- new.gp()
  add_data(p1, data.frame(x=1:10, y=21:30), y = 2:11)
  expect_error(set_mean_fun(p1,c(1,2)), regex = "length")
  expect_error(set_mean_fun(p1,function(x,c){1}), regex ="one")
  expect_error(set_mean_fun(p1,function(x) c(x,x)))
  expect_error(set_mean_fun(p1,"a"), regex = "numeric")
})





