library(shiny)
# testing manual data input -----------------------------------------------
test_that("Manual input leads to correct datatable" , {
  obj =NULL
  testServer(testing_server, {
    session$setInputs(
      input_dim = 2,
      cov = "squared_exp",
      l_sqr = 1,
      x1 = 1,
      x2 = 2,
      y = 3,
      noise = 1,
      x1_pred = 1,
      x2_pred = 2
    )
    session$setInputs(
      set_dim = 1,
      add_data = 1,
      get_pred =1
    )
    expect_equal(X_data(), data.frame(x1=1, x2=2))
    expect_equal(y_data(), 3)
  })
})


# testing work of prediction ----------------------------------------------
test_that("prediction works correct",{
  g <-new.gp() %>%
    add_data(c(1,2),3,1) %>%
    get_prediction(c(1,2))
  testServer(testing_server, {
    session$setInputs(
      input_dim = 2,
      cov = "squared_exp",
      l_sqr = 1,
      x1 = 1,
      x2 = 2,
      y = 3,
      noise = 1,
      x1_pred = 1,
      x2_pred = 2
    )
    session$setInputs(
      set_dim = 1,
      add_data = 1,
      get_pred =1
    )
    expect_equal(prediction(), g )

  })
})


# testing handling of false inputs  ---------------------------------------
test_that("False input do not cause crash", {
  testServer(testing_server, {
    session$setInputs(
      input_dim = 2,
      cov = "squared_exp",
      l_sqr = 1,
      x1 = 1,
      x2 = 2,
      y = 3,
      noise = 1,
      x1_pred = 1,
      x2_pred = 2
    )
    session$setInputs(
      set_dim = 1,
      add_data = 1,
      get_pred =1
    )
    session$setInputs(
      x1=NA
    )
    expect_no_condition(session$setInputs(add_data=1))
    session$setInputs(
      x1 = 1,
      noise = NA
    )
    expect_no_condition(session$setInputs(add_data=1))
    session$setInputs(
      noise = -1
    )
    expect_no_condition(session$setInputs(add_data=1))
    session$setInputs(
      x1 = 1,
      noise = 1,
      x2_pred = NA
    )
    expect_no_condition(session$setInputs(get_pred=1))
    session$setInputs(
      x1 = 1,
      noise = 1,
      l_sqr = -1,
      l_rat = NA,
      l_gamma = NA,
      l_exp = -3,
      alpha = -1,
      gamma = 4
    )
    expect_no_condition(session$setInputs(set_parameters=1))
    session$setInputs(
      input_dim = -12
    )
    expect_no_condition(session$setInputs(set_dim=1))
    expect_no_condition(session$setInputs(get_pred=1))
  })
})



