library(ggplot2)
library(plotly)

#plot function based on image 2.4 from Williams and Rasmussen
#' plotting a gp instance
#' @param gp gp instance
#' @param x_min start value x_axis
#' @param x_max end value x-axiy
#' @param n_points number of points generated
#' @param plotly_obj plotly object or not
#'
#'@importFrom stats dnorm quantile rnorm var integrate
#'@import plotly
#'@import ggplot2

plot_gp <- function(gp, x_min=-5, x_max=5, n_points=100, plotly_obj=FALSE){
  # Generate a grid of x values
  x_grid <- seq(x_min, x_max, length.out = n_points)
  X_learn <- unlist(gp$get_data()$X_learn)
  y_learn <- gp$get_data()$y_learn
  ranged <- x_min <= X_learn & X_learn <= x_max

  # Make predictions at each x value
  f_predict <- list()
  var_f <- list()

  for(x in x_grid){
    res <- get_prediction(gp, x)

    # Extract prediction and variance from predictions
    f_predict <- append(f_predict, res$f_predict)
    var_f <- append(var_f, res$var_f)
  }

  df <- data.frame(x = x_grid, y = unlist(f_predict),
                   var = unlist(var_f))

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    geom_ribbon(aes(ymin = y - 2*sqrt(var), ymax = y + 2*sqrt(var)), alpha = 0.2) +
    geom_point(data = data.frame(x = X_learn[ranged], y = y_learn[ranged]), aes(x = x, y = y), color = "red", size = 2) +
    xlab("Input") + ylab("Output")

  #Convert to plotly object for interactivity
  if(plotly_obj == TRUE){
    p <- ggplotly(p)
  }

  return(p)

}
#' Plotting a sample plot
#' @param gp gp instance
#' @param x_min start value x_axis
#' @param x_max end value x-axiy
#' @param n_points number of points generated
#' @param n_samples number of samples
#' @param plotly_obj plotly object or not
#'
#' @export
plot_gp_posterior <- function(gp, x_min=-5, x_max=5, n_points=100, n_samples=5, plotly_obj=FALSE){
  # Generate a grid of x values
  x_grid <- seq(x_min, x_max, length.out = n_points)
  X_learn <- unlist(gp$get_data()$X_learn)
  y_learn <- gp$get_data()$y_learn
  ranged <- x_min <= X_learn & X_learn <= x_max
  # Make predictions at each x value
  f_samples <- vector(mode = "list", length = n_points)
  var_f <- vector(mode = "list", length = n_points)

  for(i in seq_along(x_grid)){
    x <- x_grid[i]
    res <- get_prediction(gp, x)

    # Extract posterior samples and variance from predictions
    f_samples[[i]] <- rnorm(n_samples, res$f_predict, sqrt(res$var_f))
    var_f[[i]] <- res$var_f
  }

  # Create a data frame with the posterior samples for each x value
  df_samples <- data.frame(x = rep(x_grid, each = n_samples),
                           y = unlist(f_samples))

  # Create a data frame with the mean and confidence intervals of the posterior samples
  f_matrix <- do.call(rbind, f_samples)
  df_posterior <- data.frame(x = x_grid,
                             y = apply(f_matrix, 1, mean),
                             ymin = apply(f_matrix, 1, function(x) quantile(x, probs = 0.025)),
                             ymax = apply(f_matrix, 1, function(x) quantile(x, probs = 0.975)))

  # Create the plot
  p <- ggplot() +
    geom_line(data = df_posterior, aes(x = x, y = y), color = "blue") +
    geom_ribbon(data = df_posterior, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.2, fill = "green") +
    geom_point(data = data.frame(x = X_learn[ranged], y = y_learn[ranged]), aes(x = x, y = y), color = "red", size = 2) +
    geom_point(data = df_samples, aes(x = x, y = y), color = "black", alpha = 0.2, size = 0.5) +
    xlab("Input") + ylab("Output")

  #Convert to plotly object for interactivity
  if(plotly_obj == TRUE){
    p <- ggplotly(p)
  }
  return(p)
}




#example usage
#p <- new.gp(cov_fun = "squared_exp")
#add_data(p, 1:10, 1:10)
#plot_gp(p)
#plot_gp_posterior(p)
