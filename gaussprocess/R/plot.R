source("./algorithm_2.1.R")
library(ggplot2)
library(plotly)

#plot function based on image 2.4 from Williams and Rasmussen

plot_gp <- function(X_learn, y_learn, cov, noise=1, 
                    x_min=-5, x_max=5, n_points=100){
  # Generate a grid of x values
  x_grid <- seq(x_min, x_max, length.out = n_points)
  
  # Make predictions at each x value
  f_predict <- list()
  var_f <- list()
  
  for(x in x_grid){
    res <- predict_gauss(X_learn, y_learn, cov, noise=1, x)
    
    # Extract prediction and variance from predictions
    f_predict <- append(f_predict, res$f_predict)
    var_f <- append(var_f, res$var_f)
  }
  
  df <- data.frame(x = x_grid, y = unlist(f_predict), 
                   var = unlist(var_f))
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line() +
    geom_ribbon(aes(ymin = y - sqrt(var), ymax = y + sqrt(var)), alpha = 0.2) +
    geom_point(data = data.frame(x = X_learn, y = y_learn), aes(x = x, y = y), color = "red", size = 2) +
    xlab("Input") + ylab("Output")
  
  # Convert to plotly object for interactivity
  p <- ggplotly(p)
  return(p)

}

plot_gp_posterior <- function(X_learn, y_learn, cov, noise=1,
                              x_min=-5, x_max=5, n_points=100, n_samples=5){
  # Generate a grid of x values
  x_grid <- seq(x_min, x_max, length.out = n_points)
  
  # Make predictions at each x value
  f_samples <- vector(mode = "list", length = n_points)
  var_f <- vector(mode = "list", length = n_points)
  
  for(i in seq_along(x_grid)){
    x <- x_grid[i]
    res <- predict_gauss(X_learn, y_learn, cov, noise=1, x)
    
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
    geom_ribbon(data = df_posterior, aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.2, fill = "blue") +
    geom_point(data = data.frame(x = X_learn, y = y_learn), aes(x = x, y = y), color = "red", size = 2) +
    geom_point(data = df_samples, aes(x = x, y = y), color = "blue", alpha = 0.2, size = 0.5) +
    xlab("Input") + ylab("Output")
  
  # Convert to plotly object for interactivity
  p <- ggplotly(p)
  return(p)
}




#example usage
X_learn <- seq(-2, 2, length.out = 5)
y_learn <- sin(X_learn)
cov <- function(x, y) exp(-0.5 * (x - y) ^ 2)

p1 <- plot_gp(X_learn, y_learn, cov)
p1
p2 <- plot_gp_posterior(X_learn, y_learn, cov)
p2
