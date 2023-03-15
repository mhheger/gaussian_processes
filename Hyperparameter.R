# define different covariance models
linear <- function(x1,x2,theta){
  linear_cov(x1,x2,theta[["sigma"]])
}
squared_exp <- function(x1,x2, theta){
  squared_exp_cov(x1,x2,theta[["l"]])
}
constant <- function(x1,x2, theta) {
  constant_cov(x1, x2, theta[["sigma0"]])
}
exponential <- function(x1,x2, theta){
  exp_cov(x1,x2,theta[["l"]])
}
gamma_exp <- function(x1,x2, theta){
  gamma_exp_cov(x1,x2,theta[["l"]],theta[["gamma"]])
}
rational <- function(x1,x2, theta){
  rational_quadratic_cov(x1,x2,theta[["l"]],theta[["alpha"]])
}

# safe models in list 
cov_models <- list(
  
  "constant" = list(
    name = "constant",
    func = constant,
    hy_par = c("sigma0"),
    upper = NULL,
    lower = NULL
  ),
  "linear"= list(
    name = "linear",
    func = linear,
    hy_par = c("sigma"),
    upper = NULL,
    lower = NULL
  ),
  "squared_exp"= list(
    name = "squared_exp",
    func = squared_exp,
    hy_par = c("l"),
    upper = NULL,
    lower = 0.0001
  ),
  
  "exponential" = list(
    name = "exponential",
    func = exponential,
    hy_par = c("l"),
    upper = NULL,
    lower = 0.0001
  ),
  "gamma_exp" = list(
    name = "gamma_exp",
    func = gamma_exp,
    hy_par = c("l","gamma"),
    upper = c(Inf, 2),
    lower = c(0.001, 0.01)
  ),
  "rational" =list(
    name = "rational",
    func = rational,
    hy_par = c("l", "alpha"),
    upper = c(Inf, Inf),
    lower = c(0.001, 0.001)
  )
)

# define additional functions

init_cov <- function(cov_name, sigma0=0, sigma=1, l=1, gamma = 1, alpha=1){
  force(sigma0)
  force(sigma)
  force(l)
  force(gamma)
  force(alpha)
  params <- c(sigma0 = sigma0, sigma = sigma, l = l, gamma = gamma, alpha = alpha)
  #print(params)
  covariance <- cov_models[[cov_name]]
  #print(covariance)
  par <- lapply(covariance$hy_par, function(x) params[[x]])
  names(par)<- covariance$hy_par
  function(x,y) {
    covariance[["func"]](x,y, par)
  }
} 


parameter_optimization <- function(X) {
  
  # define some initial values
  min <- Inf
  solution <- list()
  best_parameters <- list(sigma0 = 1, sigma = 0, l = 1, gamma = 1, alpha = 1)
  best_method <- NULL
  initial_parameters <- X$get_parameter()
  for (i in cov_models) {
    X$set_cov(i$name)
    theta_i <- sapply(i$hy_par, function(x) initial_parameters[[x]])
    print(theta_i[[1]])
    hyper_optim <- nlminb(
      start = theta_i,
      objective = neg_llikelihood,
      X = X,
      parameter_names = i$hy_par,
      lower = i$lower,
      upper = i$upper
    )
    names(hyper_optim$par) <- i$hy_par
    lst <- list(hyper_optim$par, hyper_optim$objective)
    solution <- c(solution,list(lst))
    if(min > hyper_optim$objective){
      min <- hyper_optim$objective
      for (j in i$hy_par){
        if (j == "sigma"){
          best_parameters[["sigma"]] <- unname(hyper_optim$par)
        } else {
          best_parameters[[j]] <- hyper_optim$par[j]
        }
      }
      best_method <- i$name
    }
  }
  names(solution) <- names(cov_models)
  print(best_method)
  sigma0 <- best_parameters[["sigma0"]]
  sigma <- best_parameters[["sigma"]]
  l <- best_parameters[["l"]]
  gamma <- best_parameters[["gamma"]]
  alpha <- best_parameters[["alpha"]]
  X$set_parameter(sigma = sigma, gamma = gamma, l = l, alpha = alpha, sigma0 = sigma0)
  X$set_cov(best_method)
  return(solution)
}

# get negative log marginal likelihood, computed with algorithm 2.1 
neg_llikelihood <- function(parameter, X, parameter_names){
  sigma0 <- NULL
  sigma <- NULL
  l <- NULL
  gamma <- NULL
  alpha <- NULL
  if ("l" %in% parameter_names) l <- parameter[[1]]
  if (length(parameter_names) < 2) {
    if ("sigma0" %in% parameter_names) sigma0 <- parameter[[1]]
    if ("sigma" %in% parameter_names) sigma <- parameter[[1]]
  } else {
    if ("gamma" %in% parameter_names) gamma <- parameter[[2]]
    if ("alpha" %in% parameter_names) alpha <- parameter[[2]]
  }
  X$set_parameter(sigma0,
                  sigma,
                  l,
                  gamma,
                  alpha)
  tryCatch(
    -X$get_log_marginal_likelihood(),
    error = function(cond) return(-Inf))
}
