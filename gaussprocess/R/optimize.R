#' #' Define different covariance functions.
#' #'
#' #' @param x Input data, a numeric vector.
#' #' @param y Input data, a numeric vector.
#' #' @param par The parameters of the covariance function used, a named list.
#' #'
#' #' @return value of the selected covariance function.
#' #' @export
#' linear <- function(x,y,par){
#'   linear_cov(x,y,par$sigma)
#' }
#' squared_exp <- function(x,y, par){
#'   squared_exp_cov(x,y,par$l)
#' }
#' constant <- function(x,y, par){
#'   constant_cov(x,y,par$sigma0)
#' }
#' exponential <- function(x,y, par){
#'   exp_cov(x,y,par$l)
#' }
#' gamma_exp <- function(x,y,par){
#'   gamma_exp_cov(x,y,par$l,par$gamma)
#' }
#' rational <- function(x,y,par){
#'   rational_quadratic_cov(x,y,par$l,par$alpha)
#' }
#'
#' cov_list <- list(
#'   "linear"= list(
#'     name = "linear",
#'     fun = linear,
#'     parameter = c("sigma"),
#'     upper = NULL,
#'     lower = NULL
#'   ),
#'
#'   "squared_exp"= list(
#'     name = "squared_exp",
#'     fun = squared_exp,
#'     parameter = c("l"),
#'     upper = NULL,
#'     lower = 0.0001
#'   ),
#'
#'   "constant" = list(
#'     name = "constant",
#'     fun = constant,
#'     parameter = c("sigma0"),
#'     upper = NULL,
#'     lower = NULL
#'   ),
#'   "exponential" = list(
#'     name = "exponential",
#'     fun = exponential,
#'     parameter = c("l"),
#'     upper = NULL,
#'     lower = 0.0001
#'   ),
#'   "gamma_exp" = list(
#'     name = "gamma_exp",
#'     fun = gamma_exp,
#'     parameter = c("l","gamma"),
#'     upper = c(Inf, 2),
#'     lower = c(0.0001, 0.01)
#'   ),
#'   "rational" =list(
#'     name = "rational",
#'     fun = rational,
#'     parameter = c("l", "alpha"),
#'     upper = c(Inf, Inf),
#'     lower = c(0.0001, 0.0001)
#'   )
#' )
#'



#' init_cov
#'
#' @param covname Name of the selected covariance function
#' @param sigma A parameter
#' @param l A parameter
#' @param alpha A parameter
#' @param sigma0 A parameter
#' @param gamma A parameter
#'
#' @return The, with "covname", selected covariance function or if specified,
#' the value of the covariance function
#' @example
#' init_cov("linear", sigma = 2)(x = 1:10, y = sin(1:10))
#' t <- init_cov("rational", l = 1, alpha = 0.5)
#' t(x = 1:5, y = exp(1:5))
#' @export
init_cov <- function(covname, sigma=0, l=1, alpha=1, sigma0=1, gamma=1){
  force(sigma)
  force(l)
  force(alpha)
  force(sigma0)
  force(gamma)
  l <- list(sigma = sigma, l = l, alpha = alpha, sigma0 = sigma0, gamma = gamma)
  cov_item <- cov_list[[covname]]
  par <- lapply(cov_item$parameter, function(x) l[[x]])
  names(par)<- cov_item$parameter
  function(x,y) {
    cov_item[["fun"]](x,y, par)
  }
}
#init_cov("linear")(1:10, sin(1:10))


#' optimize_parameters
#'
#' @param X R6 object of the class "gp".
#'
#' @return List of values of parameters of each covariance function that minimize
#' the negative log marginal likelihood.
#' @example
#' T <- gp$new()
#' T$add_data(x = 1:10, y = sin(1:10))
#' T$optim_parameter()
#' @export
#'
#' @importFrom stats nlminb
optimize_parameters <- function(X){
  #------initalizing start values ----------
  min <- Inf
  best_parameters <- list(rep(1,X$get_input_dim()), l = 1, gamma = 1, alpha = 1, sigma0 = 1)
  best_method <- NULL
  solution <- list()
  start_parameters <- X$get_parameter()

  #-----optimizing parameters for the different cov-functions-----
  for(item in cov_list){
    X$set_cov(item$name)
    par_item <- sapply(item$parameter, function(x) start_parameters[[x]])
    opt <- nlminb(
      start=par_item,
      objective=help_opt,
      X=X,
      par_names = item$parameter,
      lower = item$lower,
      upper = item$upper)
    names(opt$par) <- item$parameter
    l<-list(opt$par, opt$objective)
    solution <- c(solution,list(l) )
    if(min>opt$objective){
      min<-opt$objective
      for(element in item$parameter){
        if(element == "sigma")
          best_parameters[["sigma"]] <- unname(opt$par)
        else
          best_parameters[[element]] <- opt$par[element]
      }
      best_method <- item$name
    }
  }
  names(solution)<- names(cov_list)
  print(best_method)
  #-----setting optimal parameters---------
  sigma <- best_parameters[["sigma"]]
  gamma <- best_parameters[["gamma"]]
  alpha <- best_parameters[["alpha"]]
  l <- best_parameters[["l"]]
  sigma0 <- best_parameters[["sigma0"]]
  X$set_parameter(sigma = sigma, gamma = gamma, l = l, alpha = alpha, sigma0 = sigma0)
  X$set_cov(best_method)
  return(solution)
}



#' help_opt
#' Helperfunction which return the negative log likelihood for the selected parameter.
#'
#' @param par_value Value of the parameters, named list.
#' @param X X R6 object of the class "gp".
#' @param par_names Names of the parameters, character vecotr.
#'
#' @return Negative log likelihood function.
#' @example
#' T <- gp$new()
#' T$add_data(x = 1:10, y = sin(1:10))
#' help_opt(par_values = list(sigma = 0, l = 1, gamma = 1, alpha = 1, sigma0 = 1), T, par_names = c("sigma"))
#' @export
help_opt <- function(par_value, X, par_names){
  if("sigma" %in% par_names){
    X$set_parameter(sigma = par_value)
  } else {
    names(par_value) <- par_names
    X$set_parameter(l= par_value["l"],
                    alpha = par_value["alpha"],
                    sigma0 = par_value["sigma0"],
                    gamma = par_value["gamma"])
  }
  tryCatch(
    -X$get_log_marginal_likelihood(),
    error = function(cond) return(-Inf))
}

