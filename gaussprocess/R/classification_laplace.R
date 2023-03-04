#-----definition of the different response functions---------

logit <- list(
  value = function(y,f){
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")

    (1+exp(-t(y)%*%f))^(-1)
  },
  log_p = function(y,f) {
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")
    y <- as.vector(y)
    f <- as.vector(f)
    - log(1+ exp(-t(y)%*%f))
  },
  del_f = function(y,f) {
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")
    y <- as.vector(y)
    f <- as.vector(f)
    0.5*(y+1) - (1+ exp(-1*f))^(-1)

  },
  hess_f = function(y,f){
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")
    y <- as.vector(y)
    f <- as.vector(f)
    diag(-(1+ exp(-1*f))^(-1)*(1-(1+ exp(-1*f))^(-1)),nrow=length(y))
  }
)

probit <- list(
  value = function(y,f){
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")
    pnorm(t(y)%*%f)
  },
  log_p = function(y,f) {
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")
    y <- as.vector(y)
    f <- as.vector(f)
    log(pnorm(t(y)%*%f))
  },
  del_f = function(y,f) {
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")
    y <- as.vector(y)
    f <- as.vector(f)
    (y*dnorm(f))/pnorm(y*f)

  },
  hess_f = function(y,f){
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")
    y <- as.vector(y)
    f <- as.vector(f)
    diag(- (dnorm(f))^2/(pnorm(f*y)^2)- y*f*dnorm(f)/(pnorm(y*f)),nrow=length(y))
  }
)
response_list <- list(logit =logit, probit = probit)


#' finding value of mode function using laplace-approximation
#'
#' @param K covariance-variance-matrix of dimension (n,n)
#' @param y output-vector of length n, its values are just allowed
#' to be -1 or 1
#' @param likelihood_fun name of the response function, that should be used
#' either 'probit' or 'logit'
#'
#' @return  mode: value of the mode function
#'          log_marginal_likelihood: logarithmic marginal likelihood of the
#'          result
#'
#' @examples
#' find_mode_laplace(  K = matrix( 1:4, nrow= 2),
#'                     y = c(1,-1),
#'                     likelihood_fun = "probit"
#'                     )
find_mode_laplace <- function(K, y, likelihood_fun){
  if(!is.matrix(K)) stop("K has to be a numerical matrix")
  if(! nrow(K)==ncol(K)) stop("K has to be a square matrix")
  if(!nrow(K)==length(as.vector(y))) stop("dimension of K and y do not fit")
  if(!is.numeric(K)|!is.numeric(y)) stop("input data has to be numerical")
  if(length(likelihood_fun)!=1 | typeof(likelihood_fun)!="character")
    stop("the input of likelihood_fun has to one character vector")
  if(!(likelihood_fun %in% c("logit", "probit"))) stop("the input of likelihood_fun has
                                                       to be either 'logit' or 'probit'")
  if(!all(y %in% c(-1,1))) stop("values are just allowed to be -1 or 1")
  f <- numeric(length(y))
  a <- rep(1, length(y))
  precision <- 10 ^(-5)
  response_fun <- response_list[[likelihood_fun]]
  n <- length(y)
  objective <- function(a,f){
    -0.5* t(a)%*%f + response_fun$log_p(y,f)
  }
  objective0 <- objective(a,f)
  objective1 <- 0
  for(i in 1:100) {
    W <- - response_fun$hess_f(y,f)
    L <- chol(diag(nrow = n) + sqrt(W)%*%K%*%sqrt(W))
    b <- W%*%f + response_fun$del_f(y,f)
    a <- b-backsolve(sqrt(W)%*%t(L),backsolve(L,sqrt(W)%*%K%*%b))
    f <- K%*%a
    objective1 <- objective(a,f)
    if(abs(objective0-objective1)<precision)
      break
    objective0 <- objective1
  }

  log_q <- objective(a,f) - sum(log(diag(L)))
  res <- list(mode = f, log_marginal_likelihood = log_q)
  if(i ==100)
    attr(res, "convergence") <- FALSE
  else
    attr(res, "convergence") <- TRUE
  return(res)
}


#----implementation for the prediction-algorithm-------

#' calculating the prediction probability
#'
#' @param f_mode values of the mode function, it must have the same length as y_learn
#' @param X_learn input values of the learning data, it's allowed to be a list
#' of points, a matrix or a data.frame. The number of represented points has to
#' equal the length of y_learn
#' @param y_learn output values of the learning data, it has to be a numeric vector
#' with 1 and -1 as values
#' @param cov_fun function, such that two numeric vectors x,y of the same length will
#' cause a scalar output
#' @param likelihood_fun name of the used function, a character vector, either 'probit'
#' or 'logit'
#' @param x_input a numeric vector of the same dimension as every point in X_learn,
#' point, where we want to predict the class label
#'
#' @return numeric vector, that describes the probability of the class_label 1 at point x_input
#' @export
#'
#' @examples
#' pred_laplace(f_mode = c(1,2),
#'              X_learn = c(0,2),
#'              y_learn = c(-1,1),
#'              cov_fun = function(x,y)exp(-sqrt(sum((x-y)^2))),
#'              likelihood_fun = "probit",
#'              x_input = 1)
#' -----------------------------------------------
#' pred_laplace(f_mode = c(1,2),
#'              X_learn = list(c(0,2), c(1,3)),
#'              y_learn = c(-1,1),
#'              cov_fun = function(x,y) exp(-sqrt(sum((x-y)^2))),
#'              likelihood_fun = "probit",
#'              x_input = c(1,0))
#'
pred_laplace <- function(f_mode,X_learn, y_learn, cov_fun, likelihood_fun, x_input){
  y_learn <- as.vector(y_learn)
  X_learn <- convert_to_list(X_learn, length(y_learn))

  if(length(likelihood_fun)!=1 | typeof(likelihood_fun)!="character")
    stop("the input of likelihood_fun has to one character vector")
  if(!(likelihood_fun %in% c("logit", "probit")))
    stop("the input of likelihood_fun has to be either 'logit' or 'probit'")
  if(!all(y_learn %in% c(-1,1)))
    stop("values are just allowed to be -1 or 1")
  if (!is.numeric(y_learn) | !is.numeric(X_learn[[1]]) | !is.numeric(x_input))
    stop("Cannot handle non-numeric inputs!")
  if(length(x_input)!=length(X_learn[[1]]))
    stop("Dimension of input data does not fit to dimension of training data")
  if(typeof(cov_fun)!="closure")
    stop("cov has to be a function")
  if(length(f_mode)!=length(y_learn))
    stop("dimension of f_mode doesn't fit to given learning data")
  if(!is.numeric(f_mode))
    stop("f_mode must be numeric")
  n <- length(y_learn)
  W <- -1* response_list[[likelihood_fun]]$hess_f(y_learn,f_mode)
  K <- cov_cross(X_learn, X_learn, cov_fun)
  L <- chol(diag(nrow = n) + sqrt(W)%*%K%*%sqrt(W))
  k_1 <- cov_cross(X_learn, list(x_input), cov_fun)
  f_bar <- t(k_1)%*%response_list[[likelihood_fun]]$del_f(y_learn,f_mode)
  v <- backsolve(L, sqrt(W)%*%k_1)
  var_f <- cov_cross(list(x_input), list(x_input), cov_fun) - t(v)%*%v
  integrand <- function(z){
    return(sapply(z, function(x)
      response_list[[likelihood_fun]]$value(x,1)*dnorm(x,mean=f_bar, sqrt(var_f))))
  }
  tryCatch(
    pred_class_prop <- integrate(integrand, lower = -Inf, upper = Inf)$value,
    error = function(cond) return(NA)
  )
  return(c(pred_class_prop = pred_class_prop))
}

#----overall function for solving the classification problem -----

#' using laplace-approximation for classification problems using gp
#'
#' @param X_learn input values of the learning data, it's allowed to be a list
#' of points, a matrix or a data.frame. The number of represented points has to
#' equal the length of y_learn
#' @param y_learn output values of the learning data, it has to be a numeric vector
#' with 1 and -1 as values
#' @param cov_fun function, such that two numeric vectors x,y of the same length will
#' cause a scalar output
#' @param likelihood_fun name of the used function, a character vector, either 'probit'
#' or 'logit'
#' @param x_input a numeric vector of the same dimension as every point in X_learn,
#' point, where we want to predict the class label
#'
#' @return numeric vector, that describes the probability of the class_label 1 at point x_input
#' @export
#'
#'
#' @examples
#' predict_laplace(X_learn = c(0,2),
#'                 y_learn = c(-1,1),
#'                 cov_fun = function(x,y) exp(-sqrt(sum((x-y)^2))),
#'                 likelihood_fun = "probit",
#'                 x_input = 1)
#'------------------------------------------------------------
#'predict_laplace( X_learn = list(c(0,2), c(1,3)),
#'                 y_learn = c(-1,1),
#'                 cov_fun = function(x,y) exp(-sqrt(sum((x-y)^2))),
#'                 likelihood_fun = "probit",
#'                 x_input = c(1,0))
predict_laplace <- function(X_learn, y_learn, cov_fun, likelihood_fun, x_input){
  y_learn <- as.vector(y_learn)
  X_learn <- convert_to_list(X_learn, length(y_learn))

  f_mode <- find_mode_laplace(cov_cross(X_learn, X_learn, cov_fun), y_learn, likelihood_fun)$mode
  pred <- pred_laplace(f_mode, X_learn, y_learn, cov_fun, likelihood_fun, x_input)

  if(is.na(pred)){
    warning("The value was not evaluable, so the returned value 0.5 is not trustable")
    return(0.5)
  } else{
    return(pred)
  }
}


