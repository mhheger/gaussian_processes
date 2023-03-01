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

    - log(1+ exp(-t(y)%*%f))
  },
  del_f = function(y,f) {
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")

    0.5*(y+1) - (1+ exp(-1*f))^(-1)

  },
  hess_f = function(y,f){
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")

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

    log(pnorm(t(y)%*%f))
  },
  del_f = function(y,f) {
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")

    (y*dnorm(f))/pnorm(y*f)

  },
  hess_f = function(y,f){
    if(! length(y)==length(f))
      stop("Input vectors do not have the same length!")
    if(!is.numeric(y) | !is.numeric(f))
      stop("Input vectors have to be numerical")

    diag(-(dnorm(f)/pnorm(f*y))- (y*f*dnorm(f))/pnorm(f*y),nrow=length(y))
  }
)
response_list <- list(logit =logit, probit = probit)


#----implementation of the mode-finding-algorithm------

#' finding value of mode function using laplace-approximation
#'
#' @param K covariance-matrix of dimension (n,n)
#' @param y output-vector of length n, it's values are just allowed to be -1 or 1
#' @param likelihood_fun character-vector, name of the response function,
#' either "probit" or "logit"
#'
#' @return  f_mode: value of mode function,
#'          log_marginal_likelihood: logarithmic marginal likelihood of the result
#'
#' @examples
#'
find_mode_laplace <- function(K, y, likelihood_fun){
  if(!is.matrix(K)) stop("K has to be a numerical matrix")
  if(! nrow(K)==ncol(K)) stop("K has to be a square matrix")
  if(!nrow(K)==length(as.vector(y))) stop("dimension of K and y do not fit")
  if(!is.numeric(K)|!is.numeric(y)) stop("input data has to be numerical")
  if(length(likelihood_fun)!=1 | typeof(likelihood_fun)!="character")
    stop("the input of likelihood_fun has to one character vector")
  if(!(likelihood_fun %in% c("logit", "probit"))) stop("the input of likelihood_fun has
                                                       to be either 'logit' or 'probit'")
  if(!all(y %in% c(-1,0,1))) stop("values are just allowed to be -1, 0 or 1")
  response_fct <- function(f){
    res <- -response_list[[likelihood_fun]]$log_p(y,f) +0.5*t(f)%*%solve(K,f)+0.5*(log(abs(det(K)))-length(y)*log(2*pi))
    attr(res,"gradient") <- -1*response_list[[likelihood_fun]]$del_f(y,f)+ solve(K,f)
    attr(res,"hessian") <- -1*response_list[[likelihood_fun]]$hess_f(y,f)+ solve(K)
    return(res)
  }
  f<- rep(0,length = length(y))
  f <- nlm(response_fct,f)$estimate
  tryCatch(
    log_q <- -0.5*t(f)%*%solve(K,f)+response_list[[likelihood_fun]]$log_p(y,f)-
      0.5*log(det(K)*det(solve(K)+response_list[[likelihood_fun]]$hess_f(y,f)))
  )

  return(list(mode = f, log_marginal_likelihood = log_q))
}

#----implementation for the prediction-algorithm-------
#' calculating the prediction probability
#'
#' @param f_mode values of the mode function, it must have the same length as y_learn
#' @param X_learn input values of the learning data, it's allowed to be a list
#' of points, a matrix or a data.frame. The number of represented points has to
#' equal the length of y_learn
#' @param y_learn output values of the learning data, it has to be a numeric vector
#' with 1 and -1/0 as values
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
#'              y_learn = c(0,1),
#'              cov_fun = function(x,y) exp(-abs(x*y)),
#'              likelihood_fun = "probit",
#'              x_input = 1)
#' -----------------------------------------------
#' pred_laplace(f_mode = c(1,2),
#'              X_learn = list(c(0,2), c(1,3)),
#'              y_learn = c(0,1),
#'              cov_fun = function(x,y) exp(-abs(x*y)),
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
  if(!all(y_learn %in% c(-1,0,1)))
    stop("values are just allowed to be -1, 0 or 1")
  if (!is.numeric(y_learn) | !is.numeric(X_learn[[1]]) | !is.numeric(x_input))
    stop("Cannot handle non-numeric inputs!")
  if(length(x_input)!=length(X_learn[[1]]))
    stop("Dimension of input data does not fit to dimension of training data")
  if(typeof(cov)!="closure")
    stop("cov has to be a function")
  if(length(f_mode)!=length(y_learn))
    stop("dimension of f_mode doesn't fit to given learning data")
  if(!is.numeric(f_mode))
    stop("f_mode must be numeric")

  W <- response_list[[likelihood_fun]]$hess_f(y_learn,f_mode)
  K <- cov_cross(X_learn, X_learn, cov_fun)
  k_1 <- cov_cross(X_learn, as.list(x_input), cov_fun)
  f_bar <- t(k_1)%*%response_list[[likelihood_fun]]$del_f(y_learn,f_mode)
  var_f <- cov_cross(as.list(x_input), as.list(x_input), cov_fun) -
    t(k_1)%*%solve(K+solve(W+diag(10^(-10),nrow = nrow(K))),k_1)
  n <- length(f_bar)

#problem: more-dimensional integral
  integrand <- function(z){
    return(sapply(z, function(x)
      response_list[[likelihood_fun]]$value(x,rep(1,n))*dnorm(x,mean=f_bar, sqrt(var_f))))
  }

  pred_class_prop <- integrate(integrand, lower = rep(-Inf,n), upper = rep(Inf,n))$value
  return(c(pred_class_prop = pred_class_prop))
}

#----overall function for solving the classification problem -----

#' using laplace-approximation for classification problems using gp
#'
#' @param X_learn input values of the learning data, it's allowed to be a list
#' of points, a matrix or a data.frame. The number of represented points has to
#' equal the length of y_learn
#' @param y_learn output values of the learning data, it has to be a numeric vector
#' with 1 and -1/0 as values
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
#' predict_laplace(X_learn = c(0,2),
#'                 y_learn = c(0,1),
#'                 cov_fun = function(x,y) exp(-abs(x*y)),
#'                 likelihood_fun = "probit",
#'                 x_input = 1)
predict_laplace <- function(X_learn, y_learn, cov_fun, likelihood_fun, x_input){
  y_learn <- as.vector(y_learn)
  X_learn <- convert_to_list(X_learn, length(y_learn))

  f_mode <- find_mode_laplace(cov_cross(X_learn, X_learn, cov_fun), y_learn, likelihood_fun)$mode
  pred_laplace(f_mode, X_learn, y_learn, cov_fun, likelihood_fun, x_input)

}


