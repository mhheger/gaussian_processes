
constant_cov <- function(x,y,sigma) {
  if(is.null(sigma) & (is.null(y) | is.null(x))) {
    NA
  } else if(!is.numeric(sigma)) {
    stop("parameters must be numeric!")
  } else {
    sigma ^2
  }
}



linear_cov <- function(x,y, sigma) {
  if(is.null(x) | is.null(y) | is.null(sigma)) {
    NA
  } else if(length(x) != length(y) | length(x) != length(sigma)) {
    stop("Input Vectors must have same length!")
  } else if(!is.numeric(x + y + sigma)) {
    stop("Input parameters must be numerical!")
  } else {
    sum(sigma^2 * x * y)
  }
}


squared_exp_cov <- function(x,y, l) {
  if(is.null(x) | is.null(y) | is.null(l)) {
    NA
  } else if(!is.numeric(sqrt(sum((x-y)^2)))) {
    stop("r has to be a real parameter!")
  } else if(l < 0) {
    stop("the parameter l has to be positive!") 
  } else {
    exp(- ((sqrt(sum((x-y)^2))) ^ 2 )/(2*l^2))
  } 
}



exp_cov <- function(x,y,l) {
  if(is.null(x) | is.null(y) | is.null(l)) {
    NA
  } else if(!is.numeric((sqrt(sum((x-y)^2))))) {
    stop("input vectors must be numeric!") 
  } else {
    exp(-((sqrt(sum((x-y)^2)))/l))
  }
}



gamma_exp_cov <- function(x,y, l, gamma) {
  if(is.null(x) | is.null(y) | is.null(gamma)) {
    NA
  } else if(gamma <= 0 | gamma > 2) {
    stop("Gamma has to be in the interval (0,2]!")
  } else if(l <= 0) {
    stop("l must be positive!")
  } else {
    exp(
      -((sqrt(sum((x-y)^2))) / l) ^ gamma 
    )
  }
}


rational_quadratic_cov <- function(x,y, l, alpha) {
  if(is.null(x) | is.null(y) | is.null(alpha) | is.null(l)) {
    NA
  }
  else if(l <= 0 | alpha <= 0) 
    stop("l and alpha have to be positive!")
  else {(
    1 + ((sqrt(sum((x-y)^2))) ^ 2) / (2 * alpha * l ^ 2)
  ) ^ (-alpha)
  }
}





#summarise covariance functions in list
#functions with missing parameters will be assigned NA
covariance_list <- function(l = NULL, alpha = NULL, gamma = NULL,
                            x = NULL, y = NULL, sigma = NULL, tibble = FALSE) {
  lst <- list(
    "constant" = constant_cov(x,y,sigma),
    "linear covariance" = linear_cov(x,y,sigma),
    "squared exponential covariance" = squared_exp_cov(x,y,l),
    "exponential covariance" = exp_cov(x,y,l),
    "gamma exponential covariance" = gamma_exp_cov(x,y,l,gamma),
    "rational quadratic covariance" = rational_quadratic_cov(x,y, l, alpha)
  )
  if(tibble == TRUE) {
    t(tibble::as_tibble(lst))
  } else {
    lst
  }
}

#classes for covariance functions to check if they're stationary or nondegenerate
s <- "stationary"
nd <-  "nondegenerate"
snd <- c("stationary", "nondegenerate")

class(constant_cov) <- c(class(constant_cov),s)

class(squared_exp_cov) <- c(class(squared_exp_cov), snd)
class(gamma_exp_cov) <- c(class(gamma_exp_cov), snd)
class(exp_cov) <- c(class(exp_cov), snd)
class(rational_quadratic_cov) <- c(class(rational_quadratic_cov), snd)

#check if covariance function is stationary
is.stationary <- function(f) {
  all("stationary" %in% class(f))
}
#check if covariance function is nondegenerate
is.nondegenerate <- function(f) {
  all("nondegenerate" %in% class(f))
}









