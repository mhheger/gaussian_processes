constant_cov <- function(sigma0) {
  if(is.null(sigma0)) {
    NA
  } else if(!is.numeric(sigma0)) {
    stop("parameters must be numeric!")
  } else {
    sigma0 ^2
  }
}



linear_cov <- function(x1,x2, sigma) {
  if(is.null(x1) | is.null(x2) | is.null(sigma)) {
    NA
  } else if(length(x1) != length(x2) | length(x1) != length(sigma)) {
    stop("Input Vectors must have same length!")
  } else if(!is.numeric(x1 + x2 + sigma)) {
    stop("Input parameters must be numerical!")
  } else {
    sum(sigma^2 * x1 * x2)
  }
}


squared_exp_cov <- function(r, l) {
  if(is.null(r) | is.null(l)) {
    NA
  } else if(!is.numeric(r)) {
    stop("r has to be a real parameter!")
  } else if(l < 0) {
    stop("the parameter l has to be positive!") 
  } else {
    exp(- (r ^ 2 )/(2*l^2))
  } 
}


exp_cov <- function(r,l) {
  if(is.null(r) | is.null(l)) {
    NA
  } else if(!is.numeric(r)) {
    stop("r must be numeric!") 
  } else {
 exp(-(r/l))
  }
}


gamma_exp_cov <- function(r, l, gamma) {
  if(is.null(r) | is.null(gamma)) {
    NA
  } else if(gamma <= 0 | gamma > 2) {
    stop("Gamma has to be in the interval (0,2]!")
  } else if(l <= 0) {
    stop("l must be positive!")
  } else {
          exp(
            -((r / l) ^ gamma )
          )
  }
}


rational_quadratic_cov <- function(r, l, alpha) {
  if(is.null(r) | is.null(alpha) | is.null(l)) {
    NA
  }
  else if(l <= 0 | alpha <= 0) 
      stop("l and alpha have to be positive!")
  else {(
        1 + (r ^ 2) / 2 * alpha * l ^ 2
      ) ^ (-alpha)
  }
}





#summarise covariance functions in list
#functions with missing parameters will be assigned NA
covariance_list <- function(r = NULL, l = NULL, alpha = NULL, gamma = NULL,
                            x1 = NULL, x2 = NULL, sigma = NULL, sigma0 = NULL, tibble = FALSE) {
  lst <- list(
    "constant" = constant_cov(sigma0),
    "linear covariance" = linear_cov(x1,x2,sigma),
    "squared exponential covariance" = squared_exp_cov(r,l),
    "exponential covariance" = exp_cov(r,l),
    "gamma exponential covariance" = gamma_exp_cov(r,l,gamma),
    "rational quadratic covariance" = rational_quadratic_cov(r, l, alpha)
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


















































