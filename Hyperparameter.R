R6::R6Class("Optimization",
            public =list(
              data = function(x,y, noise = 0.1) { # get data
                private$x <- x
                private$y <- y
              },
              ########## getter functions ##########
              get_theta = function(){
                theta_lst <- c("sigma", "l", "alpha", "sigma0", "gamma")
                tmp <- list()
                for (item in theta_lst){
                  lst <- list(private[[item]])
                  names(lst) <- item
                  tmp <- c(tmp, lst)
                }
                tmp
              },
              get_neg_llikelihood = function (hyper){
                sigma0 <- NULL
                sigma <- NULL
                l <- NULL
                gamma <- NULL
                alpha <- NULL
                if ("l" %in% names(hyper)) l <- hyper[[1]]
                if (length(hyper) < 2) {
                  if ("sigma0" %in% names(hyper)) sigma0 <- hyper[[1]]
                  if ("sigma" %in% names(hyper)) sigma <- hyper[[1]]
                } else {
                  if ("gamma" %in% names(hyper)) gamma <- hyper[[2]]
                  if ("alpha" %in% names(hyper)) alpha <- hyper[[2]]
                }
                self$set_hyperparameter(sigma0,
                                     sigma,
                                     l,
                                     gamma,
                                     alpha)
                n <- nrow(private$K)
                print(private$K)
                L <- chol(private$K + diag(rep(private$noise, n)))
                alpha <- solve(t(L), solve(L, private$y))
                log_det <- 2*sum(log(diag(L)))
                log_likelihood <- -0.5*private$y %*% alpha - log_det - 0.5*n*log(2*pi)
                return(-log_likelihood)
                
              },
              ########## setter functions ##########
              set_covariance = function(covariance) {
                private$cov_name <- covariance
                private$cov <- covariance_function(covariance,
                                                   sigma0 = private$sigma0,
                                                   sigma = private$sigma,
                                                   l = private$l,
                                                   gamma = private$gamma,
                                                   alpha = private$alpha)
                private$K <- covariance_matrix(private$x, private$x, private$cov)
              },
              set_hyperparameter = function(sigma0=NULL, sigma=NULL, l=NULL, gamma=NULL, alpha=NULL) {
                if(!is.null(sigma0)) private$sigma0 <- sigma0
                if(!is.null(sigma)) private$sigma <- sigma
                if(!is.null(l)) private$l <- l
                if(!is.null(gamma)) private$gamma <- gamma
                if(!is.null(alpha)) private$alpha <- alpha
                self$set_covariance(private$cov_name)
              }),
            
            
            private = list(
              x = NULL,
              y = NULL,
              noise = 0,
              sigma0 = 1,
              sigma = 1,
              l = 1,
              gamma = 1,
              alpha = 1,
              K = matrix(),
              cov_name = NULL,
              cov = NULL
            )
) -> Optim

# define additional functions
covariance_function <- function(cov_name, sigma0=0, sigma=1, l=1, gamma = 1, alpha=1){
  force(sigma0)
  force(sigma0)
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

covariance_matrix <- function(x, y, cov, ...) {
  K <- matrix(0,nrow=length(x), ncol=length(y))
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      K[i,j] <- cov(x[[i]], y[[j]], ...)
    }
  }
  return(K)
}

neg_llikelihood <- function(hyper, X) {
  X$get_neg_llikelihood(hyper)
}

# define different covariance models
linear <- function(x1,x2,theta){
  linear_cov(x1,x2,theta[["sigma"]])
}
squared_exp <- function(x1,x2, theta){
  squared_exp_cov(abs(x1-x2),theta[["l"]])
}
constant <- function(x,y, theta) {
  constant_cov(theta[["sigma0"]])
}
exponential <- function(x1,x2, theta){
  exp_cov(sqrt(sum(abs(x1-x2)^2)),theta[["l"]])
}
gamma_exp <- function(x1,x2, theta){
  gamma_exp_cov(sqrt(sum(abs(x1-x2)^2)),theta[["l"]],theta[["gamma"]])
}

rational <- function(x1,x2, theta){
  rational_quadratic_cov(sqrt(sum(abs(x1-x2)^2)),theta[["l"]],theta[["alpha"]])
}

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

# optimize Hyperparameters
parameter_optimization <- function(X) {
  
  # define some initial values
  solution <- list()
  initial_theta <- X$get_theta()
  for (i in cov_models) {
    theta_i <- sapply(i$hy_par, function(x) initial_theta[[x]])
    X$set_covariance(i$name)
    hyper_optim <- nlminb(
      start = theta_i,
      objective = neg_llikelihood,
      X = X,
      lower = i$lower,
      upper = i$upper
    )
    names(hyper_optim$par) <- i$hy_par
    lst <- list(hyper_optim$par, hyper_optim$objective)
    solution <- c(solution,list(lst))
  }
  print(solution)
}
