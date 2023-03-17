# binary laplace approximation
##definition of the different response functions---------

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

## mode finding binary laplace
#' finding value of mode function using binary laplace-approximation
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


##implementation for the binary prediction-algorithm-------

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
#'# -----------------------------------------------
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
    {pred_class_prop <- integrate(integrand, lower = -Inf, upper = Inf)$value
    return(c(pred_class_prop = pred_class_prop))
    },
    error = function(cond) return(NA)
  )


}

##overall function for solving the binary classification problem -----

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
#'#------------------------------------------------------------
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


# Multiclass Laplace ------------------------------------------------------

## Mode-finding multiclass laplace  ----------------------------------------

#' Mode-finding for multiclass laplace approximation
#'
#' @param K_list list of of variance-covariance matrices. Each class label gets
#' it's own matrix
#' @param y vector of resulting class labels, it's values are either 0 or 1.
#'
#' @return list of two elements:
#' - mode: numerical vector with mode function values
#' - log_marginal_likelihood: value of the logarithmic marginal likelihood
#' @export
#' @details
#' Assuming n points, for which the right labeling is known, and C possible
#' classes, than this vector is build up like this:
#'
#' - for each class, we consider the vector y_c, which has length n and whose
#' entries are defined as 1, if point i belongs to class c, or 0 otherwise.
#'
#'- y is created by joining all y_c - vectors.
#'
#' @examples
#' #Assuming 2 classes and 3 learning points
#' K_list <- list(
#'   diag(c(1,2,3)),
#'   diag(c(0.4,2.3,1))
#' )
#'
#' y <- c(0,0,1,1,1,0)
#'
#' find_mode_mc_laplace(K_list, y)
find_mode_mc_laplace <- function(K_list, y){     #K list of matrices, y output (0,1)
  # checks
  if(is.null(K_list) | is.null(y))
    stop("Missing arguments.")
  if(!all(y==1 | y ==0))
    stop("Values of labels have to be either 1 or 0.")
  expected_dim <- dim(K_list[[1]])
  check_dim <- sapply(K_list, dim)
  if(!all(check_dim == expected_dim))
    stop("All matrices in K_list must have the same dimension!")


  # starting calculation
  n <- nrow(K_list[[1]])
  C <- length(K_list)
  f <- numeric(n*C)
  precision <- 10 ^(-15)
  # implementing newton iteration
  for(iterations in 1:100){
    E_list <- list()
    PI <- list()
    z <- numeric(C)
    p <- sapply(seq_len(n*C), function(i){
      exp(f[i])/(sum(exp(f[seq_len(n*C)%%n == i %% n ])))
    })

    objective <- function(a,f){
      -0.5 * sum(a * f)+ sum(y * f) - sum(log(sapply(seq_len(C), function(i){
        (sum(exp(f[seq_len(n*C)%%n == i %% n ])))
      })))
    }

    objective1 <- 0

    for(c in seq_len(C)){
      p_c <- p[((c-1)*n +1) : (c*n)]
      D_c <- diag(p_c)
      PI[[c]] <- D_c
      L <- chol(diag(nrow = n)+sqrt(D_c)%*%K_list[[c]]%*%sqrt(D_c))
      E_c <- backsolve(sqrt(D_c)%*%t(L), (backsolve(L, sqrt(D_c))))
      E_list[[c]] <- E_c
      z[c] <- sum(log(diag(L)))
    }
    M <- diag(0, nrow= n)
    for(item in E_list){
      M <- M + item
    }
    #M <- chol(M)
    PI_matrix <- PI[[1]]
    for(item in 2:C){
      PI_matrix <- rbind(PI_matrix, PI[[item]])
    }
    R <- diag(nrow = n)
    for(item in 2:C){
      R <- rbind(R, diag(nrow = n))
    }

    E <- diag_block_matrix(E_list)
    K <- diag_block_matrix(K_list)
    D <- diag(p)

    b <- (D-PI_matrix %*% t(PI_matrix)) %*% f + y - p
    c <- E %*% K %*% b
    #a <- b - c + E %*% R %*%backsolve(( t(M)), backsolve(M, t(R)%*%c))
    a <- b - c + E %*% R %*%solve(M, t(R)%*%c)
    f = K%*%a

    objective0 <- objective(a,f)
    if(abs(objective0-objective1)<precision)
      break
    objective0 <- objective1

  }
  log_q <- objective(a,f) - sum(z)
  res <- list(mode = f, log_marginal_likelihood = log_q)
  if(iterations ==100)
    attr(res, "convergence") <- FALSE
  else
    attr(res, "convergence") <- TRUE
  return(res)
}


## building block matrices  ------------------------------------------------

#' Building block matrices
#'
#' @param matrix_list list of quadratic matrices
#'
#' @return one diagonal block matrix
#' @export
#'
#' @examples
#' l <- list(
#'   matrix(1:9, nrow = 3),
#'   diag(4),
#'   diag(c(1,2,3,3,3)),
#'   matrix(90:98, nrow = 3)
#' )
#'
#' diag_block_matrix(l)
#'
#'
diag_block_matrix <- function(matrix_list){
  #checking if list is non empty
  check_empty <- sapply(matrix_list, function(m) return(!is.null(m)))
  if(!all(check_empty)) stop("matrix_list is not allowed to have empty entries")
  #checking if list entries are matrices
  check_matrix <- sapply(matrix_list, is.matrix)
  if(!all(check_matrix)) stop("There are non-matrix elements in matrix_list")
  #checking if matrices are quadratic
  check_q <- sapply(matrix_list, function(m) return(nrow(m)==ncol(m)))
  if(!all(check_q)) stop("There are non-quadratic matrices in matrix_list")

  dims <- sapply(matrix_list, nrow)
  M <- diag(sum(dims))
  curr_pos <- 1
  for(i in seq_along(dims)){
    end_pos <- curr_pos-1 + dims[i]
    M[curr_pos:end_pos,curr_pos:end_pos ] <- matrix_list[[i]]
    curr_pos <- end_pos +1
  }
  return(M)
}

## internal prediction algorithm  ----------------------------------------


#' Predictions using multiclass-laplace-approximations
#'
#' @param X_learn list of input points
#' @param y_learn input class labels of learning data
#' @param K_list  list of variance-covariance-matrices
#' @param f_mode numeric vector of mode function values
#' @param cov_list list of closures, that fullfill the characteristics of a
#' covariance function
#' @param x_input numeric vector, coordinates of the input, where we want to get
#' the prediction
#' @param n_sample number of samples to approximate the probability of the class
#' labels
#'
#' @return vector describing the probabilities of the different labels
#' @export
#'
#' @examples
#' x <- 1:50
#' y <-as.integer(x<5)
#' y <- c ( y, as.integer(x>5))
#' X_learn <- as.list(x)
#' covariance_list <- list(
#'   function(x,y) exp(- 0.5*sqrt(sum((x-y)^2))),
#'   function(x,y) exp(- 0.5*sqrt(sum((x-y)^2)))
#' )
#'
#' K <- list(
#'  cov_cross(X_learn, X_learn,covariance_list[[1]]),
#'  cov_cross(X_learn, X_learn,covariance_list[[2]])
#' )
#'
#'f_mode <- find_mode_mc_laplace(K,y)$mode
#'
#'pred_mc_laplace(X_learn,y, K, f_mode, covariance_list, 20)
pred_mc_laplace <- function(X_learn, y_learn, K_list, f_mode, cov_list, x_input, n_sample = 1000){
  #checks
  if(is.null(K_list) | is.null(y_learn) | is.null(f_mode) |is.null(x_input) | is.null(cov_list))
    stop("Missing arguments.")
  if(!all(y_learn==1 | y_learn ==0))
    stop("Values of labels have to be either 1 or 0.")
  expected_dim <- dim(K_list[[1]])
  check_dim <- sapply(K_list, dim)
  if(!all(check_dim == expected_dim))
    stop("All matrices in K_list must have the same dimension!")
  if(!is.list(X_learn))
    stop("X_learn has to be a list of numeric vectors")
  check_num <- sapply(X_learn,is.numeric)
  if(!all(check_num))
    stop("X_learn has to be a list of numeric vectors")
  if(! length(f_mode)==length(y_learn))
    stop("Length of f_mode doesn't fit to the input data")
  check_clo <- sapply(cov_list, function(x) typeof(x)=="closure")
  if(!all(check_clo))
    stop("Not all items in cov_list are closures")
  if(! length(x_input) == length(X_learn[[1]]))
    stop("x_input hasn't the same length as input data")

  n <- nrow(K_list[[1]])
  C <- length(K_list)
  p <- sapply(seq_len(n*C), function(i){
    exp(f_mode[i])/(sum(exp(f_mode[seq_len(n*C)%%n == i %% n ])))
  })
  E_list <- list()
  PI <- list()
  mu <- numeric(C)
  Sigma <- diag(C)
  num_samples <- n_sample

  for(c in seq_len(C)){
    p_c <- p[((c-1)*n +1) : (c*n)]
    D_c <- diag(p_c)
    PI[[c]] <- D_c
    L <- chol(diag(nrow = n)+sqrt(D_c)%*%K_list[[c]]%*%sqrt(D_c))
    E_c <- backsolve(sqrt(D_c)%*%t(L), (backsolve(L, sqrt(D_c))))
    E_list[[c]] <- E_c
  }

  R <- diag(nrow = n)
  for(item in 2:C){
    R <- rbind(R, diag(nrow = n))
  }

  M <- diag(0, nrow= n)
  for(item in E_list){
    M <- M + item
  }
  #M <- chol(M)

  for( c in seq_len(C)){
    k_c <- cov_cross(X_learn, list(x_input), cov_list[[c]])
    mu[c] <-  t(y_learn[((c-1)*n+1):(c*n)])%*%k_c
    b <- E_list[[c]]%*% k_c
    g <- E_list[[c]]%*%(solve(M, b))

    for(d in seq_len(C)){
      k_d <- cov_cross(X_learn, list(x_input), cov_list[[d]])
      Sigma[c,d] <- t(g) %*% k_d
    }
    Sigma[c,c] <- Sigma[c,c] + cov_cross(list(x_input), list(x_input), cov_list[[c]]) +
      t(b)%*%k_c
  }

  pred <- numeric(C)
  for(i in seq_len(num_samples)){
    f <- mnormt::rmnorm(1, mean = mu, Sigma)
    pred <- pred + exp(f)/sum(exp(f))
  }

  pred <- pred / num_samples

  return(pred)
}



## Class for more-user friendly multiclass laplace approximation
R6::R6Class("gp_classification",
            public = list(
              #'# Generating a new 'gp_classification' object
              #'
              #' @param n  number of classes, that are used
              #' @param covs of length n, with names of the covariance functions, that
              #' should be used. Per default squared_exp
              #'
              #' @return a new 'gp_classification' object
              #' @export
              #'
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #'
              initialize = function(n, covs = rep("squared_exp",n)){
                private$number_categories <- n
                if(!length(covs)==n)
                  stop(stringr::str_glue("Length of covs has to equal {n}"))
                for (i in seq_len(n)){
                  private$gp_list[[i]] <- gp$new()
                  private$gp_list[[i]]$set_cov(covs[[i]])
                  private$covariances[[i]]<-private$gp_list[[i]]$get_cov()
                }
              },
              #' #Adding data to 'gp_classification'-object
              #'
              #' @param X_learn a list of points, a matrix or a data.frame of the points with
              #' known input data
              #' @param y vector of class labels (more details in details)
              #'
              #' @return modified 'gp_classification'-object
              #' @export
              #' @details
              #' Assuming n points, for which the right labeling is known, and C possible
              #' classes, than this vector is build up like this:
              #'
              #' - for each class, we consider the vector y_c, which has length n and whose
              #' entries are defined as 1, if point i belongs to class c, or 0 otherwise.
              #'
              #'- y is created by joining all y_c - vectors.
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #' y <- c(
              #'   c(0,1),       #first point belongs to the first class, second point doesn't
              #'   c(1,0)        #second point belongs to the first class, first point doesn't
              #' )
              #'
              #' X_learn <- list(
              #'   c(1,2,3),
              #'   c(-1,-2,-3)
              #' )
              #'
              #' mcgp$add_data(X_learn, y)
              #'
              add_data = function(X_learn, y){
                n <- private$number_categories
                y <- unlist(y)
                if(length(y)%%n !=0)
                  stop("Length of y is not a multiple of the number of classes!")
                if(!all(y==1 | y ==0))
                  stop("Values of labels have to be either 1 or 0.")
                k = length(y)/n
                for(i in seq_len(n)){
                  #y_i <- y[seq_len(length(y))%%n == i %% n]
                  y_i <- y[((i-1)*k +1): (i*k)]
                  private$gp_list[[i]]$add_data(X_learn, y_i)
                }
                X_learn <- convert_to_list(X_learn, length(y)/n )
                private$input_dimension <- length(X_learn[[1]])
                private$y_learn <- c(private$y_learn,y)
                private$X_learn <- c(private$X_learn,X_learn)
                private$set_mode()
              },
              #' #Getting list of covariance matrices
              #'
              #' @return list of covariance matrices, each class has it's own label
              #' @export
              #'
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #' y <- c(
              #'   c(0,1),       #first point belongs to the first class, second point doesn't
              #'   c(1,0)        #second point belongs to the first class, first point doesn't
              #' )
              #'
              #' X_learn <- list(
              #'   c(1,2,3),
              #'   c(-1,-2,-3)
              #' )
              #'
              #' mcgp$add_data(X_learn, y)
              #' mcgp$get_K_list()
              get_K_list = function(){
                n <- private$number_categories
                l <- list()
                for(i in seq_len(n)){
                  l[[i]] <- private$gp_list[[i]]$get_K()
                }
                return(l)
              },
              #' #Getting list of learning input data
              #'
              #' @return list of learning input data
              #' @export
              #'
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #' y <- c(
              #'   c(0,1),       #first point belongs to the first class, second point doesn't
              #'   c(1,0)        #second point belongs to the first class, first point doesn't
              #' )
              #'
              #' X_learn <- list(
              #'   c(1,2,3),
              #'   c(-1,-2,-3)
              #' )
              #'
              #' mcgp$add_data(X_learn, y)
              #' mcgp$get_X()
              #'
              get_X = function(){
                return(private$X_learn)
              },
              #' #Getting list of learning input labels
              #'
              #' @return list of learning input labels
              #' @export
              #'
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #' y <- c(
              #'   c(0,1),       #first point belongs to the first class, second point doesn't
              #'   c(1,0)        #second point belongs to the first class, first point doesn't
              #' )
              #'
              #' X_learn <- list(
              #'   c(1,2,3),
              #'   c(-1,-2,-3)
              #' )
              #'
              #' mcgp$add_data(X_learn, y)
              #' mcgp$get_y()
              get_y = function(){
                return(private$y_learn)
              },
              #' #Getting list of used covariance functions
              #'
              #' @return list of used covariance functions
              #' @export
              #'
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #' y <- c(
              #'   c(0,1),       #first point belongs to the first class, second point doesn't
              #'   c(1,0)        #second point belongs to the first class, first point doesn't
              #' )
              #'
              #' X_learn <- list(
              #'   c(1,2,3),
              #'   c(-1,-2,-3)
              #' )
              #'
              #' mcgp$add_data(X_learn, y)
              #' mcgp$get_covariances()
              get_covariances = function(){
                return(private$covariances)
              },
              #' #Getting class label probabilities
              #'
              #' @param x_input numeric vector of input data, where you want to get the prediction
              #' @param n_samples number of samples to get probabilities using monte-carlo method
              #'
              #' @return vector of probabiilties for each label
              #' @export
              #'
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #' y <- c(
              #'   c(0,1),       #first point belongs to the first class, second point doesn't
              #'   c(1,0)        #second point belongs to the first class, first point doesn't
              #' )
              #'
              #' X_learn <- list(
              #'   c(1,2,3),
              #'   c(-1,-2,-3)
              #' )
              #'
              #' mcgp$add_data(X_learn, y)
              #' mcgp$get_prediction(c(1,2,4))
              #'
              get_prediction = function(x_input, n_samples = 1000){
                f_mode <- private$f_mode
                K_list <- self$get_K_list()
                covs <- self$get_covariances()
                y <- self$get_y()
                X_learn <- self$get_X()
                if(is.null(f_mode)) stop("You have to add data first!")
                if(length(x_input) != private$input_dimension)
                  stop("Length of x_input does not fit to dimension of learning input")

                pred_mc_laplace(X_learn,y, K_list, f_mode, covs, x_input, n_samples)
              },

              #' #Setting parameters
              #'
              #' @param index number of the class label, whose covariance function should be
              #' modified
              #' @param cov_name character name of the covariance function that should be used.
              #' See ?set_cov to get further informations about the possible functions
              #' @param parameter_list named list of parameters, that should be used
              #'
              #' @details See ?set_cov or ?set_parameter for further information about the possible
              #' values of cov_name or parameter_list
              #' @export
              #'
              #' @examples
              #' mcgp <- gp_classification$new(n = 2)
              #' y <- c(
              #'   c(0,1),       #first point belongs to the first class, second point doesn't
              #'   c(1,0)        #second point belongs to the first class, first point doesn't
              #' )
              #'
              #' X_learn <- list(
              #'   c(1,2,3),
              #'   c(-1,-2,-3)
              #' )
              #'
              #' mcgp$add_data(X_learn, y)
              #' mcgp$set_parameters(1, "linear", list(sigma = c(1,1,1)))
              set_parameters = function(index, cov_name=NULL, parameter_list=NULL){
                if(is.null(index))
                  stop(stringr::str_glue("index has to be between 1 and {private$number_categories}"))
                if(!is.numeric(index) | is.na(index))
                  stop("index has to be a numeric vector")
                if(index <1 | index > private$number_categories)
                  stop(stringr::str_glue("index has to be between 1 and {private$number_categories}"))

                if(!is.null(cov_name))
                  private$gp_list[[index]] %>%
                  set_cov(cov_name)
                if(!is.null(parameter_list))
                  private$gp_list[[index]] %>%
                  set_parameter(parameter_list)
                private$set_mode()
              }
            ),
            private = list(
              set_mode = function(){
                K_list <- self$get_K_list()
                y <- self$get_y()
                res <- find_mode_mc_laplace(self$get_K_list(),self$get_y())
                private$f_mode <- res$mode
              },
              number_categories = NULL,
              X_learn = NULL,
              y_learn = NULL,
              gp_list = list(),
              covariances = list(),
              f_mode = NULL,
              input_dimension = NULL
            )
) -> gp_classification

