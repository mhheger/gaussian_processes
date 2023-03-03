#' predict_gauss
#'
#' @param X_learn The learning data, given as an matrix, a list of vectors or a tibble
#' @param y_learn The values of the unknown function at the positions given in X_learn
#' @param cov The covariance-function, which is used. It has to provide a function call
#' like cov(x,y) with x and y vectors of the same size and outputs a scalar.
#' @param noise The assumed value of the variance of the noise of the learning data
#' @param x_input The position we want to predict the value of the unknown function
#'
#' @return list of three elements f_predict var_f log_marginal_likelihood
#' @param f_predict The predicted value
#' @param var_f The variance
#' @export
#'
#' @examples
predict_gauss <- function(X_learn, y_learn, cov, noise, x_input){

  #do some converting to get either atomic vectors or list of vectors
  y_learn <- as.vector(y_learn)
  X_learn <- convert_to_list(X_learn, length(y_learn))

  if(length(noise)!=1) stop("length of noise is not 1")
  if(typeof(cov)!="closure") stop("cov has to be a function")

  if(noise < 0){
    stop("Noise is not allowed to be negative!")
  } else if (!is.numeric(y_learn) | !is.numeric(X_learn[[1]])
             | !is.numeric(x_input) | !is.numeric(noise)){
    stop("Cannot handle non-numeric inputs!")
  } else if(length(x_input)!=length(X_learn[[1]])) {
    stop("Dimension of input data does not fit to dimension of training data")
  }

  #need: add sanity checks for cov , handle different type
  # of input data

  #these values can be stored in a attribute of the gaussian process, they just
  #depend on the learning data
  K <- cov_cross(X_learn, X_learn, cov)
  if(det(K)==0 && noise == 0) noise = 0.01
  L <- chol(K + diag(x=noise, nrow = nrow(K)), pivot = TRUE)

  tryCatch(
    alpha <- .Internal(La_solve(K + diag(x=noise, nrow = nrow(K)),y_learn, .Machine$double.eps)),
    error = function(cond) return(NaN)
  )

  k_1 <- cov_cross(X_learn,list(x_input), cov)      # change of the vector-type of x
  v <- solve(K + diag(x=noise, nrow = nrow(K)),k_1)

  f_predict <- t(k_1)%*%alpha
  var_f <- cov(x_input,x_input)-t(k_1)%*%v
  log_marginal_likelihood <- -0.5* t(y_learn) %*% alpha - sum(log(diag(L)))-nrow(L)*0.5*2*pi

  results <- list("f_predict" = f_predict,
                  "var_f" = var_f,
                  "log_marginal_likelihood"= log_marginal_likelihood
  )
  return(results)
}



predict_gauss2 <- function(X, x_input){
  #need: add sanity checks for cov , handle different type
  # of input data

  #these values can be stored in a attribute of the gaussian process, they just
  #depend on the learning data
  learn_data <- X$get_data()
  K <- X$get_K()
  cov <- X$get_cov()
  noise <- X$get_noise()
  mean_fun <- X$get_mean_fun()

  y_learn <- learn_data$y_learn - sapply(learn_data$X_learn, mean_fun)
  #L <- chol(K + diag(x=noise, nrow = nrow(K)), pivot = TRUE)

  #tryCatch(
  # alpha <- .Internal(La_solve(K + diag(x=noise, nrow = nrow(K)),learn_data$y_learn, .Machine$double.eps)),
  #error = function(cond) return(NaN)
  #)
  alpha <- solve(K + diag(x=noise, nrow = nrow(K)),y_learn)
  k_1 <- cov_cross(learn_data$X_learn,list(x_input), cov)      # change of the vector-type of x

  v <- solve(K + diag(x=noise, nrow = nrow(K)),k_1)

  f_predict <- t(k_1)%*%alpha + mean_fun(x_input)
  var_f <- cov(x_input,x_input)-t(k_1)%*%v
  log_marginal_likelihood <- -0.5* t(learn_data$y_learn) %*% alpha - log(det(K+ diag(x=noise, nrow = nrow(K))))-nrow(K)*0.5*2*pi

  results <- list("f_predict" = f_predict,
                  "var_f" = var_f,
                  "log_marginal_likelihood"= log_marginal_likelihood
  )
  return(results)
}


#' convert_to_list
#' A package interna, that provides the converting of the input data to the
#' needed list type
#' @param x vector, matrix, data.frame, that should be converted into a list of
#' input vectors
#' @param n the dimension of each vector in the list that is returned
#'
#' @return a list of vectors of length n
#'
#' @examples
convert_to_list <- function(x, n){
  if(is.data.frame(x)) x <- as.matrix(x)
  if(length(x)%%n !=0) stop("length/dimension of input data doesn't fit")
  if(is.vector(x)) x <- matrix(x, nrow= n)
  if(length(x)-n>0 & !is.list(x)){
    if(is.matrix(x)){
      if(ncol(x)==n) x <- t(x)
      x <- apply(x,1,function(x) c(unname(x)), simplify = FALSE)
    }
    else{
      stop("The data has to be a vector, a list of vectors, a matrix or a dataframe")
    }
  }
  if(!is.list(x)) x <- as.list(x)
  return(x)
}


#' calculating the variance-covariance-matrix
#'
#' @param x first list of input vectors
#' @param y second list of input vectors
#' @param cov covariance function, that returns a scalar
#'
#' @return variance-covariance-matrix of dimension 'length(x)':'length(y)'
#'
#' @examples
#'
cov_cross <- function(x,y,cov){ #returns the variance-covariance-matrix
  if(!is.list(x) | !is.list(y) | !is.numeric(x[[1]])| !is.numeric(y[[1]])){
    stop("Input has to be a list of numerical vectors!")
  }

  m <- matrix(nrow = length(x), ncol = length(y))
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      m[i,j] <- cov(x[[i]],y[[j]])
    }
  }
  return(m)
}
