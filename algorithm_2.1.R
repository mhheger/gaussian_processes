predict_gauss <- function(X_learn, y_learn, kov, noise, x_input){
  #do some converting to get either atomic vectors or list of vectors
  X_learn <- convert_to_list(X_learn, length(y_learn))
  #stopifnot((length(x_input)==length(X_learn[[1]])) = "Function can only handle one input vector"))
  #need: add sanity checks for kov, input sizes and types, handle different type
  # of input data
  
  #these values can be stored in a attribut of the gaussian process, they just
  #depend on the learning data
  K <- kov_cross(X_learn, X_learn, kov)
  
  #Not using the cholesky-decomposition, because saving this leads to massive
  #numerical errors, even if its more stable than normal solve-function
  #L <- chol(K + diag(x=noise, nrow = nrow(K)), pivot = TRUE)
  #alpha <- solve(t(L))%*%solve(L)%*%y_learn          
  
  
  tryCatch(
    alpha <- .Internal(La_solve(K + diag(x=noise, nrow = nrow(K)),y_learn, .Machine$double.eps)),
    error = function(cond) return(NaN)
  )
  
  k_1 <- kov_cross(X_learn,list(x_input), kov)      # change of the vector-type of x 
  f_predict <- t(k_1)%*%alpha
  
  v <- solve(K + diag(x=noise, nrow = nrow(K)),k_1) 
  
  var_f <- kov(x_input,x_input)-t(k_1)%*%v
  
  log_marginal_likelihood <- -0.5* t(y_learn) %*% alpha - sum(log(diag(L)))-nrow(L)*0.5*2*pi
  
  results <- list("f_predict" = f_predict,
                  "var_f" = var_f,
                  "log_marginal_likelihood"= log_marginal_likelihood)
  return(results)
}


convert_to_list <- function(x, n){
  if(length(x)-n>0 & !is.list(x)){        
    if(is.matrix(x)){
      x <- apply(x,2,c, simplify = FALSE)
    }
    else{
      stop("The data has to be a vector, a list of vectors or a matrix")
    }
  }
  if(!is.list(x)) x <- as.list(x) 
  return(x)
}
  

kov_cross <- function(x,y,kov){ #returns the variance-kovariance-matrix 
                                #input: list of vectors
  m <- matrix(nrow = length(x), ncol = length(y))
  for (i in seq_along(x)){
    for (j in seq_along(y)){
      m[i,j] <- kov(x[[i]],y[[j]])
    }
  }
  #outer(x,y, function(x,y){mapply(kov,x,y)})
  return(m)
}
