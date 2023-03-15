#-----Definition of the gp-class------------

R6::R6Class("gp",
            public = list(
              #--------adding data---------------
              add_data = function(x,y, noise=0.1) {
                l <- convert_to_list(x, length(y))
                if(!is.null(private$input_dimension)){
                  if(length(l[[1]])!= private$input_dimension){
                    stop("new input cannot be added to the current input data.
                         You may have to initalize a new gp instance.")
                  }
                } else {
                  private$input_dimension <- length(l[[1]])
                }
                private$X_learn <- c(private$X_learn, l)
                private$y_learn <- c(private$y_learn, unlist(y))
                private$sigma <- numeric(private$input_dimension)
                self$set_noise(noise)
                if(is.null(private$covname)){
                  self$set_cov("squared_exp")
                } else {
                  self$set_cov(private$covname)
                }
                self$update_marginal_likelihood()
              },


              #-----------getter functions for the private values-----------------------------
              get_cov_name = function() private$covname,

              get_K = function() private$K,

              get_cov = function() {
                #print(self$get_cov_name())
                private$cov
              }
              ,
              get_data = function() {
                list(X_learn=private$X_learn, y_learn = private$y_learn)
              },


              get_parameter = function(){
                parameter_list <- c("sigma", "l", "alpha", "sigma0", "gamma")
                return <- list()
                for(item in parameter_list){
                  l <- list(private[[item]])
                  names(l)<- item
                  return <- c(return, l)
                }
                return
              },

              get_mean_fun = function(){
                f <- private$mean_fun
                if(is.numeric(f))
                  return(function(x) f)
                return(f)
              },

              get_noise = function() private$noise,

              get_prediction = function(input) {predict_gauss2(self, x_input = input)},

              get_log_marginal_likelihood = function() private$log_marginal_likelihood,

              get_input_dim = function() private$input_dimension,

              #-----------setter functions for the private values-----------------------------
              #argument cov: name of the covariance function
              set_cov = function(cov){
                private$covname <- cov
                private$cov <- init_cov(cov,
                                        sigma = private$sigma,
                                        l = private$l,
                                        alpha = private$alpha,
                                        sigma0 = private$sigma0,
                                        gamma = private$gamma
                )
                if(!is.null(private$input_dimension)){
                  private$K <- cov_cross(private$X_learn, private$X_learn, private$cov)
                  self$update_marginal_likelihood()
                }
              },

              set_parameter = function(sigma=NULL, l=NULL, alpha=NULL, sigma0=NULL, gamma=NULL) {
                # need to add checks in the real use functions
                if(all(!is.null(sigma))) if(all(!is.na(sigma))) private$sigma <- sigma
                if(!is.null(l)) if(!is.na(l)) private$l <- l
                if(!is.null(alpha)) if(!is.na(alpha)) private$alpha <- alpha
                if(!is.null(sigma0)) if(!is.na(sigma0)) private$sigma0 <- sigma0
                if(!is.null(gamma)) if(!is.na(gamma)) private$gamma <- gamma
                self$set_cov(private$covname)
                private$K <- cov_cross(private$X_learn, private$X_learn, private$cov)
                self$update_marginal_likelihood()
              },

              set_noise = function(noise){
                if(is.null(noise))
                  stop("Cannot handle NULL, noise has to be a positive double")
                if( length(noise)!=1)
                  stop("noise is not allowed to have a dimension higher than 1")
                if(noise<0 | !is.numeric(noise))
                  stop("noise has to be a positive double value")
                private$noise <- noise
              },

              set_mean_fun = function(f){
                if((!typeof(f)=="closure" & !is.numeric(f)) | length(f) !=1)
                  stop("mean_fun has to be either a function or a numeric vector
                       of length 1")
                if(typeof(f)=="closure" & length(as.list(args(f)))!=2)
                  stop("mean function is just allowed to work with one input variable")

                if(typeof(f)=="closure" & length(f(numeric(private$input_dimension)))!=1)
                  stop("mean function is just allowed to output numeric values of length 1")
                private$mean_fun <- f
              },

              optim_parameter = function() {optimize_parameters(self)} ,

              plot = function(x_start, x_end, name_x="x", name_y="y",mode = T,title = "") {
                if(private$input_dimension==1){
                  #getting learning points that are inside the plotted area
                  x_input <- unlist(private$X_learn)
                  ranged <- x_start <= x_input & x_input <= x_end
                  y_input <- unlist(private$y_learn)
                  input <- data.frame(xx=x_input[ranged],yy=y_input[ranged])

                  #getting predicted values inside the intervall [x_start, x_end]
                  range_x <- seq(x_start,x_end,len = 100)
                  f <-c()
                  var <-c()
                  for(x in range_x){
                    pred <- self$get_prediction(x)
                    f <- c(f, pred$f_predict)
                    var <- c(var, pred$var_f)
                  }


                  plot_data <- data.frame(x = range_x,y = f, min_f= f-var,max_f= f+var )

                  plot <- ggplot2::ggplot(plot_data, aes(x,y))+
                    geom_ribbon(aes(ymin = min_f,ymax =  max_f, fill = "variance"))+
                    geom_line(aes(x,y))+
                    labs(x = name_x, y = name_y, title = title ) +
                    theme(legend.position = "left")
                  if(mode) plot <- plot + geom_point(data=input, aes(xx,yy))
                  plot
                }
              },

              print = function(){
                if(!is.null(private$covname)){
                  parameters_name <-cov_list[[self$get_cov_name()]]$parameter
                  par_string <- ""
                  for(item in parameters_name){
                    par_string <- stringr::str_c(par_string,"    - ", item , " : ",
                                                 stringr::str_flatten(self$get_parameter()[[item]], " "), "\n")
                  }
                  s <- stringr::str_glue(
                    "------------------------------------------------------------
                  --------------------Gaussian Process-------------------------
                  type: regression \n
                  input dimension: {private$input_dimension}\n
                  size of learning data: {length(private$y_learn)} data points\n
                  used covariance function: {self$get_cov_name()} \n
                  parameters:
                  {par_string}
                  marginal log likelihood: {self$get_log_marginal_likelihood()}
                  ------------------------------------------------------------
                  ------------------------------------------------------------
                  To get further informations about the different attributes
                  use more detailed methods like get_data(<object>) or
                  get_mean_fun(<object>).
                  "
                  )
                  print(s)
                }
                else {
                  print("Empty gp object. Please add data for a more detailed output.")
                }
              },

              update_marginal_likelihood = function(){
                private$log_marginal_likelihood <-
                  unclass(self$get_prediction(numeric(private$input_dimension))$log_marginal_likelihood)
              }
            ),

            private = list(
              K = matrix(),   # variance-covariance-matrix of the input data,
              log_marginal_likelihood = Inf,
              cov = NULL , #after initializing here is a closure object
              covname = NULL,
              # values of the parameters
              sigma = 1,
              l = 1,
              alpha = 1,
              sigma0=1,
              gamma = 1,

              noise = 0,

              X_learn = list(),
              y_learn = numeric(),
              input_dimension = NULL,

              mean_fun = 0
            )

) -> gp

#------helper-functions-------------------

linear <- function(x,y,par){
  linear_cov(x,y,par$sigma)
}
squared_exp <- function(x,y, par){
  squared_exp_cov(sqrt(sum(abs(x-y)^2)),par$l)
}
constant <- function(x,y, par){
  constant_cov(par$sigma0)
}
exponential <- function(x,y, par){
  exp_cov(sqrt(sum(abs(x-y)^2)),par$l)
}
gamma_exp <- function(x,y,par){
  gamma_exp_cov(sqrt(sum(abs(x-y)^2)),par$l,par$gamma)
}
rational <- function(x,y,par){
  rational_quadratic_cov(sqrt(sum(abs(x-y)^2)),par$l,par$alpha)
}

cov_list <- list(
  "linear"= list(
    name = "linear",
    fun = linear,
    parameter = c("sigma"),
    upper = NULL,
    lower = NULL
  ),

  "squared_exp"= list(
    name = "squared_exp",
    fun = squared_exp,
    parameter = c("l"),
    upper = NULL,
    lower = 0.0001
  ),

  "constant" = list(
    name = "constant",
    fun = constant,
    parameter = c("sigma0"),
    upper = NULL,
    lower = NULL
  ),
  "exponential" = list(
    name = "exponential",
    fun = exponential,
    parameter = c("l"),
    upper = NULL,
    lower = 0.0001
  ),
  "gamma_exp" = list(
    name = "gamma_exp",
    fun = gamma_exp,
    parameter = c("l","gamma"),
    upper = c(Inf, 2),
    lower = c(0.0001, 0.01)
  ),
  "rational" =list(
    name = "rational",
    fun = rational,
    parameter = c("l", "alpha"),
    upper = c(Inf, Inf),
    lower = c(0.0001, 0.0001)
  )
)


#' Creating covariance function, that does not depend on parameters
#'
#' @param covname Character vector, that expresses the name of the used
#' covariance function. One of "constant", "linear", "squared_exp", "exponential",
#'  "rational" or "gamma"
#' @param sigma numerical vector, length has to fit to the input dimension of
#' the used data
#' @param l numerical vector, has to be a non-negative scalar
#' @param alpha numerical vector, has to be a scalar
#' @param sigma0 numerical vector, has to be a scalar
#' @param gamma numerical vector, has to be a scalar between 0 and 2
#'
#' @return covariance function of type 'closure'
#' @export
#'
#' @examples
#' cov_fun <- init_cov("linear", sigma = c(1,2,3))
#' cov_fun(c(1,1,1), c(2,3,4))
#'
init_cov <- function(covname, sigma=0, l=1, alpha=1, sigma0=1, gamma=1){
  if(!is.character(covname) | length(covname)!=1)
    stop("covname has to be a character vector of length 1")
  if(!(covname %in% names(cov_list)))
    stop(paste0(covname, " isn't a known covariance-function.
                Use 'constant', 'linear', 'squared_exp', 'exponential','gamma_exp' or 'rational"))
  force(sigma)
  force(l)
  force(alpha)
  force(sigma0)
  force(gamma)
  l <- list(sigma = sigma, l = l, alpha = alpha, sigma0 = sigma0, gamma = gamma)

  #checking if input is numeric
  for(item in names(l)){
    if(!is.numeric(l[[item]]))
      stop(paste0(item, " has to be numeric"))
  }

  cov_item <- cov_list[[covname]]
  par <- lapply(cov_item$parameter, function(x) l[[x]])
  names(par)<- cov_item$parameter

  #defining the covariance function
  function(x,y) {
    cov_item[["fun"]](x,y, par)
  }
}

#------user-interface-functions----------

#' initializing a new gp-instance
#'
#' @param mode character vector, either "regression" or "classification"
#' @param cov_fun character vector, name of the covariance function, that should
#' be used. You have the choice between:
#' "constant", "linear", "squared_exp", "exponential", "gamma_exp", "rational"
#' @param response_fun character vector, name of the response function, that should
#' be used in case of binary classification, either "logit" for logistic response
#' function or "probit" for using the cummulative normal distribution.
#'
#' @return a "gp" object
#' @export
#'
#' @examples
#' p1 <- new.gp()
#' add_data(p1, x = 1:10, y = 2:11)
#' print(p1)
#'
new.gp <- function(mode="regression", cov_fun = "squared_exp", response_fun = "probit"){
  mode_list <- c("regression", "classification")
  response_names <- names(response_list)
  cov_names <- names(cov_list)
  if(any(is.na(mode))|any(is.null(mode))|any(is.na(cov_fun))|any(is.null(cov_fun))|any(is.na(response_fun))|any(is.null(response_fun)))
    stop("Cannot handle NAs or NULL")
  if(typeof(mode)!= "character" | length(mode)!= 1)
    stop(stringr::str_glue("mode has to be a character. One of {stringr::str_flatten(mode_list, ', ')}"))
  if(typeof(cov_fun)!= "character" | length(mode)!= 1)
    stop(stringr::str_glue("cov_fun has to be a character. One of {stringr::str_flatten(cov_names, ', ')}"))
  if(typeof(mode)!= "character" | length(mode)!= 1)
    stop(stringr::str_glue("mode has to be a character. One of {stringr::str_flatten(response_names, ', ')}"))
  if(!(mode %in% mode_list))
    stop(stringr::str_glue("{mode} is not one of {stringr::str_flatten(mode_list, ', ')}"))
  if(!(cov_fun %in% cov_names))
    stop(stringr::str_glue("{cov_fun} is not one of {stringr::str_flatten(cov_names, ', ')}"))
  if(!(response_fun %in% response_names))
    stop(stringr::str_glue("{response_fun} is not one of {stringr::str_flatten(response_names, ', ')}"))

  p <- gp$new()
  p$set_cov(cov_fun)

  #need to do: handling with classification

  invisible(p)
}




#' adding data to a gp-object
#'
#' @param obj a instance of class "gp"
#' @param X_learn either a vector, a matrix, a data.frame or a list of points
#' with numerical values.
#' @param y either a numerical vector or in the case that X_learn is a data.frame
#' with named columns y is also allowed to be a character, that refers to the name
#' of one column of X_learn
#' @param noise a non-negative numerical value, that describes the standard deviation
#' of the noise
#'
#' @return a modified version of obj. It is possible to pipe this function due to
#' it's invisible output.
#' @export
#'
#' @details
#' Members of the 'gp'-class are modified via reference, such that the
#' object will be modified just by the call of 'add_data'. You do not need
#' to save the changes by yourself. If you do not want to have does changes,
#' then you might have to clone the object before passing it to 'add_data'
#'
#' @examples
#' # Generating input data
#' input_x <- 1:100
#' input_y <- sin(0.5*input_x)
#'
#' # one way to use this function
#' p <- new.gp()
#' add_data(p, input_x, input_y)
#'
#' # an equivalent way to get to the same result:
#' new.gp() %>% add_data(input_x, input_y) -> p
#'
#' # it is also possible to add data using data.frames
#'
#' df <- data.frame(x = input_x, y = input_y)
#'
#' #now we use the character "y" to refer to the column of df with name y
#' new.gp() %>% add_data(df, "y") -> p
#'

add_data <- function(obj, X_learn, y, noise = 0.1){
  if(!any(class(obj)=="gp")){
    stop("obj has to be a member of class 'gp' ")
  }
  if(is.data.frame(X_learn)){
    if(length(y)==1){
      if(is.character(y) & y %in% colnames(X_learn) ){
        name <- y
        y <- X_learn[[y]]
        X_learn[name] <- NULL
      }
      else if(is.character(y)&!(y %in% colnames(X_learn)))
        stop(stringr::str_c(y, " is not a name of a column of X_learn."))
    }
  }

  if(!is.numeric(y) | !is.numeric(X_learn[[1]]))
    stop("Input must lead to numeric data")
  obj$add_data(X_learn, y, noise)
  invisible(obj)
}



#' Getting the name of the used covariance function
#'
#' @param obj a instance of class "gp"
#'
#' @return character vector, the name of the used covariance function
#' @details For further information about usable covariance functions, you
#' might consider the documentation of "set_cov"
#' @export
#' @seealso set_cov() for setting the covariance function
#' @examples
#' p <- new.gp() %>%
#'      add_data(1:10,1:10)
#'get_cov_name(p)
#'#outputs default value: "squared_exp"
#'
#'------------------------------------------
#'p <- new.gp(cov_fun = "linear")
#'get_cov_name(p)
#'#output: "linear"

get_cov_name <- function(obj){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  return(obj$get_cov_name())
}


#' Getting the variance-covariance matrix
#'
#' @param obj a instance of class "gp"
#'
#' @return variance-covariance matrix of "obj"
#' @export
#'
#' @examples
#' #initializing a new 'gp' object
#' p <- new.gp() %>%
#'      add_data(1:3, 1:3)
#'
#' get_K(p)
#'
#'
get_K <- function(obj){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  return(obj$get_K())
}

#' Getting the used covariance function
#'
#' @param obj a instance of class "gp"
#'
#' @return a list of two elements:
#' name: the name of the covariance function
#' cov_fun: the closure used to calculate the values.
#' it's arguments are numerical vectors x and y, both of the same length like
#' the input_dimension of the used gp-object.
#' @export
#' @seealso set_cov for setting the covariance function of the gp-object
#' @examples
#' p <- new.gp() %>%
#'      add_data(1:3, 1:3)
#' f <- get_cov(p)
#' f$cov_fun(12, 34)
#'
get_cov <- function(obj){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  return(list(name = get_cov_name(obj), cov_fun = obj$get_cov()))
}


#' Getting the data of a gp
#'
#' @param obj a instance of class "gp"
#' @param df a logical value, whether the output should be a data.frame (df = T)
#' or as a listed object.
#'
#' @return either a data.frame or a list of points
#' @export
#'
#' @examples
#' p1 <- new.gp() %>% add_data(matrix(c(1:10, 2:11, 3:12),nrow=3), 1:10)
#' get_data(p1)
#' #output is a data.frame
#'
#' get_data(p1, FALSE)
#' #output is a list of listed points
#'
#'
get_data <- function(obj, df = T){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  if(df){
    data <- obj$get_data()
    X_learn <- data$X_learn
    y <- data$y_learn
    name_dim <- paste0("x", seq_len(obj$get_input_dim()))
    data_tab <- NULL
    for(item in X_learn){
      item <- as.list(item)
      names(item) <- name_dim
      if(is.null(data_tab))
        data_tab <- as.data.frame(item)
      else
        data_tab <- rbind(data_tab, as.data.frame(item))
    }
    data_tab <- cbind(data_tab, as.data.frame(y))
    return(data_tab)
  }
  else return(obj$get_data())
}

#' Getting parameters of a gp instance
#'
#' @param obj a instance of class "gp"
#' @param used logical vector, whether just the parameters used by the current
#' covariance function should be returned, FALSE per default
#'
#' @return named list of numerical vectors
#' @export
#' @details If you are interested in the parameters that are used by the different
#' available covariance functions, visit the vignette of this package
#' @examples
#' p1 <- new.gp() %>% add_data(matrix(c(1:10, 2:11, 3:12),nrow=3), 1:10)
#'
#' get_parameter(p1)
#' #outputs all parameters stored in the gp object
#'
#' get_parameter(p1, TRUE)
#' #just outputs value of the parameter l
get_parameter <- function(obj, used = F){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  pars <- obj$get_parameter()
  if(used){
    used_par <- cov_list[[get_cov_name(obj)]]$parameter
    res <- list()
    for(item in used_par)
      res[[item]] <- pars[[item]]
    return(res)
  }
  return(pars)
}

#' Getting the mean function of a gp object
#'
#' @param obj a instance of class "gp"
#'
#' @return a closure, that returns the prior mean value
#' @seealso set_mean_fun for further information about the mean function
#' @export
#'
#' @examples
#' p1 <- new.gp()
#' f <- get_mean_fun(p1)
#'
#' f(c(1,2,3))
#' # equals zero per default
#'

get_mean_fun <- function(obj){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  return(obj$get_mean_fun())
}

#' Getting the value of the noise
#'
#' @param obj  a instance of class "gp"
#'
#' @return value of the expected standard deviation of the noise according to the
#' data of the "gp" object
#'
#' @seealso set_noise() for setting this value
#' @export
#'
#' @examples
#' p <- new.gp()
#'
#' get_noise(p)
#' #default value
#'
get_noise <- function(obj){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  return(obj$get_noise())
}

#' Getting the value of the logarithmic marginal likelihood
#'
#' @param obj  a instance of class "gp"
#'
#' @return numerical value of the logarithmic marginal likelihood
#' @export
#' @details The marginal likelihood of the posteriori according to the
#' data added to the gp object. This value is highly depending on the
#' covariance function and the parameters, that are used. By maximizing
#' this value, you get one possible approach of optimzing the hyperparameters
#' of the used Gaussian Process model. If you are interested in further
#' informations about optimizing hyperparameters, see also optimize_gp
#'
#' @examples
#' p <- new.gp() %>%  add_data(1:10, 1:10)
#'
#' get_log_marginal_likelihood(p)
#'
get_log_marginal_likelihood <- function(obj){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  return(obj$get_log_marginal_likelihood())
}

#' Getting prediction using Gaussian Process Model
#'
#' @param obj a instance of class "gp"
#' @param x_input numerical vector, it's length has to equal the dimension of the
#' learning data of the gp object. It's not possible to handle more than one input
#' per function call.
#'
#' @return a list with three items:
#' f_pred: predicted value at position x_input
#' var_f : the variance of the result
#' log_marginal_likelihood: the log_marginal_likelihood of the gp object according
#' to the overall setup
#'
#' @seealso
#' predict_gauss(): if you do not want to use the class 'gp'
#' Notice: We really recommend using the 'gp' class due to efficiency in the case
#' of predicting more than one value.
#'
#' @export
#'
#' @examples
#' p <- new.gp() %>%  add_data(1:10, 1:10)
#' get_prediction(p,11)
#'
get_prediction <- function(obj, x_input){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  if(!is.numeric(x_input))
    stop("Cannot handle non-numeric inputs")
  if(!(length(x_input)==obj$get_input_dim()))
    stop(
      stringr::str_glue("Dimension of the input vector [{length(x_input)}] does not fit
                        to dimension of the data in the gp object [{obj$get_input_dim()}]"))
  return(obj$get_prediction(x_input))
}


#' Setting covariance function
#'
#' @param obj a instance of class "gp"
#' @param cov_fun a character. It has to be one of the following strings:
#' "constant", "linear", "squared_exp", "exponential","gamma_exp", "rational"
#'
#' @return invisible return of the modified object. Be aware, that all
#' setter-functions work on reference.
#' @export
#' @details
#' If you are interested in additional information about the different
#' covariance functions, you can visit the vignette of this package.
#' The main idea of those covariance functions is to describe, how the data
#' we already know is influencing the prediction of the value at some other
#' place.
#' If you are interested in optimizing the choice of the hyperparameter, try out
#' optimize_gp().
#'
#' @examples
#' p1 <- new.gp()
#' # the default value
#' get_cov(p1, "gamma_exp")
#'
#' #modifing the vallue
#' set_cov(p1,"gamma_exp")
#' get_cov(p1)
#'
#'
set_cov <- function(obj, cov_fun){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  #further checks are done in init_cov
  obj$set_cov(cov_fun)
  invisible(obj)
}

#' Setting the standard deviation of a gp instance
#'
#' @param obj a instance of class "gp"
#' @param noise a non-negative numerical value, that describes the standard deviation
#' of the noise
#'
#' @return invisible return of the modified object. Be aware, that all
#' setter-functions work on reference.
#' @export
#'
#' @examples
#' p1 <- new.gp()
#' # the default value
#' get_noise(p1)
#'
#' #modifing the vallue
#' set_noise(p1,2)
#' get_noise(p1)
#'
#'
set_noise <- function(obj, noise){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  #further checks are done in obj$set_noise(noise)
  obj$set_noise(noise)
  invisible(obj)
}

#' Setting the parameters of the used covariance function
#'
#' @param obj a instance of class "gp"
#' @param sigma either a numerical vector of length input_dim or
#' a named list of numerical vectors, such that the names of the list elements are
#' a subset of the names of the parameters settable by this function
#' @param l a numerical value greater than zero
#' @param alpha a numerical value greater than zero
#' @param sigma0 a numerical value
#' @param gamma a numerical value between 0 and 2
#'
#' @return invisible return of the modified object. Be aware, that all
#' setter-functions work on reference.
#' @export
#' @details
#' Those parameters can be used to modify the used covariance functions.
#' Here is a list, which covariance function refers to which parameter
#' constant: sigma0
#' linear: sigma
#' squared_exp: l
#' exponential: l
#' rational: l, alpha
#' gamma_exp: l, gamma
#' For further information about the influence of each parameter, check out
#' the vignette of this package.
#' If you are interested in optimizing the choice of the hyperparameter, try out
#' optimize_gp().
#'
#' @examples
#' p <- new.gp()
#' add_data(p, 1:10,1:10)
#' get_parameter(p)    #default-values
#'
#' set_parameter(p, sigma = 12, gamma = 1.5)
#' get_parameter(p)   #modified values
#'
#' #updating via list
#' par_list <- list( sigma0 = 3, alpha = 2, l = 13)
#' set_parameter(p, par_list)
#' get_parameter(p)
#'
#'
set_parameter <- function(obj, sigma=NULL, l=NULL, alpha=NULL, sigma0=NULL, gamma=NULL){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  if(!is.null(sigma)){
    if(is.list(sigma)){
      l <- sigma[["l"]]
      alpha <- sigma[["alpha"]]
      sigma0 <- sigma[["sigma0"]]
      gamma <- sigma[["gamma"]]
      sigma <- sigma[["sigma"]]
    }
  }
  if(!is.null(l)){
    if(length(l)!=1) stop("l needs length of 1")
    if(l<=0| !is.numeric(l)) stop("l has to be a positive numeric value")
  }
  if(!is.null(alpha)){
    if(length(alpha)!=1) stop("alpha needs length of 1")
    if(alpha<=0| !is.numeric(alpha)) stop("alpha has to be a positive numeric value")
  }

  if(!is.null(gamma)){
    if(length(gamma)!=1) stop("gamma needs length of 1")
    if(gamma<=0| gamma >2|!is.numeric(gamma)) stop("gamma has to be a positive numeric value between 0 and 2")
  }
  if(!is.null(sigma0)){
    if(length(sigma0)!=1) stop("sigma0 needs length of 1")
    if(!is.numeric(sigma0)) stop("sigma0 has to be a numeric value")
  }

  if(!is.null(sigma))
    if(!is.numeric(sigma)| !(length(sigma)==obj$get_input_dim()))
      stop(stringr::str_glue("sigma has to be a numeric vector with length {obj$get_input_dim()}"))

  obj$set_parameter(sigma,l, alpha, sigma0, gamma)
  invisible(obj)
}


#' Adding a mean function to a gp object
#'
#' @param obj  a instance of class "gp"
#' @param mean_fun either a numeric value or a closure, that just takes one input vector
#' and generates one-dimensional numeric output
#'
#' @return invisible return of the modified object. Be aware, that all
#' setter-functions work on reference.
#' @export
#'
#' @examples
#' p <- new.gp()
#' add_data(p, 1:10,1:10)
#' set_mean(p, function(x) x)
#'
set_mean_fun <- function(obj, mean_fun){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  #further testing is implemented in the class function
  obj$set_mean_fun(mean_fun)
  invisible(obj)
}

#' Optimizing hyperparameters of gp instance
#'
#' @param obj a instance of class "gp"
#'
#' @return invisible return of the modified object.
#'
#' @details
#' Further details about the used method.
#' Optimiziation based on the maximization of the marginal likelihood.
#' @export
#'
#' @examples
#'
#'
#'
#'
optimize_gp <- function(obj){
  if(!("gp" %in% class(obj)))
    stop("obj has to be a member of class 'gp'")
  obj$optim_parameter()
  invisible(obj)
}


