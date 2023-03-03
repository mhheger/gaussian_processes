R6::R6Class("gp",
            public = list(
              add_data = function(x,y, noise=0) {
                l <- convert_to_list(x, length(y))
                if(!is.null(private$input_dimension)){
                  if(length(l[[1]])!= private$input_dimension){
                    stop("new input cannot be added to the current input data.
                         You may have to initalize a new gp instance.")
                  }
                } else {
                  private$input_dimension <- length(l[[1]])
                }
                private$X_learn <- c(unlist(private$X_learn), l)
                private$y_learn <- c(private$y_learn, y)
                self$set_noise(noise)
                self$set_cov("squared_exp")
                self$update_marginal_likelihood()
              },
              #getter functions for the private values
              get_K = function() private$K,
              get_cov = function() private$cov,
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
              get_log_marginal_likelihood = function() self$log_marginal_likelihood,

              #setter functions for the private values
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
                #self$cov <- cov_list[cov]
                #self$parameter <- par_list[cov]
                private$K <- cov_cross(private$X_learn, private$X_learn, private$cov)
                self$update_marginal_likelihood()
              },
              set_parameter = function(sigma=NULL, l=NULL, alpha=NULL, sigma0=NULL, gamma=NULL) {
                # need to add checks in the real use functions
                if(!is.null(sigma)) if(!is.na(sigma)) private$sigma <- sigma
                if(!is.null(l)) if(!is.na(l)) private$l <- l
                if(!is.null(alpha)) if(!is.na(alpha)) private$alpha <- alpha
                if(!is.null(sigma0)) if(!is.na(sigma0)) private$sigma0 <- sigma0
                if(!is.null(gamma)) if(!is.na(gamma)) private$gamma <- gamma
                self$set_cov(private$covname)
                private$K <- cov_cross(private$X_learn, private$X_learn, private$cov)
                self$update_marginal_likelihood()
              },
              set_noise = function(noise){
                if(noise<0| length(noise)!=1 | !is.numeric(noise)) stop("noise has to be a positive double value")
                private$noise <- noise
              },
              set_mean_fun = function(f){
                if((!typeof(f)=="closure" & !is.numeric(f)) | length(f) !=1)
                  stop("mean_fun has to be either a function or a numeric vector
                       of length 1")
                if(typeof(f)=="closure" & length(args(f))>1)
                  stop("mean function is just allowed to work with one input variable")
                private$mean_fun <- f
              },
              optim_parameter = function() {optimize_parameters(self)} ,
              plot = function(x_start, x_end, name_x="x", name_y="y",mode = T,title = NULL, col= "violet") {
                if(private$input_dimension==1){
                  range_x <- seq(x_start,x_end,len = 100)
                  f <-c()
                  var <-c()
                  for(x in range_x){
                    pred <- self$get_prediction(x)
                    f <- c(f, pred$f_predict)
                    var <- c(var, pred$var_f)
                  }
                  x_input <- unlist(private$X_learn)
                  ranged <- x_start <= x_input & x_input <= x_end
                  y_input <- unlist(private$y_learn)
                  input <- data.frame(xx=x_input[ranged],yy=y_input[ranged])
                  plot_data <- data.frame(x = range_x,y = f, min_f= f-var,max_f= f+var )
                  plot <- ggplot2::ggplot(plot_data, aes(x,y))+
                    geom_ribbon(aes(ymin = min_f,ymax =  max_f, fill = "grey60"))+
                    geom_line(aes(x,y))
                  if(mode) plot <- plot + geom_point(data=input, aes(xx,yy))
                  ggplotly(plot)
                }
              },
              update_marginal_likelihood = function(){
                private$log_marginal_likelihood <- self$get_prediction(1)$log_marginal_likelihood
              }


            ),

            private = list(
              K = matrix(),   # variance-covariance-matrix of the input data,
              log_marginal_likelihood = Inf,
              cov = NULL ,
              covname = NULL,
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

plot_with_confidence_band <- function(x,y,variance, color){
  polygon(
    c(x, rev(x)),
    c(y+variance, rev(y)-rev(variance)),
    col = color,
    lty =1,
    border = NA,
    density = 100
  )
  lines(x,y, lwd = 2)
}


