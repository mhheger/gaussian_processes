library(shiny)

#' Starting Shiny App to Work with GP-Objects
#' @param obj if you want to pass a existing `gp` instance to the app
#' @description Just run gp_app() to start the Shiny-App in your RStudio or
#' browser
#' @export
#'
#' @examples
#' library(dplyr)
#' p <- new.gp() %>%
#'   add_data(data.frame(x1 = 1:10, x2 = 2:11), y = sin(3:12))
#' gp_app(p)
#'
#'@import shiny
#'@import shinyFeedback
#'@importFrom stats na.omit
#'@importFrom utils read.csv
#'@importFrom stringr str_glue
#'@import DT
gp_app <- function(obj=NULL){
  if(!is.null(obj))
    if(!("gp" %in% class(obj)))
      stop("class of object is not `gp`")

  # Defining UI -------------------------------------------------------------
  parameter_tabs <- tabsetPanel(
    id = "params",
    type = "hidden",
    tabPanel("constant",
             numericInput("sigma0", "sigma0", value = 1),
    ),
    tabPanel("linear",
             uiOutput("input_sigma")
    ),
    tabPanel("squared_exp",
             numericInput("l_sqr", "l", value = 1, min = 0.00000001),
    ),
    tabPanel("exponential",
             numericInput("l_exp", "l", value = 1, min = 0.00000001),
    ),
    tabPanel("gamma_exp",
             numericInput("l_gamma", "l", value = 1, min = 0.00000001),
             numericInput("gamma", "gamma", value = 1, min= 0.00000001, max = 2)
    ),
    tabPanel("rational",
             numericInput("l_rat", "l", value = 1, min = 0.00000001),
             numericInput("alpha", "alpha", value =1, min = 0.00000001)
    ),

  )

  ui <- fluidPage(
    shinyFeedback::useShinyFeedback(),
    titlePanel("Gaussian Processes"),
    sidebarLayout(
      sidebarPanel(
        # selectInput("modus",
        #             "Decide between regression or classification",
        #             choices = list(
        #               regression = "regression",
        #               classification = "classification"
        #             )
        # ),
        selectInput("cov",
                    "Choose for a covariance function",
                    choices = list(
                      "squared exponential" = "squared_exp",
                      "constant" = "constant",
                      "linear" = "linear",
                      "exponential" = "exponential",
                      "gamma exponential" = "gamma_exp",
                      "rational" = "rational"
                    )
        ),
        numericInput("noise",
                     "Variance of the noise:",
                     min= 0,
                     value = 1,
                     step = 0.1),
        helpText("Here you can change the value of the used parameters"),
        parameter_tabs,
        actionButton("set_parameters",
                     "Set Parameters!"
        ),
        actionButton("optimize",
                     "Optimize!"
        ),
        checkboxInput("keepcov",
                      "Keep choosen covariance function",
        )

      ),
      mainPanel(
        tabsetPanel(
          tabPanel("input data",
                   fluidRow(
                     column(width=4,
                            helpText("To handle your input data correctly,",
                                     "we need the correct input dimension."),
                            numericInput("input_dim",
                                         "What's the input dimension?",
                                         value = 1,
                                         min = 1,
                                         max = Inf,
                                         step = 1,
                            ),
                            actionButton("set_dim", "Set Dimension!"),
                     ),
                     column(width = 4,
                            offset = 2,
                            helpText("If you want to upload data via",
                                     ".csv-file, you can do this here."),
                            fileInput("input_file", "Upload .csv file:", accept = ".csv")
                     )
                   ),
                   fluidRow(
                     column(width = 4,
                            helpText("Here you can insert the coordinates of the",
                                     "points, you want to add as learning data. "),
                            uiOutput("input_data"),
                            actionButton("add_data", "Add Data!")
                     ),
                     column(width = 4,
                            offset = 2,
                            helpText("If you decide for uploading your data as a ",
                                     ".csv-file, then you can map the variables to",
                                     "to the columns of the uploaded data.frame."),
                            uiOutput("selecting_columns"),
                            actionButton("add_data_file", "Add Data!")
                     )
                   ),
                   fluidRow(
                     DT::DTOutput("data", width = "80%")
                   )
          ),
          tabPanel("get prediction",
                   fluidRow(
                     column(width = 4,
                            uiOutput("input_pred"),
                            actionButton("get_pred", "Get Prediction!")
                     ),
                     column(width = 4,
                            helpText("If you want to get a prediction for a certain point,",
                                     "just insert the coordinates of the point and click the",
                                     "button below.")
                     )
                   ),
                   verbatimTextOutput("prediction")
          ),
          tabPanel("plot",
                   plotOutput("plot",
                              dblclick = dblclickOpts("dbl_click"),
                              click = clickOpts("click_plot")),
                   fluidRow(
                     column(width = 2,  numericInput("range_min", "Min" , value = -100)),
                     column(width = 4,
                            sliderInput("xrange",
                                        "Adjust the x-range of the plot",
                                        min = -100,
                                        max = 100,
                                        value = c(0,10),
                                        step = 1
                            )
                     ),
                     column(width = 2,  numericInput("range_max", "Max" , value = 100))
                   ),
                   verbatimTextOutput("points_to_add"),
                   fluidRow(
                     column(2, actionButton("add_points_plot", "Add Point!")),
                     column(2, downloadButton("download_plot", "Download plot as")),
                     column(2, selectInput("format",label = NULL, choices = list(
                       "PDF" = ".pdf",
                       "PNG" = ".png"))
                     )
                   )

          )
        )
      )
    )
  )


  # Defining server function ------------------------------------------------
  server <- function(input,output, session){

    # global values -----------------------------------------------------------
    p1 <- gp$new()  #gp object
    X_data <- reactiveVal() #kind of stack variable of all x inputs
    y_data <- reactiveVal() #kind of stack variable of all y inputs
    X <- reactiveVal() #current x input
    y <- reactiveVal() #current y input
    prediction <- reactiveVal() #current predicted value
    update_plot <-reactiveVal(0) #variable, that causes the plot updates
    curr_plot <- reactiveVal()
    set_input_dim <- reactiveVal()
    # working with obj
    input_element <- reactiveVal()
    used_obj <- reactiveVal(FALSE)

    # handling of non-trivial input:  -----------------------------------------
    if(!is.null(obj)){
      dim <- obj$get_input_dim()
      if(!is.null(dim))
        input_element(dim)
      p1 <- obj
      data <- get_data(obj)
      if(!is.null(data)){
        X_data(subset(data, select = - c(y)))
        y_data(unlist(subset(data, select = c(y))))
      }
    }


    # render elements------------------------------------------------------------
    output$input_data <- renderUI({
      l <- c()
      for(i in seq_len(input_dim())){
        s <- paste0("x", as.character(i))
        l <- c(l,s)
      }
      l <- c(l,"y")
      res <- lapply(l,function(x) numericInput(x,x,value = 0))
    })

    output$selecting_columns <- renderUI({
      if(!is.null(file_data())){
        col_names <- colnames(file_data())
        l <- c()
        for(i in seq_len(input_dim())){
          s <- paste0("x", as.character(i))
          l <- c(l,s)
        }
        l <- c(l,"y")
        res <- lapply(l,function(x) selectInput(paste0(x,"_file"),x,choices = col_names))
      }
    })

    output$data <- DT::renderDT(
      {
        df <- cbind(X_data(), data.frame(y=y_data()))
        df
      }
    )
    output$plot <- renderPlot({
      update_plot()
      if(input_dim()==1 & !is.null(X_data())){
        note <- showNotification("Plotting...", duration = NULL, closeButton = FALSE)
        curr_plot(p1$plot(input$xrange[1],input$xrange[2]))
        on.exit(removeNotification(note), add = TRUE)
        isolate(curr_plot())
      }

    }
    )

    output$points_to_add <- renderPrint({
      if(input_dim()==1)
        print(c(X(),y()))
    }
    )

    output$input_pred <- renderUI({
      l <- c()
      for(i in seq_len(input_dim())){
        s <- paste0("x", as.character(i), "_pred")
        l <- c(l,s)
      }
      res <- lapply(l,function(x) numericInput(x,x,value = 0))
    })

    output$prediction <- renderPrint({
      if(!is.null(X_data()) & !is.null(prediction()))
        s <- str_glue("According to the setting of the Gaussian Process,
                    we get to the following prediction:

                    - predicted value: {prediction()$f_predict}

                    - variance {prediction()$var_f}

                    with log marginal likelihood:
                    - {prediction()$log_marginal_likelihood}
                    ")
      else
        s <- "No prediction without data"
      print(s)
    }
    )

    output$input_sigma <- renderUI({
      l <- c()
      for(i in seq_len(input_dim())){
        s <- paste0("sigma", as.character(i))
        l <- c(l,s)
      }
      res <- lapply(l,function(x) numericInput(x,x,value = 0))
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        paste("plot-", Sys.Date(), input$format , sep="")
      },
      content = function(file){
        ggsave(file,plot=curr_plot())
      }
    )

    #reactive events & reactives ------------------------------------------------
    input_dim <- eventReactive(
      c(input$set_dim, input_element()),
      {
        if(!is.null(input_element()) & !used_obj()){
          set_input_dim(TRUE)
          val <- input_element()
          used_obj(TRUE)
          return(val)
        }
        y_data(NULL)
        X_data(NULL)
        p1 <<- gp$new()
        if(input$input_dim>=1 & is.numeric(input$input_dim)){
          shinyFeedback::hideFeedback("input_dim")
          set_input_dim(T)
          round(input$input_dim)
        }
        else{
          shinyFeedback::feedbackDanger("input_dim", T, "Input dimension has to be
                                      a positive integer")
          1
        }
      }
    )

    test_sigma <- reactive(
      is.numeric(sigma())
    )
    sigma <- reactive({
      sigma <- c()
      for (i in seq_len(input_dim())){
        s <- paste0("sigma", as.character(i))
        sigma <- c(sigma, input[[s]])
      }
      sigma
    })

    test_l <- reactive({
      test_l <- l() > 0 & is.numeric(l())
      shinyFeedback::feedbackDanger(get_l_id(), !test_l, "l has to be positive!")
      test_l
    })

    l <- reactive({
      input[[get_l_id()]]
    }
    )
    test_sigma0 <- reactive({
      is.numeric(sigma0())
    })

    sigma0 <- reactive({
      input$sigma0
    })

    test_alpha <- reactive({
      alpha <- input$alpha
      test_alpha <- alpha > 0 & is.numeric(alpha)
      shinyFeedback::feedbackDanger("alpha", !test_alpha, "alpha has to be positive!")
      test_alpha
    })

    alpha <- reactive(
      input$alpha
    )

    test_gamma <- reactive({
      gamma <- input$gamma
      test_gamma <- 0 < gamma & gamma <= 2 & is.numeric(gamma)
      shinyFeedback::feedbackDanger("gamma", !test_gamma, "gamma has to be between 0 and 2!")
      test_gamma
    })

    gamma <- reactive(
      input$gamma
    )

    test_noise <- reactive({
      test_noise <- (is.numeric(input$noise))&(input$noise >=0)
      shinyFeedback::feedbackDanger("noise", !test_noise, "Variance of noise has to be positive!")
      test_noise
    })

    noise <- reactive(
      input$noise
    )

    file_data <- reactive({
      inFile <- input$input_file
      if (is.null(inFile)) return(NULL)
      data <- read.csv(inFile$datapath, header = TRUE)
      na.omit(data)
    }
    )

    #observing events------------------------------------------------------------

    #input of data points on the input tab via activeButton
    observeEvent(list(input$add_data),
                 {
                   if(is.null(set_input_dim())){
                     showNotification("You have to set dimension first!", type = "error")
                   } else {
                     l <- c()
                     for(i in seq_len(input_dim())){
                       s <- paste0("x", as.character(i))
                       l <- c(l,s)
                     }
                     res <- lapply(l, function(x){
                       sol <- input[[x]]
                       updateNumericInput(inputId = x, value = 0)
                       return(sol)
                     })
                     names(res) <- l
                     sol <- input$y
                     updateNumericInput(inputId = "y", value = 0)
                     if(all(!is.na(unlist(res))) & all(!is.na(sol))){
                       X(res)
                       y(sol)
                       add_data()
                     } else {
                       showNotification("Only numeric input is allowed!", type = "error")
                     }
                   }
                 })

    #input of data stored in csv_file after choosing names
    observeEvent(input$add_data_file,
                 { #add checks if input is numerical
                   if(is.null(file_data())){
                     showNotification("You have to upload a file!", type = "error")
                   } else if (is.null(set_input_dim())){
                     showNotification("You have to set dimension first!", type = "error")
                   } else if(test_noise()){
                     var_names <- paste0("x", as.character(seq_len(input_dim())))
                     selected_cols <- sapply(paste0(var_names, "_file"), function(x){
                       input[[x]]
                     })
                     if(check_cols(c(var_names, "y"))){
                       file <- file_data()
                       X_data_to_insert <- file[selected_cols]
                       y_data_to_insert <- file[input[["y_file"]]]
                       y_data_to_insert <- unname(unlist(y_data_to_insert))
                       colnames(X_data_to_insert) <- var_names
                       X_data(rbind(X_data(), X_data_to_insert))
                       y_data(c(y_data(), y_data_to_insert))
                       p1$add_data(X_data(), y_data(), noise())
                     } else {
                       showNotification("You have to select columns with numeric input",
                                        type = "error")
                     }
                   }
                 })
    #input selected data points using plot and activeButton
    observeEvent(list(input$add_points_plot, input$click),
                 {
                   add_data()
                 })



    #marking points in the plot via double clicking
    observeEvent(input$dbl_click,
                 {
                   X(list(x1=unname(input$dbl_click$x)))
                   y(c(y=unname(input$dbl_click$y)))
                 }
    )

    #dynamic UI element: just showing relevant parameters
    observeEvent(input$cov, {
      updateTabsetPanel(inputId = "params", selected = input$cov)
    })

    #setting parameters via user input
    observeEvent(input$set_parameters,{
      if(is.null(set_input_dim())){
        showNotification("You have to set dimension first!", type = "error")
      } else if (is.null(X_data())) {
        showNotification("You have to add data first!", type = "error")
      } else {
        note <- showNotification("Setting Parameters...", duration = NULL, closeButton = FALSE)
        p1$set_cov(input$cov)
        checks <- c (test_sigma(),test_sigma0(),test_alpha(), test_gamma(), test_l(), test_noise())

        if(all(checks)){
          p1$set_parameter(sigma= sigma(), l = l(), alpha = alpha(), gamma = gamma(), sigma0= sigma0())
          p1$set_noise(noise())
          update_plot(update_plot()+1)
        }
        on.exit(removeNotification(note), add = TRUE)
      }
    })

    #getting prediction for point input
    observeEvent(input$get_pred, {
      if(is.null(set_input_dim())){
        showNotification("You have to set dimension first!", type = "error")
      } else if (is.null(X_data())) {
        showNotification("You have to add data first!", type = "error")
      } else {
        note <- showNotification("Getting Prediction...", duration = NULL, closeButton = FALSE)
        input_names <- c()
        for(i in seq_len(input_dim())){
          s <- paste0("x", as.character(i), "_pred")
          input_names <- c(input_names,s)
        }
        pos <- sapply(input_names, function(x) input[[x]])
        if(all(!is.na(pos))){
          prediction(p1$get_prediction(pos))
        } else {
          showNotification("Input has to be numeric!", type = "error")
        }
        on.exit(removeNotification(note), add = TRUE)
      }
    })

    #optimizing parameters
    observeEvent(input$optimize, {
      if(is.null(set_input_dim())){
        showNotification("You have to set dimension first!", type = "error")
      } else if (is.null(X_data())) {
        showNotification("You have to add data first!", type = "error")
      } else if(test_noise()) {
        note <- showNotification("Optimizing...", duration = NULL, closeButton = FALSE)
        p1$set_noise(noise())
        para <- p1$optim_parameter()
        # if covariance function should be kept
        if(input$keepcov){
          p1$set_cov(input$cov)
          sigma <- para[[input$cov]][["sigma"]]
          alpha <- para[[input$cov]][["alpha"]]
          gamma <- para[[input$cov]][["gamma"]]
          sigma0 <- para[[input$cov]][["sigma0"]]
          l <- para[[input$cov]][["l"]]
          p1$set_parameter(sigma=sigma, alpha = alpha, gamma = gamma, l = l, sigma0 = sigma0)
        }
        updateSelectInput(inputId = "cov", selected = p1$get_cov_name())
        opt_para <- p1$get_parameter()
        for(item in names(opt_para)){
          if(item != "sigma" & item != "l"){
            updateNumericInput(inputId = item, value = unname(opt_para[[item]]))
          }
          else if(item == "l"){
            updateNumericInput(inputId = get_l_id(), value = unname(opt_para[[item]]))
          }
          else if(item == "sigma"){
            sapply(seq_len(input_dim()),function(x){
              name <- paste0("sigma", as.character(x))
              updateNumericInput(inputId = name, value = unname(opt_para$sigma[x]))})
          }
        }
        update_plot(update_plot()+1)
        on.exit(removeNotification(note), add = TRUE)
      }
    })


    #adjusting the range of the slider
    observeEvent(list(input$range_max, input$range_min),{
      if(is.numeric(input$range_max) & is.numeric(input$range_min)){
        if(input$range_max > input$range_min){
          updateSliderInput(inputId="xrange", min = input$range_min, max = input$range_max)
        }
      } else{
        showNotification("Only numeric input allowed!", type = "error")
      }
    })

    #internal functions ---------------------------------------------------------

    #adding data to the "stack" and the gp-object
    add_data <- function(){
      if(is.null(set_input_dim())){
        showNotification("You have to set dimension first!", type = "error")
      } else if(test_noise()){
        note <- showNotification("Adding data", duration = NULL, closeButton = FALSE)
        if(is.null(X_data())){
          if(!is.null(X())){
            X_data(data.frame(X()))
          }
        } else {
          X_data(rbind(X_data(),data.frame(X())))
        }
        y_data(c(y_data(), y()))
        xx <- X()
        yy <- y()
        if(!is.null(xx) & !is.null(yy)){
          p1$add_data(unlist(xx), yy, noise())
        }
        on.exit(removeNotification(note), add = TRUE)
        update_plot(update_plot()+1)
      }
    }

    #getting current l
    get_l_id <- function(){
      l_names <- c(squared_exp ="l_sqr", gamma_exp = "l_gamma", exponential = "l_exp",
                   rational = "l_rat", linear = "l_sqr", constant = "l_sqr")
      current_cov <- input$cov
      l_names[[current_cov]]
    }

    #checking if selected cols are numeric
    check_cols <- function(IDs){
      file <- file_data()
      IDs <- paste0(IDs, "_file")
      check_passed <- TRUE
      for(item in IDs){
        if(!is.numeric(unlist(file[input[[item]]]))){
          check_passed <- FALSE
          shinyFeedback::feedbackDanger(item,TRUE, "Selected column is not numeric!")
        } else{
          shinyFeedback::hideFeedback(item)
        }
      }
      return(check_passed)
    }
  }
  shinyApp(ui,server)
}


# Copied server function to run tests, ignoring obj
testing_server <- function(input,output, session){

  # global values -----------------------------------------------------------
  p1 <- gp$new()  #gp object
  X_data <- reactiveVal() #kind of stack variable of all x inputs
  y_data <- reactiveVal() #kind of stack variable of all y inputs
  X <- reactiveVal() #current x input
  y <- reactiveVal() #current y input
  prediction <- reactiveVal() #current predicted value
  update_plot <-reactiveVal(0) #variable, that causes the plot updates
  curr_plot <- reactiveVal()
  set_input_dim <- reactiveVal()
  # working with obj
  input_element <- reactiveVal()
  used_obj <- reactiveVal(FALSE)

  # # handling of non-trivial input:  -----------------------------------------
  # if(!is.null(obj)){
  #   dim <- obj$get_input_dim()
  #   if(!is.null(dim))
  #     input_element(dim)
  #   p1 <- obj
  #   data <- get_data(obj)
  #   if(!is.null(data)){
  #     X_data(subset(data, select = - c(y)))
  #     y_data(unlist(subset(data, select = c(y))))
  #   }
  # }


  # render elements------------------------------------------------------------
  output$input_data <- renderUI({
    l <- c()
    for(i in seq_len(input_dim())){
      s <- paste0("x", as.character(i))
      l <- c(l,s)
    }
    l <- c(l,"y")
    res <- lapply(l,function(x) numericInput(x,x,value = 0))
  })

  output$selecting_columns <- renderUI({
    if(!is.null(file_data())){
      col_names <- colnames(file_data())
      l <- c()
      for(i in seq_len(input_dim())){
        s <- paste0("x", as.character(i))
        l <- c(l,s)
      }
      l <- c(l,"y")
      res <- lapply(l,function(x) selectInput(paste0(x,"_file"),x,choices = col_names))
    }
  })

  output$data <- DT::renderDT(
    {
      df <- cbind(X_data(), data.frame(y=y_data()))
      df
    }
  )
  output$plot <- renderPlot({
    update_plot()
    if(input_dim()==1 & !is.null(X_data())){
      note <- showNotification("Plotting...", duration = NULL, closeButton = FALSE)
      curr_plot(p1$plot(input$xrange[1],input$xrange[2]))
      on.exit(removeNotification(note), add = TRUE)
      isolate(curr_plot())
    }

  }
  )

  output$points_to_add <- renderPrint({
    if(input_dim()==1)
      print(c(X(),y()))
  }
  )

  output$input_pred <- renderUI({
    l <- c()
    for(i in seq_len(input_dim())){
      s <- paste0("x", as.character(i), "_pred")
      l <- c(l,s)
    }
    res <- lapply(l,function(x) numericInput(x,x,value = 0))
  })

  output$prediction <- renderPrint({
    if(!is.null(X_data()) & !is.null(prediction()))
      s <- stringr::str_glue("According to the setting of the Gaussian Process,
                    we get to the following prediction:

                    - predicted value: {prediction()$f_predict}

                    - variance {prediction()$var_f}

                    with log marginal likelihood:
                    - {prediction()$log_marginal_likelihood}
                    ")
    else
      s <- "No prediction without data"
    print(s)
  }
  )

  output$input_sigma <- renderUI({
    l <- c()
    for(i in seq_len(input_dim())){
      s <- paste0("sigma", as.character(i))
      l <- c(l,s)
    }
    res <- lapply(l,function(x) numericInput(x,x,value = 0))
  })

  output$download_plot <- downloadHandler(
    filename = function() {
      paste("plot-", Sys.Date(), input$format , sep="")
    },
    content = function(file){
      ggsave(file,plot=curr_plot())
    }
  )

  #reactive events & reactives ------------------------------------------------
  input_dim <- eventReactive(
    c(input$set_dim, input_element()),
    {
      if(!is.null(input_element()) & !used_obj()){
        set_input_dim(TRUE)
        val <- input_element()
        used_obj(TRUE)
        return(val)
      }
      y_data(NULL)
      X_data(NULL)
      p1 <<- gp$new()
      if(input$input_dim>=1 & is.numeric(input$input_dim)){
        shinyFeedback::hideFeedback("input_dim")
        set_input_dim(T)
        round(input$input_dim)
      }
      else{
        shinyFeedback::feedbackDanger("input_dim", T, "Input dimension has to be
                                      a positive integer")
        1
      }
    }
  )

  test_sigma <- reactive(
    is.numeric(sigma())
  )
  sigma <- reactive({
    sigma <- c()
    for (i in seq_len(input_dim())){
      s <- paste0("sigma", as.character(i))
      sigma <- c(sigma, input[[s]])
    }
    sigma
  })

  test_l <- reactive({
    test_l <- l() > 0 & is.numeric(l())
    shinyFeedback::feedbackDanger(get_l_id(), !test_l, "l has to be positive!")
    test_l
  })

  l <- reactive({
    input[[get_l_id()]]
  }
  )
  test_sigma0 <- reactive({
    is.numeric(sigma0())
  })

  sigma0 <- reactive({
    input$sigma0
  })

  test_alpha <- reactive({
    alpha <- input$alpha
    test_alpha <- alpha > 0 & is.numeric(alpha)
    shinyFeedback::feedbackDanger("alpha", !test_alpha, "alpha has to be positive!")
    test_alpha
  })

  alpha <- reactive(
    input$alpha
  )

  test_gamma <- reactive({
    gamma <- input$gamma
    test_gamma <- 0 < gamma & gamma <= 2 & is.numeric(gamma)
    shinyFeedback::feedbackDanger("gamma", !test_gamma, "gamma has to be between 0 and 2!")
    test_gamma
  })

  gamma <- reactive(
    input$gamma
  )

  test_noise <- reactive({
    test_noise <- (is.numeric(input$noise))&(input$noise >=0)
    shinyFeedback::feedbackDanger("noise", !test_noise, "Variance of noise has to be positive!")
    test_noise
  })

  noise <- reactive(
    input$noise
  )

  file_data <- reactive({
    inFile <- input$input_file
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath, header = TRUE)
    na.omit(data)
  }
  )

  #observing events------------------------------------------------------------

  #input of data points on the input tab via activeButton
  observeEvent(list(input$add_data),
               {
                 if(is.null(set_input_dim())){
                   showNotification("You have to set dimension first!", type = "error")
                 } else {
                   l <- c()
                   for(i in seq_len(input_dim())){
                     s <- paste0("x", as.character(i))
                     l <- c(l,s)
                   }
                   res <- lapply(l, function(x){
                     sol <- input[[x]]
                     updateNumericInput(inputId = x, value = 0)
                     return(sol)
                   })
                   names(res) <- l
                   sol <- input$y
                   updateNumericInput(inputId = "y", value = 0)
                   if(all(!is.na(unlist(res))) & all(!is.na(sol))){
                     X(res)
                     y(sol)
                     add_data()
                   } else {
                     showNotification("Only numeric input is allowed!", type = "error")
                   }
                 }
               })

  #input of data stored in csv_file after choosing names
  observeEvent(input$add_data_file,
               { #add checks if input is numerical
                 if(is.null(file_data())){
                   showNotification("You have to upload a file!", type = "error")
                 } else if (is.null(set_input_dim())){
                   showNotification("You have to set dimension first!", type = "error")
                 } else if(test_noise()){
                   var_names <- paste0("x", as.character(seq_len(input_dim())))
                   selected_cols <- sapply(paste0(var_names, "_file"), function(x){
                     input[[x]]
                   })
                   if(check_cols(c(var_names, "y"))){
                     file <- file_data()
                     X_data_to_insert <- file[selected_cols]
                     y_data_to_insert <- file[input[["y_file"]]]
                     y_data_to_insert <- unname(unlist(y_data_to_insert))
                     colnames(X_data_to_insert) <- var_names
                     X_data(rbind(X_data(), X_data_to_insert))
                     y_data(c(y_data(), y_data_to_insert))
                     p1$add_data(X_data(), y_data(), noise())
                   } else {
                     showNotification("You have to select columns with numeric input",
                                      type = "error")
                   }
                 }
               })
  #input selected data points using plot and activeButton
  observeEvent(list(input$add_points_plot, input$click),
               {
                 add_data()
               })



  #marking points in the plot via double clicking
  observeEvent(input$dbl_click,
               {
                 X(list(x1=unname(input$dbl_click$x)))
                 y(c(y=unname(input$dbl_click$y)))
               }
  )

  #dynamic UI element: just showing relevant parameters
  observeEvent(input$cov, {
    updateTabsetPanel(inputId = "params", selected = input$cov)
  })

  #setting parameters via user input
  observeEvent(input$set_parameters,{
    if(is.null(set_input_dim())){
      showNotification("You have to set dimension first!", type = "error")
    } else if (is.null(X_data())) {
      showNotification("You have to add data first!", type = "error")
    } else {
      note <- showNotification("Setting Parameters...", duration = NULL, closeButton = FALSE)
      p1$set_cov(input$cov)
      checks <- c (test_sigma(),test_sigma0(),test_alpha(), test_gamma(), test_l(), test_noise())

      if(all(checks)){
        p1$set_parameter(sigma= sigma(), l = l(), alpha = alpha(), gamma = gamma(), sigma0= sigma0())
        p1$set_noise(noise())
        update_plot(update_plot()+1)
      }
      on.exit(removeNotification(note), add = TRUE)
    }
  })

  #getting prediction for point input
  observeEvent(input$get_pred, {
    if(is.null(set_input_dim())){
      showNotification("You have to set dimension first!", type = "error")
    } else if (is.null(X_data())) {
      showNotification("You have to add data first!", type = "error")
    } else {
      note <- showNotification("Getting Prediction...", duration = NULL, closeButton = FALSE)
      input_names <- c()
      for(i in seq_len(input_dim())){
        s <- paste0("x", as.character(i), "_pred")
        input_names <- c(input_names,s)
      }
      pos <- sapply(input_names, function(x) input[[x]])
      if(all(!is.na(pos))){
        prediction(p1$get_prediction(pos))
      } else {
        showNotification("Input has to be numeric!", type = "error")
      }
      on.exit(removeNotification(note), add = TRUE)
    }
  })

  #optimizing parameters
  observeEvent(input$optimize, {
    if(is.null(set_input_dim())){
      showNotification("You have to set dimension first!", type = "error")
    } else if (is.null(X_data())) {
      showNotification("You have to add data first!", type = "error")
    } else if(test_noise()) {
      note <- showNotification("Optimizing...", duration = NULL, closeButton = FALSE)
      p1$set_noise(noise())
      para <- p1$optim_parameter()
      # if covariance function should be kept
      if(input$keepcov){
        p1$set_cov(input$cov)
        sigma <- para[[input$cov]][["sigma"]]
        alpha <- para[[input$cov]][["alpha"]]
        gamma <- para[[input$cov]][["gamma"]]
        sigma0 <- para[[input$cov]][["sigma0"]]
        l <- para[[input$cov]][["l"]]
        p1$set_parameter(sigma=sigma, alpha = alpha, gamma = gamma, l = l, sigma0 = sigma0)
      }
      updateSelectInput(inputId = "cov", selected = p1$get_cov_name())
      opt_para <- p1$get_parameter()
      for(item in names(opt_para)){
        if(item != "sigma" & item != "l"){
          updateNumericInput(inputId = item, value = unname(opt_para[[item]]))
        }
        else if(item == "l"){
          updateNumericInput(inputId = get_l_id(), value = unname(opt_para[[item]]))
        }
        else if(item == "sigma"){
          sapply(seq_len(input_dim()),function(x){
            name <- paste0("sigma", as.character(x))
            updateNumericInput(inputId = name, value = unname(opt_para$sigma[x]))})
        }
      }
      update_plot(update_plot()+1)
      on.exit(removeNotification(note), add = TRUE)
    }
  })


  #adjusting the range of the slider
  observeEvent(list(input$range_max, input$range_min),{
    if(is.numeric(input$range_max) & is.numeric(input$range_min)){
      if(input$range_max > input$range_min){
        updateSliderInput(inputId="xrange", min = input$range_min, max = input$range_max)
      }
    } else{
      showNotification("Only numeric input allowed!", type = "error")
    }
  })

  #internal functions ---------------------------------------------------------

  #adding data to the "stack" and the gp-object
  add_data <- function(){
    if(is.null(set_input_dim())){
      showNotification("You have to set dimension first!", type = "error")
    } else if(test_noise()){
      note <- showNotification("Adding data", duration = NULL, closeButton = FALSE)
      if(is.null(X_data())){
        if(!is.null(X())){
          X_data(data.frame(X()))
        }
      } else {
        X_data(rbind(X_data(),data.frame(X())))
      }
      y_data(c(y_data(), y()))
      xx <- X()
      yy <- y()
      if(!is.null(xx) & !is.null(yy)){
        p1$add_data(unlist(xx), yy, noise())
      }
      on.exit(removeNotification(note), add = TRUE)
      update_plot(update_plot()+1)
    }
  }

  #getting current l
  get_l_id <- function(){
    l_names <- c(squared_exp ="l_sqr", gamma_exp = "l_gamma", exponential = "l_exp",
                 rational = "l_rat", linear = "l_sqr", constant = "l_sqr")
    current_cov <- input$cov
    l_names[[current_cov]]
  }

  #checking if selected cols are numeric
  check_cols <- function(IDs){
    file <- file_data()
    IDs <- paste0(IDs, "_file")
    check_passed <- TRUE
    for(item in IDs){
      if(!is.numeric(unlist(file[input[[item]]]))){
        check_passed <- FALSE
        shinyFeedback::feedbackDanger(item,TRUE, "Selected column is not numeric!")
      } else{
        shinyFeedback::hideFeedback(item)
      }
    }
    return(check_passed)
  }
}







