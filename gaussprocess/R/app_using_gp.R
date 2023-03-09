library(shiny)


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
           numericInput("l", "l", value = 1, min = 0.00000001),
  ),
  tabPanel("exponential",
           numericInput("l", "l", value = 1, min = 0.00000001),
  ),
  tabPanel("gamma_exp",
           numericInput("l", "l", value = 1, min = 0.00000001),
           numericInput("gamma", "gamma", value = 1, min= 0.00000001, max = 2)
  ),
  tabPanel("rational",
           numericInput("l", "l", value = 1, min = 0.00000001),
           numericInput("alpha", "alpha", value =1, min = 0.00000001)
  ),

)

ui <- fluidPage(
  titlePanel("Gaussian Processes"),
  sidebarLayout(
    sidebarPanel(
      selectInput("modus",
                  "Decide between regression or classification",
                  choices = list(
                    regression = "regression",
                    classification = "classification"
                  )
      ),
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
                   "Standard deviation of the noise:",
                   min= 0,
                   value = 1,
                   step = 0.1),
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
                 numericInput("input_dim",
                              "What's the input dimension?",
                              value = 1,
                              min = 1,
                              max = Inf,
                              step = 1,
                 ),
                 actionButton("set_dim", "Set Dimension!"),
                 fluidRow(width = 8),
                 fluidRow(
                   column(width = 4,
                          h2("Manual data input"),
                          uiOutput("input_data"),
                          actionButton("add_data", "Add Data!")
                          ),
                   column(width = 4,
                          offset = 2,
                          fileInput("input_file", "Upload .csv file:", accept = ".csv"),
                          uiOutput("selecting_columns"),
                          actionButton("add_data_file", "Add Data!")
                          )
                 ),
                 fluidRow(
                   tableOutput("data")
                 )
        ),
        tabPanel("get prediction",
                 uiOutput("input_pred"),
                 actionButton("get_pred", "Get Prediction!"),
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
server <- function(input,output){

  # global values -----------------------------------------------------------
  p1 <- gp$new()  #gp object
  X_data <- reactiveVal() #kind of stack variable of all x inputs
  y_data <- reactiveVal() #kind of stack variable of all y inputs
  X <- reactiveVal() #current x input
  y <- reactiveVal() #current y input
  prediction <- reactiveVal("") #current predicted value
  update_plot <-reactiveVal(0) #variable, that causes the plot updates
  curr_plot <- reactiveVal()

  # render elments------------------------------------------------------------
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
    col_names <- colnames(file_data())
    l <- c()
    for(i in seq_len(input_dim())){
      s <- paste0("x", as.character(i))
      l <- c(l,s)
    }
    l <- c(l,"y")
    res <- lapply(l,function(x) selectInput(paste0(x,"_file"),x,choices = col_names))
  })

  output$data <- renderTable(
    {
      df <- cbind(X_data(), data.frame(y=y_data()))
      df
    }
  )
  output$plot <- renderPlot({
    c(update_plot())
    if(input_dim()==1)
      curr_plot(p1$plot(input$xrange[1],input$xrange[2]))
      isolate(curr_plot())
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
    print("The prediction of the function value is")
    print(prediction())
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
    input$set_dim,
    {
      y_data(NULL)
      X_data(NULL)
      p1 <<- gp$new()
      input$input_dim
    }
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
                 l <- c()
                 for(i in seq_len(input$input_dim)){
                   s <- paste0("x", as.character(i))
                   l <- c(l,s)
                 }
                 res <- lapply(l, function(x){
                   sol <- input[[x]]
                   updateNumericInput(inputId = x, value = 0)
                   return(sol)
                 })
                 names(res) <- l
                 X(res)

                 sol <- input$y
                 updateNumericInput(inputId = "y", value = 0)
                 y(sol)
                 add_data()
               })

  #input of data stored in csv_file after choosing names
  observeEvent(input$add_data_file,
               { #add checks if input is numerical
                 var_names <- paste0("x", as.character(seq_len(input_dim())))
                 selected_cols <- sapply(paste0(var_names, "_file"), function(x){
                   input[[x]]
                 })
                 file <- file_data()
                 X_data_to_insert <- file[selected_cols]
                 y_data_to_insert <- file[input[["y_file"]]]
                 y_data_to_insert <- unname(unlist(y_data_to_insert))
                 colnames(X_data_to_insert) <- var_names
                 X_data(rbind(X_data(), X_data_to_insert))
                 y_data(c(y_data(), y_data_to_insert))
                 p1$add_data(X_data(), y_data(), input$noise)
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
    p1$set_cov(input$cov)
    sigma0 <- input$sigma0
    alpha <- input$alpha
    gamma <- input$gamma
    l <- input$l
    sigma <- c()
    for (i in seq_len(input$input_dim)){
      s <- paste0("sigma", as.character(i))
      sigma <- c(sigma, input[[s]])
    }
    p1$set_parameter(sigma= sigma, l = l, alpha = alpha, gamma = gamma, sigma0= sigma0)
    p1$set_noise(input$noise)
    update_plot(update_plot()+1)
  })

  #getting prediction for point input
  observeEvent(input$get_pred, {
    for(i in seq_len(input_dim())){
      s <- paste0("x", as.character(i), "_pred")
      l <- c(l,s)
    }
    pos <- sapply(s, function(x) input[[x]])
    prediction(p1$get_prediction(pos))
  })

  #optimizing parameters
  observeEvent(input$optimize, {
    p1$set_noise(input$noise)
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
    } else{
      updateSelectInput(inputId = "cov", selected = p1$get_cov_name())
    }
    updateSelectInput(inputId = "cov", selected = p1$get_cov_name())
    opt_para <- p1$get_parameter()
    for(item in names(opt_para)){
      if(item != "sigma")
        updateNumericInput(inputId = item, value = unname(opt_para[[item]]))
      else
        sapply(seq_len(input_dim()),function(x){
          name <- paste0("sigma", as.character(x))
          updateNumericInput(inputId = name, value = unname(opt_para$sigma[x]))})
    }
    update_plot(update_plot()+1)
  })


  #adjusting the range of the slider
  observeEvent(list(input$range_max, input$range_min),{
    if(input$range_max > input$range_min){
      updateSliderInput(inputId="xrange", min = input$range_min, max = input$range_max)
    }
  })

  #internal functions ---------------------------------------------------------

  #adding data to the "stack" and the gp-object
  add_data <- function(){
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
      p1$add_data(unlist(xx), yy, input$noise)
    }
    update_plot(update_plot()+1)
  }
}


shinyApp(ui,server)
