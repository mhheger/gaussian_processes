#Beispiel für die Funktionalität

#einfache Kernelfunktion 
f <- function(x,y) {      
  return (squared_exp_cov(sqrt(sum(abs(x-y)^2)),1))
  #return(0.3)
}

#einfache Plotfunktion
plot_with_confidence_band <- function(x,y,variance){
  polygon(
    c(x, rev(x)), 
    c(y+variance, rev(y)-rev(variance)),
    col = "green",
    border = NA, 
    density = 20
  )
  lines(x,y, lwd = 2)
}

#Generierung der Lerndaten 
g <- function(x) sign(x)*(0.001*x^2 +3)
x<- runif(20, -10, 10)
y <- g(x) + rnorm(20, 0 , 0.1)

#Hilfsfunktionen zur Auswertung 
eval_gauss_value <- function(z){
  predict_gauss(x,y,f, 0.1,z)$f_predict
}
eval_gauss_variance<- function(z){
  predict_gauss(x,y,f, 0.1,z)$var_f
}

#Auswertung an Testpunkten 
x_t <- seq(-10,10, 0.5)
y_test <- sapply(x_t, eval_gauss_value)
var_test <- sapply(x_t, eval_gauss_variance)

#Plotten der Ergebnisse 
plot(x,y)
lines(sort(x_t),g(x_t))
plot_with_confidence_band(x_t,y_test , var_test)




