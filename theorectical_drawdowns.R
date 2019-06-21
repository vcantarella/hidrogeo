time <- c(0.1, 1, 10.001, 10.1, 11, 20) #Time in sec
u_1 <- 1^2*0.00025/(4*(73/(24*60*60))*time)
u_2 <- 1^2*0.00025/(4*(73/(24*60*60))*(time-10))
library(pracma)


drawdown_1 <- ((150/(24*60*60))/(4*pi*(73/(24*60*60))))*expint_E1(u_1)

drawdown_2 <- ((150/(24*60*60))/(4*pi*(73/(24*60*60))))*expint_E1(u_1)+
  ((-150/(24*60*60))/(4*pi*(73/(24*60*60))))*expint_E1(u_2)

drawdown_2 <- as.numeric(drawdown_2)

drawdown <- c(drawdown_1[1:3], drawdown_2[4:6])

drawdown

library(ggplot2)
qplot(time, drawdown, geom = "smooth")

30.5*sqrt(0.015/(16.8*3.05))

integrand <- function(y) {(1/y)* exp(-y- 0.522^2/(4*y))}


(42.5/(4*pi*18.6))

a <- integrate(integrand, lower = 0.36032, upper = Inf)
hj_draw <- sapply(u, function(x){
  a <- integrate(integrand, lower = x, upper = Inf)
  a <- a[[1]]
  return(a*(42.5/(4*pi*18.6)))
})
hj_draw

c_pi <- (30.5/(4*sqrt(18.6*0.0002)))*sqrt(0.015*0.003/3.05)



h_60_draw <- sapply(u, function(x){
  integran <- function(y) {(1/y)*exp(-y)*erfc(c_pi*sqrt(x)/sqrt(y*(y-x)))}
  a <- integrate(integran, lower = x, upper = Inf)
  a <- a[[1]]
  return(a*(42.5/(4*pi*18.6)))
})
h_60_draw
