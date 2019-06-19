time <- c(0.00694, 0.04167, 0.2083, 1)
u <- 30.5^2*0.0002/(4*18.6*time)
u
library(pracma)
drawdown <- (42.5/(4*pi*18.6))*expint_E1(u)
drawdown             

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
