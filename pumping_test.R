#Pumping Test Trial:

#Pumping rate
q_l_min <- 132
q_m3_s <- 132*10^-3/60

library(optimx)
#Hantush Jacob Well-Function Optimizer:
hantush_jacob_optim <- function(Tra = 2/(24*60*60), Sa = 0.00002, Kz = 1.4e-2/(24*60*60), r, Q, b, time, drawdown){
  optimizer <- function(par){
    Tr = par[1]
    S = par[2]
    Cond = par[3]
    Alpha = r*sqrt(Kz/(Tr*b))
    u <- r^2*S/(4*Tr*time)
    hj_draw <- sapply(u, function(x){
      integrand <- function(y) {(1/y)* exp(-y- Alpha^2/(4*y))}
      a <- integrate(integrand, lower = x, upper = Inf)
      a <- a[[1]]
      return(a*(Q/(4*pi*Tr)))
    })
    RMSE_draw <- (log(hj_draw)-log(drawdown))^2
    RMSE_draw <- sum(RMSE_draw)/length(RMSE_draw)
    return(RMSE_draw)
  }
  optimx(par = c(Tra,Sa,Kz), fn = optimizer, method = "L-BFGS",lower = rep(1e-9,3), upper = rep(1e-3,3),
         control = list(maxit = 5000,all.methods = T))
  
}

time <- c(1,2,4,8,15,30,60,120,240,360,480,600,720)
time <- time*60
drawdown <- c(0.37,0.88,1.5,2.07,2.35,2.52,2.58,2.59,2.59,2.59,2.59,2.59,2.59)
hantush_jacob_optim(r = 32.6, Q = q_m3_s, b = 4.3, time = time, drawdown = drawdown)


Q <- q_m3_s
Tr = 0.0001007087
S = 2.234887e-05
Kz = 1.62037e-07
r = 32.6
Alpha = r*sqrt(Kz/4.3/Tr)
u <- r^2*S/(4*Tr*time)
hj_draw <- sapply(u, function(x){
  integrand <- function(y) {(1/y)* exp(-y- Alpha^2/(4*y))}
  a <- integrate(integrand, lower = x, upper = Inf)
  a <- a[[1]]
  return(a*(Q/(4*pi*Tr)))
})
library(ggplot2)
ggplot(data = data.frame(hj_draw, time), aes(x = log(time), y = log(hj_draw)))+geom_line()+
  geom_point(data = data.frame(drawdown, time), aes(x = log(time), y = log(drawdown)))

#Optimizer worked better with time in seconds so transmissivity would be of similar scaling than Storativity. 
#It is still hard to find a unique solution to the problem because the solution fails to predict conductance of the upper bound 
#correctly.

#Theis this time:
library(pracma)
theis_optim <- function(Tra = 2/(24*60*60), Sa = 0.00002, r, Q, time, drawdown){
  optimizer <- function(par){
    Tr = par[1]
    S = par[2]
    u <- r^2*S/(4*Tr*time)
    theis_draw <- (Q/(4*pi*Tr))*expint_E1(u)
    RMSE_draw <- (log(theis_draw)-log(drawdown))^2
    RMSE_draw <- sum(RMSE_draw)/length(RMSE_draw)
    return(RMSE_draw)
  }
  optimx(par = c(Tra,Sa), fn = optimizer, method = "L-BFGS",lower = rep(1e-9,2),
         control = list(maxit = 5000,all.methods = T))
  
}

time <- c(1,2,4,8,15,30,60,120,240)
time <- time*60
drawdown <- c(0.15,0.22,0.3,0.39,0.46,0.55,0.63,0.72,0.81)
theis_optim(r = 95, Q = 1.3/60, time = time, drawdown = drawdown)

#Perfect Match!
  
Q <- 1.3/60
Tr = 0.01387339
S = 7.631842e-05
r = 95
u <- r^2*S/(4*Tr*time)
theis_draw <- (Q/(4*pi*Tr))*expint_E1(u)
library(ggplot2)
ggplot(data = data.frame(theis_draw, time), aes(x = log(time), y = log(theis_draw)))+geom_line()+
  geom_point(data = data.frame(drawdown, time), aes(x = log(time), y = log(drawdown)))


