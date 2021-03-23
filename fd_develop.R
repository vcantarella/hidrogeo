#P 3.4
#First project to understand Finite Differences solution to partial differential equations of
#Darcy's flow.

#20 m thick confined aquifer isotropic and homogeneous
#K = 100 m/d
#S = 2e-5
#two wells pump at constant rate of 4m3/min and 12 m3/min for 30 days
#800 x 500 m grid used for this problem.

##Method 01: Using Theis solution.

#Considering the origin at the lower-left corner of the map.
#then well A is at (250,250) and B is at (550,250)
xs <- seq(0,800,by = 50)
ys <- seq(0,500, by = 50)

100*20

u_f <- function(d,time){d^2*2e-5/(4*2000*time)} #U function



library(pracma)

drawdown_Theis <- function(Q,u){(Q/(4*pi*2000))*expint_E1(u)} #Theis drawdown function

M_Theis = matrix(rep(0,length(xs)*length(ys)),ncol =length(xs))
for(i in 1:length(xs)){
  for(j in 1:length(ys)){
    d_a = sqrt((250-xs[i])^2+(250-ys[j])^2)
    d_a = ifelse(d_a == 0,0.5,d_a)
    d_b = sqrt((550-xs[i])^2+(250-ys[j])^2)
    d_b = ifelse(d_b == 0, 0.5,d_b)
    M_Theis[j,i] = drawdown_Theis(4*60*24,u_f(d_a,30))+drawdown_Theis(12*60*24,u_f(d_b,30))
    }
}

M_Theis = 100-M_Theis




library(plotly)

plot_ly(x = xs, y = ys, z = M_Theis, type = "contour",width = 400, height = 250,
               autocontour = F, contours = list(start = 80, end = 100, size = 1))




#Finite Differences model. The boundary conditions 
##are general flow boundaries defined on the outside of the border cells.
#For the model to work, I had to define one border which is wet and the other borders are no flow.




M_Theis
matrix(as.vector(M_Theis),ncol = ncol(M_Theis))


#Modelo com 4 blocos:
#(1,1) = 1,(2,1) = 2, (1,2) = 3, (2,2) = 4
#HCOF, CR,CC,0
#-CR,HCOF,0,CC
#-CC,0,HCOF,CR
#0,-CC,-CR,HCOF

#Modelo com 9 blocos:
#(1,1) = 1,(2,1) = 2, (3,1) = 3,(1,2) = 4, (2,2) = 5, (3,2) = 6,
#(1,3) = 7, (2,3) = 8, (3,3) = 9
#HCOF, CR,0,CC,0,0,0,0,0
#-CR,HCOF,0,CC
#-CC,0,HCOF,CR
#0,-CC,-CR,HCOF

#Desisto de ser inteligente, vou fazer o modelo a partir da força bruta:

#20 m thick confined aquifer isotropic and homogeneous
#K = 100 m/d
#Ss = 2e-5/20
#50 m grid

#CR = K*20*50/50
#CC = CR
#HCOF = -S*20*50*50/(1(deltat = 1 dia)*20)
#RHS = HCOF*h(t-1)
#W = 12*60*24 at W2 and 4*60*24 at w1.

#Pumping finite modeling will be calculated using the same specs of the first model.
##I'm not sure if the model is correct, however weird values were obtained as the model completely drains some cells. So I am going to expand the model further.
##Same error appeared, suggesting the problem is not on the grid size. I am making a mistake somewhere.

M_par = expand.grid(seq(1,11),seq(1,17))
names(M_par) <- c("row","column")
head(M_par)
M_par <- M_par%>% 
  mutate(CR = (100/(60*24))*20) %>% 
  mutate(CC = CR) %>% 
  mutate(h_0 = 100) %>% 
  mutate(W = 0) %>%
  mutate(bound = 0) %>% 
  mutate(HCOF = -2e-5*50*50/1440)
  

#Setting up pumping rates:
M_par$W[M_par$row == 6 & M_par$column == 6] = -4
M_par$W[M_par$row == 6 & M_par$column == 12] = -12
library(tidyverse)

m_row = max(M_par$row)
#mounting the matrix:

A = matrix(rep(0,nrow(M_par)^2),ncol = nrow(M_par))
for(i in 1:nrow(M_par)){
  r = M_par[i,]$row
  c = M_par[i,]$column
  c1 <- M_par %>% filter(row == r+1) %>% 
    filter(column ==c) %>% select(CC) %>% pull(.)
  if(is_empty(c1)) {
    c1 = 0
  }
  if(c1!=0) A[i,(m_row*(c-1)+r+1)] = c1
  c2 <- M_par %>% filter(row == r-1) %>% 
    filter(column ==c) %>% select(CC) %>% pull(.)
  if(is_empty(c2)){c2 = 0}
  if(c2!=0) A[i,(m_row*(c-1)+r-1)] = c2
  
  c3 <- M_par %>% filter(row == r) %>% 
    filter(column ==c+1) %>% select(CR) %>% pull(.)
  if(is_empty(c3)) {c3 = 0}
  if(c3!=0) A[i,(m_row*(c)+r)] = c3
  c4 <- M_par %>% filter(row == r) %>% 
    filter(column ==c-1) %>% select(CR) %>% pull(.)
  if(is_empty(c4)) {c4 = 0}
  if(c4!=0) A[i,(m_row*(c-2)+r)] = c4
  c5 <- M_par %>% filter(row == r) %>% 
    filter(column ==c) %>% select(HCOF, CR, CC) %>% 
    mutate(CR = -2*CR) %>% 
    mutate(CC = -2*CC) %>% rowSums(.)
  if(r==1 | r==m_row){
    crn <- M_par %>% filter(row == r) %>% 
      filter(column ==c) %>% select(CR) %>% 
      pull(.)
    crn <- crn - 0.0005275*50
    c5 <- c5+crn
  }
  if(c==1 | c==max(M_par$column)){
    ccn <- M_par %>% filter(row == r) %>% 
      filter(column ==c) %>% select(CC) %>% 
      pull(.)
    ccn <- ccn - 0.0005275*50
    c5 <- c5+ccn
  }
  A[i,i] = c5
}

for(d in seq(1,30,by = 1)){
  M_par <- M_par%>% 
    mutate(RHS = HCOF*h_0) %>%
    mutate(bound = 0) %>% 
    mutate(bound = ifelse(column == 1, bound - 100*0.0005275*50,bound)) %>%
    mutate(bound = ifelse(column == max(M_par$column), bound - 100*0.0005275*50,bound)) %>% 
    mutate(bound = ifelse(row == 1, bound-100*0.0005275*50,bound)) %>% 
    mutate(bound = ifelse(row == max(row), bound-100*0.0005275*50,bound))
    
  #Right hand side of the equation:
  b = M_par$RHS - M_par$W+ M_par$bound
  
  
 

  #Solve:
  h_1 = solve(A,b)
  
  M_par$h_0 = h_1
}

plot_ly(x = xs, y = ys, z = matrix(h_1,ncol = 17), type = "contour",width = 400, height = 250,
        autocontour = F, contours = list(start = 80, end = 100, size = 1))


#Full Boundary conditions:
#mutate(bound = ifelse(row == 1, bound-h_0*CR,bound)) %>% 
# mutate(bound = ifelse(row == max(row), bound-h_0*CR,bound)) %>% 
#


#P3.5: Single well pumped in steady-state homogeneous isotropic confined aquifer:

#formula: h = (Q/(2*pi*T))*log(r/r_0)+h_0

#From point (0,0) as r_0, at each point, r = r_0 + sqrt(x^2+y^2), h_0 = 6.82.
#So I need to find r_0, T and Q. And so it will be done using the known points around the square:

#I can do it for all points and then find the least sum of square errors:

t_values <- expand.grid(x = seq(0,300,100),
                       y = seq(0,300,100))

head <- c(6.82,7.56,7.99,8.29,7.19,NA,NA,8.33,7.68,NA,NA,8.41,8.04,8.18,8.36,8.53)
t_values <- t_values %>% 
  mutate(head = head)

t_points <- t_values %>% filter(!(x == 0 & y == 0)) %>% 
  filter(!is.na(head))

#head_func calculates the mean square error for points to find best approximation to parameters:
head_func <- function(pars){
  h = (pars[1]/(2*pi*pars[2]))*log((sqrt((t_points$x-pars[3])^2 + (t_points$y-pars[4])^2))/sqrt((0-pars[3])^2+(0-pars[4]^2)))+6.82
  h_real = t_points$head
  return(sum((h-h_real)^2)/nrow(t_points))
}

h = (4/(2*pi*2000))*log((3000+sqrt((0-t_points$x)^2 + (0-t_points$y^2)))/3000)+6.82


optim(par = c(4,2000,-1200,0),head_func, method = "Nelder-Mead")

t_values %>% 
  filter(is.na(head)) %>% 
  mutate(head = (455.8485253/(2*pi*68.4706368))*log(sqrt((x+99.4032385)^2 + (y-0.1741757)^2)/sqrt(99.2337795^2+0.1741757^2))+6.82)


#Write the Finite Difference approximation:

