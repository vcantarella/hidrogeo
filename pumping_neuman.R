#Equation Newman Solution:
#

u_0_y = (1-exp(-t_s*Beta*(y^2-gamma_0^2)))*cosh(gamma_0*z_D)/
  ((y^2 + (1+sigma)*gamma_0^2 - (y-gamma_0^2)^2/sigma)*cosh(gamma_0))

u_n_y = (1-exp(-t_s*Beta*(y^2-gamma_n^2)))*cosh(gamma_n*z_D)/
  ((y^2 + (1+sigma)*gamma_n^2 - (y-gamma_n^2)^2/sigma)*cosh(gamma_n))

##Find roots of equation:

sigma*gamma_0*sinh(gamma_0)-(y^2-gamma_0)^2*cosh(gamma_0) = 0 #gamma_0^2 < y^2

sigma*gamma_n+(y^2+gamma_n^2)*cos(gamma_n) = 0

(2*n-1)*(pi/2) < gamma_n < n*pi, n>1


s = (Q/(4*pi*T))*integrate(funtion(y){
  4*y*besselJ(y*sqrt(Beta),0)*(u_0_y+)
})