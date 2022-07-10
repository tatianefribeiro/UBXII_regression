####################################################################
###############         Re-parametrized UBXII       ################
#Author: Tatiane F. Ribeiro
#Date: 22th February, 2020.
#Subject: Incomplete moments
####################################################################

####################   First re-parameterization  #################
rm(list = objects())

#mgf BXII
#setwd("")
library(expint)

#\tauth-quantile
tau = 0.5

#Parameters
q = 0.01  #quantile  in (0,1)
c = .2  #shape param (prec or disp??) in (0,Infty)
d = log(1/tau)/log(1+(log(1/q))^c)  
h = 1

R = 200  

inc_mom_UBXII_rtau <- function(y){
  y^h*
    c*log(1/tau)*(log(1/y))^(c-1)/
    (y*log(1+(log(1/q))^c)*(1+(log(1/y))^c)^
       (1+log(1/tau)/log(1+(log(1/q))^c)))
}

exp_inc_mom_UBXII_rtau <- function(z){
  if(z>0 && z<exp(-1)){
    vec = vector()
    for (j in 0:(R-1)) {
      vec[j+1] = choose(-d-1,j)*h^(c*j)*gammainc(-c*(d+j),h*log(1/z))
    }
    aux = c*d*h^(c*d)*sum(vec) 
    return(aux)
  }
  if(z>=exp(-1) && z<=1){
    vec = vector()
    for (j in 0:(R-1)) {
      vec[j+1] = choose(-d-1,j)*(
        h^(-c*(j+1))*(gammainc(c*(j+1),h*log(1/z))-gammainc(c*(j+1),h))+
          h^(c*(d+j))*gammainc(-c*(d+j),h)
      )
    }
    aux = c*d*sum(vec) 
    return(aux)
  }
}

z = 0.2   
integrate(inc_mom_UBXII_rtau,0,z)  #true
exp_inc_mom_UBXII_rtau(z)

##############   Sem expansão    ##############################
int_p1 <- function(u){
  d*exp(-h*u^(1/c))*(1+u)^(-d-1)
}
integrate(int_p1,(log(1/z))^c,Inf)
z=.8
integrate(inc_mom_UBXII_rtau,0,z)  
integrate(int_p1,(log(1/z))^c,1)$value+integrate(int_p1,1,Inf)$value

