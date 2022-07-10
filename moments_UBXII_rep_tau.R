####################################################################
###############         Re-parametrized UBXII       ################
#Author: Tatiane F. Ribeiro
#Date: 20th February, 2020.
#Subject: Moments
####################################################################

####################   First re-parameterization  #################
rm(list = objects())

#mgf BXII
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/UBXIIreg2022/GitHub_UBXII_regression")
source("BXII_mgf.R")

#\tauth-quantile
tau = 0.5

#Parameters
q = 0.1  #quantile  in (0,1)
c = 1  #shape param (prec or disp??) in (0,Infty)

d = log(1/tau)/log(1+(log(1/q))^c)  #re-par.
s = 1                               #scale parameter.
h = 2
t = -h
R = 50

par <- c(c,d,s,t,R)

hthmom_UBXII_rtau <- function(y){
  y^h*
  c*log(1/tau)*(log(1/y))^(c-1)/
    (y*log(1+(log(1/q))^c)*(1+(log(1/y))^c)^
       (1+log(1/tau)/log(1+(log(1/q))^c)))
}

######################################
#E(Y^h) = M_X(-h); X~BXII(c,d,s)
integrate(hthmom_UBXII_rtau,0,1)  #analytic
BXIImgf(par)                       ##expansion

