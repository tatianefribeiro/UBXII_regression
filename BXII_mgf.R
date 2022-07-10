rm(list = ls())
c = 0.5
d = 2.1
s = 1
t = -50
R = 50
par <- c(c,d,s,t,R)

BXIImgf<-function(par){
  par[1]->c; par[2]->d; par[3]->s; par[4]->t; par[5]->R
  Expfunc<-function(n,z){
    temp <- function(t) exp(-z*t)/t^n
    return(as.numeric(integrate(temp,1,Inf)[1]))
  }
  vec <- vector()
  for(i in 0:(R-1)){
    vec[i+1]<-choose(-d-1,i)*(
      (-s*t)^(-c*(1+i))*(pgamma(-s*t,c*(1+i))*gamma(c*(1+i)))+(-s*t)^(c*(d+1))*
        Expfunc(1+c*(d+i),-s*t)
    )  
  }
  return(c*d*sum(vec))
  print(c*d*sum(vec))
}

pdf_BXII <- function(x){
  exp(t*x)*
  c*d*s^(-c)*x^(c-1)*(1+(x/s)^c)^(-d-1)
}

integrate(pdf_BXII,0,Inf)
BXIImgf(par)