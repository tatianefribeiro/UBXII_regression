########################################################
# UBXII, Kw, UW, and beta regressions - USEFUL FUNCTIONS 
########################################################
#LINK FUNCTION
lfunc <- function(beta_vec,X){  
  qi_logit_link <- 1/(1+exp(-(X%*%beta_vec))) #=logit inverse div by exp(eta_i)
  return(qi_logit_link)    #Note that here eta_i = b1+b2*x_{i2} because x2 is
  #a n-dim vector with n obs from covariate X2.
}
#################
# UBXII functions
################
r_UBXII <- function(n)   #It generates occurences of y_i ~ UBXII (q_i, c)
{
  u = runif(n)
  y_UBXII = exp(-(u^(-1/(log(1/tau)/log(1+(log(1/q_i))^c)))-1)^(1/c)) #UBXII qf
  return(y_UBXII)
}

l_UBXIIreg <- function(par){  #Log-likelihood function from UBXII regression model
  beta_vec <- par[-length(par)]  #it receives the vector par components less the last that is c
  c <- par[length(par)]           #shape parameter (c)
  q_i <- lfunc(beta_vec,X)  #It call the link function
  
  ell =  n*log(-c*log(tau))- 
    sum(log(y))+
    (c-1)*sum(log(log(1/y)))-
    sum(log(1+(log(1/y))^c))-
    sum(log(log(1+(log(1/q_i))^c)))-
    log(1/tau)*sum((log(1+(log(1/y))^c))/(log(1+(log(1/q_i))^c)))
  return(ell)
}

cdf_UBXII <- function(y,q,c)   #UBXII cdf (useful for qqplot)
{
  (1+(log(1/y))^c)^(log(tau)/log(1+(log(1/q))^c))
}

l_UBXIIreg_aug <- function(par){  #Log-likelihood function from UBXII regression model
  beta_vec <- par[-length(par)]  #it receives the vector par components less the last that is c
  c <- par[length(par)]           #shape parameter (c)
  q_i <- lfunc(beta_vec,X_aug)  #It call the link function
  
  ell =  n*log(-c*log(tau))- 
    sum(log(y))+
    (c-1)*sum(log(log(1/y)))-
    sum(log(1+(log(1/y))^c))-
    sum(log(log(1+(log(1/q_i))^c)))-
    log(1/tau)*sum((log(1+(log(1/y))^c))/(log(1+(log(1/q_i))^c)))
  return(ell)
}

###############
# KW functions
##############
r_Kw <- function(n)   #It generates occurences of y_i ~ Kw (q_i, c)~Kw (q_i, c)
{
  u = runif(n)
  y_Kw = (1-(1-u)^(log(1-q_i^(1/c))/log(.5)))^c #  Kw qf
  return(y_Kw)
}

l_Kwreg <- function(par){ # Kw Log-likelihood 
  beta_vec <- par[-length(par)]  
  c <- par[length(par)]          
  q_i <- lfunc(beta_vec,X)   
  
  ell =  sum(
    log(
      log(.5)/(c*log(1-q_i^(1/c)))*
        y^(1/c-1)*(1-y^(1/c))^
        (log(0.5)/log(1-q_i^(1/c))-1)
    )
  )
  
  return(ell)
}

# Kw log-likelihood (augmented model)
l_Kwreg_aug <- function(par){  
  beta_vec <- par[-length(par)]  
  c <- par[length(par)]          
  q_i <- lfunc(beta_vec,X_aug)   
  
  ell =  sum(
    log(
      log(.5)/(c*log(1-q_i^(1/c)))*
        y^(1/c-1)*(1-y^(1/c))^
        (log(0.5)/log(1-q_i^(1/c))-1)
    )
  )
  
  return(ell)
}

###############
# UW functions
##############
pdf_UW <- function(y,q_i,c){  # UW pdf
  c/y*(log(tau)/log(q_i))*(log(y)/log(q_i))^(c-1)*
    tau^((log(y)/log(q_i))^c)
}

cdf_UW <- function(y,q,b)   # UW cdf
{
  tau^((log(y)/log(q))^b)
}

r_UW <- function(n)    #It generates occurences of y_i ~ UW (q_i, \beta)
{
  u = runif(n)
  y_UW = exp(log(q_i)*((log(u)/log(tau))^(1/c)))
  return(y_UW)
}

l_UWreg <- function(par){   # UW log-likelihood
  beta_vec <- par[-length(par)]
  c <- par[length(par)]
  
  q_i <- lfunc(beta_vec,X)
  ell =  sum(log(pdf_UW(y,q_i,c)))
  return(ell)
}


l_UWreg_aug <- function(par){    # UW log-likelihood (augmented model)
  beta_vec <- par[-length(par)]
  b <- par[length(par)]
  
  q_i <- lfunc(beta_vec,X_aug)
  ell =sum(log(pdf_UW(y,q_i,b)))
  return(ell)
}


#################
# Beta functions
#################
pdf_beta <- function(y,q_i,c){  # UW pdf
  p = q_i*c
  q = (1-q_i)*c
  dbeta(y,p,q)
}
#integrate(pdf_beta,0,1)

r_beta <- function(n)    #It generates occurrences of y_i ~ UW (q_i, \beta)
{
  p = q_i*c
  q = (1-q_i)*c
  rbeta(n,p,q)
}

l_betareg <- function(par){   # beta log-likelihood
  beta_vec <- par[-length(par)]
  c <- par[length(par)]
  
  q_i <- lfunc(beta_vec,X)
  ell =  sum(log(pdf_beta(y,q_i,c)))
  return(ell)
}


l_betareg_aug <- function(par){    # beta log-likelihood (augmented model)
  beta_vec <- par[-length(par)]
  b <- par[length(par)]
  
  q_i <- lfunc(beta_vec,X_aug)
  ell = sum(log(pdf_beta(y,q_i,c)))
  return(ell)
}

