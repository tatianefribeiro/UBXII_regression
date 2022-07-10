###################################################################################
# Useful functions to fit unit regressions in the application study
###################################################################################

# Logit link function
lfunc <- function(beta_vec,X){  
  qi_logit_link <- 1/(1+exp(-(X%*%beta_vec)))
  return(qi_logit_link)
}

###################################################################################
############  UBXII regression model - USEFUL FUNCTIONS   #########################
###################################################################################
r_UBXII <- function(n,q_i,c)   #It generates occurences of z_i ~ UBXII (q_i, c)
{
  u = runif(n)
  z_UBXII = exp(-(u^(-1/(log(1/tau)/log(1+(log(1/q_i))^c)))-1)^(1/c)) #UBXII qf
  return(z_UBXII)
}

lfunc <- function(beta_vec,X){  #logit link function --> MATRIX FORM
  qi_logit_link <- 1/(1+exp(-(X%*%beta_vec)))
  return(qi_logit_link)
}

l_UBXIIreg <- function(par){  #Log-likelihood function from UBXII regression model
  beta_vec <- par[-length(par)]  #it receives the vector par components less the last that is c
  c <- par[length(par)]           #shape parameter (c)
  q_i <- lfunc(beta_vec,X)  #It call the link function
  
  ell =  n*log(-c*log(tau))- 
    sum(log(z))+
    (c-1)*sum(log(log(1/z)))-
    sum(log(1+(log(1/z))^c))-
    sum(log(log(1+(log(1/q_i))^c)))-
    log(1/tau)*sum((log(1+(log(1/z))^c))/(log(1+(log(1/q_i))^c)))
  return(ell)
}

cdf_UBXII <- function(z,q,c)   #UBXII cdf (useful for qqplot)
{
  (1+(log(1/z))^c)^(log(tau)/log(1+(log(1/q))^c))
}

l_UBXIIreg_aug <- function(par){  #Log-likelihood function from UBXII regression model
  beta_vec <- par[-length(par)]  #it receives the vector par components less the last that is c
  c <- par[length(par)]           #shape parameter (c)
  q_i <- lfunc(beta_vec,X_aug)  #It call the link function
  
  ell =  n*log(-c*log(tau))- 
    sum(log(z))+
    (c-1)*sum(log(log(1/z)))-
    sum(log(1+(log(1/z))^c))-
    sum(log(log(1+(log(1/q_i))^c)))-
    log(1/tau)*sum((log(1+(log(1/z))^c))/(log(1+(log(1/q_i))^c)))
  return(ell)
}


#########    NULL MODEL LOG-LIK   ##############3
l_UBXII_noreg <- function(par){ #UBXII log-likelihood function
  c <- par[1]
  q <- par[2]
  
  ell <- n*log(c*log(1/tau))+(c-1)*sum(log(log(1/z)))-
    sum(log(z))-n*log(log(1+(log(1/q))^c))-
    (1+log(1/tau)/log(1+(log(1/q))^c))*
    sum(log(1+(log(1/z))^c))
  return(ell)
}



################################################################################
############  Kw regression model - USEFUL FUNCTIONS   #########################
################################################################################
r_Kw <- function(n,w_i,dp)   #It generates occurences of Z_i ~ Kw (w_i, dp)
{
  u = runif(n)
  z_Kw = (1-(1-u)^(log(1-w_i^(1/dp))/log(.5)))^dp #  Kw qf
  return(z_Kw)
}

l_Kwreg <- function(par){ # Kw Log-likelihood 
  beta_vec <- par[-length(par)]  
  dp <- par[length(par)]          
  w_i <- lfunc(beta_vec,X)   
  
  ell =  sum(
    log(
      log(.5)/(dp*log(1-w_i^(1/dp)))*
        z^(1/dp-1)*(1-z^(1/dp))^
        (log(0.5)/log(1-w_i^(1/dp))-1)
    )
  )
  
  return(ell)
}

cdf_Kw <- function(z,w,dp)   # Kw cdf
{
  1-(1-z^(1/dp))^(log(0.5)/log(1-w^(1/dp)))
}

l_Kwreg_aug <- function(par){  # Kw log-likelihood (augmented model)
  beta_vec <- par[-length(par)]  
  dp <- par[length(par)]          
  w_i <- lfunc(beta_vec,X_aug)   
  
  ell =  sum(
    log(
      log(.5)/(dp*log(1-w_i^(1/dp)))*
        z^(1/dp-1)*(1-z^(1/dp))^
        (log(0.5)/log(1-w_i^(1/dp))-1)
    )
  )
  
  return(ell)
}

#########    NULL MODEL LOG-LIK   ##############3
l_Kw_noreg <- function(par){ 
  w <- par[1]
  dp <- par[2]
  
  ell <- sum(
    log(
      log(.5)/(dp*log(1-w^(1/dp)))*
        z^(1/dp-1)*(1-z^(1/dp))^
        (log(0.5)/log(1-w^(1/dp))-1) 
    )
  )
  return(ell)
}




################################################################################
############  UW quantile regression model - USEFUL FUNCTIONS   ################
################################################################################
pdf_UW <- function(z,q,b){  # UW pdf
  b/z*(log(tau)/log(q))*(log(z)/log(q))^(b-1)*
    tau^((log(z)/log(q))^b)
}

cdf_UW <- function(z,q,b)   # UW cdf
{
  tau^((log(z)/log(q))^b)
}

r_UW <- function(n,q,b)    #It generates occurences of Z_i ~ UW (q_i, \beta)
{
  u = runif(n)
  z_UW = exp(log(q)*((log(u)/log(tau))^(1/b)))
  return(z_UW)
}

l_UWreg <- function(par){   # UW log-likelihood
  beta_vec <- par[-length(par)]
  b <- par[length(par)]

  q_i <- lfunc(beta_vec,X)
  ell =  sum(log(pdf_UW(z,q_i,b)))
  return(ell)
}



l_UWreg_aug <- function(par){    # UW log-likelihood (augmented model)
  beta_vec <- par[-length(par)]
  b <- par[length(par)]

  q_i <- lfunc(beta_vec,X_aug)
  ell =sum(log(pdf_UW(z,q_i,b)))
  return(ell)
}

#########    NULL MODEL LOG-LIK   ##############3
l_UW_noreg <- function(par){
  q <- par[1]
  b <- par[2]

  ell <- sum(log(pdf_UW(z,q,b)))
  return(ell)
}


################################################################################
############  Beta regression model - USEFUL FUNCTIONS   ################
################################################################################

pdf_beta <- function(z,q_i,c){  # beta pdf
  p = q_i*c
  q = (1-q_i)*c
  dbeta(z,p,q)
}

cdf_beta <- function(z,q_i,c){  # beta pdf
  p = q_i*c
  q = (1-q_i)*c
  pbeta(z,p,q)
}

r_beta <- function(n)   
{
  p = q_i*c
  q = (1-q_i)*c
  rbeta(n,p,q)
}

l_betareg <- function(par){   # beta log-likelihood
  beta_vec <- par[-length(par)]
  c <- par[length(par)]
  
  q_i <- lfunc(beta_vec,X)
  ell =  sum(log(pdf_beta(z,q_i,c)))
  return(ell)
}

l_betareg_aug <- function(par){    # beta log-likelihood (augmented model)
  beta_vec <- par[-length(par)]
  c <- par[length(par)]
  
  q_i <- lfunc(beta_vec,X_aug)
  ell = sum(log(pdf_beta(z,q_i,c)))
  return(ell)
}

#########    NULL MODEL LOG-LIK   ##############3
l_beta_noreg <- function(par){
  q <- par[1]
  c <- par[2]
  
  ell <- sum(log(pdf_beta(z,q,c)))
  return(ell)
}


