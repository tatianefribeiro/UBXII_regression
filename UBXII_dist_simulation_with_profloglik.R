rm(list = ls())
set.seed(2020)
library(stargazer)
#Scenario 1    ####  CHANGE HERE!    ####################
c = 1.1
q = 0.4       #+LINE 86 (LAST LINE), i.e., is, FILE NAME .Rdata
              #Se for simular novamente mude a saída, pq essa TAB é a completa.
###################  CHANGE HERE!    ####################
tau = 0.5
guess = 1.0
param<-c(c,q)
vn <- c(25, 75, 150, 300)
R <- 10000

qf_UBXII <- function(u)
{
  exp(-(u^(-1/(log(1/tau)/log(1+(log(1/q))^c)))-1)^(1/c))
}

l_UBXII_prof <- function(theta){ 
  c = theta[1]
  
  loglik_UBXII = -n+n*log(-c*log(tau))-sum(log(y))-sum(log(1+(log(1/y))^c))+
    (c-1)*sum(log(log(1/y)))-
    n*log(-1/n*log(tau)*sum(log(1+(log(1/y))^c)))
  
  return(loglik_UBXII)
}

grr_analytic <- function(theta){
  c = theta[1]
  
  n/c+sum(log(log(1/y)))-sum(((log(1/y))^c*log(log(1/y)))/(1+(log(1/y))^c))-
    (n*sum((log(1/y))^c*log(log(1/y))/(1+(log(1/y))^c))/
       sum(log(1+(log(1/y))^c)))}

k = 1
results_fin = matrix(NA,length(vn),11)
################ ARRAY BOXPLOT #####################
arr_boxplot = array(NA, c(R,2,length(vn)))
colnames(arr_boxplot) <- c("hat_c","hat_q")
####################################################
for (n in vn) {
  time_start = Sys.time()
  i<-bug<-0
  estim<-matrix(NA,nrow = R,ncol = 2)
  
  
  while (i<R) {
    u <- runif(n)
    y <- qf_UBXII(u)
    res <- try(optim(guess,l_UBXII_prof, grr_analytic,
                     method = "BFGS",
                     control = list(maxit=100,reltol=1e-38, fnscale=-1)),
               silent = T)  
    if(class(res)=="try-error" || res$conv != 0) # a classe dos objetos que contem o erro, 
    {
      bug<-bug+1
    }else{
      i <- i+1
      mle_q <- exp(-(exp((1/n)*log(1/tau)*
                           sum(log(1+(log(1/y))^res$par)))-1)^(1/res$par))
      
      estim[i,]<-c(res$par,mle_q)
      
      #Salvando TODAS as estimativas e não só as do último n como na simu anterior:
      arr_boxplot[i,,k] <- c(res$par,mle_q)
    }
  }
  
  mean <- apply(estim[,],2,mean)
  rb <- (mean-param)/param*100
  mse <- apply(estim[,],2,var)+(mean-param)^2  #não usei o vr como estava sendo usado
  
  #Adicionei
  rmse <- sqrt(mse)
  se <-apply(estim[,],2,sd)
  
  colnames(results_fin) <- c("c","q","n","m_c","m_q","rb_c","rb_q","se_c","se_q","rmse_c","rmse_q")
  results_fin[k,] = c(param,n,mean,rb,se,rmse)
  
  print(k)
  k = k+1
  
  time_fin = Sys.time()
}
time = time_fin-time_start
results_fin
stargazer::stargazer(results_fin,digits = 4)
save.image("/home/tatiane/Desktop/master_thesis20200410/art_capts_2020/scripts/UBXII_simulation/UBXII_scen4.RData")
time
