rm(list = objects())
source("UBXII_KW_UW_beta_functions_forSc1_simu.R")

######################## FIXED QUANTITIES  ##################
(R = 5000)    # number of replications
(n = 100)   # sample size
(tau = 0.5) # tau = 0.5 indicates modeling median

(beta1 = 1.3)                 # beta_1
(beta2 = 1.4)                 # beta_2
(c = 2.0)                     # shape parameter, precision 
(beta_true = c(beta1, beta2)) # regression parameters
(par_true = c(beta_true,c)) # true parameters vector

set.seed(2020)
x2 = rnorm(n)  
X <- cbind(1,x2)         # co-variables matrix
p <- dim(X)[2]+1         # number of parameters
q_i = lfunc(beta_true,X) # it call logit link function
############################################################

# Index and results matrix initialization
ind = indM = 1 

# To compute AIC and BIC proportions 
aicUBXII = aicKw = aicUW = aicbeta = 0
bicUBXII = bicKw = bicUW = bicbeta = 0
aickb = aickw = aicwb = aicuw = aicbeb = aicbeta = 0
bickb = bickw = bicwb = bicuw = bicbeb = bicbeta = 0

# MLEs and SE matrix
MLEs_SE_mat <- matrix(NA,16,6) 
rownames(MLEs_SE_mat) <- c("UBXII", "Kw", "UW", "beta",
                           "UBXII", "Kw", "UW", "beta",
                           "UBXII", "Kw", "UW", "beta",
                           "UBXII", "Kw", "UW", "beta")

colnames(MLEs_SE_mat) <-  c("m_est_b1","m_est_b2", "m_est_c", 
                            "SE_b1","SE_b2","SE_c")
# Measures comparison matrix
measuares_comp_mat <- matrix(NA,12,4)
colnames(measuares_comp_mat) <-c("UBXII", "Kw", "UW", "beta")
rownames(measuares_comp_mat) <- c("MSE", "AIC", "BIC",
                                  "MSE", "AIC", "BIC",
                                  "MSE", "AIC", "BIC",
                                  "MSE", "AIC", "BIC")


final_res_matrix <- matrix(NA,9, 10)
colnames(final_res_matrix) <- c("UBXII", "Kw", "UW", "beta",
                                "UBXII", "Kw",
                                "UBXII", "UW",
                                "UBXII", "beta")

# Initial time of the simulation  
time_start = Sys.time() 

for (s in 1:4) {
  
  # Response for initial guess
 
  if(s==1){
    set.seed(2020)
    yr <- r_UBXII(n)
  }
  if(s==2){
    set.seed(2020)
    yr <- r_Kw(n)
  }
  if(s==3){
    set.seed(2020)
    yr <- r_UW(n)
  }
  if(s==4){
    set.seed(2020)
    yr <- r_beta(n)
  }
  
  # To avoid 0's e 1's
  for (j in 1:length(yr)) {
    if(yr[j] == 1) {yr[j] <- 0.9999}
    if(yr[j] == 0) {yr[j] <- 0.0001}
  }
  
  # Response variable for MQO as initial guess 
  z =  log(yr/(1-yr))
  
  # Initial guess - Betas
  beta_guess = lm(z~x2)  #fitted linear regression (simple in this case)
  beta1_guess <- as.numeric(beta_guess$coefficients[1]) #guess for beta1
  beta2_guess <- as.numeric(beta_guess$coefficients[2]) #guess for beta2
  c_guess = 1.0  #guess for shape parameter
  guess = c(beta1_guess,beta2_guess,c_guess)
  
  # Monte Carlo (MC) Simulation
  bug = 0  #to count the bug numbers from MC simulation
  i = 1  # index MLEs matrix
  estim = estim_Kw = estim_UW = estim_beta = matrix(NA,R,3)  #matrix for to salve the MLEs
  
  # MSE, AIC, and BIC vectors
  mse_vec_UBXII = mse_vec_Kw = mse_vec_UW = mse_vec_beta =
    aic_vec_UBXII = aic_vec_Kw = aic_vec_UW = aic_vec_beta =
    bic_vec_UBXII = bic_vec_Kw = bic_vec_UW = bic_vec_beta <-vector()
 set.seed(2020)
  # Loop MC simulation
   while (i <= R) {  
    
     # Response variable's generation
     if(s==1){
       y <- r_UBXII(n)
     }
     if(s==2){
       y <- r_Kw(n)
     }
     if(s==3){
       y <- r_UW(n)
     }
     if(s==4){
       y <- r_beta(n)
     }
    
    # To avoid 0's e 1's
    for (j in 1:length(y)) {
      if(y[j] == 1) {y[j] <- 0.9999}
      if(y[j] == 0) {y[j] <- 0.0001}
    }
     
    # Estimation by ML of UBXII, Kw, UW, and beta regression
    res <- try(optim(guess,l_UBXIIreg,
                     method = "BFGS",
                     control = list(maxit = 100, reltol = 1e-38,
                                    fnscale = -1)),
               silent = T)
    
    res_Kw <- try(optim(guess,l_Kwreg,
                        method = "BFGS",
                        control = list(maxit = 100, reltol = 1e-38,
                                       fnscale = -1)),
                  silent = T)
    
    res_UW <- try(optim(guess,l_UWreg,
                        method = "BFGS",
                        control = list(maxit = 100, reltol = 1e-38,
                                       fnscale = -1)),
                  silent = T)
    
    res_beta <- try(optim(guess,l_betareg,
                          method = "BFGS",
                          control = list(maxit = 100, reltol = 1e-38,
                                         fnscale = -1)),
                    silent = T)
    if(class(res)=="try-error" || res$conv != 0 ||
       class(res_Kw)=="try-error" || res_Kw$conv != 0 ||
       class(res_UW)=="try-error" || res_UW$conv != 0||
       class(res_beta)=="try-error" || res_beta$conv != 0)
    {
      bug <- bug+1  # counting the bugs
    }else{
      estim[i,] <- res$par     
      estim_Kw[i,] <- res_Kw$par     
      estim_UW[i,] <- res_UW$par     
      estim_beta[i,] <- res_beta$par     
      
      # MSE computation
      y_hat_UBXII <- lfunc(c(res$par[1:(dim(X)[2])]),X)
      mse_vec_UBXII[i] <- 1/n*(sum(y-y_hat_UBXII)^2)
      
      y_hat_Kw <- lfunc(c(res_Kw$par[1:(dim(X)[2])]),X)
      mse_vec_Kw[i] <- 1/n*(sum(y-y_hat_Kw)^2)
      
      
      y_hat_UW <- lfunc(c(res_UW$par[1:(dim(X)[2])]),X)
      mse_vec_UW[i] <- 1/n*(sum(y-y_hat_UW)^2)
      
      y_hat_beta <- lfunc(c(res_beta$par[1:(dim(X)[2])]),X)
      mse_vec_beta[i] <- 1/n*(sum(y-y_hat_beta)^2)
      
      
      #AIC computation
      aic_vec_UBXII[i] <- 2*(p-res$value)
      aic_vec_Kw[i] <- 2*(p-res_Kw$value)
      aic_vec_UW[i] <- 2*(p-res_UW$value)
      aic_vec_beta[i] <- 2*(p-res_beta$value)
      
      #BIC computation
      bic_vec_UBXII[i] <- p*log(n)-2*res$value
      bic_vec_Kw[i] <- p*log(n)-2*res_Kw$value
      bic_vec_UW[i] <- p*log(n)-2*res_UW$value
      bic_vec_beta[i] <- p*log(n)-2*res_beta$value
      
      #Comparisons - prop AIC e BIC
      if(s==1){
        #AIC
      if(min(aic_vec_UBXII[i], 
             aic_vec_Kw[i],
             aic_vec_UW[i], 
             aic_vec_beta[i]) == aic_vec_UBXII[i])
      { aicUBXII <- aicUBXII +1}
        
        if(min(aic_vec_UBXII[i], 
               aic_vec_Kw[i],
               aic_vec_UW[i], 
               aic_vec_beta[i]) == aic_vec_Kw[i])
        { aicKw <- aicKw +1}
        
        if(min(aic_vec_UBXII[i], 
               aic_vec_Kw[i],
               aic_vec_UW[i], 
               aic_vec_beta[i]) == aic_vec_UW[i])
        { aicUW <- aicUW +1}
        
        if(min(aic_vec_UBXII[i], 
               aic_vec_Kw[i],
               aic_vec_UW[i], 
               aic_vec_beta[i]) == aic_vec_beta[i])
        { aicbeta <- aicbeta +1}
        
      #BIC
        if(min(bic_vec_UBXII[i], 
               bic_vec_Kw[i],
               bic_vec_UW[i], 
               bic_vec_beta[i]) == bic_vec_UBXII[i])
        { bicUBXII <- bicUBXII +1}
        
        if(min(bic_vec_UBXII[i], 
               bic_vec_Kw[i],
               bic_vec_UW[i], 
               bic_vec_beta[i]) == bic_vec_Kw[i])
        { bicKw <- bicKw +1}
        
        if(min(bic_vec_UBXII[i], 
               bic_vec_Kw[i],
               bic_vec_UW[i], 
               bic_vec_beta[i]) == bic_vec_UW[i])
        { bicUW <- bicUW +1}
        
        if(min(bic_vec_UBXII[i], 
               bic_vec_Kw[i],
               bic_vec_UW[i], 
               bic_vec_beta[i]) == bic_vec_beta[i])
        { bicbeta <- bicbeta +1}

      }
      
   
      if(s==2){
        if(aic_vec_UBXII[i] < aic_vec_Kw[i]){ aickb <- aickb +1}
        else{aickw = aickw+1}
        
        if(bic_vec_UBXII[i] < bic_vec_Kw[i]){ bickb <- bickb +1}
        else{bickw = bickw+1}
      }
      
      if(s==3){
        if(aic_vec_UBXII[i] < aic_vec_UW[i]){ aicwb <- aicwb +1}
        else{aicuw = aicuw+1}
        
        if(bic_vec_UBXII[i] < bic_vec_UW[i]){ bicwb <- bicwb +1}
        else{bicuw = bicuw+1}
      }
      
      if(s==4){
        if(aic_vec_UBXII[i] < aic_vec_beta[i]){ aicbeb <- aicbeb +1}
        else{aicbeta = aicbeta+1}
        
        if(bic_vec_UBXII[i] < bic_vec_beta[i]){ bicbeb <- bicbeb +1}
        else{bicbeta = bicbeta+1}
      }
      # print number of replication
     print(i)
      i <- i+1  #it is update the replication
      
    }
  }# end loop MC
  
  # Mean of the estimates
  mean_MLEs_UBXII <- apply(estim, 2, mean)#mean of the MLEs (2 --> by columns)
  mean_MLEs_Kw <- apply(estim_Kw, 2, mean)   
  mean_MLEs_UW <- apply(estim_UW, 2, mean) 
  mean_MLEs_beta <- apply(estim_beta, 2, mean) 
  
  # Standard errors
  SE_UBXII <-  sqrt(apply(estim[,],2,var))
  SE_Kw <-  sqrt(apply(estim_Kw[,],2,var))
  SE_UW <-  sqrt(apply(estim_UW[,],2,var))
  SE_beta <-  sqrt(apply(estim_beta[,],2,var))
  
  # Results matrix -  MLEs and SE
  MLEs_SE_mat[indM,] <- c(mean_MLEs_UBXII,SE_UBXII)
  MLEs_SE_mat[indM+1,] <- c(mean_MLEs_Kw,SE_Kw)
  MLEs_SE_mat[indM+2,] <- c(mean_MLEs_UW,SE_UW)
  MLEs_SE_mat[indM+3,] <- c(mean_MLEs_beta,SE_beta)
  indM = indM+4
  
  # Results matrix - Measures comparison matrix
  measuares_comp_mat[ind,] <- c(mean(mse_vec_UBXII),
                                mean(mse_vec_Kw),
                                mean(mse_vec_UW),
                                mean(mse_vec_beta))
  measuares_comp_mat[ind+1,] <- c(mean(aic_vec_UBXII),
                                  mean(aic_vec_Kw),
                                  mean(aic_vec_UW),
                                  mean(aic_vec_beta))
  
  measuares_comp_mat[ind+2,] <- c(mean(bic_vec_UBXII),
                                  mean(bic_vec_Kw),
                                  mean(bic_vec_UW),
                                  mean(bic_vec_beta))
  ind = ind+3
  
  # AIC and BIC (%)  - simul UBXII
  propAIC <- c(aicUBXII, aicKw, aicUW, aicbeta)/R*100
  propBIC <- c(bicUBXII, bicKw, bicUW, bicbeta)/R*100
  
  # AIC and BIC (%)  - simul Kw, UW, beta
  propAIC_Kw <- c(aickb, aickw)/R*100
  propAIC_UW <- c(aicwb, aicuw)/R*100
  propAIC_beta <- c(aicbeb, aicbeta)/R*100
  
  propBIC_Kw <- c(bickb, bickw)/R*100
  propBIC_UW <- c(bicwb, bicuw)/R*100
  propBIC_beta <- c(bicbeb, bicbeta)/R*100
  

  if(s==1){
  final_res_matrix[1,1:4] <- c(mean_MLEs_UBXII[1],
                            mean_MLEs_Kw[1],
                            mean_MLEs_UW[1],
                            mean_MLEs_beta[1])
  final_res_matrix[2,1:4] <- c(SE_UBXII[1],SE_Kw[1], SE_UW[1], SE_beta[1])
  
  final_res_matrix[3,1:4] <- c(mean_MLEs_UBXII[2],
                            mean_MLEs_Kw[2],
                            mean_MLEs_UW[2],
                            mean_MLEs_beta[2])
  final_res_matrix[4,1:4] <- c(SE_UBXII[2],SE_Kw[2], SE_UW[2], SE_beta[2])
  
  final_res_matrix[5,1:4] <- c(mean_MLEs_UBXII[3],
                            mean_MLEs_Kw[3],
                            mean_MLEs_UW[3],
                            mean_MLEs_beta[3])
  final_res_matrix[6,1:4] <- c(SE_UBXII[3],SE_Kw[3], SE_UW[3], SE_beta[3])
  
  final_res_matrix[7,1:4] <- c(propAIC)
  final_res_matrix[8,1:4] <- c(propBIC)
  final_res_matrix[9,1:4] <-  c(mean(mse_vec_UBXII),
                              mean(mse_vec_Kw),
                              mean(mse_vec_UW),
                              mean(mse_vec_beta) )
  }
  
  if(s==2){
    final_res_matrix[1:6,5:6] <- rbind(
                                c(mean_MLEs_UBXII[1], mean_MLEs_Kw[1]),
                                c(SE_UBXII[1],SE_Kw[1]),
                                c(mean_MLEs_UBXII[2], mean_MLEs_Kw[2]),
                                c(SE_UBXII[2],SE_Kw[2]),
                                c(mean_MLEs_UBXII[3], mean_MLEs_Kw[3]),
                                c(SE_UBXII[3],SE_Kw[3])
                                )

    final_res_matrix[7:9,5:6] <- rbind(propAIC_Kw,
                                     propBIC_Kw,
                                     c(mean(mse_vec_UBXII),
                                       mean(mse_vec_Kw))
                                     )
  }
  
  if(s==3){
    final_res_matrix[1:6,7:8] <- rbind(c(mean_MLEs_UBXII[1], mean_MLEs_UW[1]),
                                       c(SE_UBXII[1],SE_UW[1]),
                                       c(mean_MLEs_UBXII[2], mean_MLEs_UW[2]),
                                       c(SE_UBXII[2],SE_UW[2]),
                                       c(mean_MLEs_UBXII[3], mean_MLEs_UW[3]),
                                       c(SE_UBXII[3],SE_UW[3])
    )
      
    
    final_res_matrix[7:9,7:8] <- rbind(propAIC_UW,
                                       propBIC_UW,
                                       c(mean(mse_vec_UBXII),
                                         mean(mse_vec_UW))
    )
  }
  
  if(s==4){
    final_res_matrix[1:6,9:10] <- rbind(c(mean_MLEs_UBXII[1], mean_MLEs_beta[1]),
                                       c(SE_UBXII[1],SE_beta[1]),
                                       c(mean_MLEs_UBXII[2], mean_MLEs_beta[2]),
                                       c(SE_UBXII[2],SE_beta[2]),
                                       c(mean_MLEs_UBXII[3], mean_MLEs_beta[3]),
                                       c(SE_UBXII[3],SE_beta[3])
    )
    
    
    final_res_matrix[7:9,9:10] <- rbind(propAIC_beta,
                                       propBIC_beta,
                                       c(mean(mse_vec_UBXII),
                                         mean(mse_vec_beta))
    )
  }
  

  
  print(s)
}

# Print results
MLEs_SE_mat
measuares_comp_mat
final_res_matrix

# Results for tex
stargazer::stargazer(MLEs_SE_mat,digits=4)
stargazer::stargazer(measuares_comp_mat,digits=4)
stargazer::stargazer(final_res_matrix, digits=2)
# Save the results
ttt <-
  paste0("simu_",n,
         "_R",R,
         "_beta1_", beta1,
         "_beta2_", beta2,
         "_c_", c,
         ".RData")

save.image(ttt)
