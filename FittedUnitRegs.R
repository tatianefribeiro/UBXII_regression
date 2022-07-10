###################################################################################
# PAPER: The Unit Burr XII regression: properties, simulation and application
# GOAL: Fiting the UBXII, KW, UW, BETA, and SIMPLEX regression models
# AUTHOR: Tatiane Fontana Ribeiro
# LAST UPDATE: October 31, 2020
###################################################################################
rm(list = objects())
library(tidyverse)
library(ggplot2)
set.seed(2020) 

##Directory
#Insert your directory
setwd("/media/tatiane/239fb286-e175-4578-b0d3-48c750946446/UBXIIreg2022/GitHub_UBXII_regression")
#####################              AUXILIARY FUNCTIONS          ############################### 
source("unit_regressions_fit.R") # It fits a unit reg model, provides LOOCV measure, and plots.
source("UBXII_Kw_UW_beta_functions_forApp.R")    #link function and others.

#COVID-19 data set
data.set <- read.csv("evasao_zootecnia2009.CSV")  #complete data set
head(data.set,1)

#View(data.set)

#*****************************************************************************************
#                             UBXII Regression Model
#*****************************************************************************************
#Selected covariates
#The rv + the best covariates combination found
#####################################
#new covariate:  QT_TEC_FUND_COMP_tot
#####################################
data.set.select <- data.set%>%filter(TX_EVASAO!=1)%>%
  dplyr::select(NO_IES, TX_EVASAO,
                QT_VAGAS_MATUTINO,
                QT_TEC_FUND_COMP_MASC,
                QT_TEC_FUND_COMP_FEM,
                IN_OUTRAS_FORMAS_INGRESSO,
                IN_TRADUTOR_LIBRAS,
                IN_SINTESE_VOZ   ,
                IN_INTEGRAL_CURSO,
                QT_INGRESSO_CURSO,
                IN_AJUDA_DEFICIENTE,
                IN_NOTURNO_CURSO
  )%>%mutate(QT_TEC_FUND_COMP_MASC = QT_TEC_FUND_COMP_MASC+QT_TEC_FUND_COMP_FEM)%>%
  rename(QT_TEC_FUND_COMP_tot = QT_TEC_FUND_COMP_MASC,
         IN_ACCESSIBILITY = IN_AJUDA_DEFICIENTE,
         IN_NIGHT_COURSE = IN_NOTURNO_CURSO)
######################################

data.set.select <- data.set.select%>%
  dplyr::select(NO_IES, TX_EVASAO,
                QT_VAGAS_MATUTINO,
                QT_TEC_FUND_COMP_tot,
                IN_OUTRAS_FORMAS_INGRESSO,
                IN_SINTESE_VOZ   ,
                IN_INTEGRAL_CURSO,
                IN_ACCESSIBILITY,
                IN_NIGHT_COURSE
              #  QT_INGRESSO_CURSO
                )

data.set.select <- data.set.select%>%
  dplyr::select( TX_EVASAO,
                QT_VAGAS_MATUTINO,
                IN_ACCESSIBILITY,
                IN_NIGHT_COURSE)

#Response 
z <- data.set.select$TX_EVASAO

# sample size
n <- length(z)
# Covariates matrix
X <- model.matrix(TX_EVASAO~.,data = data.set.select)  

# FITTED UBXII REGRESSION MODEL
tau = .5
UBXII_reg = UnitReg.fit(z,X,regression = "UBXII")
round(UBXII_reg$summary,4) 
UBXII_reg$diagnosticMEASURES
(res=round(LOOCV.unit_reg(z,X, regression = "UBXII"),4))


# The parameter estimates and the 95% pointwise conﬁdence intervals
# for the UBXII model and τ = 0.1, 0.2, . . . , 0.8 and 0.9
alpha = 0.05
tau_vec <- seq(0.1,0.9,0.1)
i = 1
L_beta1 = E_beta1= U_beta1 = vector()
L_beta2 = E_beta2= U_beta2 = vector()
L_beta3 = E_beta3= U_beta3 = vector()
L_beta4 = E_beta4= U_beta4 = vector()
L_c = E_c = U_c = vector()
for (tau in tau_vec) {
  
  UBXII_reg = UnitReg.fit(z,X,regression = "UBXII")
  
  L_beta1[i] <- UBXII_reg$summary[1,1]-qnorm(1-alpha/2)*UBXII_reg$summary[1,2]
  E_beta1[i] <- UBXII_reg$summary[1,1] 
  U_beta1[i] <- UBXII_reg$summary[1,1]+qnorm(1-alpha/2)*UBXII_reg$summary[1,2]
  
  L_beta2[i] <- UBXII_reg$summary[2,1]-qnorm(1-alpha/2)*UBXII_reg$summary[2,2]
  E_beta2[i] <- UBXII_reg$summary[2,1] 
  U_beta2[i] <- UBXII_reg$summary[2,1]+qnorm(1-alpha/2)*UBXII_reg$summary[2,2]
  
  L_beta3[i] <- UBXII_reg$summary[3,1]-qnorm(1-alpha/2)*UBXII_reg$summary[3,2]
  E_beta3[i] <- UBXII_reg$summary[3,1] 
  U_beta3[i] <- UBXII_reg$summary[3,1]+qnorm(1-alpha/2)*UBXII_reg$summary[3,2]
  
  L_beta4[i] <- UBXII_reg$summary[4,1]-qnorm(1-alpha/2)*UBXII_reg$summary[4,2]
  E_beta4[i] <- UBXII_reg$summary[4,1] 
  U_beta4[i] <- UBXII_reg$summary[4,1]+qnorm(1-alpha/2)*UBXII_reg$summary[4,2]
  
  L_c[i] <- UBXII_reg$summary[5,1]-qnorm(1-alpha/2)*UBXII_reg$summary[5,2]
  E_c[i] <- UBXII_reg$summary[5,1] 
  U_c[i] <- UBXII_reg$summary[5,1]+qnorm(1-alpha/2)*UBXII_reg$summary[5,2]
  
  
  i = i+1

}

# Plots sensibility - quantiles
df <- data.frame(x = tau_vec,
                 L = L_beta1,
                 E = E_beta1,
                 U =U_beta1)


plot_beta1 <- ggplot(df, aes(x = x, y = E)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))+theme_bw()+
  labs(x = "Quantile", y = expression(hat(beta[1])))+ 
  scale_x_continuous(breaks = tau_vec)

df <- data.frame(x = tau_vec,
                 L = L_beta2,
                 E = E_beta2,
                 U =U_beta2)
plot_beta2 <- ggplot(df, aes(x = x, y = E)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))+theme_bw()+
  labs(x = "Quantile", y = expression(hat(beta[2])))+ 
  scale_x_continuous(breaks = tau_vec)


df <- data.frame(x = tau_vec,
                 L = L_beta3,
                 E = E_beta3,
                 U =U_beta3)
plot_beta3 <- ggplot(df, aes(x = x, y = E)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))+theme_bw()+
  labs(x = "Quantile", y = expression(hat(beta[3])))+ 
  scale_x_continuous(breaks = tau_vec)


df <- data.frame(x = tau_vec,
                 L = L_beta4,
                 E = E_beta4,
                 U =U_beta4)
plot_beta4 <- ggplot(df, aes(x = x, y = E)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))+theme_bw()+
  labs(x = "Quantile", y = expression(hat(beta[4])))+ 
  scale_x_continuous(breaks = tau_vec)


df <- data.frame(x = tau_vec,
                 L = L_c,
                 E = E_c,
                 U =U_c)
plot_c <- ggplot(df, aes(x = x, y = E)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = U, ymin = L))+theme_bw()+
  labs(x = "Quantile", y = expression(hat(c)))+ 
  scale_x_continuous(breaks = tau_vec)


library(cowplot)
plot_grid(plot_beta1,plot_beta2, plot_beta3,plot_beta4,plot_c)

# Plots
w1<-9
h2<-7
postscript(file = "sensib_tau.eps",horizontal=F,paper="special",width = w1, height = h2,family = "Times")
{
  plot_grid(plot_beta1,plot_beta2, plot_beta3,plot_beta4,plot_c)
}
dev.off()

################################
# The generalized Cook distance 
################################
tau = 0.5
UBXII_reg = UnitReg.fit(z,X,regression = "UBXII")
theta_hat <- as.numeric(UBXII_reg$summary[,1])
J_hat <- UBXII_reg$J

UBXII_reg_i = UnitReg.fit(z[-1],X[-1,],regression = "UBXII")
theta_hat_i <- as.numeric(UBXII_reg_i$summary[,1])

t(theta_hat_i-theta_hat)%*%J_hat%*%(theta_hat_i-theta_hat)

GD <- vector()
for (i in 1:n) {
  UBXII_reg_i <- UnitReg.fit(z[-i],X[-i,], regression = "UBXII")
  theta_hat_i <- as.numeric(UBXII_reg_i$summary[,1])
  
  GD[i] <- t(theta_hat_i-theta_hat)%*%
          J_hat%*%
          (theta_hat_i-theta_hat)
}

#Plot GD
w1<-8
h2<-5
postscript(file = "gd_plot.eps",horizontal=F,paper="special",width = w1, height = h2,family = "Times")
{
  index <- 1:n
  df <- data.frame(x = index,
                   GD = GD)
  ggplot(df, aes(x=1:n, y=GD)) +    # make a plot
    geom_point() +
    geom_segment(aes(index,xend=index,0, yend=GD), data=df)+
    labs(x='Observation Index', 
         y='Generalized Cook distance')+
    theme_bw()+scale_x_continuous(breaks = seq(1,n,3))
}
dev.off()

#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "UBXII")
q_hat2 <- (unit_reg$q.fv)^2
q_hat3 <- (unit_reg$q.fv)^3
data_aug_mod <- data.frame(z,X[,-1],q_hat2,q_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "UBXII")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

# Plots
plots_quantres(UBXII_reg$residuals,UBXII_reg$q.fv,n)

#*****************************************************************************************
#                             Kw Regression Model
#*****************************************************************************************
#Response 
z <- data.set.select$TX_EVASAO
# Covariates matrix
X <- model.matrix(TX_EVASAO~.,data = data.set.select)  

Kw_reg = UnitReg.fit(z,X,regression = "Kw")
round(Kw_reg$summary,4)                       
Kw_reg$diagnosticMEASURES
(res_Kw=round(LOOCV.unit_reg(z,X, regression = "Kw"),4))


#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "Kw")
omega_hat2 <- (unit_reg$q.fv)^2
omega_hat3 <- (unit_reg$q.fv)^3
data_aug_mod <- data.frame(z,X[,-1],omega_hat2,omega_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "Kw")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

# Plots
plots_quantres(Kw_reg$residuals,Kw_reg$mu.fv,n)

#*****************************************************************************************
#                             UW Regression Model
#*****************************************************************************************
#*
#Response 
z <- data.set.select$TX_EVASAO
# Covariates matrix
X <- model.matrix(TX_EVASAO~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
UW_reg = UnitReg.fit(z,X,regression = "UW")
round(UW_reg$summary,4)                       
UW_reg$diagnosticMEASURES
(res_UW=round(LOOCV.unit_reg(z,X, regression = "UW"),4))


#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "UW")
q_hat2 <- (unit_reg$q.fv)^2
q_hat3 <- (unit_reg$q.fv)^3
data_aug_mod <- data.frame(z,X[,-1],q_hat2,q_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "UW")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

# Plots
plots_quantres(UW_reg$residuals,UW_reg$mu.fv,n)



#*****************************************************************************************
#                             BETA Regression Model
#*****************************************************************************************
#Response 
z <- data.set.select$TX_EVASAO
# Covariates matrix
X <- model.matrix(TX_EVASAO~.,data = data.set.select)  

# FITTED RUBXII REGRESSION MODEL
BETA_reg = UnitReg.fit(z,X,regression = "beta")
round(BETA_reg$summary,4)                          #!!!!!!!! NA, não sei como resolver.
round(BETA_reg$diagnosticMEASURES,4)
(res_BETA=round(LOOCV.unit_reg(z,X, regression = "beta"),4))

#RESET-type test
unit_reg <- UnitReg.fit(z,X,regression = "beta")
mu_hat2 <- (unit_reg$q.fv)^2
mu_hat3 <- (unit_reg$q.fv)^3
data_aug_mod <- data.frame(z,X[,-1],mu_hat2,mu_hat3)
X_aug <- model.matrix(z~.,data = data_aug_mod)  
X <- X_aug
mod_aug <- UnitReg.fit(z,X,regression = "beta")
ell_rest <- unit_reg$max_ell
ell_unr <- mod_aug$max_ell
(w_LR <- 2*(ell_unr-ell_rest))
(p_value_LR <- 1-pchisq(w_LR,2)) 

# Plots
plots_quantres(BETA_reg$residuals,BETA_reg$mu.fv,n)

