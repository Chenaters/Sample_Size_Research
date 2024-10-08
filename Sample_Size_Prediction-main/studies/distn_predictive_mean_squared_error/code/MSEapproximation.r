# simulation to get the approximation result for prediction square error
library(MASS)
source("Lib_CovarianceMatrix.R")

set.seed(5)
p = 1 # n observations and p predictors
MU = rep(1,p)
SIGMA =  polyCor(p, 0.3)
alpha = 2
N = 100000  
sigma2 = 1
BETA = rnorm(p,0,1)
sigma_error = 1  #under assumption


for (n in c(10,25,50,75)){
  #unconditional simulation
  MSE  = array()  
  for (i in 1:10000){
    # for unconditional, data z_1,...,z_n not given
    Z = mvrnorm(n=n, mu=MU, Sigma=SIGMA);#Matrix of covariates
    Y = alpha + Z%*%BETA + rnorm(n,0,sigma2)
    data = data.frame(Y,Z)
    covariateNames = paste("Z", 1:p, sep=""); #Covariates' names
    names(data) = c("Y",covariateNames);
    
    #new data
    Z_0 = mvrnorm(n=1, mu=MU, Sigma=SIGMA);#Matrix of covariates
    Y_0 = alpha + as.vector(Z_0)%*%BETA + rnorm(1,0,sigma2)
    newData = data.frame(t(Z_0));
    names(newData) = covariateNames;
    lmformula <- as.formula(paste("Y ~ ", paste(covariateNames[1:p], collapse= "+")));
    model <- lm(lmformula , data = data)
    
    # predicting the target variable
    prediction <- cbind(1,t(Z_0))%*%as.matrix(model$coefficients)
    MSE[i] = (Y_0- prediction)^2
  }
  hist(MSE)  #empirical distribution of (y_hat-y*)^2
  F_hat_CDF = ecdf(MSE)  #empirical CDF from 10000 observation
  plot(F_hat_CDF)
  
  range = 25
  list = seq(0,range,0.01)   # for fitted value
  F_hat_fitted = F_hat_CDF(list)
  hist(F_hat_fitted)

  MSE_true = sigma_error^2*(n+1)*(n-2)/(n*(n-p-2))# true MSE under assumption
  Var_est = (sigma_error^4*(n+1)^2*(n-2))/(n^2*(n-p-2))*((3*n-12)/(n-p-4)-(n-2)/(n-p-2))
  MSE_ep = mean(MSE)  #empirical MSE
  Var_ep = var(MSE)
  sigma_prime = sqrt(MSE_true)
  diff_MSE = abs(MSE_ep-MSE_true)
  diff_Var_SE = abs(Var_ep-Var_est)
  diff_Var_SE

  #approximation with the 1st term
  F_1_CDF = function(t) {2*pnorm(sqrt(t)/sigma_prime)-1}
  F_1_fitted = F_1_CDF(list)
  hist(F_1_fitted)
    
  #difference between F_hat and F_1
  plot(F_1_fitted, x = list, xlim = c(0,range), col = "blue", type = "l");
  lines(F_hat_fitted, x = list, xlim = c(0,range), col = "red");
  
  diff_F1 = max(abs(F_hat_fitted-F_1_fitted))
  
  
  #approximation with 1st and 2nd term 
  F_2 = function(t){
    2*pnorm(sqrt(t)/sigma_prime)-1-p/(2*sigma_error^3*(n-2)*(n-p-4))*
                            ((sqrt(t)/sigma_prime)^3-3*(sqrt(t)/sigma_prime))*dnorm(sqrt(t)/sigma_prime)
  }
  F_2_fitted = F_2(list)#fitted value
  
  plot(F_hat_fitted, x = list, col = "red", type = "l")
  lines(F_2_fitted, x = list, col = 'purple', type = "l")
  
  diff_F2 = max(abs(F_hat_fitted-F_2_fitted))
  
  #approximation with gamma
  gamma_shape = ((n-p-4)*(n-2))/(3*(n-4)*(n-p-2)-(n-2)*(n-p-4))
  gamma_scale = sigma_error^2*((n+1)/n)*((3*n-12)/(n-p-4)-(n-2)/(n-p-2))
  #validation
  gamma_est_mean = gamma_shape*gamma_scale
  gamma_est_var = gamma_shape*(gamma_scale^2)
  c(abs(gamma_est_mean-MSE_true),abs(gamma_est_mean-MSE_ep))
  c(abs(gamma_est_var-Var_est), abs(gamma_est_var-Var_ep))
  F_3 = function(t){
    pgamma(t,shape = gamma_shape,
           scale = gamma_scale)
  }
  F_3_fitted = F_3(list)#fitted value
  plot(F_hat_fitted, x = list, col = "red", type = "l")
  lines(F_3_fitted, x = list, col = 'green')
  diff_F3 = max(abs(F_hat_fitted-F_3_fitted))
  

  
  #
  F_4 = ecdf(
  (sigma_prime^2)*rchisq(N, df = 1))
  F_4_fitted = F_5(list)
  plot(F_hat_fitted, x = list, col = "red", type = "l")
  lines(F_4_fitted, x = list, col = 'blue')
  diff_F4 = max(abs(F_hat_fitted-F_5_fitted))
  
  output = as.matrix(t(c(diff_F1, diff_F2, diff_F3, 
           diff_F4, diff_MSE, diff_Var_SE)))
  View(output)

}
