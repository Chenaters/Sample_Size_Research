####################
# Simulations for varifying the prediction mean square error (PMSE) calculation 
# by \cite{narula1974predictive} "Predictive mean square error and stochastic regressor variables"
####################


#generate data
k=5
p=2
n=100



simu <- function(times){
  output = list();
  estimation = list();
  for (i in 1:times) 
  {
    sigma_para <- sample(1:10,1)
    mu_para <- sample(1:10,1)
    alpha <- matrix(rnorm(10,mu_para,sigma_para),n,1)
    beta <- matrix(rnorm(10,mu_para,sigma_para),k,1)  
    mu_z <- matrix(sample(1:10,1),k,1)
    sigma_z <- matrix(sample(1:10,1),k,k)
    Z<-matrix( rnorm(n*k,mean=mu_z,sd=sigma_z), n, k) 
    zbar <- apply(Z,2,mean)
    X <- Z-zbar
    Xbar <- apply(X,2,mean)
    sigma_ks= sample(1:20,1)
    error <- matrix(rnorm(100, 0, sd=sigma_ks),n,1)  
    Y <- alpha + Z %*% beta  + error
    data <- data.frame(Y,X)
    
    #estimation
    s_00 <- sum((Y-mean(Y))^2)/(n-1)
    s <- t(X)%*%(Y-mean(Y))/(n-1)
    S <- t(X)%*%X/(n-1)
    betahat <- as.matrix(solve(S)%*%s)
    test <- data[sample(1:100,1),]
    y_0 <- test[,1]
    x_0 <- as.matrix(test[,-1])
    y_0hat <- mean(Y) + x_0%*%betahat
    PMSE <- (y_0-y_0hat)^2
    est <- sigma_ks*(1+1/n)*(n-2)/(n-k-2)
    
    output[i] = PMSE
    estimation[i]=est
  }
  diff <- mean(sapply(output, mean,2)-sapply(estimation,mean,2))
  diff
}

simu(100)


#############################
library(MASS);

###Parameter settings
n = 10; #Sample size
k = 5; #Total number of covariates in regression
p = 2; #Number of covariates in subset/partial model

sigmak = 1; #Standard deviation of the error term in regression 

#SIGMA=diag(1, k, k); #Covariates' variance matrix. Require it is positive definite. 
SIGMA = polyCor(5, 0.3);
MU=rep(0, k); #Covariates' mean vector 
covariateNames = paste("Z", 1:k, sep=""); #Covariates' names

alpha = 1; #Intercept in regression
BETA = array(1, k); #Coefficient vector in regression

simuN = 1e5; #Number of simulations


### Simulation loops
pmse = array(NA, simuN); #prediction squared errors
pme = array(NA, simuN); #prediction errors
for (i in 1:simuN) {
  ###Estimation data
  Z = mvrnorm(n=n, mu=MU, Sigma=SIGMA); #Matrix of covariates
  EPS = rnorm(n=n, mean=0, sd=sigmak); #Vector of errors in regression.
  Y = alpha + Z%*%BETA + EPS; #Vector of response values. 
  ###New data
  Z0 = mvrnorm(n=1, mu=MU, Sigma=SIGMA); #New covariate value
  Y0 = alpha + Z0%*%BETA + rnorm(n=1, mean=0, sd=sigmak); #New response value. 
  
  ###Calculations by the lm function (seems even faster than the math formula-based implementation)
  estiData = data.frame(Y, Z); #data for estimation
  names(estiData) = c("Y", covariateNames);
  
  lmformula = as.formula(paste("Y ~ ", paste(covariateNames, collapse= "+"))); #lm model formula
  lmfit = lm(lmformula, data=estiData); #LSE by lm
  #Y.hat = fitted(lmfit); #Fitted reponse values
  
  ###Prediction
  newData = data.frame(t(Z0));#New observation for prediction
  names(newData) = covariateNames;
  Y.hat.0 = predict(lmfit, newData); #predicted response.
  
  
  # ###Calculations by the math formulas provided in the paper,
  # ## the results matchs with above lm-based results.
  X = scale(Z, scale=F); #Centralized covariates
  S <- t(X)%*%X/(n-1); #
  s = t(X)%*%(Y-mean(Y))/(n-1);
  betahat = solve(S)%*%s; #Same as the coefficients given in lmfit.
  Yhat = mean(Y) + X%*%betahat; #Same as fitted(lmfit) above.

  Zbar = apply(Z,2,mean);
  X0 = Z0 - Zbar;
  Y.hat.0 = mean(Y) + X0%*%betahat; #predicted response. Same as above Y.hat.0
  
  # ###check calculaitons
  # ones = rep(1, n);
  # (t(ones)%*% (diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)) %*% ones)^(-1)*t(ones)%*%Y - (1/n)*t(ones)%*%Z %*% solve(t(Z)%*%(diag(n)-ones%*%t(ones)/n)%*%Z) %*%t(Z)%*%ones%*%t(ones)%*%Y/n
  # mean(Y); #same as above complex formula
  # 
  # n*(t(ones)%*% (diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)) %*% ones)^(-1) - (1/n)*t(ones)%*%Z %*% solve(t(Z)%*%(diag(n)-ones%*%t(ones)/n)%*%Z) %*%t(Z)%*%ones
  # 1; #equal to above
  # 
  # t(ones)%*% (diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)) %*% ones;
  # n - (1/n)*t(ones)%*%Z %*% solve(t(Z)%*%(diag(n)-ones%*%t(ones)/n)%*%Z) %*%t(Z)%*%ones %*% t(ones)%*% (diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)) %*% ones
  # 
  # 
  # Z %*% solve(t(Z)%*%(diag(n)-ones%*%t(ones)/n)%*%Z) %*%t(Z)%*%ones %*% t(ones)%*% (diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z))/n;
  # Z%*%solve(t(Z)%*%Z)%*%t(Z); #not equal, but the summations of all components are equal. Summation 
  # 
  # 
  # t(ones)%*% Z %*% solve(t(Z)%*%(diag(n)-ones%*%t(ones)/n)%*%Z) %*%t(Z)%*%ones %*% t(ones)%*% (diag(n) - Z%*%solve(t(Z)%*%Z)%*%t(Z)) %*%ones / n;
  # t(ones)%*%  Z%*%solve(t(Z)%*%Z)%*%t(Z) %*%ones; #equal
  
  
  ###prediction errors
  pmse[i] = (Y0-Y.hat.0)^2;
  pme[i] = Y0-Y.hat.0;
}

mean(pmse); 
sigmak^2*(1+1/n)*(n-2)/(n-k-2);
mean(pme);


#Distribution of PMSE: Exponential shape. 
hist(pmse); 
median(pmse); #median is much smaller than mean. 
sd(pmse);

#Distribution of PME: bellshape mean 0, SD close to 1. 
hist(pme); #exponential distribution
median(pme); #median is much smaller than mean. 
sd(pme);


####for subset
n = 100; #Sample size
k = 5; #Total number of covariates in regression
p = 2; #Number of covariates in subset/partial model

sigmak = 1; #Standard deviation of the error term in regression 
#sigmap^2 = sigmak^2 + t(sigma)Sigma^(-1)sigma-t(sigma_1)Sigma_11^(-1)sigma_1 

#SIGMA=diag(1, k, k); #Covariates' variance matrix. Require it is positive definite. 
SIGMA = polyCor(5, 0.3);
SIGMA_p = SIGMA[1:p, 1:p];
MU=rep(0, k); #Covariates' mean vector 
MU_p=rep(0, p); 
covariateNames = paste("Z", 1:k, sep=""); #Covariates' names
covariateNames_p = paste("Z", 1:p, sep=""); 

alpha = 1; #Intercept in regression
BETA = array(1, k);#Coefficient vector in regression


simuN = 1e5; #Number of simulations


### Simulation loops
pmse_p = array(NA, simuN); #prediction squared errors
pme_p = array(NA, simuN);#prediction errors

for (i in 1:simuN) {
  ###Estimation data
  Z = mvrnorm(n=n, mu=MU, Sigma=SIGMA);#Matrix of covariates
  Z_p = Z[,1:p]
  EPS = rnorm(n=n, mean=0, sd=sigmak); #Vector of errors in regression.
  Y = alpha + Z%*%BETA + EPS; #Vector of response values. 
  sigma_00 = t(BETA)%*%SIGMA%*%BETA+sigmak^2#variance of Y
  s = SIGMA%*%as.matrix(BETA) #true covariance of  Z and Y
  s_p = s[1:2] # true cov(Z_1,Y)
  SIGMA_11 = SIGMA[1:2,1:2]
  sigmap = sigma_00-t(s_p)%*%solve(SIGMA_11)%*%s_p #true sigmap

  ###New data
  Z0 = mvrnorm(n=1, mu=MU, Sigma=SIGMA); #New covariate value
  Z0_p = Z0[1:2]
  Y0 = alpha + Z0%*%BETA + rnorm(n=1, mean=0, sd=sigmak); #New response value. 
  
  ###Calculations by the lm function (seems even faster than the math formula-based implementation)
  estiData = data.frame(Y, Z); #data for estimation
  names(estiData) = c("Y", covariateNames);
  estiData_p = estiData[,1:3]
  
  lmformula_p = as.formula(paste("Y ~ ", paste(covariateNames_p, collapse= "+"))); #lm model formula
  lmfit_p = lm(lmformula_p, data=estiData_p);#LSE by lm
  #Y.hat = fitted(lmfit); #Fitted reponse values

  
  ###Prediction
  newData_p = data.frame(t(Z0_p));#New observation for prediction
  names(newData_p) = covariateNames_p;
  Y.hat.0_p = predict(lmfit_p, newData_p); #predicted response.
  
  
  # ###Calculations by the math formulas provided in the paper,
  # ## the results matchs with above lm-based results.
  # X = scale(Z, scale=F); #Centralized covariates
  # S <- t(X)%*%X/(n-1); #
  # s = t(X)%*%(Y-mean(Y))/(n-1);
  # sigmap[i] = sigmak^2+t(s)%*%solve(SIGMA)%*%s-t(s_p)%*%solve(SIGMA_11)%*%s_p
  # betahat = solve(S)%*%s; #Same as the coefficients given in lmfit.
  # Yhat = mean(Y) + X%*%betahat; #Same as fitted(lmfit) above.
  # 
  # Zbar = apply(Z,2,mean);
  # X0 = Z0 - Zbar;
  # Y.hat.0 = mean(Y) + X0%*%betahat; #predicted response. Same as above Y.hat.0
  
  ###prediction errors
  pmse_p[i] = (Y0-Y.hat.0_p)^2;
  pme_p[i] = Y0-Y.hat.0_p;

}

sigmap
sigmap*(1+1/n)*(n-2)/(n-p-2);




hist(pmse_p); 
median(pmse_p); 
sd(pmse_p);


hist(pme_p); 
median(pme_p); 
sd(pme_p);
##########====== Functions =========##########

########################################################################
# Description: Generate polynomially decaying correlation matrix with
#             the (i,j)th element = rho + (1-rho)/(1+|i-j|)^lambda, potentially truncated
#             by band width and/or value cutoff. 
# Depends on: 
# Arguments: 
#			 p: Dimension of the correlation matrix.
#			 rho: A basis correlation value given in the formula. Require value in [0, 1].
#			 lambda: A correlation decaying coefficient given in the formula. Require value > 0.
#      truncBandWidth: Band width from the dignal, outside which the correlations are 0.
#      truncCutoff: Cutoff value to truncate smaller correlations to 0.
# Details:  
# Value: A legitimate correlation matrix   
# Note: The default parameter values give a correlation matrix with homogenous off-diagnal elements being rho. 
# Author(s): HZhang and ZWu
# References: 
# See Also: 
# Example: polyCor(5, 0.3); #All off-diagnal elements are 0.3. 
########################################################################
polyCor <- function(p, rho=0, lambda=Inf, truncBandWidth=p, truncCutoff=0){
  firstRow <- rho + 1/(1+(0:(p-1)))^lambda*(1-rho); 
  
  if (truncBandWidth < p) { #truncation by band width
    firstRow <- firstRow * c(rep(1, truncBandWidth), rep(0, p-truncBandWidth)); 
  }
  
  if (truncCutoff > 0) {    #truncation by cutoff
    firstRow[firstRow < truncCutoff] <- 0; 
  }
  
  return(toeplitz(firstRow));
}




