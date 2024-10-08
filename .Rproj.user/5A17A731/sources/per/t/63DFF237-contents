####################
# Simulations for varifying the prediction mean square error (PMSE) calculation 
# by \cite{narula1974predictive} "Predictive mean square error and stochastic regressor variables"
####################


#############################

PMSE = function(estiData, newData, p=1)
{ 
  ###Calculations by the lm function (seems even faster than the math formula-based implementation)
  lmformula = as.formula(paste("Y ~ ", paste(covariateNames[1:p], collapse= "+"))); #lm model formula
  lmfit = lm(lmformula, data=estiData);#LSE by lm
  #Y.hat = fitted(lmfit); #Fitted reponse values
  ###Prediction
  Y.hat.0 = predict(lmfit, newData); #predicted response.
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
  pmse = (Y0-Y.hat.0)^2;
  pme = Y0-Y.hat.0;
  pmse
  #print(cbind(pmse,pme))
}

PME = function(estiData, newData,p=1)
{ 
  lmformula = as.formula(paste("Y ~ ", paste(covariateNames[1:p], collapse= "+"))); #lm model formula
  lmfit = lm(lmformula, data=estiData);#LSE by lm
  Y.hat.0 = predict(lmfit, newData); #predicted response.
  pme = Y0-Y.hat.0;
  pme
}


library(MASS);
source("Lib_CovarianceMatrix.R")
####unconditional


###Parameter settings
  #n : Sample size
  #k : Total number of covariates in regression
  #p : Number of covariates in subset/partial model
  #when p=k, the partial model expands into the whole model
  
  #sigmak2 : variance of the error term in full model regression 
  #sigmap^2 = sigmak^2 + t(sigma)Sigma^(-1)sigma-t(sigma_1)Sigma_11^(-1)sigma_1, variance of error term in partial model regression 
  
  #SIGMA=diag(1, k, k); #Covariates' variance matrix. Require it is positive definite. 
  #SIGMA = polyCor(k, 0.3), covriance of total covariantes;
  #alpha : Intercept in regression
  #BETA : Coefficient vector in regression
  #MU : Covariates' mean vector 

k = 10
p = 5
n = 20
alpha = 1
SIGMA =  polyCor(k, 0.3)
MU = rep(0,k)
BETA = rep(1,k)
sigmak2 = 1
simuN = 1e3

SIGMA_p = SIGMA[1:p,1:p];
MU_p=MU[1:p]; 
covariateNames = paste("Z", 1:k, sep=""); #Covariates' names
covariateNames_p = paste("Z", 1:p, sep=""); 



#true parameters calculation
sigma_00 = t(BETA)%*%SIGMA%*%BETA+sigmak2;#variance of Y
s = SIGMA%*%as.matrix(BETA);#true covariance of  Z and Y
SIGMA_11 = SIGMA[1:p,1:p];#partial variance of Z
sigmap2 = sigma_00-t(s[1:p,1])%*%solve(SIGMA_11)%*%s[1:p,1]; #variance of partial error term


### Simulation loops
  pmse = array(NA, simuN); #prediction squared errors
  pme = array(NA, simuN);#prediction errors 
for (i in 1:simuN) {
  ###Estimation data
  Z = mvrnorm(n=n, mu=MU, Sigma=SIGMA);#Matrix of covariates
  EPS = rnorm(n=n, mean=0, sd=sigmak2); #Vector of errors in regression.
  Y = alpha + Z%*%BETA + EPS; #Vector of response values. 
  fullData = data.frame(Y, Z); #data for estimation
  names(fullData) = c("Y", covariateNames);
  estiData = fullData[,0:p+1]
  ###New data
  Z0 = mvrnorm(n=1, mu=MU, Sigma=SIGMA); #New covariate value
  Y0 = alpha + Z0%*%BETA + rnorm(n=1, mean=0, sd=sigmak2); #New response value.
  newData = data.frame(t(Z0[1:p]));#New observation for prediction
  names(newData) = covariateNames_p;
  pmse[i] = PMSE(estiData, newData,p = 5)
  pme[i] = PME(estiData, newData,p = 5)
}
c(mean(pmse),sigmap2*(1+1/n)*(n-2)/(n-p-2),mean(pme))


#hist(pmse); c(median(pmse),sd(pmse));
#hist(pme); c(median(pme),sd(pme));




### COnditional PMSE
k = 10
p = 5
n = 20
alpha = 1
SIGMA =  polyCor(k, 0.3)
MU = rep(0,k)
BETA = rep(1,k)
sigmak2 = 1
simuN = 1e3
### same parameters setting as before
SIGMA_p = SIGMA[1:p,1:p];
MU_p=MU[1:p]; 
covariateNames = paste("Z", 1:k, sep=""); #Covariates' names
covariateNames_p = paste("Z", 1:p, sep=""); 


#true parameters calculation
sigma_00 = t(BETA)%*%SIGMA%*%BETA+sigmak2;#variance of Y
s = SIGMA%*%as.matrix(BETA);#true covariance of  Z and Y
SIGMA_11 = SIGMA[1:p,1:p];#partial variance of Z
sigmap2 = sigma_00-t(s[1:p,1])%*%solve(SIGMA_11)%*%s[1:p,1]; #variance of partial error term


### New data

Z0 = mvrnorm(n=1, mu=MU, Sigma=SIGMA); #New covariate value

if(k>p)
  {PHI_1 = (as.matrix(BETA[1:p]))+solve(SIGMA[1:p,1:p])%*%(SIGMA[1:p,(p+1):k])%*%(as.matrix(BETA[(p+1):k])); # expected value of partial covariate coefficient
pmse_sc = sigmak2 + sigmap2/n + sigmap2*(t(Z0[1:p]-MU[1:p])%*%solve(SIGMA[1:p,1:p])%*%(Z0[1:p]-MU[1:p])+p/n)/(n-p-2) + ((t(Z0-MU))%*%BETA - (t(Z0[1:p]-MU[1:p]))%*%PHI_1)^2}#expected partial conditional model pmse 
{pmse_fc=sigmak2*(1+1/n)+sigmak2*((Z0-MU)%*%solve(SIGMA)%*%(Z0-MU)+k/n)/(n-k-2)}
#expected full model conditional pmse as in the paper

# loop 
pmse_c = array(NA,simuN)
for (i in 1:simuN) {
  ###Estimation data
  Z = mvrnorm(n=n, mu=MU, Sigma=SIGMA);#Matrix of covariates
  EPS = rnorm(n=n, mean=0, sd=sigmak2); #Vector of errors in regression.
  Y = alpha + Z%*%BETA + EPS; #Vector of response values.  
  Y0 = alpha + Z0%*%BETA + rnorm(n=1, mean=0, sd=sigmak2); #New response value.
  fullData = data.frame(Y, Z); #data for estimation
  names(fullData) = c("Y", covariateNames);
  estiData = estiData[,0:p+1]  
  newData = data.frame(t(Z0[1:p]));#New observation for prediction
  names(newData) = covariateNames_p;
  pmse_c[i] = PMSE(estiData, newData)
} 

if(k>=p+1)
  {c(mean(pmse_c),pmse_sc)}
  {c(mean(pmse_c),pmse_fc)}










