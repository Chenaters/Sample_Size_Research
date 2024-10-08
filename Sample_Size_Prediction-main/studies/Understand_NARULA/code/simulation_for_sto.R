####################
# Simulations for varifying the prediction mean square error (PMSE) calculation 
# by \cite{narula1974predictive} "Predictive mean square error and stochastic regressor variables"
####################


#############################

PMSE = function(estiData, newData, p)
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

PME = function(estiData, newData,p)
{ 
  lmformula = as.formula(paste("Y ~ ", paste(covariateNames[1:p], collapse= "+"))); #lm model formula
  lmfit = lm(lmformula, data=estiData);#LSE by lm
  Y.hat.0 = predict(lmfit, newData); #predicted response.
  pme = Y0-Y.hat.0;
  pme
}


library(MASS);
source("Lib_CovarianceMatrix.R")
install.packages("combinat")                # Install combinat package
library("combinat")

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

k = 4 
n = 83 
alpha = 1
SIGMA =  polyCor(k, 0.3)
MU = rep(0,k)
BETA = rep(1,k)
sigmak2 = 1
simuN = 1e4
Z = mvrnorm(n=n, mu=MU, Sigma=SIGMA);#Matrix of covariates
EPS = rnorm(n=n, mean=0, sd=sigmak2); #Vector of errors in regression.
Y = alpha + Z%*%BETA + EPS; #Vector of response values. 
fullData = data.frame(Y, Z); #data for estimation
covariateNames = paste("Z", 1:k, sep=""); #Covariates' names
names(fullData) = c("Y", covariateNames);

comb = unlist(lapply(seq_along(covariateNames),combn,x=covariateNames,simplify=FALSE),recursive=FALSE)
all_formula = list(NA)
for ( i in 1:15) {
  all_formula[[i]] = as.formula(paste("Y ~ ", paste(comb[[i]], collapse= "+")))
}
all_formula

# 25 are selected at random and are used to estimate the parameters
estiData = fullData[sample(nrow(fullData), 25), ]
#the remaining 58 observations are predicted
predData = list(Y = fullData[-sample(nrow(fullData), 25), 1], 
                Z = fullData[-sample(nrow(fullData), 25), 2:5])


#fit model
lm.est = list(NA)
for (i in 1:15) {
  lm.est[[i]] = lm(all_formula[[i]],data = estiData)
}  ##we got 15 subset models

s_00 = t(Y-mean(Y))%*%(Y-mean(Y))/(n-1) #sample variance of Y
s = t(t(Y-mean(Y))%*%Z/(n-1))
S = t(Z)%*%Z/(n-1)  ####
var(Y);cov(Y,Z);var(Z)

#full model
PMSE_full = sum((predData$Y-predict.lm(lm.est[[15]],newdata = predData$Z))^2)

#choose the subset model with smallest PMSE for each observation

PMSE = list(NA)
for (n in 1:58) {
  for (i in 1:15) {
  pred = predict.lm(lm.est[[i]],newdata = predData$Z)
  PMSE[[i]] = sum((predData$Y-pred)^2)
  print(PMSE[[i]])
}
min(unlist(PMSE))
}


sigmak2 = sigma_error
sigmap2 = sigma_00-t(s[1:p,1])%*%solve(SIGMA_11)%*%s[1:p,1]; 
#variance of partial error term
cri = sigmap2(1+1/n)*(n-2)/(n-p-2)













