#calculation for "new" predictors based on 
#Baker TA, Buchanan NT, Corson N. Factorsinfluencing chronic pain intensity in older black women: examining depression, locus of control, and physical health.


library(MASS)
library(Matrix)
source("Lib_CovarianceMatrix.R")
set.seed(1)

#The correlation matrix refers to the data in the paper
CORR = read.table("SIGMA.txt")

#parameter settings
k = 12 #total k predictors in full regression
p = 7 #p predictors in reduced regression
n = 500  # sample size
alpha = 1 # intercept
CORR = as.matrix(read.table("SIGMA.txt"))#The correlation matrix refers to the data in the paper
SIGMA =  polyCor(k, 0.3)
MU = rep(0,k)  # mean of covariates
BETA_star = c(-0.2, -0.03, -0.02, -0.04, 0.12, 0.18, 0.26, 0.25, -0.01, 0.08, -0.26, 0.21)   #effect size
sigmak = 1  #standard deviation of error term in full regression
sigmak2 = sigmak^2  #variance of error term in full regression
lambda = (n-p-2)/(n-k-2)   #inflation factor
simuN = 1e4


sigmak2_hat = array(NA)
sigmap2_hat = array(NA)
pmse_res = array(NA)
pmse_1_res = array(NA)
pPMSEr = array(NA)
n_star = array(NA)
coh_f2 = array(NA)

for (i in 1:simuN) 
  {
#generate data
Z = mvrnorm(n=n, mu=MU, Sigma=SIGMA); #Matrix of total covariates
Y = alpha + Z%*%BETA_star + rnorm(n,0,sigmak)
data = data.frame(Y,Z)
covariateNames = paste("Z", 1:k, sep=""); #Covariates' names
names(data) = c("Y",covariateNames);
# patrition to "basic" predictors and "new" predictors
Z_1 = Z[,1:p]; Z_2 = Z[,-(1:p)] 

#Calculation
mu_0 = mean(Y)  #mu* = (mu_0, MU')'
MU_hat = mean(Z)
sigma_00 = var(Y)  
s = t(cov(Y,Z))
S = var(Z)    #estimation of covariance matrix S* = (sigma_00, s'//s, S)
s_1 = s[1:p,1];s_2 = s[-(1:p),1];  #partition the estimation of covariance vector s = (s_1',s_2')'
S_11 = S[1:p,1:p]; S_12 = S[1:p,-(1:p)]; 
S_21 = S[-(1:p),1:p]; S_22 = S[-(1:p),-(1:p)]

sigmak2_hat[i] = sigma_00 - t(s)%*%solve(S)%*%s
sigmap2_hat[i] = sigma_00-t(s_1)%*%solve(S_11)%*%s_1;


# generate new data
Z0 = mvrnorm(n=1, mu=MU, Sigma=SIGMA); 
Y0 = alpha + Z0%*%BETA_star + rnorm(n=1, mean=0, sd=sigmak); #new true response
newData = data.frame(Y0,t(Z0)); 
names(newData) = c("Y",covariateNames);


#full regression
lmformula = as.formula(paste("Y ~ ", paste(covariateNames, collapse= "+"))); 
lmfit = lm(lmformula, data=data);
Y.hat.0 = predict(lmfit, newData); 
pse_full = (Y0-Y.hat.0)^2;
pmse_full_est = sigmak2_hat[i]*(n+1)*(n-2)/(n*(n-k-2))
pmse_res[i] = pse_full-pmse_full_est


# reduced regression
data_1 = data[,1:(p+1)]
newData_1 = newData[1:(p+1)]
lmformula_1 = as.formula(paste("Y ~ ", paste(covariateNames[1:p], collapse= "+"))); 
lmfit_1 = lm(lmformula_1, data=data_1);
Y.hat.0_1 = predict(lmfit_1, newData_1); 
pse_red = (Y0-Y.hat.0_1)^2;
pmse_red_est = sigmap2_hat[i]*(n+1)*(n-2)/(n*(n-p-2))
pmse_1_res[i] = pse_red-pmse_red_est  #difference between real pmse and estimated pmse

#percentage of PMSE reduction
pPMSEr[i] = (1-sigmak2_hat[i]*(n-p-2)/(sigmap2_hat[i]*(n-k-2)))
#error variance ratio
EVR = sigmak2_hat[i]/sigmap2_hat[i]
#Efficient sample size
eff = 0.1 #alpha the efficiency
lambda_star = 1+ eff*(1/EVR - 1)
n_star[i] = p+2+(k-p)*lambda_star/(lambda_star-1)


#Effect sizes
coh_f2[i] = (1-EVR)/EVR  #sigmak2/sigmap2 = 1/(coh_f2+1)

}

output = c(n,mean(sigmak2_hat), mean(sigmap2_hat) ,mean(pmse_res), mean(pmse_1_res),mean(pPMSEr), mean(n_star),mean(coh_f2))
output

  output_mat[10,] = output
  
  output_mat[,2]/output_mat[,3]
  
  