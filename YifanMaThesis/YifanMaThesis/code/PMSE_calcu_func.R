#Related functions for prediction accuracy calculation with "new predictors"


# Description: Convert correlation matrix into covariance matrix for some data may have only correlation matrix.  
# Arguments: CORR: Correlation matrix
#            SD:standard deviation vector
# Output: Covariance matrix
Corr_to_Cov = function(CORR, SD)
  {
  return(sweep(sweep(CORR, 1, SD, "*"), 2, SD, "*"))
}


#Description: The inflation factor is related to estimation accuracy for error terms for full regression
#Arguments: n: sample size, 
#           p: number of original predictors in reduced model, 
#           k: number of predictors in full model, where k-p is the number of "new" predictors
#Output: The Inflation factor
#reference: Sawyer, R. (1982). Sample size and the accuracy of predictions made from multiple regression equations. Journal of Educational Statistics, 7(2), 91-104.
inflation = function(n,p,k){
  return((n-p-2)/(n-k-2))
}


#Variance of error term in full and reduced regression: sigmak2 & sigmap2
#for full regression
#Description: Calculate the error term variance based on covariance matrix and response variance for full model and reduced model
#Arguments: COV: covariance matrix of the response and all of the predictors
#           sig: covariance vector between response and predictors 
#           sigma_00: response variance
#Output: the error term variance sigmak2

sigmak2 = function(COV){
  sigma_00 = COV[1,1]       #response variance
  sig = as.matrix(COV[-1,1])   #covariance vector between Y and Z_i
  SIG = COV[-1,-1]  #covariate variance matrix
  return(sigma_00 - t(sig)%*%solve(SIG)%*%sig)
}

#for reduced regression
#Description: Calculate the error term variance based on covariance matrix and response variance for reduced model
#Arguments: COV: covariance matrix of the response and all of the predictors
#           p: number of original predictors in reduced model, 
#           k: number of predictors in full model, where k-p is the number of "new" predictors
#           
#           sig: covariance vector between response and predictors 
#           sigma_00: response variance 
#Output: the error term variance sigmap2 in reduced regression
sigmap2 = function(COV,p,k){
  sigma_00 = COV[1,1]
  sig = as.matrix(COV[-1,1])   #covariance vector between Y and Z_i
  SIG = COV[-1,-1]  #covariate variance matrix
  sig_1 = sig[1:p,1]; #partitioned covariance vector
  SIG_11 = SIG[1:p,1:p]; 
  return(sigma_00 - t(sig_1)%*%solve(SIG_11)%*%sig_1)
}


#pPMSEr
#Description: Evaluate how much better the performance could be after we introduce 
#             "new" predictors by Calculating reduced PMSE percentage between 
#             reduced regression and full regression (i.e. the regression with "new" predictors)
#Argument: n: sample size, 
#          p: number of original predictors in reduced model, 
#          k: number of predictors in full regression
#          sigmak2: error variance in full model
#          sigmap2: error variance in reduced model
#Output: reduced PMSE percentage, range fron 0 to 1
pPMSEr = function(sigmak2,sigmap2,n,k,p)
{
  PMSE = sigmak2*(n+1)*(n-2)/(n*(n-k-2)) #full regression with k predictors
  PMSE1 = sigmap2*(n+1)*(n-2)/(n*(n-p-2)) #reduced regression with p predictors
  return((PMSE1-PMSE)/PMSE1)
}


#effect size
#Description: calculate the effect size for predictors by LSE and partitioned covariance matrix
#Argument: COV: covariance matrix, including all predictors and response
#           p: number of original predictors in reduced model, 
#          k: number of predictors in full regression
#reference: Narula, S. C. (1974). Predictive mean square error and stochastic regressor variables. Journal of the Royal Statistical Society: Series C (Applied Statistics), 23 11–17.
#OutputL: calculated effect size
BETA_calcu = function(COV,p,k){
  sigma_00 = COV[1,1]
  sig = as.matrix(COV[-1,1])   #covariance vector between Y and Z_i
  SIG = COV[-1,-1]  #covariate variance matrix
  if (p == k)
    return(solve(SIG)%*%sig)
  else
  sig_1 = sig[1:p,1]; sig_2 = sig[-(1:p),1] #partitioned covariance vector
  SIG_11 = SIG[1:p,1:p]; SIG_12 = SIG[1:p,-(1:p)];
  SIG_21 = SIG[-(1:p),1:p]; SIG_22 = SIG[-(1:p),-(1:p)];#formulas refer to Narula
  BETA_calc1 = solve(SIG_11-SIG_12%*%solve(SIG_22)%*%SIG_21)%*%(sig_1-SIG_12%*%solve(SIG_22)%*%sig_2) 
  BETA_calc2 = solve(SIG_22-SIG_21%*%solve(SIG_11)%*%SIG_12)%*%(sig_2-SIG_21%*%solve(SIG_11)%*%sig_1)
return(rbind(BETA_calc1,BETA_calc2))
}

#efficient sample size
#Description: calculate how much sample size is large enough to make sure adding
#             new predictors do not worsen prediction accuracy.
#Argument: sigmak2: error variance in full model
#          sigmap2: error variance in reduced model
#          p: number of original predictors in reduced model, 
#          k: number of predictors in full regression
#Output: efficient sample size
eff_n = function(sigmak2,sigmap2,k,p){
  EVR = sigmak2/sigmap2
  lambda_star = 1+0.1*(1/EVR-1)
  return(p+2+(k-p)*lambda_star/(lambda_star-1))
}



