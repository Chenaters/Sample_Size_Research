####################
# Funcitons for calculating the prediction mean square error (PMSE) and 
# related measures
####################


#############################




library(MASS);
####unconditional




########################################################################
# Description: Calculate the percentage PMSE reduction (pPMSEr) based on regression model
#    and related influential factors. 
# Depends on: 
# Arguments: 
#	    n: Sample size. Require: n > p1+p2+ (LSE is meaningful s.t. the calcu is correct). 
#     p1: Number of basic predictors in basic/reduced model
#     p2: Number of new predictors in the extended/full model
#     method: calculation methods: 
#             "f2": by Cohen’s f^2 of the new predictors. Require fnew2.
#             "beta": by effect sizes according to marginal coefficients. 
#                     Require: Beta.marg, SIGMA, and R2 or sigmak2 (ignore sigmak2 if R2 given). 
#     fnew2: Cohen’s f^2 for new predictors's effects conditional on the known predictors.
#     Beta.marg: Vector of marginal coefficients/effects of all predictors (i.e., the  
#                 marginal change of response per unit increase of each predictor). 
#                 1-p1 entries: the basic predictors; p1+1, ..., p1+p2 entries: new predictors
#     SIGMA: Covariance matrix among predictors. 
#            Entry ordering is consistent with Beta.marg.
#     R2: Proportion of response’s variance accounted for by all predictors
#     sigmak2: Error variance after accounting for all predictors.
# Details: The basic formula is pPMSEr= (1 - sigma_k^2/sigma_p^2)*inflation, where sigma_k^2 and 
#         sigma_p^2 are the error variances of the extended and and basic models, respectively. 
#         inflation = (n-p1-2)/(n-k-2). 
# Value:    
# Note:  
# Author(s): ZWu
# References: 
# See Also: 
# Example: 
#     pPMSEr(n=30, p1=2, p2=5, method="f2", fnew2=0.1); #Negative value indicates 
#     pPMSEr(n=100, p1=2, p2=5, method="f2", fnew2=0.1);
########################################################################
pPMSEr <- function(n, p1, p2, method="f2", fnew2, Beta.marg, SIGMA, R2=NA, sigmak2){
  k = p1+p2; #Total number of predictors
  inflation = (n-p1-2)/(n-k-2); #The inflation factor
  
  
  if (method=="f2") return(pMSEr = 1 - inflation/(fnew2+1));
  if(method=="beta"){
    Sigma = diag(SIGMA)*Beta.marg; #Vector of marginal covariance between response and each predictors.
    
    #Beta.extd = solve(SIGMA)%*%Sigma; #Joint coeff of the extended model 
    #Beta.bas = solve(SIGMA)%*%Sigma; #Joint coeff of the basic model 
    
    var.extd = t(Sigma)%*%solve(SIGMA)%*%Sigma; #variance explained by the extended model
    var.bas = t(Sigma[1:p1])%*%solve(SIGMA[1:p1, 1:p1])%*%Sigma[1:p1]; #variance explained by the basic model
    
    sigma00 = ifelse(R2!=NA, var.extd/R2, var.extd+sigmak2);
    
    sigmak2 = sigma00 - var.extd; #Error variance of the extended model
    sigmap2 = sigma00 - var.bas; #Error variance of the basic model
    
    return(pMSEr = 1 - (sigmak2/sigmap2)*inflation);
  }
}





########################################################################
# Description: Calculate the “efficient sample size": The smallest n that makes the 
#     pPMSEr to exceed a percentage of the pPMSEr at n = Infy. 
# Depends on: 
# Arguments: 
#     p1: Number of basic predictors in basic/reduced model
#     p2: Number of new predictors in the extended/full model
#     evr: The error variance ratio: (error var of full model)/(error var of reduced model)
#     effi: Efficiency: pPMSEr(nstar)/pPMSEr(Infty) >= effi. Require: value in (0, 1). Default is 90%. 
# Details: 
# Value:  nstar: the “efficient sample size". 
#         lambda: the inflation at nstar. 
#         pPMSEr.nstar: the pMSEr at nstar.
#         pPMSEr.infty: the pMSEr at infty sample size (i.e., inflation=1)
# Note:  
# Author(s): ZWu
# References: 
# See Also: 
# Example: effi.n(p1=2, p2=5, evr=0.8); 
#          effi.n(p1=2, p2=5, evr=0.8, effi=.95); 
#          f2=0.3; #Cohen's f-square for the extra predictors in full model over reduced model. 
#          effi.n(p1=2, p2=5, evr=1/(f2+1));
########################################################################
effi.n <- function(p1, p2, evr, effi=0.9){
  alpha = 1-effi; 
  lambda = 1 + alpha*(1/evr -1);
  nstar = ceiling(p1 + 2 + p2*lambda/(lambda-1));
  pPMSEr.nstar = 1-evr*lambda;
  pPMSEr.infty = 1-evr;
  return(c(nstar=nstar, lambda=lambda, pPMSEr.nstar=pPMSEr.nstar, pPMSEr.infty=pPMSEr.infty)); 
}


