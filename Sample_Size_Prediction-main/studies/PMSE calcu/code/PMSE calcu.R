#calculation for "new" predictors based on 
#Baker TA, Buchanan NT, Corson N. Factorsinfluencing chronic pain intensity in older black women: examining depression, locus of control, and physical health.


library(MASS)
library(Matrix)
source("Lib_CovarianceMatrix.R")
set.seed(1)
#convert correlation matrix to covariance matrix
#The correlation matrix refers to the data in the paper
CORR = Matrix::forceSymmetric(as.matrix(read.table("SIGMA.txt")),uplo = "U")
SD = as.vector(rep(1,13)) #standard deviation vector
COV <-  sweep(sweep(CORR, 1, SD, "*"), 2, SD, "*")  


#parameter settings
k = 12 #total k predictors in full regression
p = 7 #p predictors in reduced regression
n = 181  # sample size
inflation = (n-p-2)/(n-k-2) #the inflation factor
alpha = 1 # intercept
SIGMA =  COV
MU = rep(0,k+1)  # mean of covariates
BETA_star = c(-0.2, -0.03, -0.02, -0.04, 0.12, 0.18, 0.26, 0.25, -0.01, 0.08, -0.26, 0.21)   #effect size
inflation = (n-p-2)/(n-k-2)   #inflation factor
simuN = 1e4
 
#calculation from covariance matrix
sigma_00 = COV[1,1]
sig = as.matrix(COV[-1,1])   #covariance vector between Y and Z_i
SIG = COV[-1,-1]  #covariate variance matrix
sig_1 = sig[1:p,1]; sig_2 = sig[(p+1):k,1] #partitioned covariance vector
SIG_11 = SIG[1:p,1:p]; SIG_12 = SIG[1:p,-(1:p)];
SIG_21 = SIG[-(1:p),1:p]; SIG_22 = SIG[-(1:p),-(1:p)];
BETA_true = SIG%*%sig #full model effects
BETA_true1 = BETA_true[1:p]; BETA_true2 = BETA_true[-(1:p)]; 
#not same as BETA_true?
BETA_calc1 = solve(SIG_11-SIG_12%*%solve(SIG_22)%*%SIG_21)%*%(sig_1-SIG_12%*%solve(SIG_22)%*%sig_2) 
BETA_calc2 = solve(SIG_22-SIG_21%*%solve(SIG_11)%*%SIG_12)%*%(sig_2-SIG_21%*%solve(SIG_11)%*%sig_1)
#error variance
sigmak2 = sigma_00 - t(sig)%*%solve(SIG)%*%sig  #variance of error term in full regression
sigmap2 = sigma_00 - t(sig_1)%*%solve(SIG_11)%*%sig_1  #variance of error term in full regression
EVR = sigmak2/sigmap2 #error variance ratio

BETA_redu1 = solve(SIG_11)%*%sig_1#reduced model effects
#pPMSEr
PMSE_calc = sigmak2*(n+1)*(n-2)/(n*(n-k-2)) #full regression with k predictors
PMSE1_calc = sigmap2*(n+1)*(n-2)/(n*(n-p-2)) #reduced regression with p predictors
pPMSEr_calc = (PMSE1_calc-PMSE_calc)/PMSE1_calc 
#efficient sample size
lambda_star = 1+0.1*(1/EVR-1)
n_calc = p+2+(k-p)*lambda_star/(lambda_star-1) #sample size in the paper is 181
#cohen's f^2
R2 = 0.31; R2_1 = 0.06; 
f2_calcu = R2-R2_1/(1-R2)

