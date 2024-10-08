#calculation for "new" predictors based on 
#Baker TA, Buchanan NT, Corson N. Factorsinfluencing chronic pain intensity in older black women: examining depression, locus of control, and physical health.


library(MASS)
library(Matrix)
source("PMSE_calcu_func.R")

set.seed(1)
#convert correlation matrix to covariance matrix
#The correlation matrix refers to the data in the paper
CORR = Matrix::forceSymmetric(as.matrix(read.table("SIGMA.txt")),uplo = "U")
SD = as.vector(rep(1,13)) #standard deviation vector
COV <-  Corr_to_Cov(CORR,SD)


#parameter settings
k = 12 #total k predictors in full regression
p = 7 #p predictors in reduced regression
n = 181  # sample size
inflation = inflation(n,p,k) #the inflation factor
alpha = 1 # intercept
SIGMA =  COV[-1,-1]
MU = rep(0,k+1)  # mean of covariates
BETA_star = c(-0.2, -0.03, -0.02, -0.04, 0.12, 0.18, 0.26, 0.25, -0.01, 0.08, -0.26, 0.21)   #effect size
simuN = 1e4
 
#calculation from covariance matrix
BETA_true = BETA_calcu(COV,k,k) #full model effects
#not same as BETA_true?
BETA_cal = BETA_calcu(COV=COV,p=p,k=k)
#error variance
sigmak2 = sigmak2(COV=COV)  #variance of error term in full regression
sigmap2 = sigmap2(COV=COV,p=p,k=k)  #variance of error term in  reduced regression
EVR = sigmak2/sigmap2 #error variance ratio

BETA_redu1 = solve(SIG_11)%*%sig_1#reduced model effects
#pPMSEr
pPMSEr_calc = pPMSEr(n = n, p = p, k = k, sigmak2 = sigmak2, sigmap2 = sigmap2)
#efficient sample size
n_star = eff_n(sigmak2 = sigmak2, sigmap2 = sigmap2,k = k,p = p) #sample size in the paper is 181
#cohen's f^2
R2 = 0.31; R2_1 = 0.06; 
f2_calcu = R2-R2_1/(1-R2)

ggplot(BETA_true, type = "l")
lines(BETA_redu1, col = "red")
lines(BETA_star, col = "blue")
