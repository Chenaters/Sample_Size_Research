
####
# Sample size by simulations of predictive models of binary responses
# For HEAL proposal 2023
####

library(MASS)
library(nlme);
source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Prediction/Lib_Prediction.R"); 	
source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Simulations/Simulate_Data/Lib_simu_genetic_data.R"); 
source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Prediction/generate_response_formula.R"); 	




###Parameter setting

###Parameters on the effects / coefficients
beta0 = -0.5; #For binary trait, the risk at X=0 is 1/(1+exp(-beta0)). 
#The average risk (prevalence) at given coefficients and predictors can be calculated using the 
#code from function get.Y.logit: First set beta0=0, then run the code to get Xmatrix, 
#then set XData=Xmatrix, then use the code in get.Y.logit to get Xbeta, 
#then choose beta0's value such that mean(Y) reach the wanted prevalence. 
# Y = array(NA, length(Xbeta));
#  for (i in 1:length(Y)) {
#    Py1x = 1/(1 + exp(-Xbeta[i]) - beta0); # P(y=1 | x) based on logistic model
#    Y[i] = ifelse(runif(1) < Py1x, 1, 0);
#  }
#  mean(Y); #prevalance of Y.

#beta0 = -7; # for binary trait. Risk when X=0 is 1/(1+exp(7))=0.0009110512. Also, at the given coefficients and factors, beta0 = -7 makes the prevalence of Y roughly 0.5. 
#beta0 = -3; # for binary trait. Risk when X=0 is 1/(1+exp(3))=0.04742587. Also, at the given coefficients and factors, beta0 = -3 makes the prevalence of Y, i.e., mean(Y), roughly 0.7. 


raceN = 3; #Number of ethnic groups
sigma.race = 0.05; #pair-wise covariance among individuals in the same ethnic group
beta.race = 1; #treat the "coefficient" of the mix-effect b_i be 1. 


## for binary outcome
beta.sex = log(1.1); #1;
beta.age = log(1.25); #-1;
beta.edu = log(1.26); #1;
beta.len = log(1.13); #1;
beta.bas = log(0.80);#1;
beta.conc = log(0.77); #1;


###Parameters on prediction process
isRandomCV=T; #Random cross-validation in prediction 
nfold=5; #The number of folds in cross-validation
nrepeat=2; #number of repeats of cross-validation

###Parameters on simulations
simuN = 100; #The number of simulations.

###Data simulation and prediction outcomes
predProp = c(0, 0.25, 0.5, 0.75, 1); #Proportion of true predictors besides names.basic that are included in prediction model. 
models = vector(mode = "list", length(predProp)); #prediction models .
outputs = vector(mode = "list", length(predProp)); #prediction outputs .

###Parameters on data
groupSampleSizes = seq(20, 120, by=20);
#groupSampleSizes = seq(500, 2000, by=100);
#groupSampleSizes = seq(150, 500, by=50);
AUC = array(NA, dim=c(length(groupSampleSizes), length(predProp))); #Store AUC over sample sizes and proportions of predictors used. 

for(gi in 1:length(groupSampleSizes)){
  groupSampleSize = groupSampleSizes[gi]; #sample size for each of the three groups: control, mbsr, and acupuncture. 
  
  ###Looping through simulations
  for(i in 1:simuN) {
    ###Generate data
    x0 = rep(1, groupSampleSize*3);
    
    #The mixed-effect term for racial group 
    xrace = sample(1:raceN, size=groupSampleSize*3, replace=T, prob=rep(1/raceN, raceN));
    #Assume equal chance for each racial group to be sampled. 
    b.xrace = array(NA, dim=groupSampleSize*3); #b.xrace is the vector of b_i values. 
    for (racei in 1:raceN){ b.xrace[which(xrace==racei)] = rnorm(1, sd=sigma.race); }
    #assign the same b_i value for the all in the ith racial group. 
    
    #The "basic" factors
    xsex = rbinom(n=groupSampleSize*3, size=1, prob=0.5); #50% recruited are males??
    xage = rnorm(n=groupSampleSize*3, mean=0, sd=1); #standardized age. 
    xedu = rbinom(n=groupSampleSize*3, size=1, prob=0.5);
    xlen = rbinom(n=groupSampleSize*3, size=1, prob=0.5);
    xbas = rnorm(n=groupSampleSize*3, mean=0, sd=1);
    xconc = rbinom(n=groupSampleSize*3, size=1, prob=0.5);
    
    #PainMarker data
    xbiom = data.frame(matrix(rnorm(n=groupSampleSize*3*markerN), ncol=markerN)); #Assume biomarker values are N(0,1);
    names(xbiom) = names.xbiom; 
    
    Xmatrix = cbind(x0, b.xrace, xsex, xage, xedu, xlen, xbas, xconc, xbiom);
    
    #### Generate response 
    names.mainEff  = c("x0", "b.xrace", "xsex", "xage", "xedu", "xlen", "xbas", "xconc", names.xbiom);
    coeffs.mainEff = c(beta0, beta.race, beta.sex, beta.age, beta.edu, beta.len, beta.bas, beta.conc, beta.biom);
    
    ####----Binary Response-----
    resp = get.Y.logit(XData=Xmatrix, names.mainEff, coeffs.mainEff);
    
    ####Combine data for analysis 
    xrace = as.factor(xrace); #Convert race into factor variable, which is used in data analysis.
    dat = data.frame(Y=resp$Y, Xmatrix, xrace);
    
    
    ####Predictive analysis
    names.control = c("xrace", "xsex", "xage", "xedu", "xlen", "xbas", "xconc"); #Controlling factors in prediction models
    for (mi in 1:length(predProp)){
      ##Create model formula based on proportion of predictors used. 
      xbiom.used    = round(length(names.xbiom)*predProp[mi]);
      names.main=c(names.control, names.xbiom[0:xbiom.used]);
      models[[mi]]= formula.f.r(names.main=names.main, names.random=c("1"), names.group=c("xrace"));
      
      ### ---- Predict binary outcome ----
      out = predEvaluCV.glm(formula=models[[mi]]$fixed, dat=dat, nfold=nfold, nrepeat=nrepeat, isRandomCV=isRandomCV);
      #outputs[[mi]] = rbind(outputs[[mi]], t(c(prob=out[1], sensi=out[2], speci=out[3], AUC=out[4])));
      outputs[[mi]] = rbind(outputs[[mi]], t(unlist(out)));
    }
  }
  
  #Prediction accuracies
  for (mi in 1:length(predProp)){
    ##Binary trait
    AUC[gi, mi] = apply(outputs[[mi]], 2, mean)[4];
  }
  print(c(sampleSize=groupSampleSize*3));
} 

# #True underlying model
c(prevalence=1/(1+exp(-beta0)), raceN=raceN, sigma.race = sigma.race,beta.sex = beta.sex, beta.age = beta.age, beta.edu=beta.edu, 
        beta.len=beta.len, beta.bas=beta.bas, beta.conc=beta.conc, markerN = markerN, beta.biom.v = beta.biom.v);

cbind(groupSampleSizes*3, AUC); 

