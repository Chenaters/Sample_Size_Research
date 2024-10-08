#simulation final version for thesis

library(MASS)
source("NARULA/sto_simu _functions.R")
k = 12
p = 3
predNum = seq(3,12);
#sigmak2.full = sigmak2(COV);
b.fixed = as.matrix(c(b.age, b.edu, b.income, b.como, b.location, b.medic, b.phys, b.depress,
                      b.satisf, b.loc.chance, b.loc.power, b.loc.internal));

names.pred = c("Age", "Education", "Income", "Comorbidities",
               "PainLocations", "Medications", "PhysicalFunctioning",
               "DepressiveSymptoms", "LifeSatisfaction",
               "LOCchance", "LOCpowerful", "LOCinternal")
N.pred = length(names.pred)
#suppose predictors data are standardized
MU = rep(0,(N.pred+1))


### residual variance for continous response
sigmak2 <- 0.4687399; # from calculation 


###Parameters on prediction process
simuN = 5000; #The number of simulations.
samplesizeN = seq(30, 600, by=30); #sample size levels

predNum = factor(1:9); # number of 'non basic' predictors over 3 basic predictors
outputs = vector(mode = "list", length(predNum)); #prediction outputs using CORR generate matrix.
basic.output = vector(mode = "list", 1); #basic model output


#Store correlation and MSE over sample sizes and proportions of predictors using CORR generate matrix. 
PMSE = array(NA, dim=c(length(samplesizeN), length(predNum))); 

#Basic model output
Basic.PMSE = array(NA,dim=c(length(samplesizeN), 1));

PMSE.simu.matrix = matrix(nrow = length(samplesizeN), ncol = length(predNum))

for(ni in 1:length(samplesizeN)){
  subjN = samplesizeN[ni]; #sample size 
  
  ###Looping through simulations
  for(si in 1:simuN) {
    ###Generate data  1. generate all predictors and response using the covariance matrix
    dat = as.data.frame(mvrnorm(n= subjN, mu = MU,Sigma = COV));
    names(dat) = c("y",names.pred)
    newData = as.data.frame(t(mvrnorm(n=1, mu=MU, Sigma=COV))); 
    names(newData) = c("y",names.pred);
    
    #Basic model setting
    BasicTerm = paste("y ~", paste(names.pred[1:3], collapse = " + "));
    BasicTerm.formula = as.formula(paste("y ~", paste(names.pred[1:3], collapse = " + ")));
    ###Full model definition for lm
    full.model.formula.lm = as.formula(paste("y ~ ", paste(names.pred, collapse= " + ")));
    
    for (pi in 1:length(predNum)) {
      lm.formula = as.formula(paste("y ~", paste(names.pred[1:(3+pi)], collapse = " + ")));
      lmfit = lm(lm.formula, data=dat);
      Y.hat.0 = predict(lmfit, newData); 
      outputs[[pi]] = rbind(outputs[[pi]], (newData$y-Y.hat.0)^2);
    }
    lmfit.basic = lm(BasicTerm.formula, data=dat);
    Y.hat.0 = predict(lmfit.basic, newData); 
    basic.output[[1]] = rbind(basic.output[[1]], (newData$y-Y.hat.0)^2);
  }
  
  #Prediction accuracies
  for (pi in 1:length(predNum))
  {
    PMSE[ni, pi] = apply(outputs[[pi]], 2, mean)[1];
    }
  Basic.PMSE[ni,1] = apply(basic.output[[1]], 2, mean)[1];
  
  print(subjN);
} 

###Final output
#data.frame( b.age, b.edu, b.income, b.como, b.location, b.medic, b.phys, b.depress,
           # b.satisf, b.loc.chance, b.loc.power, b.loc.internal);

#output  and MSE for multinormal y
MSE.output = as.data.frame(cbind(samplesizeN, Basic.PMSE, PMSE));
names(MSE.output) = c("Sample Size", "Basic Model", nonbasic.names);
write.csv(MSE.output, file = "PMSEApr22.csv")
