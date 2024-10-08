#### Simulation based on 
#### #Baker TA, Buchanan NT, Corson N. Factorsinfluencing chronic pain intensity in older black women: examining depression, locus of control, and physical health.


library(MASS)
library(nlme);
library(simr);
source("~/Documents/GitHub/Sample_Size_Prediction/code/Lib_Prediction.R")

CORR = Matrix::forceSymmetric(as.matrix(read.table("SIGMA.txt")),uplo = "U");
SD = as.vector(rep(1,13)); #standard deviation vector
COV <-  sweep(sweep(CORR, 1, SD, "*"), 2, SD, "*");  

#parameters settings from Baker Table2
# Not using intercept 
#for basic predictors
b.age = -0.13; #*
b.edu = 0.07;
b.income = 0.1;
#for set 2 
b.como = 0.5;
b.location = 0.23;
b.medic = -0.15;
b.phys = 0.18; #**
#set 3
b.depress = -0.05; #*
b.satisf = -0.47;
b.loc.chance = 0.43;
b.loc.power = 0.14; #*
b.loc.internal = 0.21; #*

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
isRandomCV=F; #Random cross-validation in prediction 
nfold=5; #The number of folds in cross-validation
nrepeat=2; #number of repeats of cross-validation
simuN = 1000; #The number of simulations.
samplesizeN = seq(30, 600, by=30); #sample size levels


predNum = factor(1:9); # number of 'non basic' predictors over 3 basic predictors
models = vector(mode = "list", length(predNum)); #prediction models .
outputs = vector(mode = "list", length(predNum)); #prediction outputs using CORR generate matrix.
outputs_1 = vector(mode = "list", length(predNum)); #for contraction with y = XB+error
basic.output = vector(mode = "list", 1); #basic model output
basic.output_1 = vector(mode = "list", 1); #basic model output


#Store correlation and MSE over sample sizes and proportions of predictors using CORR generate matrix. 
CORR.mat = array(NA, dim=c(length(samplesizeN), length(predNum))); 
MSE = array(NA, dim=c(length(samplesizeN), length(predNum))); 

#Contraction with y = XB+error
CORR.mat_1 = array(NA, dim=c(length(samplesizeN), length(predNum))); 
MSE_1 = array(NA, dim=c(length(samplesizeN), length(predNum))); 

#Basic model output
Basic.CORR = array(NA,dim=c(length(samplesizeN), 1));
Basic.MSE = array(NA,dim=c(length(samplesizeN), 1));
Basic.CORR.1 = array(NA,dim=c(length(samplesizeN), 1));
Basic.MSE.1 = array(NA,dim=c(length(samplesizeN), 1));



#for added predictors
for(ni in 1:length(samplesizeN)){
  subjN = samplesizeN[ni]; #sample size 

  ###Looping through simulations
  for(si in 1:simuN) {
    
    ###Generate data  1. generate all predictors and response using the covariance matrix
    dat = as.data.frame(mvrnorm(n= subjN, mu = MU,Sigma = COV));
    names(dat) = c("y",names.pred)
    
    #####check the consistent 2. Using y = XB+ error to generate y
    pred.full = as.data.frame(mvrnorm(n= subjN, mu = MU[1:12],Sigma = COV[-1,-1]));
    names(pred.full) = names.pred;
    dat_1 = data.frame(y = as.matrix(pred.full)%*%b.fixed+rnorm(subjN, mean = 0, sd = sqrt(sigmak2)),
                       pred.full);
     
    #Basic model setting
    BasicTerm = paste("y ~", paste(names.pred[1:3], collapse = " + "));
    BasicTerm.formula = as.formula(paste("y ~", paste(names.pred[1:3], collapse = " + ")));
    ###Full model definition for lm
    full.model.formula.lm = as.formula(paste("y ~ ", paste(names.pred, collapse= " + ")));
                         
                     
    for (pi in 1:length(predNum)){
      ##Create model formula based on number of non-basic predictors used. 
      nonbasic.used.N   = as.numeric(predNum[pi]);
      nonbasic.names = names.pred[4:12]
      Partial.model.formula = as.formula(paste(BasicTerm, "+", 
                          paste(nonbasic.names[1:predNum[pi]], collapse= " + ")));

      
      #########################
      #Prediciton by lm
      #Using CORR generated y.
      out = meanPredEvaluCV.lme(fixed=Partial.model.formula, dat=dat,model_R='lm', 
                                randomf=NULL, loopn=nrepeat, cvNumber=nfold); 
            outputs[[pi]] = rbind(outputs[[pi]], 
                            t(c(MSE=out[1], L2normRatio=out[2], 
                                L1normRatio=out[3], correlation=out[4], 
                                MSEoverObsVar=out[5])));
      #Using y = XB+ error 
      out_1 = meanPredEvaluCV.lme(fixed=Partial.model.formula, dat=dat_1,model_R='lm', 
                                randomf=NULL, loopn=nrepeat, cvNumber=nfold);
      outputs_1[[pi]] = rbind(outputs_1[[pi]], 
                            t(c(MSE=out_1[1], L2normRatio=out_1[2], 
                                L1normRatio=out_1[3], correlation=out_1[4], 
                                MSEoverObsVar=out_1[5])));
    }
    #basic model Using CORR generated y
    basic.out = meanPredEvaluCV.lme(fixed=BasicTerm.formula, dat=dat,model_R='lm', 
                                    randomf=NULL, loopn=nrepeat, cvNumber=nfold); 
    basic.output[[1]] = rbind(basic.output[[1]], 
                       t(c(MSE=basic.out[1], L2normRatio=basic.out[2], 
                           L1normRatio=basic.out[3], correlation=basic.out[4], 
                           MSEoverObsVar=basic.out[5])));
    
    #basic model Using y = XB +error
    basic.out_1 = meanPredEvaluCV.lme(fixed=BasicTerm.formula, dat=dat_1,model_R='lm', 
                                    randomf=NULL, loopn=nrepeat, cvNumber=nfold); 
    basic.output_1[[1]] = rbind(basic.output_1[[1]], 
                              t(c(MSE=basic.out_1[1], L2normRatio=basic.out_1[2], 
                                  L1normRatio=basic.out_1[3], correlation=basic.out_1[4], 
                                  MSEoverObsVar=basic.out_1[5])));
  }
  
  #Prediction accuracies
  for (pi in 1:length(predNum))
    {
      CORR.mat[ni, pi] = apply(outputs[[pi]], 2, mean)[4];
      MSE[ni, pi] = apply(outputs[[pi]], 2, mean)[1];
      
      CORR.mat_1[ni, pi] = apply(outputs_1[[pi]], 2, mean)[4];
      MSE_1[ni, pi] = apply(outputs_1[[pi]], 2, mean)[1];
  }
  
  Basic.CORR[ni,1] = apply(basic.output[[1]], 2, mean)[4];
  Basic.MSE[ni,1] = apply(basic.output[[1]], 2, mean)[1];
  
  Basic.CORR.1[ni,1] = apply(basic.output_1[[1]], 2, mean)[4];
  Basic.MSE.1[ni,1] = apply(basic.output_1[[1]], 2, mean)[1];
  
  print(subjN);
} 

###Final output
data.frame( b.age, b.edu, b.income, b.como, b.location, b.medic, b.phys, b.depress,
           b.satisf, b.loc.chance, b.loc.power, b.loc.internal);

#output CORR and MSE for multinormal y
CORR.output = as.data.frame(cbind(samplesizeN, Basic.CORR, CORR.mat));
names(CORR.output) = c("Sample Size", "Basic Model", nonbasic.names);
MSE.output = as.data.frame(cbind(samplesizeN, Basic.MSE, MSE));
names(MSE.output) = c("Sample Size", "Basic Model", nonbasic.names);

#output CORR and MSE table for y = XB+error
CORR.output.1 = as.data.frame(cbind(samplesizeN, Basic.CORR.1, CORR.mat_1));
names(CORR.output.1) = c("Sample Size", "Basic Model", nonbasic.names);
MSE.output.1 = as.data.frame(cbind(samplesizeN, Basic.MSE.1, MSE_1));
names(MSE.output.1) = c("Sample Size", "Basic Model", nonbasic.names);

#output csv files
write.csv(CORR.output, "PMSEsimu_CORR.csv")
write.csv(MSE.output, "PMSEsimu_MSE.csv")
write.csv(CORR.output.1, "PMSEsimu_CORR_1.csv")
write.csv(MSE.output.1, "PMSEsimu_MSE_1.csv")



1. compare part from computation and simulation 
2. corr_1 and MSE_1 observation file and code !
3. function code description and arguement. !
4. keep updating
5. upload presentation !

  
1. add PMSE calcu table
2. change the comparison file
3. analysis the difference between pmse.calcu, PMSE.cov, pmse.linear, 