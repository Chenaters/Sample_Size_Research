
####
# Sample size by simulations of predictive models of quantitative responses
# For HEAL proposal 2023
####

library(MASS)
library(nlme);
library(simr);

source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Prediction/Lib_Prediction.R"); 	
source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Simulations/Simulate_Data/Lib_simu_genetic_data.R"); 
source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Prediction/generate_response_formula.R"); 	

###Model settings and parameters
  ###Type of responses and predictors
  respoType = "contiRespo"; #continuous response
  #respoType = "binaryRespo"; #binary response

  ###Subjects and time variable
  time = c(0, 1, 2, 6); #baseline, 1, 2, and 6 months; consider it quantitative (not factor)
  
  ###Intercept and slopes
  intercept = 4.7; 
  b.predi = 0.2;
  b.time = 0.2;
  b.predi.time = 0.1; 
  
  biomN = 20;                   #The # of biomarkers
  b.biom =1;                #The value of the coeff of the biomarkers
  B.biom = rep(b.biom, biomN);  #Coeff of the biomarkers
  names.xbiom = paste("xbiom", 1:biomN, sep=""); #variable names for biomarkers
  
  fixed <- c(intercept, b.predi, b.time, b.predi.time, B.biom);
  
  ### Random intercept to account for the intra-subject correlation over time points 
  cor.subj = 0.1; #intra-subject correlation
  rand <- list(cor.subj);
  
  ### residual variance for continous response
  res <- 2;


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
  subjNs = seq(100, 300, by=50);
  
  CORR = array(NA, dim=c(length(subjNs), length(predProp))); #Store correlation over sample sizes and proportions of predictors used. 
  AUC = array(NA, dim=c(length(subjNs), length(predProp))); #Store AUC over sample sizes and proportions of predictors used. 

for(li in 1:length(subjNs)){
  subjN = subjNs[li]; #sample size for each of the three groups: control, mbsr, and acupuncture. 
  subj <- factor(1:subjN); #subjects are factors 
  
  ###Looping through simulations
  for(si in 1:simuN) {
      ###Generate data
      subj_full <- rep(subj, length(time)); #Data of subjects
      time_full <- rep(time, each=subjN);   #Data of times
      predi_full = rnorm(subjN*length(time)); #Data of predictor
      
      biom_full = as.data.frame(matrix(rnorm(subjN*length(time)*biomN), ncol=biomN, dimnames = list(NULL, names.xbiom)));
    
      ###Combine into a data fame
      covars <- data.frame(id=subj_full, time=time_full, predictor=predi_full, biom_full);
    
      ###Model definition for lmer
      full.model.formula.lmer = as.formula(paste("y ~ predictor*time +", paste(names.xbiom, collapse= " + "), "+ (1|id)")); 

      
      ### Generate response 
      if (respoType == "binaryRespo") {
        model <- makeGlmer(full.model.formula.lmer, family="binomial", fixef=fixed, VarCorr=rand, data=covars);
      }else {
        model <- makeLmer(full.model.formula.lmer, sigma=res, fixef=fixed, VarCorr=rand, data=covars);
      }
      dat = getData(model); #get the simulated data. 

      for (pi in 1:length(predProp)){
        ##Create model formula based on proportion of predictors used. 
        fixedTerm = paste("y ~ predictor*time");
        xbiom.used.N   = round(length(names.xbiom)*predProp[pi]);
        if(xbiom.used.N > 0){
          fixedTerm = paste(fixedTerm, "+", paste(names.xbiom[0:xbiom.used.N], collapse= " + "));
        }
        fixedTerm  = as.formula(fixedTerm);
        randomTerm = as.formula("~1 | id"); 
        
        if (respoType == "binaryRespo") {
          out = predEvaluCV.glm(formula=fixedTerm, dat=dat, nfold=nfold, nrepeat=nrepeat, isRandomCV=isRandomCV);
          outputs[[pi]] = rbind(outputs[[pi]], t(unlist(out)));
        }else{
          out = meanPredEvaluCV.lme(fixed=fixedTerm, dat=dat, randomf=randomTerm, model_R='lme', predLevel=1, loopn=nrepeat, cvNumber=nfold); #Prediction by lme.
          #out = meanPredEvaluCV.lme(fixed=models[[pi]]$fixed, dat=dat, randomf=NULL, model_R='lm', loopn=nrepeat, cvNumber=nfold); #Prediciton by lm.
          outputs[[pi]] = rbind(outputs[[pi]], t(c(MSE=out[1], L2normRatio=out[2], L1normRatio=out[3], correlation=out[4], MSEoverObsVar=out[5])));
        }
      }
  }

  #Prediction accuracies
  for (pi in 1:length(predProp)){
    if (respoType == "binaryRespo") {
      AUC[li, pi] = apply(outputs[[pi]], 2, mean)[4];
    }else{
      CORR[li, pi] = apply(outputs[[pi]], 2, mean)[4];
    }
  }
  
  print(subjN);
} 

###Final output
  data.frame(respoType, intercept, b.predi, b.time, b.predi.time, biomN, b.biom);
  if (respoType == "binaryRespo") {
    print("percentage of response:");
    print(mean(dat$y)); #print the frequency of binary response.
    print(cbind(subjNs, AUC));
  }else{
    print("the baseline mean and sd of continuous response:");
    print(c(mean(dat$y[which(dat$time==0)]), sd(dat$y[which(dat$time==0)]))); #print the baseline mean and sd of continuous response
    print(cbind(subjNs, CORR));
  }

