####
# Minimal sample size calculation for prediction based on pmsamplesize.
# Reference: ?pmsampsize
#            Descriptions and examples given in paper \cite{riley2020calculating} Calculating the sample size required for developing a clinical prediction model
####

library(pmsampsize);
library(rms);
source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Prediction/mimimum_samplesize_prediction.R"); 	


############## 
### Predict quantitative response by regression model.
## Response: PEG range is 0-10. Consider mean 5.85 and sd 2.43: 
  respMean = 5.85; 
  respSD = 2.43; 
  R2s = seq(0.4, 0.9, by=0.1); #Coefficient of determination 
  predN = seq(10, 50, by = 5); #Number of predictors
  pNs = predN + 2; #The number of parameters in predictive model (# of predictors + intercept and error sd)
  Ns = array(NA, dim=c(length(R2s), length(pNs))); #The sample sizes
  for (ri in 1:length(R2s)){
    for (pi in 1:length(pNs)){
      Ns[ri, pi] = pmsampsize(type = "c", rsquared =R2s[ri], parameters = pNs[pi], intercept = respMean, sd = respSD)$sample_size; 
    }
  }
  respMean; respSD; 
  R2s; 
  predN;
  Ns;



############## 
### Predict binary response by logit model.
## Response: Pain relief Consider outcome prevalence is
# anticipated to be 0.3, and C-statistic (i.e., AUC)
# Note: Cox and Snell R2 is R2_C&S = 1 – (L0 / LM)2/n https://statisticalhorizons.com/r2logistic
#       We can either use R2_C&S or C-statistic
  prevalence=0.7; #The prevalence of outcome.
  shrinkage = 0.85; 
  delta.NagR2 = 0.1; #Small absolute difference in the model's apparent and adjusted Nagelkerke's R-squared value
  delta.aveRisk = 0.1;
  Cstats = seq(0.7, 0.95, by=0.05); # C-statistic (i.e., AUC)
  predN = seq(10, 50, by = 5); # The # of predictors
  pNs = predN + 1;        #The number of parameters in predictive model= # of predictors + intercept
  Ns = array(NA, dim=c(length(Cstats), length(pNs)));
  for (ci in 1:length(Cstats)){
    for (pi in 1:length(pNs)){
      #size by criteria 1
      n1=pmsampsize(type = "b", cstatistic =Cstats[ci], parameters = pNs[pi], prevalence=prevalence, shrinkage=shrinkage)$results_table[1, 1]; 
      #size by criteria 2
      n2 = round(minSampleSize.NagR2(P=pNs[pi], Cstat=Cstats[ci], prevalence=prevalence, delta=delta.NagR2));
      #size by criteria 3
      n3 = round(minSampleSize.aveRisk(prevalence=prevalence, delta=delta.aveRisk));
      Ns[ci, pi] = max(n1, n2, n3);
    }
  }
  c(prevalence=prevalence, shrinkage = shrinkage, delta.NagR2 = delta.NagR2, delta.aveRisk=delta.aveRisk); #settings for min sample size calculation
  Cstats;
  predN
  Ns;





