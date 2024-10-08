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
## Response: Pain impact range is 8-50. Consider mean 29 and sd 20: 
  R2s = seq(0.4, 0.9, by=0.1); #Coefficient of determination 
  pNs = seq(10, 50, by = 5); #The number of parameters in predictive model
  Ns = array(NA, dim=c(length(R2s), length(pNs)));
  for (ri in 1:length(R2s)){
    for (pi in 1:length(pNs)){
      Ns[ri, pi] = pmsampsize(type = "c", rsquared =R2s[ri], parameters = pNs[pi], intercept = 29, sd = 20)$sample_size; 
    }
  }
  R2s;
  pNs;
  Ns;



############## 
### Predict binary response by logit model.
## Response: Pain response. Consider outcome prevalence is
# anticipated to be 0.5, and C-statistic (i.e., AUC)
# Note: Cox and Snell R2 is R2_C&S = 1 – (L0 / LM)2/n https://statisticalhorizons.com/r2logistic
#       We can either use R2_C&S or C-statistic
  prevalence=0.5; #The prevalence of outcome.
  shrinkage = 0.9; 
  delta.NagR2 = 0.05; #Small absolute difference in the model's apparent and adjusted Nagelkerke's R-squared value
  delta.aveRisk = 0.05;
  Cstats = seq(0.7, 0.95, by=0.05); # C-statistic (i.e., AUC)
  pNs = seq(10, 50, by = 5);        #The number of parameters in predictive model
  Ns = array(NA, dim=c(length(R2s), length(pNs)));
  for (ci in 1:length(Cstats)){
    for (pi in 1:length(pNs)){
      #size by criteria 1
      n1=pmsampsize(type = "b", cstatistic =Cstats[ci], parameters = pNs[pi], prevalence=prevalence, shrinkage=shrinkage)$results_table[1, 1]; 
      #size by criteria 2
      n2 = round(minSampleSize.NagR2(P=pNs[pi], Cstat=Cstats[ci], prevalence=prevalence, delta=delta.NagR2));
      #size by criteria 3
      n3 = round(minSampleSize.aveRisk(prevalence=prevalence, delta=delta));
      Ns[ci, pi] = max(n1, n2, n3);
    }
  }
  c(prevalence=prevalence, shrinkage = shrinkage, delta.NagR2 = delta.NagR2, delta.aveRisk=delta.aveRisk); #settings for min sample size calculation
  Cstats;
  pNs
  Ns;





