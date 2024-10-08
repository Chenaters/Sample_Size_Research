####-------------------------------------------------
#   Sample size analysis based on GLMM by SIMR  
# Reference: 
#   https://humburg.github.io/Power-Analysis/simr_power_analysis.html
####-------------------------------------------------


library(simr);

### Bonferroni control for Type I error rate in multiple-hypothesis testing
  prediN=100; #number of candidate predictors to be tested simultaneously
  alphaLevel = 0.05/prediN; 
  
  subjN =400; #Restriction for the number of subjects
  breaksN = seq(150, subjN, by=50); #Numbers of subjects at which power values are obtained
  #breaksN = seq(50, 100, by=10);
  

###Model settings and parameters
  ###Type of responses and predictors
  #respoType = "contiRespo"; #continuous response
  respoType = "binaryRespo"; #binary response
  
  prediType = "contiPredi"; #assume normal distn with mean=0 and sd=1
  #prediType = "binaryPredi"; #assume balanced distribution along subjects
  
  ###Subjects and time variable
  time = c(0, 1, 2, 6); #baseline, 1, 2, and 6 months; consider it quantitative (not factor)
  
  ###Intercept and slopes
  intercept = -0.5; 
  b.predi = 0.1;
  b.time = 0.2;
  b.predi.time = 0.1; 
  fixed <- c(intercept, b.predi, b.time, b.predi.time);
  
  ### Random intercept to account for the intra-subject correlation over time points 
  subj <- factor(1:subjN); #subjects are factors 
  cor.subj = 0.1; #intra-subject correlation
  rand <- list(cor.subj);
  
  ### residual variance for continous response
  res <- 2.4;
  
###Generate data
  subj_full <- rep(subj, length(time)); #Data of subjects
  time_full <- rep(time, each=subjN);   #Data of times
  
  #Data of predictor
  if (prediType == "binaryPredi") {
    predi = 0:1;
    predi_full = rep(rep(predi, each=subjN/2),  length(time)); #assume balanced distribution
  }else {
    predi_full = rnorm(subjN*length(time));
  }
  
  ###Combine into a data fame
  covars <- data.frame(id=subj_full, time=time_full, predictor=predi_full);


###Generate model
  if (respoType == "binaryRespo") {
    model <- makeGlmer(y ~ predictor*time + (1|id), family="binomial", fixef=fixed, VarCorr=rand, data=covars);
    dat = getData(model); #get the simulated data. 
    print(mean(dat$y)); #print the frequency of binary response.
    #print(mean(simY[which(covars$time==0)])); #print the frequency of binary response.
  }else {
    model <- makeLmer(y ~ predictor*time + (1|id), sigma=res, fixef=fixed, VarCorr=rand, data=covars);
    dat = getData(model); #get the simulated data. 
    print(c(mean(dat$y[which(dat$time==0)]), sd(dat$y[which(dat$time==0)]))); #print the baseline mean and sd of continuous response
  }

  

###Power analysis. The null model is that predictor has no effects: y~time. 
 
  ###Power at the given data
  #sim_predictor <- powerSim(model, nsim=100, test = fcompare(y~time), alpha=alphaLevel); 
  #sim_predictor;
  
  ###Power over a series of sample size.
  p_curve_predictor <- powerCurve(model, test=fcompare(y~time), along="id", breaks=breaksN, alpha=alphaLevel, progress=F);
  
  ###Output results
  model; #recall the model info
  outPath = "/Users/zheyangwu/Grants/Power_Pred_Studies/Studies/King_HEAL_2023/output/";
  fileName = paste("Pow_", respoType, "_", prediType, "_b0_", intercept, "_b.predi_", b.predi, "_b.time_", b.time, "_b.predi.time_", 
                   b.predi.time, sep="");
  simY = doSim(model); #simulate response, each run is a random generation (without setting seeds)
  if (respoType == "binaryRespo") {
    propY = round(mean(dat$y), 1); #proportion of Y (positive response)
    print(propY);
    fileName = paste(fileName, "_propYb_", propY, "_prediN_", prediN,  ".png", sep="");
  }else {
    meanYb = round(mean(dat$y[which(dat$time==0)]), 1); #mean and sd of Y at baseline 
    sdYb = round(sd(dat$y[which(dat$time==0)]), 1); 
    print (c(meanYb, sdYb));
    fileName = paste(fileName, "_meanYb_", meanYb, "_sdYb_", sdYb, "_prediN_", prediN, ".png", sep="");
  }
  
  ##the power curve
  png(paste(outPath, fileName, sep=""));
  plot(p_curve_predictor);
  dev.off();
  


  #---- SIMR allows for changing effect size or number of subjects based on model.
# ###Changing the effect size. This process is equivalent to above process with changed slope for the interaction term.
#   model2 <- model;
#   fixef(model2)['predictor:time'] <- 0.1;
#   
#   sim_predictor2 <- powerSim(model2, nsim=100, test = fcompare(y~time), alpha=alphaLevel, progress=F);
#   sim_predictor2;
# 
# 
# ### Changing the number of participants per race
#   #The increased rows simply replicate the original y values. 
#   #Not good b/c the observations are not independent anymore? 
#   #Although, the power output is the similar to that by changing the subjN. 
#   
#   model_ext_subj <- extend(model, along="id", n=20) #The increased rows simply replicate the original y values. Not good b/c the observations
#   model_ext_subj;
#   nrow(getData(model_ext_subj));
#   covars_ext = getData(model_ext_subj);
#   
#   sim_predictor_subj <- powerSim(model_ext_subj, nsim=100, test = fcompare(y~time), alpha=alphaLevel);
#   sim_predictor_subj; #This power is 
#   
#   
# ###Draw power curve  
#   p_curve_predictor <- powerCurve(model_ext_subj, test=fcompare(y~time), along="id", breaks=c(20), alpha=alphaLevel, progress=F);
#   plot(p_curve_predictor);
#   
