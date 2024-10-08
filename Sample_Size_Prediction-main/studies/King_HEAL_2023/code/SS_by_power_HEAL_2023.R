####-------------------------------------------------
#   Description: Based on assumed models 
#                this code get the statistical power and sample size (# of patients). # 
####-------------------------------------------------


 
	library(pwr);
	library(powerMediation);
	library(SIMR);



#T test based on Cohen’s d - difference between the means divided by the pooled standard deviation
meandiff = 1;
SD = 2.5;
d = meandiff/SD;
power=0.85;
alpha=0.05
pwr.t.test(n = NULL, d = d, sig.level = alpha, power = power, type = "two.sample", alternative="two.sided");


#Simple regression 
meandiff = 1;
SD = 2.5;
d = meandiff/SD;
power=0.85;
alpha=0.05
pwr.t.test(n = NULL, d = d, sig.level = alpha, power = power, type = "two.sample", alternative="two.sided");



###

####-------------------------------------------------
# Power calculation for regression model.
# Model assumption: Y = b0 + b1*X + epsilon, where X~N(0,sigma.x^2), epsilon ~ N(0, sigma^2)
####-------------------------------------------------
b1s = c(0.32, 0.5); #seq(.05, 0.32, by=.01);
	sigma.x =1;
	sigma = 1;
	sampleSizes = c(50, 100); #seq(1500, 3000, by=500); 
	
	R2=f2=corYX=array(NA, length(b1s));
	powerEst = array(NA, dim=c(length(sampleSizes), length(b1s)));
	
	for (i in 1:length(b1s)) {
	  b1=b1s[i]; 
	  R2[i] = b1^2*sigma.x^2 / (b1^2*sigma.x^2 + sigma^2); #Coefficient of determination = 1 - SS_resid/SS_total \approx 1 - sigma^2/Var(Y).
	  f2[i] = R2[i]/(1-R2[i]); #so called "effect size" =  b1^2*sigma.x^2 / sigma^2;
	  #b1 = sqrt(f2 * sigma^2 / sigma.x^2); #get b1 based on f2, sigma.x, and sigma.
	  
	  corYX[i] = b1*sigma.x^2 / (sigma.x*sqrt(b1^2*sigma.x^2 + sigma^2)); #Correlation coeff between X and Y. 

	  for (j in 1:length(sampleSizes)) {
      sampleSize = sampleSizes[j];
      #Power calculation by Cohen: regression model
      #powerEst[j, i] = pwr.f2.test(u=1, v=sampleSize-1-1, f2=f2[i], sig.level=0.1, power=NULL)$power;
      powerEst[j, i] = pwr.f2.test(u=1, v=sampleSize-1-1, f2=f2[i], sig.level=0.1/100, power=NULL)$power; #Bonferroni correction for 100 simultaneous tests
	  }
	}
	round(rbind(b1s, R2, f2, corYX, powerEst), 4);

	###Other calculation methods: 
  #Power calculation by Cohen: correlation model
  pwr.r.test(n=sampleSize, r=corYX, sig.level=0.05, power=NULL, alternative="two.sided");
   	
  #Use library "powerMediation"
  powerMediation.VSMc(n=sampleSize, b2=b1, sigma.m=sigma.x, sigma.e=sigma, corr.xm=0, alpha=0.05);
   	
   	
  ####-------------------------------------------------
  # Power calculation for logistic regression model.
  # Model assumption: logit(P(Y=1|X)) = b0 + b1*X, where X~N(0,sigma.x^2)
  ####-------------------------------------------------
  
  #preva = 4.17/10000; #The prevalence in population, i.e., P(Y=1) = E(P(Y=1|X)) = E(1/(1+exp(-b0-b1*X)));
	preva = 10.4/10000; #The prevalence in population, i.e., P(Y=1) = E(P(Y=1|X)) = E(1/(1+exp(-b0-b1*X)));
	
	#b0 = ; #For given prevalence, b1, and distn of X, we can get the corresponding b0.
	b1s = seq(.5, 3, by=.5);
	sigma.x =1;
	sampleSizes = seq(2000, 5000, by=500); 
	
	
	OR=array(NA, length(b1s));
	powerEst = array(NA, dim=c(length(sampleSizes), length(b1s)));
	
	for (i in 1:length(b1s)) {
	  b1=b1s[i]; 
	  OR[i] = exp(b1); #Odds ratio per unit increase of X. 

	  for (j in 1:length(sampleSizes)) {
	    sampleSize = sampleSizes[j];
	    #Power calculation by Cohen: regression model
 	    powerEst[j, i] = powerMediation.VSMc.logistic(n=sampleSize, b2=b1, sigma.m = sigma.x, p=preva, corr.xm=0,
                                                     alpha=0.1, verbose=F)$power;
	    #powerEst[j, i] = powerLogisticCon(n=sampleSize, OR=exp(b1), p1=preva, alpha=0.1); #Same as above
	  }
	}
	preva;
	round(rbind(b1s, OR, powerEst), 4);
	