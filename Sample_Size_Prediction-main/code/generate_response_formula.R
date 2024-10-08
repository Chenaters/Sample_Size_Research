####
# Functions for supporting sample size study 
# Generate responses by regression and logit that could contain interactions.
# Generate formulas for regression and logit that could contain interactions. 
# ZWu
####



##################################
# Function for generating the fixed and random formula objects for main and interaction terms 
#  (interactions are based on treatment and modifiers). 
#  The formula objects are consistent with those in the lme function.
##################################
# Arguments: 
# names.main: Vector of variable names for "main terms"
# name.resp:  Variable name of response. 
# names.trt: Vector of variable names of the treatments. If NULL, then no interaction terms. 
# names.mdf: List of variable names of the modifiers. Require: Each element in the list is a verctor of modifier names corresponding (in order) to the treatment varible in names.trt.
# names.random: Vector of varaible names for the random effects. 
# names.group: Vector of variable names for grouping structure. 
# Output: formula objects of fixed effects and the random effects. 
# Author: ZWu
# Example: 
# name.resp = "Y";
# names.main = c("xrace", "xacup", "xmbsr", "xsex", "xage");
# names.trt = c("xacup", "xmbsr"); #variable names of the treatments
# names.mdf.acup = paste("zmdfAcup", 1:3, sep="");
# names.mdf.mbsr = paste("zmdfMBSR", 1:3, sep="");
# names.mdf = list(names.mdf.acup, names.mdf.mbsr);
# names.random = c("1"); #c("x1", "x2");
# names.group =  c("xrace"); #NULL; #c("g1", "g2");
# formula.f.r(names.main=names.main, names.trt=names.trt, names.mdf=names.mdf, names.random=names.random, names.group=names.group);
#
formula.f.r = function(names.main, name.resp = "Y", names.trt=NULL, names.mdf=NULL, names.random=c("1"), names.group=NULL) {
  
  fixed = paste(paste(name.resp, "~"), paste(names.main, collapse= "+")); 
  if (!is.null(names.trt)) {
    for (trti in 1:length(names.trt)){
      xtrt = names.trt[trti];
      xmdf = names.mdf[[trti]];
      fixed = paste(fixed, "+", paste(paste(xtrt, "*", sep=""), xmdf, collapse="+"));
    }
  }
  
  random = paste("~ ", paste(names.random, collapse="+"), "|"); 
  if (!is.null(names.group)) random = paste(random, paste(names.group, collapse="/"));
  
  return(list(fixed=as.formula(fixed), random=as.formula(random)));
}



##################################
# Function for generating quantitative response based on regression model that contains main effects and interaction effects. 
#  The interaction terms are based on treatment variables and their modifiers. 
# Arguments:
#     XData: Data frame of independent variables with column names being the names of the corresponding variables.
#     names.mainEff: Vector of variable names of the main-effect variables.
#     coeffs.mainEff: Vector of coefficients of the main-effect variables.
#     names.trt: Vector of variable names of the treatments. No interactions if names.trt is NULL. 
#     names.mdf: List of variable names of the modifiers corresponding to the treatments. Require: Each list element is a verctor of modifier names corresponding to the treatment in names.trt.
#     coeffs.interaction = list(gama.acup, gama.mbsr); #coefficients of the treatment-modifier interaction terms. Require: Consistent with names.trt and names.mdf.
#     errSD: Standard deviation of the error term in the model. 
# Output: Response vector Y. 
#     R2: coefficient of determination of all predictors, i.e., the proportion of variation explained by all predictors
#     f2: Cohen's f^2. Consistent with R2. 
# Author: ZWu
# Example: 
#   groupSampleSize = 100; #sample size for each of the three groups: control, mbsr, and acupuncture. 
# 
#   ###==== Effects / coefficients ====
#   beta0 = 0; #intercept effect
# 
#   raceN = 3; #Number of racial groups
#   sigma.race = 1; #pair-wise covariance among individuals in the same racial group
#   beta.race = 1; #treat the "coefficient" of the mix-effect b_i be 1. 
# 
#   beta.acup = 5; #Coeff of acupuncture
#   beta.mbsr = 5; #Coeff of mbsr
#   beta.sex = 1; #positive value assumes males are easier to reduce pain than females??
#   beta.age = -1; #negative value assume the olders are hard to reduce pain ??
# 
#   markerN = 10; #The # of PainMarkers
#   beta.biom.v = 1; #The value of the coeff of the biomarkers
#   beta.biom = rep(beta.biom.v, markerN); #Coeff of the biomarkers
#   names.xbiom = paste("xbiom", 1:markerN, sep=""); #variable names for biomarkers
# 
#   mdfN.acup = 5; #The # of modifiers for acupuncture
#   beta.mdf.acup.v=1; #The value of the coeff of accupunctur's modifiers
#   beta.mdf.acup = rep(beta.mdf.acup.v, mdfN.acup); #Vector of coeffs of accupunctur's modifiers 
#   gama.acup.v = 1; #The value of the coeff of the modifier*acupucture interaction terms.
#   gama.acup = rep(gama.acup.v, mdfN.acup); #Vector of coeffs of the modifier*acupucture interaction terms.
# 
#   names.mdf.acup = paste("zmdfAcup", 1:mdfN.acup, sep=""); #variable names of the modifiers for acupuncture
# 
#   mdfN.mbsr = 5; #The # of modifiers for MBSR
#   beta.mdf.mbsr.v=1; #The value of the coeff of MBSR's modifiers
#   beta.mdf.mbsr = rep(beta.mdf.mbsr.v, mdfN.mbsr); #Vector of coeffs of accupunctur's modifiers 
#   gama.mbsr.v = 1; #The value of the coeff of the modifier*MBSR interaction terms.
#   gama.mbsr = rep(gama.mbsr.v, mdfN.mbsr); #Vector of coeffs of the modifier*MBSR interaction terms.
# 
#   names.mdf.mbsr = paste("zmdfMBSR", 1:mdfN.mbsr, sep=""); #variable names of the modifiers for MBSR
# 
# ###==== Generate data ====
#     x0 = rep(1, groupSampleSize*3); #Intercept
#     
#     #The mixed-effect term for racial group 
#     xrace = sample(1:raceN, size=groupSampleSize*3, replace=T, prob=rep(1/raceN, raceN));
#       #Assume equal chance for each racial group to be sampled. 
#     b.xrace = array(NA, dim=groupSampleSize*3); #b.xrace is the vector of b_i values. 
#     for (racei in 1:raceN){ b.xrace[which(xrace==racei)] = rnorm(1, sd=sigma.race); }
#       #assign the same b_i value for the all in the ith racial group. 
#     
#     xacup = c(rep(0, groupSampleSize*2), rep(1, groupSampleSize)); #acupuncture group indicator
#     xmbsr = c(rep(0, groupSampleSize), rep(1, groupSampleSize), rep(0, groupSampleSize)); #MBSR group indicator
#     xsex = rbinom(n=groupSampleSize*3, size=1, prob=0.5); #50% recruited are males??
#     xage = rnorm(n=groupSampleSize*3, mean=0, sd=1); #standardized age. 
#     
#     #PainMarker data
#     xbiom = data.frame(matrix(rnorm(n=groupSampleSize*3*markerN), ncol=markerN)); #Assume biomarker values are N(0,1);
#     names(xbiom) = names.xbiom; 
#     
#     #acupuncture-modifier data
#     zmdfAcup = data.frame(matrix(rnorm(n=groupSampleSize*3*mdfN.acup), ncol=mdfN.acup)); #Assume acupuncture-modifier values are N(0,1);
#     names(zmdfAcup) = names.mdf.acup;
# 
#     #mbsr-modifier data
#     zmdfMBSR = data.frame(matrix(rnorm(n=groupSampleSize*3*mdfN.mbsr), ncol=mdfN.mbsr)); #Assume acupuncture-modifier values are N(0,1);
#     names(zmdfMBSR) = names.mdf.mbsr;
# 
#     Xmatrix = cbind(x0, b.xrace, xacup, xmbsr, xsex, xage, xbiom, zmdfAcup, zmdfMBSR);
#   
#     ###====Specify true model =====
#     names.mainEff = c("x0", "b.xrace", "xacup", "xmbsr", "xsex", "xage", names.xbiom, names.mdf.acup, names.mdf.mbsr);
#     coeffs.mainEff = c(beta0, beta.race, beta.acup, beta.mbsr, beta.sex, beta.age, beta.biom, beta.mdf.acup, beta.mdf.mbsr);
# 
#     names.trt = c("xacup", "xmbsr"); #variable names of the treatments
#     names.mdf = list(names.mdf.acup, names.mdf.mbsr); #variable names of the modifiers corresponding to the treatments. Require: Each element is a verctor of modifier names corresponding to the treatment in names.trt.
#     coeffs.interaction = list(gama.acup, gama.mbsr); #coefficients of the treatment-modifier interaction terms. Require: Consistent with names.trt and names.mdf.
#     
#     ###===Generte the response ====
#     out = get.Y.reg(XData=Xmatrix, names.mainEff, coeffs.mainEff, names.trt, names.mdf, coeffs.interaction);
#   
#################################
get.Y.reg = function (XData, names.mainEff, coeffs.mainEff, names.trt=NULL, names.mdf=NULL, coeffs.interaction=NULL, errSD=1){
  #Main-effects   
  Xbeta = as.matrix(XData[, names.mainEff])%*%coeffs.mainEff;
  
  #Interaction-effects
  if (!is.null(names.trt)) {
    for (trti in 1:length(names.trt)){
      xtrt = names.trt[trti];
      xmdf = names.mdf[[trti]];
      coeffInter = coeffs.interaction[[trti]];
      for (mdfj in 1:length(xmdf)){
        Xbeta = Xbeta + XData[, xtrt]*XData[, xmdf[mdfj]]*coeffInter[mdfj];
      }
    }
  }
  
  Epsilon = rnorm(dim(XData)[1], mean=0, sd=errSD);
  Y = Xbeta + Epsilon;
  
  f2 = var(Xbeta)/errSD^2; #Estimated Cohen's f^2
  R2 = f2/(1+f2); #coefficient of determination. 
  return(list(Y=Y, f2=f2, R2=R2));
}



##################################
# Function for generating binary response based on logit model that contains main effects and interaction effects. 
#  The interaction terms are based on treatment variables and their modifiers. 
# Arguments:
#     XData: Data frame of independent variables with column names being the names of the corresponding variables.
#     names.mainEff: Vector of variable names of the main-effect variables.
#     coeffs.mainEff: Vector of coefficients of the main-effect variables.
#     names.trt: Vector of variable names of the treatments. No interactions if names.trt is NULL. 
#     names.mdf: List of variable names of the modifiers corresponding to the treatments. Require: Each list element is a verctor of modifier names corresponding to the treatment in names.trt.
#     coeffs.interaction = list(gama.acup, gama.mbsr); #coefficients of the treatment-modifier interaction terms. Require: Consistent with names.trt and names.mdf.
# Output: Response vector Y. 
# Author: ZWu
# Example: 
#################################
get.Y.logit = function (XData, names.mainEff, coeffs.mainEff, names.trt=NULL, names.mdf=NULL, coeffs.interaction=NULL){
  #Main-effects   
  Xbeta = as.matrix(XData[, names.mainEff])%*%coeffs.mainEff;
  
  #Interaction-effects
  if (!is.null(names.trt)) {
    for (trti in 1:length(names.trt)){
      xtrt = names.trt[trti];
      xmdf = names.mdf[[trti]];
      coeffInter = coeffs.interaction[[trti]];
      for (mdfj in 1:length(xmdf)){
        Xbeta = Xbeta + XData[, xtrt]*XData[, xmdf[mdfj]]*coeffInter[mdfj];
      }
    }
  }
  
  Y = array(NA, length(Xbeta));
  for (i in 1:length(Y)) {
    Py1x = 1/(1 + exp(-Xbeta[i])); # P(y=1 | x) based on logistic model
    Y[i] = ifelse(runif(1) < Py1x, 1, 0);
  }
  #mean(Y); #prevalance of Y. 
  
  return(list(Y=Y));
}