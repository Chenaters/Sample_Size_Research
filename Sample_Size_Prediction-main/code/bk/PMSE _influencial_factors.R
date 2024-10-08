####################
# Get the percentage PMSE reduction (pPMSEr) based on influencial factor: 
#     sample size n; 
#     numbers of predictors in basic/reduced and full model
#     effect sizes
#     covariance among predictors
# Draw the functional curve of pPMSEr over sample size at given other influencial factors (or "patterns" of them). 
#     Is there a "shoulder point", where increasing sample size does not improve pPMSEr much? 
####################

source("PMSE_functions.R");
source("/Users/zheyangwu/ResearchDoc/Computation/CodingLibraries/myRLibrary/Visualizations/Func_Figs_Power_Size.R"); #For plotting


###
#Study how sample size change the percentage PMSE reduction (pPMSEr), given the number of new predictors 
#   and their effects in terms of Cohen's f2 .
###

N = seq(50, 300, by=10); #Sample sizes
F2 = seq(0.1, 0.5, by=0.2); #f2 values 
P2 = seq(2, 10, by=4); #Number of new predictors
p1 = 2; #Number of basic predictors


for (f2 in F2){
  pPMSErs = NULL; #storing pPMSEr
  for (p2 in P2){
    pPMSErs = rbind(pPMSErs,  
                    unlist(lapply(N, function(x){pPMSEr(n=x, p1=p1, p2=p2, fnew2=f2)})));
  }
  rownames(pPMSErs) = paste0("p2_", P2);
  trendCurves(dat=pPMSErs, trendValue=N, 
              outpath="/Users/zheyangwu/ResearchProj/SampleSize/Sample_Size_Prediction/reports/PMSE improvement by new predictors/figures/",
              fileName=paste("Cohenf2_p1_", p1, "_f2_", f2, ".png", sep=""),
              xlab = "sample size", ylab="pPMSEr", main=paste0("f2 ", f2),ylim=c(-0.3, 0.5));
}


for (p2 in P2){
  pPMSErs = NULL; #storing pPMSEr
  for (f2 in F2){
    pPMSErs = rbind(pPMSErs,  
                    unlist(lapply(N, function(x){pPMSEr(n=x, p1=p1, p2=p2, fnew2=f2)})));
  }
  rownames(pPMSErs) = paste0("f2_", F2);
  # plotTrendCurves(dat=pPMSErs, trendValue=N,  
  #                 xlab = "sample size", ylab="pPMSEr", main=paste0("p2=", p2), 
  #                 ylim=c(-1, 1));
  trendCurves(dat=pPMSErs, trendValue=N, 
              outpath="/Users/zheyangwu/ResearchProj/SampleSize/Sample_Size_Prediction/reports/PMSE improvement by new predictors/figures/",
              fileName=paste("Cohenf2_p1_", p1, "_p2_", p2, ".png", sep=""),
              xlab = "sample size", ylab="pPMSEr", main=paste0(p2, " new predictors"),ylim=c(-0.3, 0.5));
}

