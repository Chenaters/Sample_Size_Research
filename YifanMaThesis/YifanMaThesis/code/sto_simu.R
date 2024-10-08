####################
# Simulations for varifying the prediction mean square error (PMSE) calculation 
# by \cite{narula1974predictive} "Predictive mean square error and stochastic regressor variables"
####################


#############################




library(MASS);
source("....")
####unconditional


ConPMSE(100,10,10,sigmak2 = 1,MU = rep(0,10),SIGMA = polyCor(10,0.3),alpha = 1, BETA = array(1,10)) #full model
ConPMSE(100,10,5,sigmak2 = 1,MU = rep(0,10),SIGMA = polyCor(10,0.3),alpha = 1, BETA = array(1,10)) #partial model
