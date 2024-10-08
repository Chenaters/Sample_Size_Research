	#library(randomForest);

	library(MASS)

	
#### ------------ Functions for model fitting and prediction --------------
		
	#######################################################
	#Description: Function splits a vector almost evenly into a list of nfold parts of the vector
	#Depends on: 
	#Arguments: vect: vector need to be splited. Require: length(vect) >= nfold.
	#			nfold: number of folds the vector need to be splitted 
	#Details: 
	#Value: Output a list of nfold of parts of the vector
	#Note:  
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: splitVector(1:10/10, 5);
	#         splitVector(1:10, 3);
	#######################################################
	splitVector <- function(vect, nfold) {
		N = length(vect);
        mi <- trunc(N/nfold);  # number of elments in the ith group, where i=2..nfold
        m1 <- N - mi*(nfold-1); # number of elements in the 1st group
        index=c(rep(1, m1), sapply(2:nfold, function(x) rep(x, mi)));
        return(split(vect, index));
	}



	#######################################################
	#Description: Function splits a vector s.t. values within clusters are distributed 
	#             almost evenly into a list of nfold parts. 
	#Depends on: splitVector()
	#Arguments: vect: vector need to be splited. Require: length(vect) >= nfold.
	#			nfold: number of folds the vector need to be splitted 
	#			vectClusters: vector of clusters. Require: Its order corresponds to vect
	#Details: 
	#Value: Output a list of nfold of parts of the vector
	#Note:  
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: splitVectorInCluster(1:10/10, 5, c(rep(1, 5), rep(2, 5)));
	#######################################################
	splitVectorInCluster <- function(vect, nfold, vectClusters) {
     	groups = vector("list", nfold);
        clusters = unique(vectClusters);
       	for (i in 1:length(clusters)) {
        	clustIndx = which(vectClusters ==clusters[i]); #find data index corresp to current clusterN
        	if (length(clustIndx) < nfold) { #if cannot fully split, then split as much as possible
       			cvGroups = splitVector(vect[clustIndx], length(clustIndx));
       			for (j in 1:length(clustIndx)) groups[[j]] = c(groups[[j]], cvGroups[[j]]);      		
       		} else { # if cannot fully split, then split it into two groups
        		cvGroups = splitVector(vect[clustIndx], nfold);
        		for (j in 1:nfold) groups[[j]] = c(groups[[j]], cvGroups[[j]]);
        	}
       	}
       	return(groups);
	}


	#######################################################
	#Description: Function of model selection based on prediction
	#Depends on: 
	#Arguments: Y: response variable
	#		Xs: data frame of predictor variables with columne names. Require X has at least two columnes
	#		candiVars: vector of candicate variables to be selected from. Require: either included in the colnames(Xs), or
	#				   by default "ALL" indicates all columne variables of Xs. Require: 
	#		strategy: model search strategies: "forward", "backward", "stepwise", "exhaustive" 
	#		criterion: criterion of prediction: "corr" (correlation b/w observed and predicted responses), "mse"
	#		trace: 0: no trace; 1: print current selected model; 2: print each evaluated model and its evaluation.
	#		steps: max number of steps to stop searching for strategies EXCEPT "exhaustive". 
	#		bestModelN: the max number of best models to be returned (or the # of models traced, if bestModelN is too large).
	#		evalFunc: objective evaluation function for prediction criteria
	#		...: options to be passed to evalFunc.
	#Details: If the model to be selected is Linear Mixed Effect model, then only the fixed terms are selected, and
	#	   random effect term needs to be provided explicitly through ... for evalFunc
	#Value: The best model with the best prediction result according to the criterion.
	#Note: This function is originally used for gene expression prediction project with Weng. 
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	modelSelByPred <- function(Xs, Y, candiVars="ALL", strategy="forward", criterion="corr", trace=0, steps = 1000, bestModelN=20, evalFunc, ...) {
		#Y=Y; Xs=dataX[, 1:2]; candiVars="ALL"; strategy="stepwise"; criterion="corr"; evalFunc=predEvaluCV.lme; nfold=10; trace=1; steps = 100; bestModelN=10;
		
		if (candiVars=="ALL") candiVars = colnames(Xs); 
				
		if (strategy == "stepwise") {
			model  = array(F, length(candiVars)); 	#current model indexing candiVars
			modelC = !model; 						#complementary of current model
			modelEval = -1; 						#initial value for evaluation of current model
			mi = 0; 								#indx for stored models
			flag.DirectionSwitch = F; 				#indicate the direction is just switched.
			models=c(); modelEvals=c(); 			#store models and their evaluations
			
			while(mi < steps) { #model selection
				if (sum(model)==0) {selDirection = "forward"; }
				if (sum(modelC)==0) {selDirection = "backward"; }

				if(selDirection == "forward") {
					if (trace>=1) print(selDirection);
					evals = array(-1, sum(modelC)); #store current evaluations
					modelCIndx = which(modelC);
					for (i in 1:length(modelCIndx)) {
						m = model; m[modelCIndx[i]] = T; #model to be evaluated
						fmla = as.formula(paste("Y ~", paste(candiVars[m], collapse= "+")));
						datm = data.frame(Y, Xs[, candiVars[m]]); names(datm) = c("Y", candiVars[m]);
						out = evalFunc(fmla, dat=datm, ...);
						if (criterion == "corr") { evals[i] = cor(out[[1]], out[[2]]); }
						if (criterion == "mse") {  }
						if (trace>=2) print(list(candiVars[m], evals[i]));
					}
					if (max(evals) > modelEval) {
						model[modelCIndx[which.max(evals)]]=T; 
						modelC[modelCIndx[which.max(evals)]]=F;
						modelEval = max(evals);
						models=rbind(models, model); modelEvals=c(modelEvals, modelEval);
						mi = mi+1;
						if (trace>=1) { print("Currently selected model:"); print(candiVars[model]); print(modelEval); }
						flag.DirectionSwitch = F;
					}else{
						selDirection = "backward";
						if (flag.DirectionSwitch) break; #break the while loop if direction just switched and current direction is not working
						flag.DirectionSwitch = T;
					}
				}
				if (selDirection == "backward") {
					if (trace>=1) print(selDirection);
					evals = array(-1, sum(model)); #store current evaluations
					modelIndx = which(model);
					for (i in 1:length(modelIndx)) {
						m = model; m[modelIndx[i]] = F; #model to be evaluated
						fmla = as.formula(paste("Y ~", paste(candiVars[m], collapse= "+")));
						datm = data.frame(Y, Xs[, candiVars[m]]); names(datm) = c("Y", candiVars[m]);
						out = evalFunc(fmla, dat=datm, ...);
						if (criterion == "corr") { evals[i] = cor(out[[1]], out[[2]]); }
						if (criterion == "mse") {  }
						if (trace>=2) print(list(candiVars[m], evals[i]));
					}
					if (max(evals) > modelEval) {
						model[modelIndx[which.max(evals)]]=F; 
						modelC[modelIndx[which.max(evals)]]=T;
						modelEval = max(evals);
						models=rbind(models, model); modelEvals=c(modelEvals, modelEval);
						mi = mi+1;
						if (trace>=1) { print("Currently selected model:"); print(candiVars[model]); print(modelEval); }
						flag.DirectionSwitch = F;
					}else{
						selDirection = "forward";
						if (flag.DirectionSwitch) break; #break the while loop if direction just switched and current direction is not working.
						flag.DirectionSwitch = T;
					}
				}
			}
		}
		
		if (strategy == "forward") {
			model  = array(F, length(candiVars)); 	#current model indexing candiVars
			modelC = !model; 						#complementary of current model
			modelEval = -1; 						#initial value for evaluation of current model
			mi = 0; 								#indx for stored models
			models=c(); modelEvals=c(); 			#store models and their evaluations
			
			while(mi < steps) { #model selection
				evals = array(-1, sum(modelC)); #store current evaluations
				modelCIndx = which(modelC);
				for (i in 1:length(modelCIndx)) {
					m = model; m[modelCIndx[i]] = T; #model to be evaluated
					fmla = as.formula(paste("Y ~", paste(candiVars[m], collapse= "+")));
					datm = data.frame(Y, Xs[, candiVars[m]]); names(datm) = c("Y", candiVars[m]);
					out = evalFunc(fmla, dat=datm, ...);
					if (criterion == "corr") { evals[i] = cor(out[[1]], out[[2]]); }
					if (criterion == "mse") {  }
					if (trace>=2) print(list(candiVars[m], evals[i]));
				}
				if (max(evals) > modelEval) {
					model[modelCIndx[which.max(evals)]]=T; 
					modelC[modelCIndx[which.max(evals)]]=F;
					modelEval = max(evals);
					models=rbind(models, model); modelEvals=c(modelEvals, modelEval);
					mi = mi+1;
					if (trace>=1) { print("Currently selected model:"); print(candiVars[model]); print(modelEval); }
					if (sum(modelC) == 0) break;
				}else {	break; }
			}			
		}

		if (strategy == "backward") {
			model = array(T, length(candiVars)); 	#current model indexing candiVars
			modelC = !model; 						#complementary of current model
			modelEval = -1; 						#initial value for evaluation of current model
			mi = 0; 								#indx for stored models
			models=c(); modelEvals=c(); 			#store models and their evaluations
			
			while(mi < steps) { #model selection
				evals = array(-1, sum(model)); #store current evaluations
				modelIndx = which(model);
				for (i in 1:length(modelIndx)) {
					m = model; m[modelIndx[i]] = F; #model to be evaluated
					fmla = as.formula(paste("Y ~", paste(candiVars[m], collapse= "+")));
					datm = data.frame(Y, Xs[, candiVars[m]]); names(datm) = c("Y", candiVars[m]);
					out = evalFunc(fmla, dat=datm, ...);
					if (criterion == "corr") { evals[i] = cor(out[[1]], out[[2]]); }
					if (criterion == "mse") {  }
					if (trace>=2) print(list(candiVars[m], evals[i]));
				}
				if (max(evals) > modelEval) {
					model[modelIndx[which.max(evals)]]=F; 
					modelC[modelIndx[which.max(evals)]]=T;
					modelEval = max(evals);
					models=rbind(models, model); modelEvals=c(modelEvals, modelEval);
					mi = mi+1;
					if (trace>=1) { print("Currently selected model:"); print(candiVars[model]); print(modelEval); }
					if (sum(model) == 1) break;
				}else{ break; }
			}
		}
		
		if (strategy == "exhaustive") {
			models = (expand.grid(rep(list(0:1), length(candiVars)))==1)[-1, ]; #Each row indicates predictor locations of a model
			modelN = nrow(models);
			modelEvals = array(-1, modelN);
			for (i in 1:modelN) { #consider all models
				fmla = as.formula(paste("Y ~", paste(candiVars[models[i, ]], collapse= "+")));
				datm = data.frame(Y, Xs[, candiVars[models[i, ]]]); 
				names(datm) = c("Y", candiVars[models[i, ]]);
				out = evalFunc(fmla, dat=datm, ...);
				if (criterion == "corr") { modelEvals[i] = cor(out[[1]], out[[2]]); }
				if (criterion == "mse") {  }
				if (trace>=2) print(list(candiVars[models[i, ]], modelEvals[i]));
			}
		}

		bestMIndx = order(modelEvals, decreasing=T)[1:min(length(modelEvals), bestModelN)];	
		if (length(bestMIndx)==1) bestModels = data.frame(t(models[bestMIndx, ])); 
		if (length(bestMIndx) > 1) bestModels= models[bestMIndx, ];
		bestModelEvals=modelEvals[bestMIndx];
		colnames(bestModels)=candiVars; rownames(bestModels)=1:nrow(bestModels);
		return( list(criterion=criterion, bestModels=bestModels, bestModelEvals=bestModelEvals));
	}
	
	

	#######################################################
	#Description: Exhaustively study all models for prediction, based on the given response(s) and a set of candidate predictors 
	#Depends on: predEvaluCV.glmmPQL
	#Arguments: respVar: Name of response variable. Require: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#		predVars: Vector of the names of candidate predictors for the formula the model
	#		random: a formula describing the random effect. 
	#		dat: data contains columns of respVar and predVars. For response variable: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#		predLevel: integer increasing from outermost to innermost clustering level. #0: population predictions. For model_R="lme".
	#		nfold: integer number of fold. Require nfold >= 2 to validate prediction-based evaluation
	#		nrepeat: number of times (>=1) to repeat random nfold-CV to stablize the results. Useful only if isRandomCV is TRUE. 
	#		isRandomCV: If TRUE, randomize the order of observations (but keep the association) before CV.
	#		outFile: the path and file name to output model and prediciton information. Default is NULL indicating no file output.
	#		graphPath: the path for outputing the graph file. Default is NULL indicating no graph output. 
	#		correlation: the correlation argument for lme function (used in glmmPQL function). Default is NULL indicating no within-group correlations. 
	#		bestModelN: the max number of best models to be returned (or the # of models traced, if bestModelN is too large).
	#Details:
	#Value: ???The best model with the best prediction result according to the criterion.
	#Note: 
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	allModelPred.glmmPQL <- function(respVar, predVars, random, dat, predLevel=1, nfold=5, nrepeat=1, isRandomCV=T, outFile=NULL, graphPath=NULL, correlation=NULL, bestModelN=20) {
	
		models = (expand.grid(rep(list(0:1), length(predVars)))==1)[-1, ]; #Each row indicates predictor locations of a model
		modelN = nrow(models);
		modelEvals = array(NA, modelN); #store the evaluation values for models.
		fixedTerms = array(NA, modelN); #store the fixed term for models.
		
		for (i in 1:modelN) { #consider all models
			fixedTerms[i] = paste(paste(respVar, "~"), paste(predVars[models[i, ]], collapse= "+"));
			fixed = as.formula(fixedTerms[i]); 
			print(fixed);
			
			if(is.null(graphPath)) {
				graphFile = NULL;
			}else {
				graphFile = paste(graphPath, respVar, "_", paste(predVars[models[i, ]], collapse= "_"), "_", nfold, "fold_", nrepeat, "repeat", sep="");
			}
										
			out = predEvaluCV.glmmPQL(fixed=fixed, random=random, dat=dat, predLevel=predLevel, nfold=nfold, nrepeat=nrepeat, isRandomCV= isRandomCV, graphFile=graphFile, correlation=correlation);

			write.table(t(c(fixedTerms[i], out[[1]], out[[2]], out[[3]], out[[2]]+out[[3]], out[[4]])), 
								 file = outFile, row.names=F, col.names=F, quote=F, append=T);
			
			#modelEvals[i] = out[[2]]+out[[3]]; #model evaluation based on sensitivity+specificity
			modelEvals[i] = out[[4]]; #model evaluation based on AUC
		}

		#record the best models based on the largest values in modelEvals
		bestMIndx = order(modelEvals, decreasing=T)[1:min(length(modelEvals), bestModelN)];	
		if (length(bestMIndx)==1) bestModels = data.frame(t(models[bestMIndx, ])); 
		if (length(bestMIndx) > 1) bestModels= models[bestMIndx, ];
		bestModelEvals=modelEvals[bestMIndx];
		colnames(bestModels)=predVars; rownames(bestModels)=1:nrow(bestModels);
		return( list(fixedTerms=fixedTerms[bestMIndx], bestModels=bestModels, bestModelEvals=bestModelEvals));
	}
	
	
	
	
	
	#######################################################
	#Description: Function of evaluating prediction by cross-validation through glmmPQL 
	#Depends on: glmmPQL
	#Arguments: fixed: a two-sided linear formula giving fixed-effects part of the model.
	#			random: a formula or list of formulae describing the random effect. 
	#			dat: dataframe of response and predictors. Require to include all vars in options fixed and random. For response variable: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#			predLevel: integer increasing from outermost to innermost clustering level. #0: population predictions. For model_R="lme".
	#			nfold: integer number of fold. Require nfold >= 2 to validate prediction-based evaluation
	#			nrepeat: number of times (>=1) to repeat random nfold-CV to stablize the results. Useful only if isRandomCV is TRUE. 
	#			probCutoffs: Probability cutoffs to calculate sensitivity and specificity. Either "ALL" (default), which considers all possible cutoffs, i.e., all the predicted probabilities, or a given vector of cutoffs. 
	#			isRandomCV: If TRUE, randomize the order of observations (but keep the association) before CV.
	#			isBestCutoff: If TRUE, output result corresponding to the best probability cutoff for prediction, where the summation of sensitivity+specificity is maximal. Otherwise, output all probThresholds and the corresponding sensitivities and specificities. 
	#			cvInClusterVar: NULL or a vector of variable names that indicate clusters from which we randomly select training/testing
	#			graphFile: the path and file name for the outputing graph on probability cutoffs. Default is NULL indicating no graph output. 
	#			correlation: the correlation argument for lme function (used in glmmPQL function). Default is NULL indicating no within-group correlations. 
	#			...: other arguments for glmmPQL()
	#Details: 1) If cvInClusterVar is defined, training/testing data sets for cross-validation can be selected within each of the clusters. 
	#Value: Output a list of probThresholds and the corresponding sensitivities and specificities.
	#Note: This function is originally used for gene expression prediction project with Weng. 
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	predEvaluCV.glmmPQL <- function(fixed, random, dat, predLevel=0, nfold=5, nrepeat=1, probCutoffs="ALL", isRandomCV=T, isBestCutoff=T, cvInClusterVar=NULL, graphFile=NULL, correlation=NULL, ...){		
 		require(MASS);
        Y_observed=c(); prob_predicted=c(); #To store the observed value and predicted probabilities
	    respVar = as.character(fixed[[2]]); #name of the response variable
	
		for (ir in 1:nrepeat) {
			#suffle observations for random cross-validation
	   		if (isRandomCV) { dat=dat[sample(1:nrow(dat)), ]; } 

	        # divide the dataset into n-fold
	        if (!is.null(cvInClusterVar)) { 
	        	groups=	splitVectorInCluster(1:nrow(dat), nfold, dat[, cvInClusterVar]);
	        } else { groups = splitVector(1:nrow(dat), nfold); }
	
	        # take the kth group in turn to get predicted values for the kth group
	        for(k in 1:nfold) {
	        	ind = groups[[k]];
	            testset = dat[ind, ]; 
	            trainset = dat[-ind, ]; 
	
				#gm = glmmPQL(fixed=fixed,  random=random, family=binomial, data=trainset, ...);	
				gm = glmmPQL(fixed=fixed,  random=random, family=binomial, data=trainset, correlation=correlation, niter=10, verbose=T);	
				predProb= predict(gm, testset, type="response", level= predLevel);
	           
	            Y_observed = c(Y_observed, testset[, respVar]);
	            prob_predicted = c(prob_predicted, predProb);
	        }
        }


        # ##Implmentation that can take arbitrary cutoffs to calculate sensitivities and specificities
        # if(probCutoffs=="ALL"){ probCutoffs = sort(prob_predicted); }
		# sensitivities = specificities = array(NA, length(probCutoffs));
		# for(k in 1:length(probCutoffs)) {
			# predOut = as.integer(prob_predicted >= probCutoffs[k]);
			# sensitivities[k] = sum(Y_observed==1 & predOut==1) / sum(Y_observed==1);
			# specificities[k] = sum(Y_observed==0 & predOut==0) / sum(Y_observed==0);
		# }
		# senNspe = sensitivities+ specificities;
		
		##Default of using ordered prob_predicted as cutoffs to calculate  sensitivities, specificities, and AUC
		results = getROC_AUC(Y_observed, prob_predicted);				
		probCutoffs = results$cutoffs;
		sensitivities = results$truePositiveRates;
		specificities = 1 - results$falsePositiveRates;		
		senNspe = sensitivities+ specificities;
			
		#output graph for the relationship between probCutoffs and senNspe
		if(!is.null(graphFile)) {
			pdf(paste(graphFile, "_sNp.pdf", sep="")); 
			plot(probCutoffs, senNspe, type="l");
			dev.off();
			
			pdf(paste(graphFile, "_ROC.pdf", sep="")); 
			drawROC(results$falsePositiveRates, results$truePositiveRates, results$auc);
			dev.off();
		}
			
		#the best results that has maximal sensitivity+specificity.
		if (isBestCutoff) {			
			index.best = which(senNspe==max(senNspe))[1];
			probCutoffs = probCutoffs[index.best];
			sensitivities = sensitivities[index.best]; 
			specificities = specificities[index.best]; 			
		}

		return( list(probCutoffs=probCutoffs, sensitivities=sensitivities, specificities=specificities, auc=results$auc));		
	}



	#######################################################
	#Description: Function of evaluating prediction by cross-validation through (classification-aided) regression models
	#Depends on: lm or lme (either could be aided by randomForest)
	#Arguments: fixedf: formula of the fixed effects
	#			dat: dataframe of response and predictors. Require to include all vars in options fixed and random.
	#			randomf: formula of the random effect. Required if model_R="lme". 
	#			model_R: model for regression: "lm", "lme"
	#			predLevel: integer increasing from outermost to innermost clustering level. #0: population predictions. For model_R="lme".
	#			nfold: integer number of fold. Require nfold >= 2 to validate prediction-based evaluation
	#			isRandomCV: If TRUE, randomize the order of observations (but keep the association) before CV.
	#			model_C: model for classification: NULL (no calssification), or "randomForest". Only for fixed effect covariates.
	#			cvInClusterVar: NULL or a vector of variable names that indicate clusters from which we randomly select training/testing
	#			zero: a positive cutoff to indicate a nonzero response. 
	#Details: 1) If cvInClusterVar is defined, training/testing data sets for cross-validation can be 
	#            selected within each of the clusters. 
	#		  2) If classification method (e.g., randomForest) is defined, a small value defined by "zero" 
	#        is used to categorize the training responses to be 1 or 0. Also, the regression-predicted 
	#        values are replaced by this small value defined by "zero" whenever the corresponding 
	#        classification-predicted values are 0. 
	#Value: Output the observed and predicted response values. 
	#Note: This function is originally used for gene expression prediction project with Weng. 
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	predEvaluCV.lme <- function(fixedf, dat, randomf=NULL, model_R='lm', predLevel=0, nfold=2, isRandomCV=FALSE, model_C='none', cvInClusterVar=NULL, zero=0.01){
		#fixedf=as.formula("Y~X"); dat=dat; randomf=as.formula("~ 1|gIndx"); model_C='none'; model_R='lme';predLevel=1; isRandomCV=T; cvInClusterVar="gIndx"; nfold=20; zero=0;
		#fixedf=as.formula("Y~X+gIndx"); dat=dat; randomf=NULL; model_C='randomForest'; model_R='lm';predLevel=0; isRandomCV=F; cvInClusterVar=NULL; nfold=2; zero=0;

		#fixedf=fixedf; dat=dataX; randomf=as.formula("~ 1|cluster"); model_R='lme'; predLevel=1; nfold=5; isRandomCV=T; model_C='none'; cvInClusterVar=NULL; zero=0;
	
   		if (isRandomCV) { dat=dat[sample(1:nrow(dat)), ]; } #rownames(dat)=1:nrow(dat); }
    	respVar = as.character(fixedf[[2]]); #name of the response variable

        Y_observed=c(); Y_predicted=c(); Y_prediction_C=c();

        # divide the dataset into n-fold, randomly
        if (!is.null(cvInClusterVar)) { 
        	groups = vector("list", nfold);
         	clusters = unique(dat[, cvInClusterVar]);
        	for (i in 1:length(clusters)) {
        		clustIndx = which(dat[, cvInClusterVar] ==clusters[i]); #find data index corresp to current clusterN
        		if (length(clustIndx) < nfold) { #if cannot fully split, then split as much as possible
         			cvGroups = splitVector(clustIndx, length(clustIndx));
        			for (j in 1:length(clustIndx)) groups[[j]] = c(groups[[j]], cvGroups[[j]]);      		
        		} else { # if cannot fully split, then split it into two groups
         			cvGroups = splitVector(clustIndx, nfold);
        			for (j in 1:nfold) groups[[j]] = c(groups[[j]], cvGroups[[j]]);
        		}
        	}
        } else { groups = splitVector(1:nrow(dat), nfold); }

        # take the kth group in turn to get predicted values for the kth group
        for(k in 1:nfold) {
        	ind = groups[[k]];
          testset = dat[ind, ]; 
          trainset = dat[-ind, ]; 

          # classification
          prediction_C = 1;
          if (model_C == 'randomForest') {
          	require(randomForest);
           	trainsetC = trainset;
           	trainsetC[, respVar] = as.factor(ifelse(abs(trainsetC[, respVar]) > zero, 1, 0));

         		FOREST_model = randomForest(fixedf, data=trainsetC, ntree=100, importance = TRUE);
     				prediction_C = predict(FOREST_model, testset, type='response');  
            if(is.factor(prediction_C)) prediction_C = as.numeric(as.character(prediction_C));  # convert factor to numeric
			    }
            
            # regression
            if (model_R=="lm") {
            	LINEAR_model <- lm(fixedf, data=trainset, na.action = na.omit);
            	prediction_R <- as.vector(predict(LINEAR_model, type="response", testset, na.action=na.omit));
            }
            if (model_R=="lme") {
				      #require(nlme);
           		LINEAR_model <- lme(fixed=fixedf, random = randomf, data = trainset);
           		LINEAR_model$call$fixed <- fixedf; #To avoid R complaining that object "fixed" not found. See https://stackoverflow.com/questions/63680096/using-formula-to-predict-inside-of-r-function-generates-object-not-found-error
            	prediction_R = as.vector(predict(LINEAR_model, newdata=testset, level=predLevel, na.action=na.omit));
          		if (any(is.na(prediction_R))) { #if obs cannot be pred at high predLevel, replace them by level=0.
           			naIndx = which(is.na(prediction_R));
           			prediction_R0 = as.vector(predict(LINEAR_model, newdata=testset, level=0, na.action=na.omit));
          			prediction_R[naIndx] = prediction_R0[naIndx];
          		}
            }
           
            # final prediction combining classification and regression
            Y_observed = c(Y_observed, testset[, respVar]);
            Y_predicted = c(Y_predicted, prediction_C * prediction_R);
            Y_prediction_C=c(Y_prediction_C, prediction_C);
        }
        Y_predicted[Y_prediction_C==0]=zero;
        return(list(Y_observed, Y_predicted));
	}


	
	#######################################################
	#Description: Function of getting mean prediction measures based on several loops of Cross-Validation
	#Depends on: lme function in library nlme.
	#Arguments: 
  #   	fixed: formula of the fixed effects
	#			dat: dataframe of response and predictors. Require to include all vars in options fixed and random.
	#			randomf: formula of the random effect. Required if model_R="lme". 
	#			model_R: model for regression: "lm", "lme"
	#			predLevel: integer increasing from outermost to innermost clustering level. #0: population predictions. For model_R="lme".
	#			loopn: number of loops of CV
	#			cvNumber: number of CV folds
	#Details:
	#Values: prediction accuracy measures based on multiple CVs and loops
	#   average MSE, 
	#   L2normRatio (or Relative Absolute Error: https://www.statisticshowto.com/relative-absolute-error/) 
	#   L1normRatio 
	#   correlation b/w the observed and the predicted values. 
	#   MSEoverObsVar: prediction MSE over variance of observed value (sample mean is a naive prediction)
	#Note:  
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	meanPredEvaluCV.lme <- function (fixed, dat, randomf=NULL, model_R="lm", predLevel=0, loopn=1, cvNumber=10) {
	  
	  respVar = as.character(fixed[[2]]); #name of the response variable
	  subIndices = 1:dim(dat)[1]; #indices of the subjects
	  
		MSEs=L1normRatios=L2normRatios=correlations= MSEoverObsVar=array(0, c(loopn, cvNumber)); #Store predi criteria.
		for(i3 in 1:loopn) {
			testIndices = matrix(sample(subIndices, size=(length(subIndices) %/% cvNumber)*cvNumber, replace=FALSE), ncol= cvNumber); 
				#Randomly group subjects according to cvNumber. At most cvNumber genes are not in test set
			for (i4 in 1: cvNumber) {
				dat.test =  dat[testIndices[, i4], ]; #Split data into testing and training sets
				dat.train = dat[setdiff(subIndices, testIndices[, i4]), ];
			
				###Observed response values
				observed = dat.test[, respVar]; 

				###Predicted response values
				if(model_R=="lm"){
				  gm = lm(fixed, data = dat.train, na.action = na.omit); #model fitting using training data
				  predOut = predict(gm, newdata=dat.test, na.action=na.omit); #predict outcomes
				} else if(model_R=="lme"){
				  #gm <- do.call("lme", list(fixed=fixed, random=randomf, data=quote(dat.train)));
				  gm <- lme(fixed=fixed, random = randomf, data = dat.train);
				  gm$call$fixed <- fixed; #To avoid R complaining that object "fixed" not found. See https://stackoverflow.com/questions/63680096/using-formula-to-predict-inside-of-r-function-generates-object-not-found-error
				  predOut = as.vector(predict(gm, newdata=dat.test, level=predLevel, na.action=na.omit));

				  if (any(is.na(predOut))) { #if obs cannot be pred at high predLevel, replace them by level=0.
				    naIndx = which(is.na(predOut));
				    predOut0 = as.vector(predict(gm, newdata=dat.test, level=0, na.action=na.omit));
				    predOut[naIndx] = predOut0[naIndx];
				  }
				} else return("model_R not supported")

		    ###Prediction accuracy measures
				MSEs[i3, i4]         = mean((observed - predOut)^2, na.rm=TRUE); 
				L2normRatios[i3, i4] = sqrt(sum((observed - predOut)^2, na.rm=TRUE)/sum(observed^2, na.rm=TRUE));				
				L1normRatios[i3, i4] = sum(abs(observed - predOut), na.rm=TRUE)/sum(abs(observed), na.rm=TRUE);
				correlations[i3, i4] = cor(observed, predOut, use="pairwise.complete.obs");
				MSEoverObsVar[i3, i4]= MSEs[i3, i4]/var(observed, na.rm=TRUE);
				#MSEoverObsVar[i3, i4]= MSEs[i3, i4]/mean((observed-mean(observed, na.rm=T))^2, na.rm=TRUE);
			}
		}
		MSE = mean(MSEs); 
		L2normRatio = mean(L2normRatios); 
		L1normRatio=mean(L1normRatios); 
		correlation=mean(correlations); 
		MSEoverObsVar=mean(MSEoverObsVar);
		return(c(MSE, L2normRatio, L1normRatio, correlation, MSEoverObsVar));
	}
	

	
	
	
	
	##########===========================
	
	
	
	
	
	
	
	#######################################################
	#Description: Exhaustively study all glm models for prediction, based on the given response(s) and a set of candidate predictors 
	#Depends on: predEvaluCV.glm
	#Arguments: respVar: Name of response variable. Require: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#		predVars: Vector of the names of candidate predictors 
	#		dat: data contains columns of respVar and predVars. For response variable: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#		nfold: integer number of fold. Require nfold >= 2 to validate prediction-based evaluation
	#		nrepeat: number of times (>=1) to repeat random nfold-CV to stablize the results. Useful only if isRandomCV is TRUE. 
	#		isRandomCV: If TRUE, randomize the order of observations (but keep the association) before CV.
	#		outFile: the path and file name to output model and prediciton information. Default is NULL indicating no file output.
	#		graphPath: the path for outputing the graph file. Default is NULL indicating no graph output. 
	#		bestModelN: the max number of best models to be returned (or the # of models traced, if bestModelN is too large).
	#		modelSize: only search models with the size <= modelSize. Default is the # of covariates in predVars.
	#		confidenceLevel: A value between 0 and 1 for the confidence interval of prediction values, default is NULL meaning no confidence calculation
	#		simuN: Number of simulations to get the confidence interval. Default is 100. 
	#		isBootstrap: Whether we use bootstrap to estimate the confidence interval. Default is FALSE. 
	#Details:
	#Value: 1) If outFile is not NULL, then output a corresponding file, where each row is model with: 
	#       model fomula, the probability probCutoff that which sensitivity+specificity is maximized, sensitivity, 
	#       specificity, AUC. And, if confidenceLevel is not NULL, the corresponding AUC mean and confidence interval
	#       obtained by repeating the original data or the bootstrapped data. 
	#       2) If graphPath is not NULL, then for each model output graphes of sensitivity+specificity curve and ROC curve 
	#       over a series of probability cutoffs.
	#       3) Return bestModelN of models with the highest possible sensitivity+specificity. 
	#Note: This function is similar as allModelPred.glmmPQL, but for independent samples. 
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	allModelPred.glm <- function(respVar, predVars, dat, nfold=5, nrepeat=1, isRandomCV=T, outFile=NULL, 
	                             graphPath=NULL, bestModelN=20, modelSize=length(predVars), 
	                             confidenceLevel=NULL, simuN=100, isBootstrap=FALSE) {
	
		models = (expand.grid(rep(list(0:1), length(predVars)))==1)[-1, ]; #Each row indicates predictor locations of a model
		models = models[apply(models, 1, sum)<= modelSize, ]; #Restrict searching to model size. 
		
		modelN = nrow(models);
		modelEvals = array(NA, modelN); #store the evaluation values for models.
		formulaTerms = array(NA, modelN); #store the formula term for models.
		
		for (i in 1:modelN) { #consider all models
			formulaTerms[i] = paste(paste(respVar, "~"), paste(predVars[models[i, ]], collapse= "+"));
			formula = as.formula(formulaTerms[i]); 
			print(formula);
			
			if(is.null(graphPath)) {
				graphFile = NULL;
			}else {
				graphFile = paste(graphPath, respVar, "_", paste(predVars[models[i, ]], collapse= "_"), "_", nfold, "fold_", nrepeat, "repeat", sep="");
			}
			
			#Calculate sensitivity, specificity, and AUC
			out = predEvaluCV.glm(formula=formula, dat=dat, nfold=nfold, nrepeat=nrepeat, isRandomCV= isRandomCV, graphFile=graphFile);
			
			if (!is.null(confidenceLevel)){ #Calculate mean and C.I. for AUC. 
				auc = predEvaluCV.glm.CI(formula=formula, dat=dat, respVar=respVar, nfold=nfold, nrepeat=nrepeat, isRandomCV=T, 
				                         cvInClusterVar= respVar, confidenceLevel=confidenceLevel, simuN=simuN, isBootstrap=isBootstrap);
				outputValue = t(c(formulaTerms[i], out[[1]], out[[2]], out[[3]], out[[2]]+out[[3]], out[[4]], auc$meanAUC, auc$CI));
			}else {
				outputValue = t(c(formulaTerms[i], out[[1]], out[[2]], out[[3]], out[[2]]+out[[3]], out[[4]]));
			}

			write.table(outputValue, file = outFile, row.names=F, col.names=F, quote=F, append=T);
			
			modelEvals[i] = out[[2]]+out[[3]];
		}

		#record the best models based on the largest values in modelEvals, i.e., sensitivity+specificity
		bestMIndx = order(modelEvals, decreasing=T)[1:min(length(modelEvals), bestModelN)];	
		if (length(bestMIndx)==1) bestModels = data.frame(t(models[bestMIndx, ])); 
		if (length(bestMIndx) > 1) bestModels= models[bestMIndx, ];
		bestModelEvals=modelEvals[bestMIndx];
		colnames(bestModels)=predVars; rownames(bestModels)=1:nrow(bestModels);
		return( list(formulaTerms=formulaTerms[bestMIndx], bestModels=bestModels, bestModelEvals=bestModelEvals));
	}
	

	#######################################################
	#Description: Function of evaluating prediction by cross-validation through glm. 
	#             Give confidence interval of the prediction accuracy
	#Depends on: glm
	#Arguments: formula: a two-sided linear formula of the model.
	#			dat: dataframe of response and predictors. Require to include all vars in formula. For response variable: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#			respVar: Name of response variable. Require: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#			nfold: integer number of fold. Require nfold >= 2 to validate prediction-based evaluation
	#			nrepeat: number of times (>=1) to repeat random nfold-CV to stablize the results. Useful only if isRandomCV is TRUE. 
	#			isRandomCV: If TRUE, randomize the order of observations (but keep the association) before CV.
	#			cvInClusterVar: NULL or a vector of variable names that indicate clusters from which we randomly select training/testing, usually it could be respVar. 
	#			confidenceLevel: A value between 0 and 1 for the confidence interval of AUC. Default 0.95. 
	#			simuN: Number of simulations in estimating C.I. Default 100. 
	#			isBootstrap: TRUE of FALSE, whether we simulate data based on empirical distributions (by bootstrap) of the two groups, respectively. Default is FALSE, meaning it repeats predictions using original data. 
	#			...: other arguments. 
	#Details: 1) If cvInClusterVar is defined, training/testing data sets for cross-validation can be selected within each of the clusters. 
	#Value: Output the mean and C.I. for AUC.
	#Note: .  
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	predEvaluCV.glm.CI <- function(formula, dat, respVar, nfold=5, nrepeat=1, isRandomCV=T, cvInClusterVar=NULL, 
	                               confidenceLevel=0.95, simuN=100, isBootstrap=FALSE, ...){		 
    
	  Y_observed=c(); prob_predicted=c(); #To store the observed value and predicted probabilities	
         
		AUCs = array(NA, simuN); 
		for (si in 1:simuN) {
			
			### Simulate data: Based on empirical distributions (by bootstrap) of the two groups, respectively
			if (isBootstrap) {
				caseIndx = which(dat[, respVar]==1);
				controlIndx = which(dat[, respVar]==0);
				caseSamp = sample(caseIndx, length(caseIndx), replace=TRUE);
				controlSamp = sample(controlIndx, length(controlIndx), replace=TRUE);
				simuDat = dat[c(caseSamp, controlSamp), ];	
			}else {
				simuDat = dat; 
			}		

			### Prediction by cross-validation
			for (ir in 1:nrepeat) {
				#suffle observations for random cross-validation
		   		if (isRandomCV) { simuDat=simuDat[sample(1:nrow(simuDat)), ]; } 
	
		        # divide the dataset into n-fold
		        if (!is.null(cvInClusterVar)) { 
		        	groups=	splitVectorInCluster(1:nrow(simuDat), nfold, simuDat[, cvInClusterVar]);
		        } else { groups = splitVector(1:nrow(simuDat), nfold); }
		
		        # take the kth group in turn to get predicted values for the kth group
		        for(k in 1:nfold) {
		        	ind = groups[[k]];
		            testset = simuDat[ind, ]; 
		            trainset = simuDat[-ind, ]; 
		
					gm = glm(formula= formula, family=binomial, data=trainset);	
					predProb= predict(gm, testset, type="response");
		           
		            Y_observed = c(Y_observed, testset[, respVar]);
		            prob_predicted = c(prob_predicted, predProb);
		        }
	        }
			
			results = getROC_AUC(Y_observed, prob_predicted);	
	        AUCs[si] = results$auc; 
		}
					
		return(list(meanAUC=mean(AUCs), CI=quantile(AUCs, c((1-confidenceLevel)/2, (1+confidenceLevel)/2))));
	}

	
	
	
	
	#######################################################
	#Description: Function of evaluating prediction by cross-validation through glm 
	#Depends on: glm
	#Arguments: formula: a two-sided linear formula of the model.
	#			dat: dataframe of response and predictors. Require to include all vars in formula. For response variable: 1 indicates positive (for sensitivity); 0 indicates negative (for specificity). 
	#			nfold: integer number of fold. Require nfold >= 2 to validate prediction-based evaluation
	#			nrepeat: number of times (>=1) to repeat random nfold-CV to stablize the results. Useful only if isRandomCV is TRUE. 
	#			isRandomCV: If TRUE, randomize the order of observations (but keep the association) before CV.
	#			isBestCutoff: If TRUE, output result corresponding to the best probability cutoff for prediction, where the summation of sensitivity+specificity is maximal. Otherwise, output all probThresholds and the corresponding sensitivities and specificities. 
	#			cvInClusterVar: NULL or a vector of variable names that indicate clusters from which we randomly select training/testing
	#			graphFile: the path and file name for the outputing graph on probability cutoffs. Default is NULL indicating no graph output. 
	#			...: other arguments for glmmPQL()
	#Details: 1) If cvInClusterVar is defined, training/testing data sets for cross-validation can be selected within each of the clusters. 
	#Value: A list of probCutoffs, sensitivities, specificities, AUC. 
	#Note: This function is similar as predEvaluCV.glm but for independent samples.  
	#Author(s): ZWu
	#References:
	#See Also:
	#Example: 
	#######################################################
	predEvaluCV.glm <- function(formula, dat, nfold=5, nrepeat=1, isRandomCV=T, isBestCutoff=T, cvInClusterVar=NULL, 
	                            graphFile=NULL, ...){		
 
    Y_observed=c(); prob_predicted=c(); #To store the observed value and predicted probabilities
	  respVar = as.character(formula[[2]]); #name of the response variable
	    
	  cvInClusterVar= respVar;
	
		for (ir in 1:nrepeat) {
			#suffle observations for random cross-validation
	   	if (isRandomCV) { dat=dat[sample(1:nrow(dat)), ]; } 

	    # divide the dataset into n-fold
	    if (!is.null(cvInClusterVar)) { 
	        	groups=	splitVectorInCluster(1:nrow(dat), nfold, dat[, cvInClusterVar]);
	        } else { groups = splitVector(1:nrow(dat), nfold); }
	
	    # take the kth group in turn to get predicted values for the kth group
	    for(k in 1:nfold) {
	      ind = groups[[k]];
	      testset = dat[ind, ]; 
	      trainset = dat[-ind, ]; 
	      
				gm = glm(formula= formula, family=binomial, data=trainset);	
				predProb= predict(gm, testset, type="response");
	           
	      Y_observed = c(Y_observed, testset[, respVar]);
	      prob_predicted = c(prob_predicted, predProb);
      }
    }
        
		results = getROC_AUC(Y_observed, prob_predicted);	
			
		probCutoffs = results$cutoffs;
		sensitivities = results$truePositiveRates;
		specificities = 1 - results$falsePositiveRates;		
		senNspe = sensitivities+ specificities;
			
		#output graph for the relationship between probCutoffs and senNspe
		if(!is.null(graphFile)) {
			pdf(paste(graphFile, "_sNp.pdf", sep="")); 
			plot(probCutoffs, senNspe); #, type="l");
			dev.off();
			
			pdf(paste(graphFile, "_ROC.pdf", sep="")); 
			drawROC(results$falsePositiveRates, results$truePositiveRates, results$auc);
			dev.off();
		}
			
		#the best results that has maximal sensitivity+specificity.
		if (isBestCutoff) {			
			index.best = which(senNspe==max(senNspe))[1];
			probCutoffs = probCutoffs[index.best];
			sensitivities = sensitivities[index.best]; 
			specificities = specificities[index.best]; 			
		}

		return( list(probCutoffs=probCutoffs, sensitivities=sensitivities, specificities=specificities, auc=results$auc));
	}


	#######################################################
	#Description: Function of calculating true positive rate and false postive rate based on true outcomes and the corresponding probabilities 
	#Depends on:
	#Arguments: trueOutcomes: Vector of true outcomes. Require: 1 indicates true positive outcomes; 0 indicate true negative outcomes.
	#			probs: For each true outcome, the probability of the true positive outcome. 
	#Details:  
	#Value: probability cutoffs; false positive rates; true positive rates; area under ROC curve.
	#Note:  
	#Author(s): ZWu (based on http://stackoverflow.com/questions/4903092/calculate-auc-in-r)
	#References:
	#See Also:
	#Example: 
		# true_Y = c(1,1,1,1,0,1,0,1,0,0);
	   # probs = c(1,0.999,0.999,0.973,0.568,0.421,0.382,0.377,0.146,0.11);
		# getROC_AUC(true_Y, probs);
	#######################################################
	getROC_AUC <- function(trueOutcomes, probs){
	    probsSort = sort(probs, decreasing = TRUE, index.return = TRUE);
	    cutoffs = unlist(probsSort$x); #cutoffs are the sorted probabilies. 
	    idx = unlist(probsSort$ix);
	
	    roc_y = trueOutcomes[idx]; #true outcomes sorted by cutoffs. 
	    truePositiveRates = cumsum(roc_y == 1)/sum(roc_y == 1);    #Above a specific cutoff, 1s are true positives. 
	    falsePositiveRates = cumsum(roc_y == 0)/sum(roc_y == 0); #Above a specific cutoff, 0s are false positives. 
	
	    auc = sum((falsePositiveRates[2:length(roc_y)]-falsePositiveRates[1:(length(roc_y)-1)])*truePositiveRates[2:length(roc_y)]);
	    
	    return(list(cutoffs=cutoffs, falsePositiveRates=falsePositiveRates, truePositiveRates=truePositiveRates, auc=auc))
	}
	
	
	
	#######################################################
	#Description: Function of drawing ROC curve 
	#Depends on: 
	#Arguments: falsePositiveRates: The x-axis values (i.e., 1 - specificities). 
	#			truePositiveRates: The y-axis values (i.e., sensitivities). 
	#			auc: Area under the ROC curve. 
	#Details:  
	#Value: 
	#Note:  
	#Author(s): ZWu (based on http://stackoverflow.com/questions/4903092/calculate-auc-in-r)
	#References:
	#See Also:
	#Example: 		
		# true_Y = c(1,1,1,1,0,1,0,1,0,0);
	   # probs = c(1,0.999,0.999,0.973,0.568,0.421,0.382,0.377,0.146,0.11);
		# results=getROC_AUC(true_Y, probs);
		# drawROC(results$falsePositiveRates, results$truePositiveRates, results$auc);
	#######################################################
	drawROC <- function(falsePositiveRates, truePositiveRates, auc=NULL){
		plot(falsePositiveRates, truePositiveRates, type = "l", col = "blue", xlab = "False Positive Rate", ylab = "True Positive Rate", main = "ROC");
		axis(1, seq(0.0,1.0,0.1));
		axis(2, seq(0.0,1.0,0.1));
		abline(h=seq(0.0,1.0,0.1), v=seq(0.0,1.0,0.1), col="gray", lty=3);
		if (!is.null(auc)) {
			legend(0.7, 0.3, sprintf("%3.3f",auc), lty=c(1,1), lwd=c(2.5,2.5), col="blue", title = "AUC");
		}
	}

	 