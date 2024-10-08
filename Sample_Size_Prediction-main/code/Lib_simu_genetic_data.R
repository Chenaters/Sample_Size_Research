######################
# This libarary collect functions that generate genetic data, including:
#     Genotype data
#     Response traits
######################

###-------Dependent libraries-------

#library(mvtBinaryEP); #package ‘mvtBinaryEP’ is not available (for R version 3.3.3)
library(Rlab);

###-------Useful functions from external sources-------


###-------Main functions-------

  
  # Problem: mvtBinaryEP is not available for the lastest R version. 
  #############################################################################
  # Description: To generate a matrix of correlated genotypes containing values 0, 1, 2 (copy number of minor allele.);
  # Depends on: function ep in libarary mvtBinaryEP, and self-implemented function genIndpBX
  # Arguments: mafs: A vector of minor allele frequences
  #		         Corr: The correlation coefficient matrix among the genotypes
  #		         numb: Number of subjects (i.e., sample size) 
  # Detail:  We first generate iid random vectors of binary values under the given correlation coefficients. 
  #          Genotype are obtained by adding two adjacent independent random vector observations, which gives 
  #          the same correlation coefficients. 
  # Value: A matrix of genotype data, where rows correspond to subjects; columns correspond to markers. 	
  # Note: It inherits all the limits from function ep. E.g., some extreme correlation and MAF cannot be reliably 
  #       realized (refer to the ep function). 
  # Author(s): SHe
  # Reference: Wu et al. 2013 AOAS
  # See also:
  # Example: 
  #############################################################################
  #gen_bx <- function(mafs, Corr, numb) { #original name
  gen_corr_genotype <- function(mafs, Corr, numb) {
    if( sum(abs(Corr-diag(nrow(Corr))))!=0 ){	# if the correlation matrix is not a diagonal matrix
      da <- ep(mu=mafs,R=Corr,nRep= (numb*2),seed=NULL)$y
      x  <- da[seq(2,(numb*2),by=2),]+da[seq(1,(2*numb-1),by=2),]
    } else {
      x <- genIndpBX(mafs, numb)
    }	
    return (x)
  }
  
  gen_corr_genotype <- function(mafs, Corr, numb, normal=FALSE) {
    if( sum(abs(Corr-diag(nrow(Corr))))!=0 ){	# if the correlation matrix is not a diagonal matrix
      if(normal){
        da <- rmvbin(n= (numb*2),margprob=mafs,sigma =Corr)
      }else{
        da <- rmvbin(n= (numb*2),margprob=mafs,bincorr =Corr)
      }
      
      x  <- da[seq(2,(numb*2),by=2),]+da[seq(1,(2*numb-1),by=2),]
    } else {
      x <- gen_indp_genotype(mafs, numb)
    }	
    return (x)
  }
  
  
  #############################################################################
  # Description: To generate a matrix of independent genotypes containing values 0, 1, 2 (copy number of minor allele. );
  # Depends on: function rbern in package Rlab 
  # Arguments: mafs: A vector of minor allele frequences
  #		     numb: Number of subjects (i.e., sample size) 
  # Detail:  We first generate iid random vectors of binary values. Genotype are obtained by adding two adjacent 
  #          independent random vector observations. 
  # Value: A matrix of genotype data, where rows correspond to subjects; columns correspond to markers. 	
  # Note: This is much faster than function gen_bx for generating independent genotypes. 
  # Author(s): SHe
  # Reference: Wu et al. 2013 AOAS
  # See also:
  # Example: 
  #############################################################################
  #genIndpBX<- function(mafs, numb){
  gen_indp_genotype <- function(mafs, numb){
    da <- matrix(unlist(lapply(mafs, rbern, n=numb*2)), nrow=numb*2, ncol=length(mafs), byrow= F)
    x  <- da[seq(2,(numb*2),by=2),]+da[seq(1,(2*numb-1),by=2),]
    return (x)
  }

  
  #############################################################################
  # Description: Generate a vector of response variable y based on Gaussian regression model. 
  # Depends on: 
  # Arguments: 
  #   X         : a matrix of data set where row represents number of subjects and column represents predictors.
  #		coefs     : an array of predictors'coefficients. For the true signals the coefficients are non-zero, otherwise 0.
  #		errorStd  : standard deviation of the random error. 
  # Details: 
  # Value: A random vector of quantitative responses (containing the Gaussian random error term). 
  # Note: Here errorStd should be the same as the argument errorStd used in function genBeta_WS, whenever coefs is generated from fuction genBeta_WS.  
  # Author: SHe and Wu
  # Reference(s): Wu et al. AOAS 2013
  # See also: 
  # Example: 
  #############################################################################
  #genYMarks <- function(X, coefs, errorStd=1) { #orginal name
  gen_Y_reg <- function(X, coefs, errorStd=1) {
    epsilon <- rnorm(dim(X)[1], sd=errorStd); #random errors with the same number of subjects
    Y	     <- X%*%coefs;
    Y	     <- Y + epsilon;
    return (Y);
  }
  
  
  
  
  
  #############################################################################
  # Description: Generate a vector of binary response y based on logistic regression model. 
  # Depends on: 
  # Arguments: X   : a matrix of data set where row represents number of subjects and column represents predictors.
  #		         beta: an vector of predictors'coefficients. For the true signals the coefficients are non-zero, 
  #                  otherwise 0. Require: The length of beta equals the ncol(X). 
  # Details:
  # Value: A random vector of binary responses, each element corresponding to the conditional probability given 
  #        each row of x, based on logistic model. 
  # Note: If the first column of X contains 1s, then the first element of beta is the intercept term. 
  # Author: ZWu
  # Reference(s): Wu et al. AOAS 2013
  # See also: 
  # Example: 
  #############################################################################
  #genBinY.logit <- function(X, beta) { #original name
  gen_Y_logit <- function(X, beta) {
    Y = apply(X, 1, function(x, beta){
      Py1x = 1/(1 + exp(-sum(x*beta))); # P(y=1 | x) based on logistic model
      return( ifelse(runif(1) < Py1x, 1, 0) ); #Generate 1 based on P(y=1 | x)
    },
    beta=beta);
    
    return (Y);
  }  
  
  
  #############################################################################
  # Description: Generate case-control data, including genotype and binary response values, 
  #              by retrieving cases and controls from a simulated population of genotypes and diseases. 
  # Depends on: self-implemented functions gen_bx and genBinY.logit
  # Arguments: 
  #   mafs    : Vector of minor allele frequence for genotypes.
  #	  GenoCorr: Correlation coefficient matrix for the genotypes.
  #		caseN   : Number of cases.
  #	  controlN: Number of controls.
  #		beta    : An vector of predictors'coefficients. For the true signals the coefficients 
  #             are non-zero, otherwise 0. Require: The length of beta equals ncol(X)+1, which 
  #             the first element being the intercept. 
  #		popnSize: The population size, for which the population data is going to be generated. 
  #             By "default", it is the 10 times of the caseN. 
  # Details: Disease model is defined by logistic regression. We generate large population, 
  #          and draw cases and controls from the populatio data.
  # Value: y vector of status (0 is control, 1 is case), corresponding x matrix of genotypes. 
  # Note:  You need to make popnSize large enough to make sure there always more cases in the 
  #        population data than caseN to be drawn. 
  # Author: SHe and ZWu
  # Reference(s): Wu et al. AOAS 2013
  # See also: 
  # Example: 
  #############################################################################
  genCaseControlData.fromPopn <- function(mafs, GenoCorr, caseN, controlN, beta, popnSize=caseN*10) {
    
    #Generate population data
    xs <- gen_bx(mafs, GenoCorr, popnSize);
    Ys = genBinY.logit(cbind(array(1, ncol(xs)), xs), beta);
    
    SampleCases	<- sample(which(Ys==1), caseN, replace = F);
    SampleContl	<- sample(which(Ys==0), controlN, replace = F);
    
    #random sample from population data
    x<- xs[c(SampleContl,SampleCases), ]; 
    y<- Ys[c(SampleContl,SampleCases)];			
    
    return(list(y=y, x=x));
  }  
  
  
  
  
  
  #############################################################################
  # Description: Generate case-control data, including genotype and binary response values, 
  #              based on marginal genetic effects. 
  # Depends on: 
  # Arguments: 
  #   caseN   : Number of cases.
  #		controlN: Number of controls.
  #  	maf     : Minor allele frequence (same for all genotypes).
  #		pD      : Disease prevalence, P(Disease), in population.
  # 	grr     : Marginal genetic relative risk. 
  #		NtrueSNP: Number of true SNPs. 
  #		NfalseSNP: Number of false SNPs. Require at least one of NtrueSNP and NfalseSNP be positive integer. 
  #		aEffType: Type of marginal allelic effect: "multiplicative", "dominant", "recessive", or "additive". 
  #		shuffleX: Whether we shuffle the column of X such that the true and false SNPs are located randomly. 
  # Details: Genotypes are indepdendent. Each genotype is randomly generated based on detectance 
  #          (i.e., P(genotype | Disease)). 
  # Value: y vector of status (0 is control, 1 is case). The corresponding x matrix of genotypes. 
  #        The first NtrueSNP columns are the genotypes of the true SNPs. A genotype value is the number of 
  #        copies of the minor allele. 
  # Note:   
  # Author:  ZWu
  # Reference(s): 
  # See also: 
  # Example: 
  #############################################################################
  genCaseControlData.margi <- function(caseN, controlN, maf, pD, grr, NtrueSNP=1, NfalseSNP=0, 
                                       aEffType="multiplicative", shuffleX=F) {
    
    paa = (1-maf)^2; pAa = 2*maf*(1-maf); pAA = maf^2; #genotype freq. Disease allele A is the minor allele
    pX = c(paa, pAa, pAA); #Marginal genotype frequency in population.
    genotype = c(0, 1, 2);
    
    ###penetrance P(D|aa), P(D|Aa), P(D|AA):
    if (aEffType =="multiplicative") { 
      pD.aa = pD / (paa + grr*pAa + grr^2*pAA); 
      pD.Aa = grr*pD.aa; pD.AA = grr^2*pD.aa; 
    } 
    if (aEffType =="dominant") { 
      pD.aa = pD / (paa + grr*pAa + grr*pAA); 
      pD.Aa = grr*pD.aa; pD.AA = grr*pD.aa; 
    } 
    if (aEffType =="recessive") { 
      pD.aa = pD / (paa + pAa + grr*pAA); 
      pD.Aa = pD.aa; pD.AA = grr*pD.aa; 
    }
    if (aEffType =="additive") 	{ 
      pD.aa = pD / (paa + grr*pAa + 2*grr*pAA); 
      pD.Aa = grr*pD.aa; pD.AA = 2*grr*pD.aa; 
    }
    if (pD.Aa > 1) stop("Error: P(D|Aa) > 1"); 
    if (pD.AA > 1) stop("Error: P(D|AA) > 1");
    
    ###Detectance
    pX.D = c(pD.aa, pD.Aa, pD.AA)*c(paa, pAa, pAA)/pD;
    pX.C = (1-c(pD.aa, pD.Aa, pD.AA))*c(paa, pAa, pAA) / (1-pD) ;
    
    if(NtrueSNP > 0) {
      x = matrix(array(NA, (caseN+controlN)*NtrueSNP), ncol=NtrueSNP);
      for (i in 1:NtrueSNP) {
        x[, i] = c(sample(genotype, caseN, replace=T, prob=pX.D), 
                   sample(genotype, controlN, replace=T, prob=pX.C));
      }
    }
    if (NfalseSNP > 0) 	{
      if (NtrueSNP > 0) {
        x = cbind(x, matrix(sample(genotype, (caseN+controlN)*NfalseSNP, replace=T, prob=pX), ncol= NfalseSNP));
      }else{
        x = matrix(sample(genotype, (caseN+controlN)*NfalseSNP, replace=T, prob=pX), ncol= NfalseSNP);
      }
    }
    
    y = c(rep(1, caseN), rep(0, controlN));
    
    if (shuffleX) {
      x = x[, sample(1:ncol(x))];
    }
    
    return(list(y=y, x=x));
  }  
  
  
  
