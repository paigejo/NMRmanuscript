# truth: realised value
# my.est: estimator
# my.var: variance
# lower: lower CI limit
# upper: upper CI limit
library(PearsonDS)
# library(gsl)
# library(moments)
logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  res = exp(x)/(1+exp(x))
  res[x > 100] = 1
  res[x < -100] = 0
  res
}

# this function computes all scoring rules, with and without binomial variation
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# numChildren: a vector containing the number of children with the same length as truth
# logitEst: the predictive mean on the logit scale, of the same length as truth
# logitVar: the predictive variance on the logit scale, of the same length as truth
# est: the estimate on the probability scale, of the same length as truth. By default calculated from probMat
# var: the estimate variance on the probability scale, of the same length as truth. By default calculated from probMat
# logitL: the lower end of the credible interval on the logit scale
# logitU: the upper end of the credible interval on the logit scale
# probMat: a matrix of joint draws of probability values, with number of rows equal to the length of truth, a number of 
#          columns equal to the number of draws. If not included, a gaussian distribution is assumed on the logit scale.
# significance: the significance level of the credible interval. By default 80%
# nsim: the number of sample quantiles from the gaussian distribution when approximating the proportion scale 
#       distributions from the logit scale
# NOTE: Discrete, count level credible intervals are estimated based on the input probMat along with coverage and CRPS
getScores = function(truth, numChildren, logitEst, logitVar, est=NULL, var=NULL, probMat=NULL, significance=.8, nsim=10, 
                     do1row=TRUE, includeBVarResults=FALSE, getRelativeBias=FALSE) {
  # first calculate bias (the same with and without binomial variation). Since on empirical proportion scale, set n to 1
  thisBias = bias(truth, logitEst, logit=FALSE, logitVar, n=1)
  
  if(getRelativeBias)
    thisRelativeBias = relativeBias(truth, logitEst, logit=FALSE, logitVar, n=1)
  
  # calculate average predictive variance on the logit scale
  # thisLogitVar = mean(logitVar)
  
  ##### The below code was commented out because it is done automatically in the coverage function:
  # # if necessary, calculate credible interval boundaries on the logit scale without binomial variation
  # if(is.null(logitL))
  #   logitL = qnorm((1 - significance) / 2, logitEst, sqrt(logitVar))
  # if(is.null(logitU))
  #   logitU = qnorm(1 - (1 - significance) / 2, logitEst, sqrt(logitVar))
  # 
  # # if necessary, do the same with binomial variation based on the joint simulations in the probability matrix
  # if(is.null(probL))
  #   probL = apply(probMat, 1, function(x) {quantile(x, prob=(1 - significance) / 2)})
  # if(is.null(probU))
  #   probU = apply(probMat, 1, function(x) {quantile(x, prob=1 - (1 - significance) / 2)})
  # 
  # # calculate credible interval width on the proportion scale with and without binomial variation
  # thisWidthNoBinom = mean(expit(logitU) - expit(logitL))
  # thisWidthBinom = mean(probU - probL)
  
  # compute central estimates if probMat is not null
  if(!is.null(probMat)) {
    if(is.null(est))
      est = rowMeans(probMat)
  }
  else {
    quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(logitEst)), ncol=nsim, byrow=TRUE)
    logitP = matrix(qnorm(quantiles, logitEst, sd=sqrt(logitVar)), ncol=nsim)
    probMat = expit(logitP)
    # if(is.null(var))
    #   var = apply(pMat, 1, sd)^2
    if(is.null(est))
      est = rowMeans(probMat)
  }
  
  # first calculate bias (the same with and without binomial variation). Since on empirical proportion scale, set n to 1
  out = mse(truth, logit(est), logit=FALSE, nsim=nsim, decompose = TRUE)
  thisMSE = out$MSE
  thisBias = out$bias
  thisVar = out$var
  
  # calculate coverage and credible interval width with and without binomial variation
  coverageNoBinom = coverage(truth, doLogit=FALSE, bVar=FALSE, numChildren=numChildren, logitEst=logitEst, logitVar=logitVar, 
                             lowerRejectProb=lowerRejectProb, upperRejectProb=upperRejectProb, probMat=probMat, 
                             significance=significance, returnIntervalWidth=TRUE, nsim=nsim)
  thisCoverageNoBinom = coverageNoBinom[1]
  thisWidthNoBinom = coverageNoBinom[2]
  if(includeBVarResults) {
    coverageBinom = coverage(truth, doLogit=FALSE, bVar=TRUE, numChildren=numChildren, logitEst=logitEst, logitVar=logitVar, 
                             lowerRejectProb=lowerRejectProb, upperRejectProb=upperRejectProb, probMat=probMat, 
                             significance=significance, returnIntervalWidth=TRUE, returnPMFs=TRUE, nsim=nsim)
    thisCoverageBinom = coverageBinom$intervalInfo[1]
    thisWidthBinom = coverageBinom$intervalInfo[2]
    
    # also get the probability mass functions for the binomial predictive distributions for the CRPS calculations
    pmfs = coverageBinom$pmfs
  }
  else {
    thisCoverageBinom = NA
    thisWidthBinom = NA
    pmfs = NULL
  }
  
  # calculate CRPS on the proportion scale with and without binomial variation
  thisCRPSNoBinom = crpsNormal(truth, logitEst, logitVar, logit=FALSE, n=numChildren, bVar=FALSE, probMat=probMat, nsim=nsim)
  if(includeBVarResults) {
    thisCRPSBinom = crpsNormal(truth, logitEst, logitVar, logit=FALSE, n=numChildren, bVar=TRUE, probMat=probMat, nsim=nsim, pmfs=pmfs)
  }
  else {
    thisCRPSBinom = NA
  }
  
  # return the results
  if(!getRelativeBias) {
    if(do1row) {
      results = matrix(c(thisBias, thisVar, thisMSE, thisCRPSNoBinom, thisCRPSBinom, thisCoverageNoBinom, thisCoverageBinom, 
                         thisWidthNoBinom, thisWidthBinom), nrow=1)
      colnames(results) = c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")
      
      as.data.frame(results)
    } else {
      results = matrix(c(thisBias, thisVar, thisMSE, thisCRPSNoBinom, thisCoverageNoBinom, thisWidthNoBinom, 
                         thisBias, thisVar, thisMSE, thisCRPSBinom, thisCoverageBinom, thisWidthBinom), 
                       nrow=2, byrow = TRUE)
      colnames(results) = c("bias", "var", "mse", "crps", "coverage", "length")
      as.data.frame(results)
    }
  }
  else {
    if(do1row) {
      results = matrix(c(thisBias, thisRelativeBias, thisVar, thisMSE, thisCRPSNoBinom, thisCRPSBinom, thisCoverageNoBinom, thisCoverageBinom, 
                         thisWidthNoBinom, thisWidthBinom), nrow=1)
      colnames(results) = c("bias", "relativeBias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")
      
      as.data.frame(results)
    } else {
      results = matrix(c(thisBias, thisRelativeBias, thisVar, thisMSE, thisCRPSNoBinom, thisCoverageNoBinom, thisWidthNoBinom, 
                         thisBias, thisRelativeBias, thisVar, thisMSE, thisCRPSBinom, thisCoverageBinom, thisWidthBinom), 
                       nrow=2, byrow = TRUE)
      colnames(results) = c("bias", "relativeBias", "var", "mse", "crps", "coverage", "length")
      as.data.frame(results)
    }
  }
}

# SPDE model scores are calculated differently, because it has three levels of spatial aggregation:
#  - inexact integration: 
#    - continuously integrate predictions with respect to population density, ignoring binomial variation
#  - exact integration without binomial variation: 
#    - integrate predictions at discrete enumeration area locations, ignoring binomial variation
#  - exact integration with binomial variation: 
#    - integrate predictions at discrete enumeration area locations, accounting for binomial variation
# The first two cases use the probMat input but don't have pmfs, whereas the latter uses pmfs but not probMat. 
# Estimates are already assumed to be debiased and provided at the logit scale in all cases.
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# numChildren: a vector containing the number of children with the same length as truth
# logitEst: the predictive mean on the logit scale, of the same length as truth
# logitVar: the predictive variance on the logit scale, of the same length as truth
# est: the estimate on the probability scale, of the same length as truth. By default calculated from probMat
# var: the estimate variance on the probability scale, of the same length as truth. By default calculated from probMat
# logitL: the lower end of the credible interval on the logit scale
# logitU: the upper end of the credible interval on the logit scale
# probMat: a matrix of joint draws of probability values, with number of rows equal to the length of truth, a number of 
#          columns equal to the number of draws. If not included, a gaussian distribution is assumed on the logit scale.
# significance: the significance level of the credible interval. By default 80%
# nsim: the number of sample quantiles from the gaussian distribution when approximating the proportion scale 
#       distributions from the logit scale
# bVar: whether or not to include binomial variation
# pmfs: a list of probability mass functions with the same length as truth
# marginals: list of INLA marginals on the logit probability scale
# NOTE: Discrete, count level credible intervals are estimated based on the input probMat along with coverage and CRPS
getScoresSPDE = function(truth, numChildren, logitEst, est=NULL, var=NULL, probMat=NULL, significance=.8, nsim=10, 
                         bVar, pmfs=NULL, logitVar=NULL, logitL=NULL, logitU=NULL, getRelativeBias=FALSE) {
  
  # calculate predictive estimates at the proportion scale
  if(!is.null(pmfs)) {
    if(!bVar && is.null(probMat))
      stop("bVar FALSE, but only pmfs included")
    
    # var = sapply(1:length(pmfs), function(i) {sum(((0:numChildren[i]) / numChildren[i])^2 * pmfs[[i]]) - (sum(((0:numChildren[i]) / numChildren[i]) * pmfs[[i]]))^2})
    est = sapply(1:length(pmfs), function(i) {sum((0:length(pmfs[[i]])) * pmfs[[i]]) / length(pmfs[[i]])})
  }
  else if(!is.null(probMat)) {
    # if(is.null(var))
    # var = apply(probMat, 1, sd)^2
    
    if(is.null(var))
      est = rowMeans(probMat)
  }
  else {
    stop("For the SPDE model, either marginals or probMat must be provided")
  }
  
  # calculate mean squared error and bias at the proportion scale
  out = mse(truth, est, logit=TRUE, nsim=nsim, n=1, decompose=TRUE)
  thisMSE = out$MSE
  thisVar = out$var
  thisBias = out$bias
  
  # get relative bias (relative to the truth)
  if(getRelativeBias) {
    thisRelativeBias = relativeBias(truth, logitEst, logit=FALSE, logitVar, n=1)
  }
  
  # if CIs are included, but binomial variation is not, calculate credible intervals based on them
  if(!(any(is.null(logitL)) || any(is.null(logitU))) && !bVar) {
    lower = logitL
    upper = logitU
  } else {
    lower = NULL
    upper = NULL
  }
  
  # calculate coverage and credible interval width with and without binomial variation
  cvg = coverage(truth, doLogit=FALSE, bVar=bVar, numChildren=numChildren, logitEst=logitEst, logitVar=logitVar, 
                 lowerRejectProb=lowerRejectProb, upperRejectProb=upperRejectProb, probMat=probMat, 
                 significance=significance, returnIntervalWidth=TRUE, returnPMFs=TRUE, pmfs=pmfs, nsim=nsim, 
                 lower=lower, upper=upper)
  
  thisCoverage = cvg$intervalInfo[1]
  thisWidth = cvg$intervalInfo[2]
  
  # also get the probability mass functions for the binomial predictive distributions for the CRPS calculations
  if(bVar)
    pmfs = cvg$pmfs
  
  # calculate CRPS on the proportion scale with and without binomial variation
  thisCRPS = crpsNormal(truth, logitEst, logitVar, logit=FALSE, n=numChildren, bVar=bVar, probMat=probMat, 
                        nsim=nsim, pmfs=pmfs)
  
  # return the results
  if(!getRelativeBias) {
    results = matrix(c(thisBias, thisVar, thisMSE, thisCRPS, thisCoverage, thisWidth), nrow=1, byrow = TRUE)
    colnames(results) = c("bias", "var", "mse", "crps", "coverage", "length")
  }
  else {
    results = matrix(c(thisBias, thisRelativeBias, thisVar, thisMSE, thisCRPS, thisCoverage, thisWidth), nrow=1, byrow = TRUE)
    colnames(results) = c("bias", "relativeBias", "var", "mse", "crps", "coverage", "length")
  }
  as.data.frame(results)
}

##### the following function computed scoring rules at the county level after leaving out county data. We decided 
##### to only calculate scoring rules at the cluster level
# # a few resources for scoring rules:
# # http://webpages.math.luc.edu/~ebalderama/myfiles/modelchecking101_pres.pdf
# # http://www.math.chalmers.se/~bodavid/GMRF2015/Lectures/Flab4.pdf
# # http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3409851/
# # https://www.researchgate.net/profile/Havard_Rue/publication/226713867_Posterior_and_Cross-validatory_Predictive_Checks_A_Comparison_of_MCMC_and_INLA/links/551faa8c0cf29dcabb088bc6.pdf
# # original definition of cpo in Geisser's discussion in:
# # https://www.jstor.org/stable/pdf/2982063.pdf?refreqid=excelsior%3Aff026889bf010b97c9d0c7ae7bed6935
# # NOTE: for cluster level validation scores set the direct estimate values to be the truth and the direct estimate variance to 0
# # doLogit: If it is possible to take the logit of the direct estimates, then set this to TRUE 
# getValidationScores = function(logitDirectEsts, logitDirectVars, logitEsts=NULL, logitVars=NULL, ests=NULL, 
#                                directEsts=NULL, probMat=NULL, weights=rep(1, length(logitDirectEsts)), 
#                                logitWeights=rep(1, length(logitDirectEsts)), usePearson=TRUE, nsim=10, 
#                                n=25, bVar=FALSE) {
#   clusterValidation = all(logitDirectVars == 0)
#   
#   # weight each county prediction proportional to the inverse direct estimate variance
#   weights = weights / sum(weights)
#   logitWeights = logitWeights / sum(logitWeights)
#   
#   # calculate logit scale predictive distribution if necessary
#   if(is.null(logitEsts))
#     logitEsts = rowMeans(logit(probMat))
#   if(is.null(logitVars))
#     logitVars = apply(logit(probMat), 1, var)
#   if(is.null(directEsts))
#     directEsts = logitNormMean(cbind(logitDirectEsts, logitDirectVars))
#   
#   # first calculate mean squared error on probability and logit scales, decomposed into variance and bias squared
#   if(is.null(ests))
#     out = mse(directEsts, logitEsts, logit=FALSE, logitVars, nsim=nsim, weights=weights, decompose=TRUE)
#   else
#     out = mse(directEsts, ests, logit=TRUE, decompose=TRUE)
#   MSE = out$MSE
#   Var = out$var
#   Biassq = out$biassq
#   
#   if(!clusterValidation) {
#     outLogit = mse(logitDirectEsts, logitEsts, logit=TRUE, logitVars, nsim=nsim, weights=logitWeights, decompose=TRUE)
#     MSELogit = outLogit$MSE
#     VarLogit = outLogit$var
#     BiassqLogit = outLogit$biassq
#   }
#   else {
#     MSELogit = NULL
#     VarLogit = NULL
#     BiassqLogit = NULL
#   }
#   
#   # the other scoring rules require us to know the full predictive distribution. Either 
#   # use a gaussian or pearson approximation
#   if(!bVar) {
#     if(usePearson && !is.null(probMat)) {
#       # for the Pearson distribution to the probabilities
#       skewnessVals = apply(probMat, 1, skewness)
#       kurtosisVals = apply(probMat, 1, kurtosis)
#       
#       # using Pearson type IV approximation
#       momentMat = cbind(logitEsts, logitVars, skewnessVals, kurtosisVals)
#       # pearsonParMat = apply(momentMat, 1, pearsonFitM)
#       pearsonParMat = mapply(pearsonFitM, momentMat[,1], momentMat[,2], momentMat[,3], momentMat[,4])
#       pearsonParMat = split(t(pearsonParMat), 1:ncol(pearsonParMat))
#       
#       # calculate CRPS
#       
#       # calculate CPO, log score as defined by LPML in http://webpages.math.luc.edu/~ebalderama/myfiles/modelchecking101_pres.pdf
#       cpos = mapply(dpearson, x=logitDirectEsts, params=pearsonParMat)
#       cpo = mean(cpos)
#       logScore = mean(log(cpos))
#       
#       crps = crpsPearsonDirect(logitDirectEsts, logitDirectVars, pearsonParMat)
#       crps = out[1]
#       
#       # calculate PIT
#       pits = mapply(ppearson, q=logitDirectEsts, params=pearsonParMat)
#     }
#     else {
#       # log score as defined by LPML in http://webpages.math.luc.edu/~ebalderama/myfiles/modelchecking101_pres.pdf
#       crps = crpsNormalDirect(logitDirectEsts, logitDirectVars, logitEsts, logitVars)
#       cpos = dnorm(logitDirectEsts, logitEsts, sqrt(logitVars))
#       cpo = mean(cpos)
#       logScore = mean(log(cpos))
#       pits = pnorm(logitDirectEsts, logitEsts, sqrt(logitVars))
#     }
#   }
#   else {
#     # in this case, include binomial variation
#     out = getCPO(directEsts, logitEsts, logitVars, n=n, bVar = TRUE)
#     cpo = out$CPO
#     logScore = out$logScore
#     pits = out$PITs
#     
#     crps = crpsNormal(directEsts, logitEsts, logitVars, logit=FALSE, bVar = TRUE)
#   }
#   
#   # return the results
#   if(!clusterValidation) {
#     results = matrix(c(MSE, Var, Biassq, MSELogit, VarLogit, BiassqLogit, cpo, crps, logScore), nrow=1, byrow = TRUE)
#     colnames(results) = c("MSE", "Var", "Bias^2", "Logit MSE", "Logit Var", "Logit Bias^2", "CPO", "CRPS", "logScore")
#   }
#   else {
#     results = matrix(c(MSE, Var, Biassq, cpo, crps, logScore), nrow=1, byrow = TRUE)
#     colnames(results) = c("MSE", "Var", "Bias^2", "CPO", "CRPS", "logScore")
#   }
#   list(scores=as.data.frame(results), pit=pits, weights=weights, logitWeights=logitWeights)
# }

# a few resources for scoring rules:
# http://webpages.math.luc.edu/~ebalderama/myfiles/modelchecking101_pres.pdf
# http://www.math.chalmers.se/~bodavid/GMRF2015/Lectures/Flab4.pdf
# http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3409851/
# https://www.researchgate.net/profile/Havard_Rue/publication/226713867_Posterior_and_Cross-validatory_Predictive_Checks_A_Comparison_of_MCMC_and_INLA/links/551faa8c0cf29dcabb088bc6.pdf
# original definition of cpo in Geisser's discussion in:
# https://www.jstor.org/stable/pdf/2982063.pdf?refreqid=excelsior%3Aff026889bf010b97c9d0c7ae7bed6935
# INPUTS
# truth: true empirical proportions in the clusters
# urbanVec: a logical vector with the same length as logitDirectEsts. Used to average scoring rules 
#           over urban and rural clusters separately
# filterType: although all scoring rules are calculated, this determines an additional filtered set of final 
#             scoring rules to set aside for display purposes
# NOTE: binomial variation is always included. Do not use this function for the smoothed direct model
getValidationScores = function(truth, logitEsts=NULL, logitVars=NULL, ests=NULL, 
                               probMat=NULL, weights=rep(1, length(truth)), 
                               usePearson=FALSE, nsim=10, n=25, urbanVec=NULL, 
                               filterType=c("inSample", "leftOutCounty")) {
  
  filterType = match.arg(filterType)
  
  # weight each county prediction proportional to the inverse direct estimate variance
  weights = weights / sum(weights)
  
  # calculate logit scale predictive distribution if necessary
  if(is.null(logitEsts))
    logitEsts = rowMeans(logit(probMat))
  if(is.null(logitVars))
    logitVars = apply(logit(probMat), 1, var)
  
  # first calculate mean squared error on probability and logit scales, decomposed into variance and bias squared
  if(is.null(ests))
    out = mse(truth, logitEsts, logit=FALSE, logitVars, nsim=nsim, weights=weights, decompose=TRUE)
  else
    out = mse(truth, ests, logit=TRUE, decompose=TRUE)
  MSE = out$MSE
  Var = out$var
  Bias = out$bias
  
  # calculate the rest of the scoring rules
  out = getCPO(truth, logitEsts, logitVars, n=n, probMat=probMat, bVar = TRUE)
  cpo = out$CPO
  logScore = out$logScore
  pits = out$PITs
  
  crps = crpsNormal(truth, logitEsts, logitVars, logit=FALSE, bVar = TRUE)
  
  allResults = matrix(c(MSE, Var, Bias, cpo, crps, logScore), nrow=1, byrow = TRUE)
  colnames(allResults) = c("MSE", "Var", "Bias", "CPO", "CRPS", "logScore")
  
  if(filterType == "inSample")
    filterI = c(1:3, 5)
  else if(filterType == "leftOutCounty")
    filterI = c(1:5)
  else
    stop(paste0("unrecognized value of filterType ", filterType))
  
  outUrban = NULL
  outRural = NULL
  if(!is.null(urbanVec)) {
    urbanProbMat = NULL
    ruralProbMat = NULL
    if(!is.null(probMat)) {
      urbanProbMat = probMat[urbanVec,]
      ruralProbMat = probMat[!urbanVec,]
    }
    outUrban = getValidationScores(truth[urbanVec], logitEsts[urbanVec], logitVars[urbanVec], ests[urbanVec], urbanProbMat, 
                                   weights[urbanVec], usePearson, nsim, n[urbanVec], NULL, filterType)
    outRural = getValidationScores(truth[!urbanVec], logitEsts[!urbanVec], logitVars[!urbanVec], ests[!urbanVec], ruralProbMat, 
                                   weights[!urbanVec], usePearson, nsim, n[!urbanVec], NULL, filterType)
    
    allResults = matrix(c(MSE, Var, Bias, cpo, crps, logScore, outUrban$allResults, outRural$allResults),
                        nrow=1, byrow = TRUE)
    allNames = c("MSE", "Var", "Bias", "CPO", "CRPS", "logScore")
    temp1 = paste0("Urban ", allNames)
    temp2 = paste0("Rural ", allNames)
    allNames = c(allNames, temp1, temp2)
    colnames(allResults) = allNames
    
    filterI = c(filterI, filterI + 6, filterI + 2 * 6)
  }
  
  # filter out unwanted results for the final results table
  scores = data.frame(as.list(allResults[,filterI]))
  allResults = data.frame(as.list(allResults))
  
  # return the results
  list(scores=scores, pit=pits, weights=weights, 
       allResults=allResults, resultsUrban=outUrban, resultsRural=outRural)
}

# based off of calculations in:
# http://www.stat.columbia.edu/~gelman/research/published/waic_understand3.pdf
# http://www.r-inla.org/faq#TOC-How-are-the-Devicance-Information-Criteria-DIC-and-The-Watanabe-Akaike-information-criterion-WAIC-computed-
getDICWAIC = function(logitDirectEsts, logitDirectVars, logitEsts=NULL, logitVars=NULL, 
                      probMat=NULL, weights, usePearson=TRUE, nsim=10) {
  
  # weight each county prediction proportional to the inverse direct estimate variance
  weights = weights / sum(weights)
  
  # calculate logit scale predictive distribution if necessary
  if(is.null(logitEsts))
    logitEsts = rowMeans(logit(probMat))
  if(is.null(logitVars))
    logitVars = apply(logit(probMat), 1, var)
  
  # the other scoring rules require us to know the full predictive distribution. Either 
  # use a gaussian or pearson approximation
  if(usePearson && !is.null(probMat)) {
    # for the Pearson distribution to the probabilities
    skewnessVals = apply(probMat, 1, skewness)
    kurtosisVals = apply(probMat, 1, kurtosis)
    
    # using Pearson type IV approximation
    momentMat = cbind(logitEsts, logitVars, skewnessVals, kurtosisVals)
    # pearsonParMat = apply(momentMat, 1, pearsonFitM)
    pearsonParMat = mapply(pearsonFitM, momentMat[,1], momentMat[,2], momentMat[,3], momentMat[,4])
    pearsonParMat = split(t(pearsonParMat), 1:ncol(pearsonParMat))
    
    # calculate DIC = -2log p(y | \hat{theta})
    cpos = mapply(dpearson, x=logitDirectEsts, params=pearsonParMat)
    cpo = sum(cpos * weights)
    
    # calculate CRPS, log score as defined by LPML in http://webpages.math.luc.edu/~ebalderama/myfiles/modelchecking101_pres.pdf
    out = crpsPearsonDirect(logitDirectEsts, logitDirectVars, pearsonParMat)
    crps = out[1]
    logScore = out[2]
    
    # calculate PIT
    pits = mapply(ppearson, q=logitDirectEsts, params=pearsonParMat)
  }
  else {
    # calculate CRPS using Gaussian approximation
    cpos = dnorm(logitDirectEsts, logitEsts, logitVars)
    cpo = sum(cpos * weights)
    
    # calculate CRPS, log score as defined by LPML in http://webpages.math.luc.edu/~ebalderama/myfiles/modelchecking101_pres.pdf
    out = crpsNormalDirect(logitDirectEsts, logitDirectVars, logitEsts, logitVars)
    crps = out[1]
    logScore = out[2]
    
    # calculate PIT
    pits = mapply(pnorm, q=logitDirectEsts, logitEsts, logitVars)
  }
  
  # return the results
  results = matrix(c(MSE, cpo, crps, logScore), nrow=1, byrow = TRUE)
  colnames(results) = c("MSE", "CPO", "CRPS", "logScore")
  list(scores=as.data.frame(results), pit=pits, weights=weights)
}

mse <- function(truth, my.est, logit=TRUE, my.var=NULL, nsim=10, n=1, weights=NULL, decompose=FALSE, urbanVec=NULL){
  if(!is.null(weights))
    weights = weights / sum(weights)
  
  if(!logit){
    # this is the MSE of the median prediction
    # truth = expit(truth)
    # my.est = expit(my.est)
    
    # calculate the MSE of the mean prediction:
    # first draw different values of p from the predictive distribution
    if(is.null(my.var)) {
      my.est = expit(my.est)
    } else {
      quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(my.est)), ncol=nsim, byrow=TRUE)
      logitP = matrix(qnorm(quantiles, my.est, sd=sqrt(my.var)), ncol=nsim)
      my.est = rowMeans(expit(logitP))
    }
  }
  res = my.est - truth
  
  if(decompose) {
    thisBias = my.est - truth
  }
  
  if(!is.null(weights)) {
    MSE = sum(res^2 * weights)
    if(decompose) {
      bias=sum(res * weights)
      thisVar = (res - sum(res*weights))^2
      thisVar = sum(thisVar * weights)
      out = list(MSE=MSE, bias=bias, var=thisVar)
    }
    else
      out = MSE
  }
  else {
    MSE = mean(res^2)
    if(decompose) {
      bias=mean(res)
      thisVar = (res - mean(res))^2
      thisVar = mean(thisVar)
      out = list(MSE=MSE, bias=bias, var=thisVar)
    }
    else
      out = MSE
  }
  
  if(!is.null(urbanVec) && decompose) {
    # compute MSE and related values for urban and rural areas if necessary
    urbanN = n[urbanVec]
    urbanN = urbanN[!is.na(urbanN)]
    ruralN = n[ruralVec]
    ruralN = ruralN[!is.na(ruralN)]
    outUrban = mse(truth[urbanVec], my.est[urbanVec], logit, my.var[urbanVec], nsim, urbanN, weights[urbanVec], decompose)
    outRural = mse(truth[!urbanVec], my.est[!urbanVec], logit, my.var[!urbanVec], nsim, ruralN, weights[!urbanVec], decompose)
    
    # combine results
    out = list(MSE=MSE, bias=bias, var=thisVar, MSEUrban=outUrban$MSE, biasUrban=outUrban$bias, varUrban=outUrban$var, MSERural=outRural$MSE, biasRural=outRural$bias, varRural=outRural$var)
  }
  
  out
}

# logit: if TRUE use the raw estimates, otherwise assume estimates are on the logit scale
bias <- function(truth, my.est, logit=TRUE, my.var=NULL, n=1, weights=NULL){
  if(!is.null(weights))
    weights = weights / sum(weights)
  
  if(!logit) {
    # this is the bias of the median prediction, not the mean
    # truth = expit(truth)
    # my.est = expit(my.est)
    
    # debias the estimate if necessary
    if(is.null(my.var))
      warning("logit variance not included for bias calculation, so assuming estimate does not need to be debiased")
    else
      my.est = logitNormMean(cbind(my.est, sqrt(my.var)))
  }
  res <- my.est - truth
  
  if(!logit) {
    if(!is.null(weights))
      return(n * sum(res * weights))
    else
      return(mean(n * res))
  }
  
  if(is.null(weights))
    mean(res)
  else
    sum(res * weights)
}

# logit: if TRUE use the raw estimates, otherwise assume estimates are on the logit scale
relativeBias <- function(truth, my.est, logit=TRUE, my.var=NULL, n=1, weights=NULL){
  if(!is.null(weights))
    weights = weights / sum(weights)
  
  if(!logit) {
    # this is the bias of the median prediction, not the mean
    # truth = expit(truth)
    # my.est = expit(my.est)
    
    # debias the estimate if necessary
    if(is.null(my.var))
      warning("logit variance not included for bias calculation, so assuming estimate does not need to be debiased")
    else
      my.est = logitNormMean(cbind(my.est, sqrt(my.var)))
  }
  res <- (my.est - truth) / truth
  
  if(!logit) {
    if(!is.null(weights))
      return(n * sum(res * weights))
    else
      return(mean(n * res))
  }
  
  if(is.null(weights))
    mean(res)
  else
    sum(res * weights)
}

# truth: a vector of observations on the desired scale
# my.est: a vector of logit-scale predictions of the same length as truth 
# my.var: a vector of logit-scale predictive variances of the same length as truth
# logit: whether CRPS should be computed on the logit scale or the empirical proportion scale
# nsim: if binomial variation is included, the number of probabilities over which to integrate in the pmf approximation
# n: the binomial sample size. Can be a vector of the same length as truth
# bVar: whether or not binomial variation should be included in the predictive distribution
# probMat: if included, use these probability draws in the pmf integration/approximation. If the estimates are normal 
#          on the logit scale, do not include this since just using the estimates and variance will be accurate.
# pmfs: a list of probability mass functions of the same length as truth
# marginals: a list of inla marginals of the same length as truth
crpsNormal <- function(truth, my.est, my.var=NULL, logit=TRUE, nsim=10, n=25, bVar=TRUE, probMat=NULL, pmfs=NULL, marginals=NULL){
  # if(resultType == "EA" || resultType == "pixel")
  #   logit=FALSE
  
  if(logit) {
    sig = sqrt(my.var)
    x0 <- (truth - my.est) / sig
    res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
    
    ## sign as in Held (2008)
    res <- -res
  }
  else {
    if(bVar) {
      # first draw different values of p from the predictive distribution if necessary
      if(is.null(probMat) && is.null(pmfs)) {
        quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(my.est)), ncol=nsim, byrow=TRUE)
        logitP = matrix(qnorm(quantiles, my.est, sd=sqrt(my.var)), ncol=nsim)
        probMat = expit(logitP)
      }
      
      # for each value of p, calculate the crps and return the mean
      crpsVals = crpsBinomial(truth, probMat, n=n, pmfs=pmfs)
      res = mean(crpsVals)
    } else {
      # no binomial variation is included. Integrate numerically on the proportion scale
      
      if(!is.null(marginals) && is.null(probMat) && is.null(my.var)) {
        ## too slow
        # CRPSMarginal = function(i) {
        #   integrate(function(p) {(inla.pmarginal(logit(p), marginals[[i]]) - (p >= truth[i]))^2}, 0, 1)$value
        # }
        # crpsVals = sapply(1:length(marginals), CRPSMarginal)
        
        ## also too slow
        # quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(my.est)), ncol=nsim, byrow=TRUE)
        # quantilesMarginal = function(i) {
        #   inla.qmarginal(quantiles, marginals[[i]])
        # }
        # logitP = sapply(1:length(marginals), quantilesMarginal)
        # probMat = expit(logitP)
        
        quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(my.est)), ncol=nsim, byrow=TRUE)
        logitP = matrix(qnorm(quantiles, my.est, sd=sqrt(my.var)), ncol=nsim)
        probMat = expit(logitP)
      }
      
      # # compute the crps for this row of truth
      # crpsRow = function(rowI) {
      #   thisTruth = truth[rowI]
      #   
      #   # either build the predictive cdf assuming normality on the logit scale or from the 
      #   # empirical distribution given by probMat if the user supplies it
      #   if(is.null(probMat)) {
      #     thisEst = my.est[rowI]
      #     thisVar = my.var[rowI]
      #     thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
      #   } else {
      #     thisCdf = ecdf(probMat[rowI,])
      #   }
      #   
      #   intFun = function(ws) {
      #     (thisCdf(ws) - (ws >= thisTruth))^2
      #   }
      #   
      #   if(is.null(probMat)) {
      #     # when integrating we will set bounds on the integral to be reasonable to avoid 
      #     # faulty "divergent integral" error. The bounds will be 20 standard errors out 
      #     # of the estimate, making sure to include the truth.
      #     lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
      #     upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
      #     integrate(intFun, lowerBound, upperBound)$value
      #   }
      #   else {
      #     # since we are using the empirical distribution, there is a closed form for the integral
      #     ps = probMat[rowI,]
      #     allPoints = sort(c(ps, thisTruth))
      #     deltas = diff(allPoints)
      #     sum(deltas * intFun(allPoints[1:length(ps)]))
      #   }
      # }
      # 
      # crpsRow2 = function(rowI) {
      #   thisTruth = truth[rowI]
      #   
      #   # either build the predictive cdf assuming normality on the logit scale or from the 
      #   # empirical distribution given by probMat if the user supplies it
      #   if(is.null(probMat)) {
      #     thisEst = my.est[rowI]
      #     thisVar = my.var[rowI]
      #     thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
      #   } else {
      #     # thisCdf = ecdf(probMat[rowI,])
      #     sorted = sort(probMat[rowI,])
      #     thisCdf = approxfun(sorted, (1:length(sorted))/length(sorted), 
      #                         method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
      #   }
      #   
      #   intFun = function(ws) {
      #     (thisCdf(ws) - (ws >= thisTruth))^2
      #   }
      #   
      #   if(is.null(probMat)) {
      #     # when integrating we will set bounds on the integral to be reasonable to avoid 
      #     # faulty "divergent integral" error. The bounds will be 20 standard errors out 
      #     # of the estimate, making sure to include the truth.
      #     lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
      #     upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
      #     integrate(intFun, lowerBound, upperBound)$value
      #   }
      #   else {
      #     # since we are using the empirical distribution, there is a closed form for the integral
      #     allPoints = sort(c(sorted, thisTruth))
      #     firstGreater = match(TRUE, sorted >= thisTruth)
      #     if(is.na(firstGreater))
      #       allPoints = c(sorted, thisTruth)
      #     else if(firstGreater == 1)
      #       allPoints = c(thisTruth, sorted)
      #     else
      #       allPoints = c(sorted[1:(firstGreater - 1)], thisTruth, sorted[firstGreater:length(sorted)])
      #     
      #     deltas = diff(allPoints)
      #     sum(deltas * intFun(allPoints[1:length(sorted)]))
      #   }
      # }
      
      crpsRow = function(rowI) {
        thisTruth = truth[rowI]
        
        # either build the predictive cdf assuming normality on the logit scale or from the 
        # empirical distribution given by probMat if the user supplies it
        if(is.null(probMat)) {
          thisEst = my.est[rowI]
          thisVar = my.var[rowI]
          thisCdf = function(ws) {pnorm(logit(ws), thisEst, sqrt(thisVar))}
          
          intFun = function(ws) {
            (thisCdf(ws) - (ws >= thisTruth))^2
          }
        } else {
          # thisCdf = ecdf(probMat[rowI,])
          sorted = probMat[rowI,] # already sorted
          
          # since we are using the empirical distribution, there is a closed form for the integral
          allPoints = sort(c(sorted, thisTruth))
          deltas = diff(allPoints)
          firstGreater = match(TRUE, sorted >= thisTruth)
          vals = (1:length(sorted))/length(sorted)
          if(is.na(firstGreater))
            return(sum((vals)^2 * deltas))
          else if(firstGreater == 1)
            return(deltas[1] + sum((1-vals[1:(length(sorted)-1)])^2 * deltas[2:length(deltas)]))
          else {
            left = sum(vals[1:(firstGreater-1)]^2 * deltas[1:(firstGreater-1)])
            mid = sum((1 - vals[firstGreater-1])^2 * deltas[firstGreater])
            right = ifelse(firstGreater == length(vals), 0, sum((1 - vals[firstGreater:(length(vals)-1)])^2 * deltas[(firstGreater+1):length(deltas)]))
            return(left+mid+right)
          }
          
          # intFun = function(ws) {
          #   (thisCdf(ws) - (ws >= thisTruth))^2
          # }
        }
        
        if(is.null(probMat)) {
          # when integrating we will set bounds on the integral to be reasonable to avoid 
          # faulty "divergent integral" error. The bounds will be 20 standard errors out 
          # of the estimate, making sure to include the truth.
          lowerBound = max(0, min(thisTruth - .01, expit(thisEst - 20 * sqrt(thisVar))))
          upperBound = min(1, max(thisTruth + .01, expit(thisEst + 20 * sqrt(thisVar))))
          integrate(intFun, lowerBound, upperBound)$value
        }
        else {
          stop("should have already returned result")
        }
      }
      
      if(!is.null(probMat))
        probMat = t(apply(probMat, 1, sort))
      crpsVals = sapply(1:length(truth), crpsRow)
      
      res = mean(crpsVals)
    }
  }
  
  res
}

crpsNormalDirect <- function(directLogitEsts, directLogitVars, logitEsts, logitVars) {
  
  integrand = function(p, i) {
    directSD = sqrt(directLogitVars[i])
    directMu <- directLogitEsts[i]
    estSD = sqrt(logitVars[i])
    estMu = logitEsts[i]
    
    diff = pnorm(logit(p), directMu, directSD) - pnorm(logit(p), estMu, estSD)
    diff^2
  }
  getOneCRPS = function(i) {
    minVal = min(directLogitEsts[i] - 5 * sqrt(directLogitVars[i]), logitEsts[i] - 5 * sqrt(logitVars[i]))
    maxVal = max(directLogitEsts[i] + 5 * sqrt(directLogitVars[i]), logitEsts[i] + 5 * sqrt(logitVars[i]))
    integrate(integrand, expit(minVal), expit(maxVal), i=i)$value
  }
  
  res = sapply(1:length(directLogitEsts), getOneCRPS)
  mean(res)
}

crpsPearsonDirect <- function(directLogitEsts, directLogitVars, pearsonParameters) {
  
  integrand = function(p, i) {
    directSD = sqrt(directLogitVars[i])
    directMu <- directLogitEsts[i]
    theseParameters = unlist(pearsonParameters[i])
    
    diff = pnorm(logit(p), directMu, directSD) - ppearson(logit(p), params=theseParameters)
    diff^2
  }
  getOneCRPS = function(i) {
    theseParameters = unlist(pearsonParameters[i])
    minVal = expit(min(directLogitEsts[i] - 5 * sqrt(directLogitVars[i]), qpearson(pnorm(-5), params=theseParameters)))
    maxVal = expit(max(directLogitEsts[i] + 5 * sqrt(directLogitVars[i]), qpearson(pnorm(5), params=theseParameters)))
    integrate(integrand, minVal, maxVal, i=i)$value
  }
  
  res = sapply(1:length(directLogitEsts), getOneCRPS)
  mean(res)
}

# for a set of independent fixed binomials, determine the CRPS.
# trueProportion: can be a vector
# p.est: can be a matrix, with number of rows equal to length(trueProportion). In that case, the 
#        number of columns can be thought of as the number of draws of p in the case of a random p.
# n: can be a vector, with length equal to length(trueProportion)
# empiricalProportion: if true, calculate the CRPS of the empirical proportion distribution on [0, 1] 
#                      rather than the binomial distribution on [0, n]
# pmfs: a list of probability mass functions with the same length as trueProportion
crpsBinomial <- function(trueProportion, p.est=NULL, n=25, parClust=NULL, empiricalProportion=TRUE, pmfs=NULL) {
  # get observations at the count scale, make sure n is a vector
  trueCounts = trueProportion * n
  if(length(n) == 1)
    n = rep(n, length(trueProportion))
  
  if(is.null(pmfs)) {
    p.est = matrix(p.est, nrow=length(trueProportion))
  }
  else
    cdfs = lapply(pmfs, function(pmf) {cumsum(pmf)})
  
  rowFun = function(rowI) {
    thisCount = eval(trueCounts[rowI])
    toN = eval(0:n[rowI])
    
    if(is.null(pmfs)) {
      ## in this case we must compute the probability mass function for this region ourselves
      thisN = eval(n[rowI])
      
      crpsBinomialSingle = function(colI) {
        thisP = p.est[rowI, colI]
        sum((pbinom(toN, thisN, thisP) - (thisCount <= toN))^2)
      }
      
      # return the region CRPS, taking the expectation over the distribution of p
      mean(unlist(sapply(1:ncol(p.est), crpsBinomialSingle)))
    } else
      sum((cdfs[[rowI]] - (thisCount <= toN))^2)
  }
  
  if(is.null(parClust))
    crpsByRegion = lapply(1:length(trueProportion), function(i) {rowFun(i)})
  else
    crpsByRegion = parLapply(parClust, 1:length(trueProportion), function(i) {rowFun(i)})
  
  if(empiricalProportion)
    mean(unlist(crpsByRegion) / n) # multiply by delta p for the discrete integral
  else
    mean(unlist(crpsByRegion))
}

# for a discrete random variable, determine the CRPS. 
# probs: a vector of probabilities of length n + 1
crpsCounts <- function(trueCount, probs, parClust=NULL, empiricalProportion=TRUE) {
  n = length(probs) - 1
  toN = eval(0:n)
  
  # get cdf
  cdf = cumsum(probs)
  
  if(is.null(parClust))
    res = sum((cdf - (trueCount <= toN))^2)
  else {
    res = parSapply(parClust, toN, function(x) {trueCount <= x})
    res = parSapply(parClust, 1:(n + 1), function(x) {(cdf[x] - res[x])^2})
    res = sum(res)
  }
  
  if(empiricalProportion)
    res / n
  else
    res
}

# for a discrete random variable, determine the CRPS. 
# probs: a vector of probabilities of length n + 1
crpsCountsNoBVar <- function(trueProportion, logitEst, logitVar, parClust=NULL, empiricalProportion=TRUE) {
  n = length(probs) - 1
  toN = eval(0:n)
  
  # get cdf
  cdf = cumsum(probs)
  
  if(is.null(parClust))
    res = sum((cdf - (trueCount <= toN))^2)
  else {
    res = parSapply(parClust, toN, function(x) {trueCount <= x})
    res = parSapply(parClust, 1:(n + 1), function(x) {(cdf[x] - res[x])^2})
    res = sum(res)
  }
  
  if(empiricalProportion)
    res / n
  else
    res
}
# truth: a vector of observations on the desired scale
# my.est: a vector of logit-scale predictions of the same length as truth 
# my.var: a vector of logit-scale predictive variances of the same length as truth
# logit: whether CRPS should be computed on the logit scale or the empirical proportion scale
# nsim: if binomial variation is included, the number of probabilities over which to integrate in the pmf approximation
# n: the binomial sample size. Can be a vector of the same length as truth
# bVar: whether or not binomial variation should be included in the predictive distribution-
getCPO <- function(truth, my.est, my.var=NULL, nsim=10, n=25, probMat=NULL, bVar=TRUE){
  if(bVar) {
    # first draw different values of p from the predictive distribution if necessary
    if(is.null(probMat)) {
      quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(my.est)), ncol=nsim, byrow=TRUE)
      logitP = matrix(qnorm(quantiles, my.est, sd=sqrt(my.var)), ncol=nsim)
      probMat = expit(logitP)
    }
    
    # for each value of p, calculate the crps and return the mean
    out = cpoBinomial(truth, probMat, n=n)
    cpoVals = out$CPO
    logScore = out$logScore
    PIT = out$PITs
  }
  else {
    stop("getCPO does not yet support not including binomial variation.  bVar must be set to TRUE")
  }
  
  list(CPO=cpoVals, logScore=logScore, PITs=PIT)
}

# for a set of independent fixed binomials, determine the CRPS.
# trueProportion: can be a vector
# p.est: can be a matrix, with number of rows equal to length(trueProportion). In that case, the 
#        number of columns can be thought of as the number of draws of p in the case of a random p.
# n: can be a vector, with length equal to length(trueProportion)
# empiricalProportion: if true, calculate the CRPS of the empirical proportion distribution on [0, 1] 
#                      rather than the binomial distribution on [0, n]
cpoBinomial = function(trueProportion, p.est=NULL, n=25, parClust=NULL, empiricalProportion=TRUE) {
  # get observations at the count scale, make sure n is a vector
  trueCounts = trueProportion * n
  if(length(n) == 1)
    n = rep(n, length(trueProportion))
  
  p.est = matrix(p.est, nrow=length(trueProportion))
  
  getPMF = function(rowI) {
    ps = p.est[rowI,]
    thisN = n[rowI]
    
    rowMeans(sapply(ps, function(p) {dbinom(0:thisN, thisN, p)}))
  }
  pmfs = lapply(1:nrow(p.est), getPMF)
  cdfs = lapply(pmfs, cumsum)
  PITs = lapply(1:length(cdfs), function(i) {
    cdf = cdfs[[i]]
    trueCount = trueCounts[i]
    cdf[trueCount+1]
  })
  
  rowFunCPO = function(rowI) {
    # thisCount = eval(trueCounts[rowI])
    # 
    # ## in this case we must compute the probability mass function for this region ourselves
    # thisN = eval(n[rowI])
    # 
    # cpoBinomialSingle = function(colI) {
    #   thisP = p.est[rowI, colI]
    #   dbinom(thisCount, thisN, thisP)
    # }
    # 
    # # return the region CRPS, taking the expectation over the distribution of p
    # mean(unlist(sapply(1:ncol(p.est), cpoBinomialSingle)))
    
    pmfs[[rowI]][trueCounts[rowI] + 1]
  }
  
  if(is.null(parClust))
    cpoByRegion = lapply(1:length(trueProportion), function(i) {rowFunCPO(i)})
  else
    cpoByRegion = parLapply(parClust, 1:length(trueProportion), function(i) {rowFunCPO(i)})
  
  cpoByRegion= unlist(cpoByRegion)
  
  logScore = mean(log(cpoByRegion))
  
  list(CPO = mean(cpoByRegion), logScore=logScore, PITs=PITs)
}

dss = function(truth, my.est, my.var){
  
  my.sd = sqrt(my.var)
  # proportional to logarithmic score assuming normality
  dss_score <- (2*log(my.sd) + ((truth - my.est)/my.sd)^2)
  
  return(mean(dss_score))
}

# either include both lower and upper, or include either: 
#    - the probability draw matrix
#    - estimates and variances on the logit scale
# nsim: the number of draws from the distribution of p in the case that the logit estimates and variance are provided
#       to estimate the probability mass function
# truth: the true empirical proportions of mortality rates within the regions or enumeration areas of interest
# lower: the logit scale lower end of the credible interval
# upper: the logit scale upper end of the credible interval
# numChildren: a vector containing the number of children with the same length as truth
# logitEst: the predictive mean on the logit scale, of the same length as truth
# logitVar: the predictive variance on the logit scale, of the same length as truth
# probMat: a matrix of joint draws of probability values, with number of rows equal to the length of truth, a number of 
#          columns equal to the number of draws. If not included, a gaussian distribution is assumed on the logit scale.
# significance: the significance level of the credible interval. By default 80%
# pmfs: a list of pmfs (vectors of probabilities) with length the same as truth
coverage = function(truth, lower=NULL, upper=NULL, lowerRejectProb=NULL, upperRejectProb=NULL, doLogit = TRUE, 
                    bVar=FALSE, numChildren=NULL, probMat=NULL, logitEst=NULL, logitVar=NULL, pmfs=NULL, 
                    nsim=10, significance=.8, returnIntervalWidth=FALSE, returnPMFs=FALSE){
  
  if(any(is.null(lower)) || any(is.null(upper))) {
    # if the user did not supply their own credible intervals, we must get them ourselves given the other information
    
    if(is.null(probMat) && (is.null(logitEst) || is.null(logitVar)) && is.null(pmfs))
      stop("either include both lower and upper, probMat, pmfs, or both logitEst and logitVar")
    
    if(!is.null(logitEst) && !is.null(logitVar) && bVar) {
      # in this case, we must generate the probabilities over which we will later integrate to get the 
      # pmfs and credible intervals
      
      ## first we estimate the probability mass function
      # get the probabilities over which we will integrate to get the probability mass function
      # quantiles = matrix(rep(seq(0, 1, l=nsim + 2)[c(-1, -(nsim + 2))], length(logitEst)), ncol=nsim, byrow=TRUE)
      quantiles = matrix(rep(seq(0, 1 - 1 / nsim, by = 1/nsim) + 1 / (2 * nsim), length(logitEst)), ncol=nsim, byrow=TRUE)
      logitP = matrix(qnorm(quantiles, logitEst, sd=sqrt(logitVar)), ncol=nsim)
      probMat = expit(logitP)
    }
    else if(!is.null(logitEst) && !is.null(logitVar) && !bVar) {
      # We don't want to include binomial variation, and we have the gaussian distribution parameters on the logit scale. 
      # Use those parameters to generate the credible intervals
      
      lower = qnorm((1 - significance) / 2, logitEst, sqrt(logitVar))
      upper = qnorm(1 - (1 - significance) / 2, logitEst, sqrt(logitVar))
    }
    else {
      # we don't have information about the predictive distribution on the logit scale, and don't assume normality. 
      # Instead, use the user supplied to probability matrix probMat
      
      if(!bVar) {
        # we aren't accounting for binomial variation, so just take the quantiles of the probability draws
        CIs = logit(apply(probMat, 1, function(ps) {quantile(ps, probs=c((1 - significance) / 2, 1 - (1 - significance) / 2))}))
        lower = CIs[1,]
        upper = CIs[2,]
      }
      else {
        # don't need to do anything here, since it will be already done to probMat in the following if statement
      }
    }
    
    if(bVar) {
      # for each region, take the average binomial probability mass, averaging over that drawn probabilities (integrate)
      # do this for all possible count values to get the probability mass for each count
      if(is.null(pmfs)) {
        getPMF = function(rowI) {
          ps = probMat[rowI,]
          n = numChildren[rowI]
          
          rowMeans(sapply(ps, function(p) {dbinom(0:n, n, p)}))
        }
        pmfs = lapply(1:nrow(probMat), getPMF)
      }
      
      ## now that we have the pmfs, we can calculate the credible intervals (on the logit rather than count scale)
      # NOTE: it is fine if we have a 0 or 1 interval boundary, since if we take the logit and get infinite results, 
      #       comparisons still work
      # intervals = sapply(pmfs, getHDI, significance=significance)
      intervals = sapply(pmfs, getQuantileInterval, significance=significance)
      lower = logit(intervals[1,] / numChildren)
      upper = logit(intervals[2,] / numChildren)
      lowerRejectProb = intervals[3,]
      upperRejectProb = intervals[4,]
    }
  }
  
  # if(!doLogit){
  lower = expit(lower)
  upper = expit(upper)
  # truth = expit(truth)
  # }
  
  if(any(lower >= upper)) {
    warning("lower >= upper, reordering")
    tmp = lower
    wrongOrder = lower >= upper
    lower[wrongOrder] = upper[wrongOrder]
    upper[wrongOrder] = tmp[wrongOrder]
  }
  
  res = lower <= truth & upper >= truth
  
  if(returnIntervalWidth)
    width = upper - lower
  
  if(bVar) {
    if(any(is.null(lowerRejectProb)) || any(is.null(upperRejectProb)))
      stop("binomial variation included, but no coverage rejection probabilities included")
    
    # in the case of discrete credible intervals, reject if the truth is at the boundary with some probability
    atLower = truth == lower
    atUpper = truth == upper
    rejectLower = runif(sum(atLower)) < lowerRejectProb[atLower]
    rejectUpper = runif(sum(atUpper)) < upperRejectProb[atUpper]
    res[atLower] = res[atLower] - rejectLower
    res[atUpper] = res[atUpper] - rejectUpper
    
    # adjust interval width for random rejection if necessary
    if(returnIntervalWidth) {
      width[atLower] = width[atLower] - rejectLower / numChildren[atLower]
      width[atUpper] = width[atUpper] - rejectUpper / numChildren[atUpper]
    }
  }
  
  allResults = c(coverage=mean(res))
  if(returnIntervalWidth)
    allResults = c(allResults, width=mean(width))
  if(returnPMFs)
    allResults = list(intervalInfo=allResults, pmfs=pmfs)
  
  allResults
}

intervalWidth = function(lower, upper, logit = TRUE) {
  if(!logit){
    lower = expit(lower)
    upper = expit(upper)
  }
  
  res = mean(abs(upper - lower))
  if(any(lower > upper))
    warning("some interval widths are negative, taking absolute value")
  res
}

# taken from logitnorm package.  Calculates the mean of a distribution whose 
# logit is Gaussian. Each row of muSigmaMat has a mean and standard deviation 
# on the logit scale
logitNormMean = function(muSigmaMat, parClust=NULL, ...) {
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust))
      apply(muSigmaMat, 1, logitNormMean, ...)
    else
      parApply(parClust, muSigmaMat, 1, logitNormMean, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    if(sigma == 0)
      expit(mu)
    else {
      fExp <- function(x) exp(plogis(x, log.p=TRUE) + dnorm(x, mean = mu, sd = sigma, log=TRUE))
      integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
    }
  }
}

# calculate expectation and variance of the logit-normal distribution
# code based on logitnorm package
logitNormMoments = function(muSigmaMat, parClust=NULL, ...) {
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust))
      t(apply(muSigmaMat, 1, logitNormMoments, ...))
    else
      t(parApply(parClust, muSigmaMat, 1, logitNormMoments, ...))
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    fExp <- function(x) exp(plogis(x, log.p=TRUE) + dnorm(x, mean = mu, sd = sigma, log=TRUE))
    expectation = integrate(fExp, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
    
    fVar = function(x) (plogis(x) - expectation)^2 * dnorm(x, mean = mu, sd = sigma)
    variance = integrate(fVar, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
    c(mean = expectation, var = variance)
  }
}

# calculates the E[UV] for U and V bernoulli random variables with the given distribution
# of p on the logit scale
logitNormCrossExpectation = function(muSigmaMat=matrix(c(0, 1), ncol=2), ...) {
  
  # now calculate the expectation of the product of two of the random variables
  if(length(muSigmaMat) > 2) {
    if(is.null(parClust))
      apply(muSigmaMat, 1, logitNormMean, ...)
    else
      parApply(parClust, muSigmaMat, 1, logitNormMean, ...)
  }
  else {
    mu = muSigmaMat[1]
    sigma = muSigmaMat[2]
    f11 <- function(x) exp(2 * plogis(x, log.p = TRUE) + dnorm(x, mean = mu, sd = sigma, log = TRUE))
    p11 = integrate(f11, mu-10*sigma, mu+10*sigma, abs.tol = 0, ...)$value
    
    p11
  }
}

# calculates the cor(U, V) for U and V bernoulli random variables with the given distribution
# of p on the logit scale
logitNormCor = function(muSigmaMat=matrix(c(0, 1), ncol=2)) {
  # first calculate marginal expectations
  expectedProbs = logitNormMean(muSigmaMat)
  
  # calculate the cross expectations
  crossExpectations = logitNormCrossExpectation(muSigmaMat)
  
  # now calculate the correlations
  numerator = crossExpectations - expectedProbs^2
  denominator = expectedProbs * (1 - expectedProbs)
  numerator / denominator
}

# generate a binomial confidence interval with at least the significance requested. The
# interval is generated by starting at the maximum probability count and expanding 
# outward to the highest probability counts.
# NOTE: the return confidence interval is on counts, not the index of the probs vector
generateBinomialInterval = function(probs, significance = .80) {
  if(any(is.na(probs)))
    return(list(lower=NA, upper=NA, leftRejectProb=NA, rightRejectProb=NA))
  # Find the maximum probability, end expand outward from there
  maxI = which.max(probs)
  leftI = maxI
  rightI = maxI
  
  expandInterval = function(leftI, rightI) {
    if(length(leftI) > 1 || length(rightI) > 1) {
      print("woah there")
    }
    currentProb = sum(probs[leftI:rightI])
    
    if(currentProb >= significance)
      return(list(leftI = leftI, rightI = rightI))
    
    if(leftI == 1)
      leftProb = 0
    else
      leftProb = probs[leftI - 1]
    if(rightI == length(probs))
      rightProb = 0
    else
      rightProb = probs[rightI + 1]
    
    if(rightProb == leftProb) {
      # in this case, expand symmetrically to the mode, or randomly if distances from mode are the same
      modeI = which.max(probs)
      leftDist = modeI - leftI
      rightDist = rightI - modeI
      if(leftDist == rightDist) {
        if(runif(1) >= .5)
          return(list(leftI = leftI, rightI = rightI + 1))
        else
          return(list(leftI = leftI - 1, rightI = rightI))
      } else if(leftDist > rightDist)
        return(list(leftI = leftI, rightI = rightI + 1))
      else
        return(list(leftI = leftI - 1, rightI = rightI))
    } else if(rightProb >= leftProb)
      return(list(leftI = leftI, rightI = rightI + 1))
    else
      return(list(leftI = leftI - 1, rightI = rightI))
  }
  
  # keep expanding the interval until were done (i doesn't matter here)
  for(i in 1:length(probs)) {
    out = expandInterval(leftI, rightI)
    if(leftI == out$leftI && rightI == out$rightI)
      break
    
    leftI = out$leftI
    rightI = out$rightI
  }
  
  # determine probability of rejection at boundary
  totalProb = sum(probs[leftI:rightI])
  extraProb = totalProb - significance
  leftProb = probs[leftI]
  rightProb = probs[rightI]
  if(leftProb <= rightProb) {
    leftRejectProb = extraProb
    rightRejectProb = 0
  }
  else {
    rightRejectProb = extraProb
    leftRejectProb = 0
  }
  
  # return results (subtract one to the index to get the count)
  c(lower=leftI - 1, upper=rightI - 1, leftRejectProb=leftRejectProb, rightRejectProb=rightRejectProb)
}

# generate a highest density interval with at least the significance requested
# NOTE: the returned confidence interval is on counts, not the index of the probs vector
getHDI = function(probs, significance = .80) {
  if(any(is.na(probs)))
    return(list(lower=NA, upper=NA, leftRejectProb=NA, rightRejectProb=NA))
  
  # Find the maximum probability, and expand outward from there
  maxI = which.max(probs)
  n = length(probs) - 1
  cumulativeProbs = cumsum(probs)
  
  getIntervalWidth = function(leftI) {
    # totalProbs = probs[leftI] + cumulativeProbs[leftI:(n + 1)] - cumulativeProbs[leftI]
    # width = binarySearchMatch(1, totalProbs >= significance)
    width = binarySearchMatch(1, cumulativeProbs[leftI:(n + 1)], function(x, i) {(x[i] - x[1] + probs[leftI]) >= significance})
    rightI = leftI - 1 + width
    c(width = width, leftI = leftI, rightI = rightI, prob=cumulativeProbs[rightI] - cumulativeProbs[leftI] + probs[leftI])
  }
  
  # calculate interval widths for relevant starting points
  maxRelevantI = binarySearchMatch(1, cumulativeProbs >= 1 - significance)
  intervalMatrix = sapply(1:maxRelevantI, getIntervalWidth)
  
  # determine which interval is the highest density interval, and get its information
  # NOTE: there may be multiple highest density intervals, since the distribution is discrete. 
  #       Taking any of them is fine.
  minI = which.min(intervalMatrix[1,])
  leftI = intervalMatrix[2, minI]
  rightI = intervalMatrix[3, minI]
  totalProb = intervalMatrix[4, minI]
  
  # determine probability of rejection at boundary
  extraProb = totalProb - significance
  leftProb = probs[leftI]
  rightProb = probs[rightI]
  if(leftProb <= rightProb) {
    leftRejectProb = extraProb
    rightRejectProb = 0
  }
  else {
    rightRejectProb = extraProb
    leftRejectProb = 0
  }
  
  # return results (subtract one to the index to get the count)
  c(lower=leftI - 1, upper=rightI - 1, leftRejectProb=leftRejectProb, rightRejectProb=rightRejectProb)
}

# generate a credible interval with at least the significance requested by including utmost (1 - significance) / 2
# probability mass on each side
# NOTE: the returned confidence interval is on counts, not the index of the probs vector
getQuantileInterval = function(probs, significance = .80) {
  if(any(is.na(probs)))
    return(list(lower=NA, upper=NA, leftRejectProb=NA, rightRejectProb=NA))
  
  # get the cdf
  n = length(probs) - 1
  cumulativeProbs = cumsum(probs)
  
  # get left and right endpoints
  leftI = max(binarySearchMatch(1, cumulativeProbs >= (1 - significance) / 2) - 1, 1)
  probLeft = cumulativeProbs[leftI] - probs[leftI]
  extraProbLeft = (1 - significance) / 2 - probLeft
  rightI = max(binarySearchMatch(1, cumulativeProbs >= 1 - (1 - significance) / 2), 1)
  probRight = 1 - cumulativeProbs[leftI]
  extraProbRight = (1 - significance) / 2 - probRight
  
  # determine probability of rejection at boundaries
  leftRejectProb = extraProbLeft
  rightRejectProb = extraProbRight
  
  # return results (subtract one to the index to get the count)
  c(lower=leftI - 1, upper=rightI - 1, leftRejectProb=leftRejectProb, rightRejectProb=rightRejectProb)
}

# approximate the sum of binomial distributions using the method from:
# https://stackoverflow.com/questions/15926448/approximate-the-distribution-of-a-sum-of-binomial-random-variables-in-r
# http://www.dtic.mil/dtic/tr/fulltext/u2/a266969.pdf
# use the first 4 moments of the sum, fits with a pearson distribution
dSumBinom = function(k, ns=25, ps=.5) {
  if(length(ns) == 1)
    ns = rep(ns, length(ps))
  if(length(ps) == 1)
    ps = rep(ps, length(ns))
  
  # get the first 4 cumulants
  k.1<-sum(ns*ps)
  k.2<-sum(ns*ps*(1-ps))
  k.3<-sum(ns*ps*(1-ps)*(1-2*ps))
  k.4<-sum(ns*ps*(1-ps)*(1-6*ps*(1-ps)))
  
  # obtain skewness and excess kurtosis
  beta.1<-k.3^2/k.2^3
  beta.2<-k.4/k.2^2
  
  # get the moments, and return the probability mass
  # (by integrating the density of the continuous pearson distribution)
  moments <- c(mean=k.1,variance=k.2,skewness=sqrt(beta.1),kurtosis=beta.2 + 3)
  if(!identical(k, 0:sum(ns))) {
    results = k
    zeros = k == 0
    results[zeros] = ppearson(.5,moments=moments)
    results[!zeros] = ppearson(k[!zeros] + .5,moments=moments) - ppearson(k[!zeros] - .5,moments=moments)
    results
  } else {
    # the commented out code is technically correct, but the code below is a good approximation that saves time
    # results = ppearson(k + .5,moments=moments)
    # results[2:length(results)] = diff(results)
    # results
    results = dpearson(k, moments=moments)
    results * (1 / sum(results))
  }
}

# approximate the sum of binomial distributions using the method from:
# https://stackoverflow.com/questions/15926448/approximate-the-distribution-of-a-sum-of-binomial-random-variables-in-r
# http://www.dtic.mil/dtic/tr/fulltext/u2/a266969.pdf
# use the first 4 moments of the sum, fits with a pearson distribution
# normalize: if true, normalize the returned densities by their sum (use this when approximating full pmf)
dPearsonMoments = function(k, moments, normalize=TRUE) {
  # moments <- c(mean=k.1,variance=k.2,skewness=sqrt(beta.1),kurtosis=beta.2 + 3)
  
  results = dpearson(k, moments=moments)
  if(normalize)
    results * (1 / sum(results))
  results
}

# a previous version of the following function that is slightly less efficient for computing full pmfs
dSumBinomRandom = function(k, ns=25, pMat, parClust=NULL) {
  if(is.null(parClust) || length(ns) <= 100 || ncol(pMat) <= 30)
    dSumBinoms = matrix(apply(pMat, 2, dSumBinom, k=k, ns=ns), ncol=ncol(pMat))
  else
    dSumBinoms = matrix(parApply(parClust, pMat, 2, dSumBinom, k=k, ns=ns), ncol=ncol(pMat))
  rowMeans(dSumBinoms)
}

# approximate the sum of binomials where the probabilities or sampled from a random distribution
# using a Pearson distribution
# pMat: a matrix of probabilities, the number of rows equaling the length of ns, and the number 
#       of columns corresponding to the number of samples of p for each discrete distribution
dSumBinomRandom2 = function(k, ns=25, pMat, parClust=NULL, normalize=TRUE) {
  # get the first 4 cumulants
  temp = sweep(pMat, 1, ns, "*")
  k.1s = colSums(temp)
  temp = temp*(1 - pMat)
  k.2s = colSums(temp)
  temp2 = temp * (1 - 2 * pMat)
  k.3s = colSums(temp2)
  temp = temp * (1 - 6 * pMat * (1 - pMat))
  k.4s = colSums(temp)
  # k.1<-sum(ns*ps)
  # k.2<-sum(ns*ps*(1-ps))
  # k.3<-sum(ns*ps*(1-ps)*(1-2*ps))
  # k.4<-sum(ns*ps*(1-ps)*(1-6*ps*(1-ps)))
  
  # obtain skewness and excess kurtosis
  beta.1s = k.3s^2 / k.2s^3
  beta.2s = k.4s / k.2s^2
  # beta.1<-k.3^2/k.2^3
  # beta.2<-k.4/k.2^2
  
  # get the moments, and return the probability mass
  # (by integrating the density of the continuous pearson distribution)
  moments <- rbind(mean=k.1s,variance=k.2s,skewness=sqrt(beta.1s),kurtosis=beta.2s + 3)
  # moments <- c(mean=k.1,variance=k.2,skewness=sqrt(beta.1),kurtosis=beta.2 + 3)
  
  # if(is.null(parClust) || length(ns) <= 100 || ncol(pMat) <= 30)
  #   dSumBinoms = matrix(apply(pMat, 2, dSumBinom, k=k, ns=ns), ncol=ncol(pMat))
  # else
  #   dSumBinoms = matrix(parApply(parClust, pMat, 2, dSumBinom, k=k, ns=ns), ncol=ncol(pMat))
  
  dSumBinoms = apply(moments, 2, dPearsonMoments, k=k, normalize=normalize)
  rowMeans(dSumBinoms)
}

# approximate the sum of binomials where the probabilities or sampled from a random distribution 
# using simulations from binomial distributions
# pMat: a matrix of probabilities, the number of rows equaling the length of ns, and the number 
#       of columns corresponding to the number of samples of p for each discrete distribution
dSumBinomRandomSim = function(k, ns=25, pMat, nSim = 1000) {
  dSumBinoms = matrix(apply(pMat, 2, dSumBinomSim, k=k, ns=ns, nSim=nSim), ncol=ncol(pMat))
  rowMeans(dSumBinoms)
}

# same as dSumBinom except we used simulations from binomial distributions conditional on p 
# to approximate the distribution rather than using the pearson distribution approximation
dSumBinomSim = function(k, ns=25, ps=.5, nSim=1000) {
  if(length(ns) == 1)
    ns = rep(ns, length(ps))
  if(length(ps) == 1)
    ps = rep(ps, length(ns))
  
  countSims = matrix(rbinom(length(ns) * nSim, ns, ps), nrow=length(ns), ncol=nSim)
  cdf = ecdf(colSums(countSims))
  cdf(k) - cdf(k - 1)
}

# function for adding on binomial variation to simulation draws
# logitProbMat: matrix of simulated logit probability draws with nrows=length(nsList)=number of regions 
#               and ncols=length(clustVars)=number of simulations
# nsList: If ms is NULL, a list of vectors of length equal to the number of regions, where the length of each vector 
#         corresponds to the number of EAs in that region, and the values within each vector are the 
#         number of children in each EA. If ms is not NULL, this is a vector of assumed number of children
#         per EA
# clustVars: a vector of length equal to ncol(logitProbMat), which is the number of probability draws 
#            from the predictive distribution. Each element is the cluster variance from that predictive 
#            distributions draw
# 
simBinVar = function(logitProbMat, nsList, clustVars=rep(0, nrow(logitProbMat)), ms=NULL) {
  ns = matrix(rep(ns, ncol(logitProbMat)), ncol=ncol(logitProbMat))
  clustVars = matrix(rep(clustVars, nrow(logitProbMat)), nrow=nrow(logitProbMat))
  mapply(simBinVarHelper, logitProb=logitProbMat, ns=ns, clustVar=clustVars)
}

# function for adding on binomial variation to simulation draws
simBinVarHelper = function(logitProb, ns, clustVar=0) {
  if(clustVar == 0)
    logit(rbinom(1, expit(logitProb), sum(ns)) / sum(ns))
  else
    logit(sum(rbinom(length(ns), expit(logitProb + rnorm(length(ns), sd=sqrt(clustVar))), ns)) / sum(ns))
}







