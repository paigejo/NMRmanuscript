##################################################################################
# designBased.R                                                                  #
#    Produce county level estimates using pixel-wise spatial model with          #
#    aggregation using true population and BYM spatial model                     #
##################################################################################

# same as runBYM, except fits the BYM2 reparameterized model (iid component in the BYM2 is that the county level)
runBYM2 = function(tausq=0.1^2, test=FALSE, includeUrbanRural=TRUE, includeCluster=TRUE, maxDataSets=NULL, 
                   aggregateByPopulation=FALSE, margVar=0.15^2, gamma=-1, plotPriorPost=FALSE, strictPrior=FALSE, 
                   seed=NULL, effRange=150) {
  if(!is.null(seed))
    set.seed(seed)
  
  # load and relevant data
  rangeText = ifelse(effRange == 150, "", "Range50")
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0", rangeText, ".RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0Test", rangeText, ".RData"))
  
  # Get true ratios of urban/rural
  urbRatio = vector('numeric', length = 47)
  counties = sort(unique(as.character(SRSDat$eaDat$admin1)))
  if(!aggregateByPopulation) {
    for(i in 1:47){
      # compute the exact proportion of children that are in urban areas
      idx = which(as.numeric(factor(SRSDat$eaDat$admin1)) == i)
      urbanI = idx[SRSDat$eaDat$urban[idx]]
      urbRatio[i] = sum(SRSDat$eaDat$numChildren[urbanI])/sum(SRSDat$eaDat$numChildren[idx])
    }
    urbRatioAdjusted = urbRatio
  } else {
    
    # now compute the proportion of the population that is urban (possibly adjusting for 
    # the number of children per unit population being different in urban and rural areas)
    getUrbanRatios = function(aggregationWeightTable) {
      urbRatio = aggregationWeightTable$popUrb / aggregationWeightTable$popTotal
      sortI = matchMultiple(counties, aggregationWeightTable$County)
      urbRatio = urbRatio[sortI]
    }
    urbRatio = getUrbanRatios(poppc)
    urbRatioAdjusted = getUrbanRatios(adjustPopulationPerCountyTable())
  }
  
  sortI = matchMultiple(counties, easpc$County)
  clustersPerUrban = easpc$EAUrb[sortI]
  clustersPerRural = easpc$EARur[sortI]
  clustersPerCounty = rowSums(cbind(clustersPerUrban, clustersPerRural))
  
  # Define formula
  if(!strictPrior)
    hyperList = list(param=c(1, 0.01), prior="pc.prec")
  else
    hyperList = list(param=c(.15, 0.01), prior="pc.prec")
  if(!strictPrior)
    clusterList = list(param=c(1, 0.01), prior="pc.prec")
  else
    clusterList = list(param=c(.15, 0.01), prior="pc.prec")
  if(includeUrbanRural) {
    if(includeCluster) {
      formula = y ~ urban +
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=hyperList, phi=list(param=c(0.5, 2/3), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = clusterList))
    } else {
      formula = y ~ urban + 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=hyperList, phi=list(param=c(0.5, 2/3), prior="pc")))
    }
  } else {
    if(includeCluster) {
      formula = y ~ f(idx, model="bym2",
                      graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
                      hyper=list(prec=hyperList, 
                                 phi=list(param=c(0.5, 2/3), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = clusterList))
    } else {
      formula = y ~ 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=clusterList, phi=list(param=c(0.5, 2/3), prior="pc")))
    }
  }
  
  
  # Number of simulations for producing results
  Nsim = 10000
  
  # Help functions
  logit = function(x){
    return(log(x/(1-x)))
  }
  expit = function(x){
    return(1/(1+exp(-x)))
  }
  
  if(includeCluster)
    parNames = c("Intercept", "Urban", "Cluster Var", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var", "Cluster SD", "BYM2 Tot. SD", "BYM2 Spatial SD", "BYM2 iid SD")
  else
    parNames = c("Intercept", "Urban", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var", "BYM2 Tot. SD", "BYM2 Spatial SD", "BYM2 iid SD")
  includeI = c(1, rep(2, includeUrbanRural), 3:length(parNames))
  parNames = parNames[includeI]
  
  # set the maximum number of data sets to analyze
  if(is.null(maxDataSets)) {
    maxDataSets = length(SRSDat$clustDat)
  }
  
  # Go through datasets for SRSDat
  sampCountySRSDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountySRSDatAdjusted = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountySRSDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  sampClusterSRSDatUrban = array(NA, dim = c(47, Nsim, maxDataSets))
  sampClusterSRSDatRural = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelSRSDatUrban = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelSRSDatRural = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelSRSDatUrbanMod = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelSRSDatRuralMod = array(NA, dim = c(47, Nsim, maxDataSets))
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountySRSDatPar = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountySRSDatSD = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountySRSDat10 = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountySRSDat50 = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountySRSDat90 = matrix(NA, length(parNames), ncol=maxDataSets)
  rownames(sampCountySRSDatPar) = parNames
  rownames(sampCountySRSDatSD) = parNames
  rownames(sampCountySRSDat10) = parNames
  rownames(sampCountySRSDat50) = parNames
  rownames(sampCountySRSDat90) = parNames
  for(i in 1:maxDataSets){
    # Extract data
    currData = SRSDat$clustDat[[i]]
    currData$admin1 = factor(currData$admin1)
    
    # INLA data
    dat = list(y = currData$died,
               Ntrials = currData$numChildren,
               urban = currData$urban,
               idx = as.numeric(currData$admin1),
               idxEps = 1:length(currData$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
    dat$urban = c(rep(c(1,0), each = 47), dat$urban)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
    
    # Run model
    print(paste0("fitting SRS model for dataset ", i, "/", maxDataSets))
    result = inla(formula = formula, 
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat, 
                  control.compute = list(config = TRUE), 
                  quantiles=c(0.1, 0.5, 0.9))
    
    if(plotPriorPost) {
      # maxX = max(result$marginals.hyperpar[[1]][,1]) * 1.1
      maxX = inla.qmarginal(0.8, result$marginals.hyperpar[[1]])
      xs = exp(seq(-5, log(maxX), l=1500))
      ys = inla.pc.dprec(xs, 1, 0.01)
      maxY = max(max(ys), max(result$marginals.hyperpar[[1]][,2]))
      plot(xs, ys, type="l", col="blue", xlab="BYM2 Precision", ylab="Density", xlim=c(0, maxX), 
           ylim=c(0, maxY), main="BYM2 Precision Prior vs Posterior")
      lines(result$marginals.hyperpar[[1]], col="red")
      legend("topright", c("Prior", "Posterior"), col=c("blue", "red"), lty=1)
    }
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountySRSDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountySRSDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountySRSDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountySRSDat50[1:(1 + includeUrbanRural), i] = result$summary.fixed[,4]
    sampCountySRSDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountySRSDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountySRSDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountySRSDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountySRSDat50[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountySRSDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    # BYM2 hyperparameter phi
    sampCountySRSDatPar[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,1]
    sampCountySRSDatSD[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,2]
    sampCountySRSDat10[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,3]
    sampCountySRSDat50[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,3]
    sampCountySRSDat90[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
    transformFunction = function(x) {c(1/x[3], x[2], 1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[3]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
    if(!includeCluster)
      transformFunction = function(x) {c(1/x[1], x[2], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
    transformedOut = apply(out, 1, transformFunction)
    
    # now calculate the summary statistics of the transformed BYM2 hyperparameters
    sampCountySRSDatPar[(2 + includeUrbanRural):length(parNames), i] = rowMeans(transformedOut)
    sampCountySRSDatSD[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, sd)
    sampCountySRSDat10[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, quantile, probs=.1)
    sampCountySRSDat50[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, quantile, probs=.5)
    sampCountySRSDat90[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, quantile, probs=.9)
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
      
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          # muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
          # muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
          # sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
          # sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
          sampRuralMod[, j] = sapply(1:length(clustersPerRural), function(i) {mean(expit(sampRural[i, j] + rnorm(clustersPerRural[i], sd=clusterSigma)))})
          sampRuralMod[clustersPerRural == 0,] = 0 # we just need to set these values to something other than NaN
          sampUrbanMod[, j] = sapply(1:length(clustersPerUrban), function(i) {mean(expit(sampUrban[i, j] + rnorm(clustersPerUrban[i], sd=clusterSigma)))})
          sampClusterSRSDatRural[,j,i] = sampRural[, j] + rnorm(47, sd=clusterSigma)
          sampClusterSRSDatUrban[,j,i] = sampUrban[, j] + rnorm(47, sd=clusterSigma)
        }
      }
      sampCountySRSDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
      sampCountySRSDatAdjusted[,,i] = logit(expit(sampUrban)*urbRatioAdjusted + expit(sampRural)*(1-urbRatioAdjusted))
      sampCountySRSDatMod[,,i] = logit(sampUrbanMod*urbRatioAdjusted + sampRuralMod*(1-urbRatioAdjusted))
      if(!includeCluster) {
        sampClusterSRSDatRural[,,i] = sampRural
        sampClusterSRSDatUrban[,,i] = sampUrban
      } else {
        sampPixelSRSDatUrbanMod[,,i] = sampUrbanMod
        sampPixelSRSDatRuralMod[,,i] = sampRuralMod
      }
      sampPixelSRSDatUrban[,,i] = sampUrban
      sampPixelSRSDatRural[,,i] = sampRural
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
      sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampCounty[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          # muSigmaMat = cbind(sampCounty[, j], clusterSigma)
          # sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat
          sampCountyMod[, j] = sapply(1:length(clustersPerCounty), function(i) {mean(expit(sampCounty[i, j] + rnorm(clustersPerCounty[i], sd=clusterSigma)))})
          sampClusterSRSDatRural[,j,i] = expit(sampCounty[, j] + rnorm(47, sd=clusterSigma))
          sampClusterSRSDatUrban[,j,i] = sampClusterSRSDatRural[,j,i]
        }
      }
      sampCountySRSDat[,,i] = sampCounty
      sampCountySRSDatAdjusted[,,i] = sampCounty
      sampCountySRSDatMod[,,i] = logit(sampCountyMod)
      sampPixelSRSDatUrban[,,i] = sampCounty
      sampPixelSRSDatRural[,,i] = sampCounty
      if(!includeCluster) {
        sampClusterSRSDatRural[,,i] = sampCounty
        sampClusterSRSDatUrban[,,i] = sampCounty
      } else {
        sampPixelSRSDatUrbanMod[,,i] = sampCountyMod
        sampPixelSRSDatRuralMod[,,i] = sampCountyMod
      }
    }
  }
  
  # Go through datasets for overSampDat
  sampCountyOverSampDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountyOverSampDatAdjusted = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountyOverSampDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  sampClusterOverSampDatUrban = array(NA, dim = c(47, Nsim, maxDataSets))
  sampClusterOverSampDatRural = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelOverSampDatUrban = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelOverSampDatRural = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelOverSampDatUrbanMod = array(NA, dim = c(47, Nsim, maxDataSets))
  sampPixelOverSampDatRuralMod = array(NA, dim = c(47, Nsim, maxDataSets))
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountyOverSampDatPar = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountyOverSampDatSD = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountyOverSampDat10 = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountyOverSampDat50 = matrix(NA, length(parNames), ncol=maxDataSets)
  sampCountyOverSampDat90 = matrix(NA, length(parNames), ncol=maxDataSets)
  rownames(sampCountyOverSampDatPar) = parNames
  rownames(sampCountyOverSampDatSD) = parNames
  rownames(sampCountyOverSampDat10) = parNames
  rownames(sampCountyOverSampDat50) = parNames
  rownames(sampCountyOverSampDat90) = parNames
  for(i in 1:maxDataSets) {
    # Extract data
    currData = overSampDat$clustDat[[i]]
    currData$admin1 = factor(currData$admin1)
    
    # INLA data
    dat = list(y = currData$died,
               Ntrials = currData$numChildren,
               urban = currData$urban,
               idx = as.numeric(currData$admin1),
               idx2 = as.numeric(currData$admin1),
               idxEps = 1:length(currData$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
    dat$urban = c(rep(c(1,0), each = 47), dat$urban)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idx2 = c(rep(1:47, 2), dat$idx2)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
    
    # Run model
    print(paste0("fitting urban oversampled model for dataset ", i, "/", maxDataSets))
    result = inla(formula = formula,
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat, quantiles = c(0.1, 0.5, 0.9), 
                  control.compute = list(config = TRUE))
    
    # fixed effects
    sampCountyOverSampDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountyOverSampDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountyOverSampDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountyOverSampDat50[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountyOverSampDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    ## include parameter estimates in the table
    # fixed effects
    sampCountyOverSampDatPar[1:(1 + includeUrbanRural), i] = result$summary.fixed[,1]
    sampCountyOverSampDatSD[1:(1 + includeUrbanRural), i] = result$summary.fixed[,2]
    sampCountyOverSampDat10[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountyOverSampDat50[1:(1 + includeUrbanRural), i] = result$summary.fixed[,3]
    sampCountyOverSampDat90[1:(1 + includeUrbanRural), i] = result$summary.fixed[,5]
    
    # BYM2 hyperparameter phi
    sampCountyOverSampDatPar[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,1]
    sampCountyOverSampDatSD[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,2]
    sampCountyOverSampDat10[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,3]
    sampCountyOverSampDat50[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,3]
    sampCountyOverSampDat90[(2 + includeCluster + includeUrbanRural), i] = result$summary.hyperpar[2,5]
    
    ## transformed hyperparameters
    # sample the hyperparameters, using the marginals to improve the sampling
    out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
    transformFunction = function(x) {c(1/x[3], x[2], 1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[3]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
    if(!includeCluster)
      transformFunction = function(x) {c(1/x[1], x[2], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
    transformedOut = apply(out, 1, transformFunction)
    
    # now calculate the summary statistics of the transformed BYM2 hyperparameters
    sampCountyOverSampDatPar[(2 + includeUrbanRural):length(parNames), i] = rowMeans(transformedOut)
    sampCountyOverSampDatSD[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, sd)
    sampCountyOverSampDat10[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, quantile, probs=.1)
    sampCountyOverSampDat50[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, quantile, probs=.5)
    sampCountyOverSampDat90[(2 + includeUrbanRural):length(parNames), i] = apply(transformedOut, 1, quantile, probs=.9)
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
        
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          # muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
          # muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
          # sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
          # sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
          sampRuralMod[, j] = sapply(1:length(clustersPerRural), function(i) {mean(expit(sampRural[i, j] + rnorm(clustersPerRural[i], sd=clusterSigma)))})
          sampRuralMod[clustersPerRural == 0,] = 0 # we just need to set these values to something other than NaN
          sampUrbanMod[, j] = sapply(1:length(clustersPerUrban), function(i) {mean(expit(sampUrban[i, j] + rnorm(clustersPerUrban[i], sd=clusterSigma)))})
          sampClusterOverSampDatRural[,j,i] = sampRural[, j] + rnorm(47, sd=clusterSigma)
          sampClusterOverSampDatUrban[,j,i] = sampUrban[, j] + rnorm(47, sd=clusterSigma)
        }
      }
      sampCountyOverSampDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
      sampCountyOverSampDatAdjusted[,,i] = logit(expit(sampUrban)*urbRatioAdjusted + expit(sampRural)*(1-urbRatioAdjusted))
      sampCountyOverSampDatMod[,,i] = logit(sampUrbanMod*urbRatioAdjusted + sampRuralMod*(1-urbRatioAdjusted))
      if(!includeCluster) {
        sampClusterOverSampDatRural[,,i] = sampRural
        sampClusterOverSampDatUrban[,,i] = sampUrban
      } else {
        sampPixelOverSampDatUrbanMod[,,i] = sampUrbanMod
        sampPixelOverSampDatRuralMod[,,i] = sampRuralMod
      }
      sampPixelOverSampDatUrban[,,i] = sampUrban
      sampPixelOverSampDatRural[,,i] = sampRural
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
      sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampCounty[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          # muSigmaMat = cbind(sampCounty[, j], clusterSigma)
          # sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat
          sampCountyMod[, j] = sapply(1:length(clustersPerCounty), function(i) {mean(expit(sampCounty[i, j] + rnorm(clustersPerCounty[i], sd=clusterSigma)))})
          sampClusterOverSampDatRural[,j,i] = sampCounty[, j] + rnorm(47, sd=clusterSigma)
          sampClusterOverSampDatUrban[,j,i] = sampClusterOverSampDatRural[,j,i]
        }
      }
      sampCountyOverSampDat[,,i] = sampCounty
      sampCountyOverSampDatAdjusted[,,i] = sampCounty
      sampCountyOverSampDatMod[,,i] = logit(sampCountyMod)
      sampPixelOverSampDatUrban[,,i] = sampCounty
      sampPixelOverSampDatRural[,,i] = sampCounty
      if(!includeCluster) {
        sampClusterOverSampDatRural[,,i] = sampCounty
        sampClusterOverSampDatUrban[,,i] = sampCounty
      } else {
        sampPixelOverSampDatUrbanMod[,,i] = sampCountyMod
        sampPixelOverSampDatRuralMod[,,i] = sampCountyMod
      }
    }
  }
  
  processSamples = function(samp){
    # 80% credible intervals
    CI = t(apply(X = samp, MARGIN = 1, FUN = quantile, probs = c(0.1, 0.5, 0.9)))
    mm = rowMeans(samp)
    ss = apply(X = samp, MARGIN = 1, FUN = sd)
    
    return(list(logit = list(CI = CI,
                             mean = mm,
                             stddev = ss),
                prob = list(CI = exp(CI))))
  }
  
  ## SRSdata
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  for(i in 1:dim(sampCountySRSDat)[3]){
    tmp = processSamples(sampCountySRSDat[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resSRSdat = list(Q10 = Q10,
                   Q50 = Q50,
                   Q90 = Q90,
                   mean = mm,
                   stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatAdjusted)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatAdjusted)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatAdjusted)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatAdjusted)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatAdjusted)[3])
  for(i in 1:dim(sampCountySRSDatAdjusted)[3]){
    tmp = processSamples(sampCountySRSDatAdjusted[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resSRSdatAdjusted = list(Q10 = Q10,
                           Q50 = Q50,
                           Q90 = Q90,
                           mean = mm,
                           stddev = ss)
  
  # calculate summary statistics for cluster and pixel predictions
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatRural)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatRural)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatRural)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatRural)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatRural)[3])
  for(i in 1:dim(sampClusterSRSDatRural)[3]) {
    tmp = processSamples(sampClusterSRSDatRural[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatSRSClusterRural = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatUrban)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatUrban)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatUrban)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatUrban)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampClusterSRSDatUrban)[3])
  for(i in 1:dim(sampClusterSRSDatUrban)[3]) {
    tmp = processSamples(sampClusterSRSDatUrban[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatSRSClusterUrban = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrban)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrban)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrban)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrban)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrban)[3])
  for(i in 1:dim(sampPixelSRSDatUrban)[3]) {
    tmp = processSamples(sampPixelSRSDatUrban[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatSRSPixelUrban = list(Q10 = Q10,
                             Q50 = Q50,
                             Q90 = Q90,
                             mean = mm,
                             stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRural)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRural)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRural)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRural)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRural)[3])
  for(i in 1:dim(sampPixelSRSDatRural)[3]) {
    tmp = processSamples(sampPixelSRSDatRural[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatSRSPixelRural = list(Q10 = Q10,
                             Q50 = Q50,
                             Q90 = Q90,
                             mean = mm,
                             stddev = ss)
  
  if(includeCluster) {
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    for(i in 1:dim(sampCountySRSDatMod)[3]){
      tmp = processSamples(sampCountySRSDatMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resSRSdatMod = list(Q10 = Q10,
                        Q50 = Q50,
                        Q90 = Q90,
                        mean = mm,
                        stddev = ss)
    
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrbanMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrbanMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrbanMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrbanMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatUrbanMod)[3])
    for(i in 1:dim(sampPixelSRSDatUrbanMod)[3]) {
      tmp = processSamples(sampPixelSRSDatUrbanMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resDatSRSPixelUrbanMod = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
    
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRuralMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRuralMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRuralMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRuralMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampPixelSRSDatRuralMod)[3])
    for(i in 1:dim(sampPixelSRSDatRuralMod)[3]) {
      tmp = processSamples(sampPixelSRSDatRuralMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resDatSRSPixelRuralMod = list(Q10 = Q10,
                                  Q50 = Q50,
                                  Q90 = Q90,
                                  mean = mm,
                                  stddev = ss)
  } else {
    resSRSdatMod = NULL
    resDatSRSPixelUrbanMod = NULL
    resDatSRSPixelRuralMod = NULL
  }
  
  ## now collect the parameters
  # make rural parameter the urban parameter
  # if(includeUrbanRural) {
  #   sampCountySRSDatPar[2,] = -sampCountySRSDatPar[2,]
  #   sampCountySRSDat10[2,] = -sampCountySRSDat10[2,]
  #   sampCountySRSDat50[2,] = -sampCountySRSDat50[2,]
  #   sampCountySRSDat90[2,] = -sampCountySRSDat90[2,]
  # }
  mm = rowMeans(sampCountySRSDatPar)
  ss = rowMeans(sampCountySRSDatSD)
  Q10 = rowMeans(sampCountySRSDat10)
  Q50 = rowMeans(sampCountySRSDat50)
  Q90 = rowMeans(sampCountySRSDat90)
  resSRSdatPar = data.frame(list(Est = mm,
                                 SD = ss, 
                                 Q10 = Q10,
                                 Q50 = Q50,
                                 Q90 = Q90))
  
  ## overSampDat
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  for(i in 1:dim(sampCountyOverSampDat)[3]){
    tmp = processSamples(sampCountyOverSampDat[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatOverSamp = list(Q10 = Q10,
                        Q50 = Q50,
                        Q90 = Q90,
                        mean = mm,
                        stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatAdjusted)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatAdjusted)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatAdjusted)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatAdjusted)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatAdjusted)[3])
  for(i in 1:dim(sampCountyOverSampDatAdjusted)[3]){
    tmp = processSamples(sampCountyOverSampDatAdjusted[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatOverSampAdjusted = list(Q10 = Q10,
                        Q50 = Q50,
                        Q90 = Q90,
                        mean = mm,
                        stddev = ss)
  
  # calculate summary statistics for cluster and pixel predictions
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatRural)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatRural)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatRural)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatRural)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatRural)[3])
  for(i in 1:dim(sampClusterOverSampDatRural)[3]) {
    tmp = processSamples(sampClusterOverSampDatRural[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatOverSampClusterRural = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatUrban)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatUrban)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatUrban)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatUrban)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampClusterOverSampDatUrban)[3])
  for(i in 1:dim(sampClusterOverSampDatUrban)[3]) {
    tmp = processSamples(sampClusterOverSampDatUrban[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatOverSampClusterUrban = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrban)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrban)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrban)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrban)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrban)[3])
  for(i in 1:dim(sampPixelOverSampDatUrban)[3]) {
    tmp = processSamples(sampPixelOverSampDatUrban[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatOverSampPixelUrban = list(Q10 = Q10,
                             Q50 = Q50,
                             Q90 = Q90,
                             mean = mm,
                             stddev = ss)
  
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRural)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRural)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRural)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRural)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRural)[3])
  for(i in 1:dim(sampPixelOverSampDatRural)[3]) {
    tmp = processSamples(sampPixelOverSampDatRural[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatOverSampPixelRural = list(Q10 = Q10,
                             Q50 = Q50,
                             Q90 = Q90,
                             mean = mm,
                             stddev = ss)
  
  if(includeCluster) {
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    for(i in 1:dim(sampCountyOverSampDatMod)[3]){
      tmp = processSamples(sampCountyOverSampDatMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resDatOverSampMod = list(Q10 = Q10,
                             Q50 = Q50,
                             Q90 = Q90,
                             mean = mm,
                             stddev = ss)
    
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrbanMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrbanMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrbanMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrbanMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatUrbanMod)[3])
    for(i in 1:dim(sampPixelOverSampDatUrbanMod)[3]) {
      tmp = processSamples(sampPixelOverSampDatUrbanMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resDatOverSampPixelUrbanMod = list(Q10 = Q10,
                                  Q50 = Q50,
                                  Q90 = Q90,
                                  mean = mm,
                                  stddev = ss)
    
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRuralMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRuralMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRuralMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRuralMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampPixelOverSampDatRuralMod)[3])
    for(i in 1:dim(sampPixelOverSampDatRuralMod)[3]) {
      tmp = processSamples(sampPixelOverSampDatRuralMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resDatOverSampPixelRuralMod = list(Q10 = Q10,
                                  Q50 = Q50,
                                  Q90 = Q90,
                                  mean = mm,
                                  stddev = ss)
  } else {
    resDatOverSampMod = NULL
    resDatOverSampPixelUrbanMod = NULL
    resDatOverSampPixelRuralMod = NULL
  }
  
  ## now collect the parameters
  # make rural parameter the urban parameter
  # if(includeUrbanRural) {
  #   sampCountyOverSampDatPar[2,] = -sampCountyOverSampDatPar[2,]
  #   sampCountyOverSampDat10[2,] = -sampCountyOverSampDat10[2,]
  #   sampCountyOverSampDat50[2,] = -sampCountyOverSampDat50[2,]
  #   sampCountyOverSampDat90[2,] = -sampCountyOverSampDat90[2,]
  # }
  mm = rowMeans(sampCountyOverSampDatPar)
  ss = rowMeans(sampCountyOverSampDatSD)
  Q10 = rowMeans(sampCountyOverSampDat10)
  Q50 = rowMeans(sampCountyOverSampDat50)
  Q90 = rowMeans(sampCountyOverSampDat90)
  resDatOverSampPar = data.frame(list(Est = mm,
                                      SD = ss, 
                                      Q10 = Q10,
                                      Q50 = Q50,
                                      Q90 = Q90))
  
  # Full result
  designRes = list(SRSdat = resSRSdatAdjusted,
                   overSampDat = resDatOverSampAdjusted, 
                   SRSdatUnadjusted = resSRSdat,
                   overSampDatUnadjusted = resDatOverSamp, 
                   SRSdatPixelRural = resDatSRSPixelRural,
                   overSampDatPixelRural = resDatOverSampPixelRural, 
                   SRSdatPixelUrban = resDatSRSPixelUrban,
                   overSampDatPixelUrban = resDatOverSampPixelUrban, 
                   SRSdatClusterRural = resDatSRSClusterRural,
                   overSampDatClusterRural = resDatOverSampClusterRural, 
                   SRSdatClusterUrban = resDatSRSClusterUrban,
                   overSampDatClusterUrban = resDatOverSampClusterUrban, 
                   overSampDatPar = resDatOverSampPar, 
                   SRSdatPar = resSRSdatPar)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      includeUrbanRural, '.RData'), designRes = designRes)
  
  testText = ifelse(test, "Test", "")
  strictPriorText = ifelse(strictPrior, "strictPrior", "")
  save(file = paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                     includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 
                     maxDataSets, strictPriorText, testText, rangeText, '.RData'), 
       designRes = designRes)
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    designRes = list(SRSdat = resSRSdatMod,
                     overSampDat = resDatOverSampMod, 
                     SRSdatAdjusted = resSRSdatAdjusted,
                     overSampDatAdjusted = resDatOverSampAdjusted, 
                     SRSdatUnadjusted = resSRSdat,
                     overSampDatUnadjusted = resDatOverSamp, 
                     SRSdatPixelRural = resDatSRSPixelRuralMod,
                     overSampDatPixelRural = resDatOverSampPixelRuralMod, 
                     SRSdatPixelUrban = resDatSRSPixelRuralMod,
                     overSampDatPixelUrban = resDatOverSampPixelRuralMod, 
                     SRSdatClusterRural = resDatSRSClusterRural,
                     overSampDatClusterRural = resDatOverSampClusterRural, 
                     SRSdatClusterUrban = resDatSRSClusterUrban,
                     overSampDatClusterUrban = resDatOverSampClusterUrban, 
                     overSampDatPar = resDatOverSampPar, 
                     SRSdatPar = resSRSdatPar)
    
    save(file = paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                       includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, 'debiasedMaxDataSets', 
                       maxDataSets, strictPriorText, testText, rangeText, '.RData'), 
         designRes = designRes)
  }
  
  invisible(designRes)
}

# same as runBYM2, except fits a single data set (the ed global data frame)
# doPredsAtPostMean: if TRUE, fix all model hyperparameters at the posterior mean 
# getPosteriorDensity: EXPERIMENTAL: evaluate the posterior density of direct estimates using multivariate normal approximation
runBYM2Dat = function(dat=ed, includeUrbanRural=TRUE, includeCluster=TRUE, saveResults=TRUE, fileNameRoot="Ed", 
                      previousResult=NULL, doValidation=FALSE, predCountyI=NULL, counties=sort(unique(poppc$County)), 
                      doPredsAtPostMean=FALSE, directLogitEsts=NULL, fixedParameters=FALSE, getPosteriorDensity=FALSE) {
  
  # remove data from the given county for validation if necessary
  if(!is.null(predCountyI))
    dat$y[as.character(dat$admin1)==counties[predCountyI]] = NA
  
  includeUrban = includeUrbanRural
  
  # Get true ratios of urban/rural
  if(fileNameRoot == "Ed")
    thisPopTable = adjustPopulationPerCountyTable("women")
  else
    thisPopTable = adjustPopulationPerCountyTable("children")
  urbRatio = vector('numeric', length = 47)
  counties = sort(unique(as.character(dat$admin1)))
  urbRatio = thisPopTable$popUrb / thisPopTable$popTotal
  sortI = matchMultiple(counties, thisPopTable$County)
  urbRatio = urbRatio[sortI]
  
  # Define formula
  if(includeUrbanRural) {
    if(includeCluster) {
      formula = y ~ urban +
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 2/3), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(1,0.01))))
    } else {
      formula = y ~ urban + 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 2/3), prior="pc")))
    }
  } else {
    if(includeCluster) {
      formula = y ~ f(idx, model="bym2",
                      graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
                      hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), 
                                 phi=list(param=c(0.5, 2/3), prior="pc"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = "pc.prec", param = c(1,0.01))))
    } else {
      formula = y ~ 
        f(idx, model="bym2",
          graph="Kenyaadm1.graph", scale.model=TRUE, constr=TRUE, 
          hyper=list(prec=list(param=c(1, 0.01), prior="pc.prec"), phi=list(param=c(0.5, 2/3), prior="pc")))
    }
  }
  
  # Number of simulations for producing results
  Nsim = 1000
  
  # Help functions
  logit = function(x){
    return(log(x/(1-x)))
  }
  expit = function(x){
    return(1/(1+exp(-x)))
  }
  
  if(includeCluster)
    parNames = c("Intercept", "Urban", "Cluster Var", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var", "Cluster SD", "BYM2 Tot. SD", "BYM2 Spatial SD", "BYM2 iid SD")
  else
    parNames = c("Intercept", "Urban", "BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var", "BYM2 Tot. SD", "BYM2 Spatial SD", "BYM2 iid SD")
  includeI = c(1, rep(2, includeUrban), 3:length(parNames))
  parNames = parNames[includeI]
  
  # Go through education data set
  sampCountyDat = matrix(NA, nrow=47, ncol=Nsim)
  sampClusterDatUrban = matrix(NA, nrow=47, ncol=Nsim)
  sampClusterDatRural = matrix(NA, nrow=47, ncol=Nsim)
  sampPixelDatUrban = matrix(NA, nrow=47, ncol=Nsim)
  sampPixelDatRural = matrix(NA, nrow=47, ncol=Nsim)
  
  # for the final parameters to store, 2 fixed effects, 2 + includeCluster estimated 
  # hyperparameters, and 2 hyperparameters we will get via transformation
  sampCountyDatPar = numeric(length(parNames))
  sampCountyDatSD = numeric(length(parNames))
  sampCountyDat10 = numeric(length(parNames))
  sampCountyDat50 = numeric(length(parNames))
  sampCountyDat90 = numeric(length(parNames))
  
  names(sampCountyDatPar) = parNames
  names(sampCountyDatSD) = parNames
  names(sampCountyDat10) = parNames
  names(sampCountyDat50) = parNames
  names(sampCountyDat90) = parNames
  
  # Extract data
  dat$admin1 = factor(dat$admin1)
  
  # INLA data
  dat = list(y = dat$y,
             Ntrials = dat$n,
             urban = dat$urban,
             idx = as.numeric(dat$admin1),
             idxEps = 1:length(dat$y))
  
  # Add unobserved data to make sampling easier
  dat$y = c(rep(NA, 47*2), dat$y)
  dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
  dat$urban = c(rep(c(1,0), each = 47), dat$urban)
  dat$idx = c(rep(1:47, 2), dat$idx)
  dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
  
  # set posterior approximation for calculating validation quantities
  
  if(!doPredsAtPostMean) {
    # based on http://www.r-inla.org/faq#TOC-How-can-I-compute-cross-validation-or-predictive-measures-of-fit-
    # when calculating log posterior density, we use fewer integration points to get better estimates of the 
    # covariance matrices at each integration point
    if(!getPosteriorDensity)
      control.inla = list(strategy="laplace", int.strategy="grid", diff.logdens=4, npoints=21)
    else
      control.inla = list(strategy="laplace", int.strategy="grid", diff.logdens=4, npoints=9)
  }
  else {
    # we need to generate predictions in this case based on a fixed set of hyperparameters
    control.inla = list(strategy="laplace", int.strategy="eb") 
  }
  
  # initialize the fitting process based on a previous optimum if necessary
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousResult)) {
    # initialize the fitting process based on a previous optimum
    # modeControl$result = previousResult
    modeControl$theta = previousResult$mode$theta
    modeControl$x = previousResult$mode$x
    modeControl$restart = !fixedParameters
    modeControl$fixed = fixedParameters
  }
  
  # Run model
  print("fitting BYM model...")
  result = inla(formula = formula, 
                family="binomial",
                Ntrials = Ntrials,
                data=dat, 
                control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
                quantiles=c(0.1, 0.5, 0.9), 
                control.mode=modeControl)
  
  ## include parameter estimates in the table
  # fixed effects
  sampCountyDatPar[1:(1 + includeUrban)] = result$summary.fixed[,1]
  sampCountyDatSD[1:(1 + includeUrban)] = result$summary.fixed[,2]
  sampCountyDat10[1:(1 + includeUrban)] = result$summary.fixed[,3]
  sampCountyDat50[1:(1 + includeUrban)] = result$summary.fixed[,4]
  sampCountyDat90[1:(1 + includeUrban)] = result$summary.fixed[,5]
  
  # BYM2 hyperparameter phi
  sampCountyDatPar[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,1]
  sampCountyDatSD[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,2]
  sampCountyDat10[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,3]
  sampCountyDat50[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,4]
  sampCountyDat90[(2 + includeCluster + includeUrban)] = result$summary.hyperpar[2,5]
  
  ## transformed hyperparameters
  # sample the hyperparameters, using the marginals to improve the sampling
  out = inla.hyperpar.sample(1000, result, improve.marginals=TRUE)
  transformFunction = function(x) {c(1/x[3], x[2], 1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[3]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
  if(!includeCluster)
    transformFunction = function(x) {c(1/x[1], x[2], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
  transformedOut = apply(out, 1, transformFunction)
  
  # now calculate the summary statistics of the transformed BYM2 hyperparameters
  sampCountyDatPar[(2 + includeUrban):length(parNames)] = rowMeans(transformedOut)
  sampCountyDatSD[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, sd)
  sampCountyDat10[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, quantile, probs=.1)
  sampCountyDat50[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, quantile, probs=.5)
  sampCountyDat90[(2 + includeUrban):length(parNames)] = apply(transformedOut, 1, quantile, probs=.9)
  
  if(includeCluster)
    sampClusterSigmaSRS = 1:Nsim
  
  sortI = matchMultiple(counties, easpc$County)
  clustersPerUrban = easpc$EAUrb[sortI]
  clustersPerRural = easpc$EARur[sortI]
  
  # Simulate from posterior
  if(includeUrbanRural) {
    samp = inla.posterior.sample(n = Nsim, result = result)
    sampRural = matrix(NA, nrow = 47, ncol = Nsim)
    sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
    sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatUrban = matrix(NA, nrow = 47, ncol = Nsim)
    sampPixelDatRural = matrix(NA, nrow = 47, ncol = Nsim)
    
    for(j in 1:Nsim){
      sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
      sampUrban[, j] = samp[[j]]$latent[1:47]
      if(includeCluster) {
        # if cluster effect is included, must debias predictions in each modeled strata
        clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
        # muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
        # muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
        # sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
        # sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
        sampRuralMod[, j] = sapply(1:length(clustersPerRural), function(i) {mean(expit(sampRural[i, j] + rnorm(clustersPerRural[i], sd=clusterSigma)))})
        sampRuralMod[clustersPerRural == 0,] = 0 # we just need to set these values to something other than NaN
        sampUrbanMod[, j] = sapply(1:length(clustersPerUrban), function(i) {mean(expit(sampUrban[i, j] + rnorm(clustersPerUrban[i], sd=clusterSigma)))})
        sampClusterDatRural[, j] = sampRural[, j] + rnorm(47, sd=clusterSigma)
        sampClusterDatUrban[, j] = sampUrban[, j] + rnorm(47, sd=clusterSigma)
      }
    }
    sampCountyDat = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
    sampCountyDatMod = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
    if(!includeCluster) {
      sampClusterDatRural = sampRural
      sampClusterDatUrban = sampUrban
    } else {
      sampPixelDatUrbanMod = sampUrbanMod
      sampPixelDatRuralMod = sampRuralMod
    }
    sampPixelDatUrban = sampUrban
    sampPixelDatRural = sampRural
  } else {
    samp = inla.posterior.sample(n = Nsim, result = result)
    sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
    sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
    clustersPerCounty = rowSums(cbind(clustersPerUrban, clustersPerRural))
    for(j in 1:Nsim){
      sampCounty[, j] = samp[[j]]$latent[1:47]
      if(includeCluster) {
        # if cluster effect is included, must debias predictions in each modeled strata
        clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
        # muSigmaMat = cbind(sampCounty[, j], clusterSigma)
        # sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat
        sampCountyMod[, j] = sapply(1:length(clustersPerCounty), function(i) {mean(expit(sampCounty[i, j] + rnorm(clustersPerCounty[i], sd=clusterSigma)))})
        sampClusterDatRural[, j] = sampCounty[, j] + rnorm(47, sd=clusterSigma)
        sampClusterDatUrban[, j] = sampClusterDatRural[, j]
      }
    }
    sampCountyDat = sampCounty
    sampCountyDatMod = logit(sampCountyMod)
    sampPixelDatUrban = sampCounty
    sampPixelDatRural = sampCounty
    if(!includeCluster) {
      sampClusterDatRural = sampCounty
      sampClusterDatUrban = sampCounty
    } else {
      sampPixelDatUrbanMod = sampCountyMod
      sampPixelDatRuralMod = sampCountyMod
    }
  }
  
  processSamples = function(samp){
    # 80% credible intervals
    CI = t(apply(X = samp, MARGIN = 1, FUN = quantile, probs = c(0.1, 0.5, 0.9)))
    mm = rowMeans(samp)
    ss = apply(X = samp, MARGIN = 1, FUN = sd)
    
    return(list(logit = list(CI = CI,
                             mean = mm,
                             stddev = ss),
                prob = list(CI = exp(CI))))
  }
  
  ## generate predictions
  Q10 = numeric(47)
  Q50 = numeric(47)
  Q90 = numeric(47)
  mm = numeric(47)
  ss = numeric(47)
  
  tmp = processSamples(sampCountyDat)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDat = list(Q10 = Q10,
                 Q50 = Q50,
                 Q90 = Q90,
                 mean = mm,
                 stddev = ss)
  
  if(!is.null(predCountyI))
    resDat = data.frame(resDat)[predCountyI,]
  
  # calculate summary statistics for cluster and pixel predictions
  tmp = processSamples(sampClusterDatRural)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatClusterRural = list(Q10 = Q10,
                            Q50 = Q50,
                            Q90 = Q90,
                            mean = mm,
                            stddev = ss)
  
  if(!is.null(predCountyI))
    resDatClusterRural = data.frame(resDatClusterRural)[predCountyI,]
  
  tmp = processSamples(sampClusterDatUrban)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatClusterUrban = list(Q10 = Q10,
                            Q50 = Q50,
                            Q90 = Q90,
                            mean = mm,
                            stddev = ss)
  
  if(!is.null(predCountyI))
    resDatClusterUrban = data.frame(resDatClusterUrban)[predCountyI,]
  
  tmp = processSamples(sampPixelDatUrban)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatPixelUrban = list(Q10 = Q10,
                          Q50 = Q50,
                          Q90 = Q90,
                          mean = mm,
                          stddev = ss)
  
  if(!is.null(predCountyI))
    resDatPixelUrban = data.frame(resDatPixelUrban)[predCountyI,]
  
  tmp = processSamples(sampPixelDatRural)
  Q10 = tmp$logit$CI[,1]
  Q50 = tmp$logit$CI[,2]
  Q90 = tmp$logit$CI[,3]
  mm = tmp$logit$mean
  ss = tmp$logit$stddev
  
  resDatPixelRural = list(Q10 = Q10,
                          Q50 = Q50,
                          Q90 = Q90,
                          mean = mm,
                          stddev = ss)
  
  if(!is.null(predCountyI))
    resDatPixelRural = data.frame(resDatPixelRural)[predCountyI,]
  
  if(includeCluster) {
    tmp = processSamples(sampCountyDatMod)
    Q10 = tmp$logit$CI[,1]
    Q50 = tmp$logit$CI[,2]
    Q90 = tmp$logit$CI[,3]
    mm = tmp$logit$mean
    ss = tmp$logit$stddev
    resDatMod = list(Q10 = Q10,
                      Q50 = Q50,
                      Q90 = Q90,
                      mean = mm,
                      stddev = ss)
    
    if(!is.null(predCountyI))
      resDatMod = data.frame(resDatMod)[predCountyI,]
    
    tmp = processSamples(sampPixelDatUrbanMod)
    Q10 = tmp$logit$CI[,1]
    Q50 = tmp$logit$CI[,2]
    Q90 = tmp$logit$CI[,3]
    mm = tmp$logit$mean
    ss = tmp$logit$stddev
    resDatPixelUrbanMod = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
    
    if(!is.null(predCountyI))
      resDatPixelUrbanMod = data.frame(resDatPixelUrbanMod)[predCountyI,]
    
    tmp = processSamples(sampPixelDatRuralMod)
    Q10 = tmp$logit$CI[,1]
    Q50 = tmp$logit$CI[,2]
    Q90 = tmp$logit$CI[,3]
    mm = tmp$logit$mean
    ss = tmp$logit$stddev
    resDatPixelRuralMod = list(Q10 = Q10,
                               Q50 = Q50,
                               Q90 = Q90,
                               mean = mm,
                               stddev = ss)
    
    if(!is.null(predCountyI))
      resDatPixelRuralMod = data.frame(resDatPixelRuralMod)[predCountyI,]
    
  } else {
    resDatMod = NULL
    resDatPixelUrbanMod = NULL
    resDatPixelRuralMod = NULL
  }
  
  
  
  ## now collect the parameters
  # make rural parameter the urban parameter
  # if(includeUrban) {
  #   sampCountyDatPar[2] = -sampCountyDatPar[2]
  #   sampCountyDat10[2] = -sampCountyDat10[2]
  #   sampCountyDat50[2] = -sampCountyDat50[2]
  #   sampCountyDat90[2] = -sampCountyDat90[2]
  # }
  mm = sampCountyDatPar
  ss = sampCountyDatSD
  Q10 = sampCountyDat10
  Q50 = sampCountyDat50
  Q90 = sampCountyDat90
  resDatPar = data.frame(list(mean = mm,
                              stddev = ss, 
                              Q10 = Q10,
                              Q50 = Q50,
                              Q90 = Q90))
  
  # Full result
  designRes = list(predictions = resDat,
                   parameters = resDatPar, 
                   predictionsPixelUrban = resDatPixelUrban, 
                   predictionsPixelRural = resDatPixelRural, 
                   predictionsClusterUrban = resDatClusterUrban, 
                   predictionsClusterRural = resDatClusterRural)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      includeUrbanRural, '.RData'), designRes = designRes)
  validationText = ""
  if(doValidation)
    validationText = "ValidationFull"
  else if(!is.null(predCountyI))
    validationText = paste0("Validation", predCountyI) # really we are setting saveResults to FALSE in this case
  if(saveResults) {
    save(file = paste0('bym2', fileNameRoot, validationText, 'UrbRur',includeUrbanRural, 'Cluster', includeCluster, '.RData'), 
         designRes = designRes)
  }
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    temp = designRes
    designRes = list(predictions = resDatMod,
                     parameters = resDatPar, 
                     predictionsPixelUrban = resDatPixelUrbanMod, 
                     predictionsPixelRural = resDatPixelRuralMod, 
                     predictionsClusterUrban = resDatClusterUrban, 
                     predictionsClusterRural = resDatClusterRural)
    
    if(saveResults) {
      save(file = paste0('bym2', fileNameRoot, validationText, 'UrbRur',includeUrbanRural, 'Cluster', includeCluster, 'debiased.RData'), 
           designRes = designRes)
    }
    
    designResMod = designRes
    designRes = temp
  }
  
  # in order to compute DIC or WAIC it is necessary to generate predictive distribution at the 
  # posterior mean for each county. We do it here on a logit scale:
  if(getPosteriorDensity) {
    if(is.null(directLogitEsts))
        stop("Must supply directLogitEsts if getPosteriorDensity == TRUE")
    
    # approximate the posterior county samples using a multivariate gaussian on the logit scale
    theseCountySamples = sampCountyDat
    mu = resDat$mean
    theseCountyResiduals = sweep(theseCountySamples, 1, mu, "-")
    Sigma = (1 / (Nsim - 1)) * theseCountyResiduals %*% t(theseCountyResiduals)
    if(includeCluster) {
      theseCountySamplesMod = sampCountyDatMod
      muMod = resDatMod$mean
      theseCountyResidualsMod = sweep(theseCountySamplesMod, 1, muMod, "-")
      SigmaMod = (1 / (Nsim - 1)) * theseCountyResidualsMod %*% t(theseCountyResidualsMod)
    }
    
    # calculate posterior density of the direct estimates
    logLik = logLikGP(directLogitEsts - mu, chol(Sigma))
    if(includeCluster)
      logLikMod = logLikGP(directLogitEsts - muMod, chol(SigmaMod))
    else
      logLikMod = NULL
  }
  
  # compute cluster predictive standard deviation for each strata (square root of the the strata variance plus the cluster variance)
  if(includeCluster) {
    clustPredSD = sqrt(designRes$predictions$stddev^2 + designRes$parameters$stddev[2 + includeUrban]^2)
    designRes$predictions$clustPredSD = clustPredSD
    designResMod$predictions$clustPredSD = clustPredSD
  }
  
  # compile results
  if(!doValidation)
    out = designRes
  else if(!getPosteriorDensity)
    out = list(designRes=designRes, mod=result)
  else
    out = list(designRes=designRes, mod=result, directLogLik=logLik, directLogLikMod=logLikMod)
  if(includeCluster)
    out$designResMod = designResMod
  
  out
}

runBYM = function(tausq=0.1^2, test=FALSE, includeUrbanRural=TRUE, includeCluster=TRUE) {
  
  # load and relevant data
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar0.0225tausq", round(tausq, 4), "gamma-1HHoldVar0urbanOverSamplefrac0.RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar0.0225tausq", round(tausq, 4), "gamma-1HHoldVar0urbanOverSamplefrac0Test.RData"))
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
  
  # Get true ratios of urban/rural
  urbRatio = vector('numeric', length = 47)
  for(i in 1:47){
    idx = which(as.numeric(factor(SRSDat$eaDat$admin1)) == i)
    urbanI = idx[SRSDat$eaDat$urban[idx]]
    urbRatio[i] = sum(SRSDat$eaDat$numChildren[urbanI])/sum(SRSDat$eaDat$numChildren[idx])
  }
  
  # Define formula
  if(includeUrbanRural) {
    if(includeCluster) {
      formula = y ~ rural +
        f(idx, model="iid", 
          hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag",
          graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = 'loggamma', param = c(1,0.01))))
    } else {
      formula = y ~ rural + 
        f(idx, model="iid", 
          hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag",
          graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma")))
    }
  } else {
    if(includeCluster) {
      formula = y ~ f(idx, model="iid", 
                      hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag", graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma"))) +
        f(idxEps, model = "iid",
          hyper = list(prec = list(prior = 'loggamma', param = c(1,0.01))))
    } else {
      formula = y ~ 
        f(idx, model="iid", 
          hyper=list(prec=list(param=c(0.5, 0.001488), prior="loggamma"))) + 
        f(idx2, model="besag", graph="Kenyaadm1.graph", 
          hyper=list(prec=list(param=c(0.5, 0.00360), prior="loggamma")))
    }
  }
  
  
  # Number of simulations for producing results
  Nsim = 1000
  
  # Help functions
  logit = function(x){
    return(log(x/(1-x)))
  }
  expit = function(x){
    return(1/(1+exp(-x)))
  }
  
  # Go through datasets for SRSDat
  sampCountySRSDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountySRSDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  for(i in 1:maxDataSets){
    # Extract data
    currData = SRSDat$clustDat[[i]]
    currData$admin1 = factor(currData$admin1)
    
    # INLA data
    dat = list(y = currData$died,
               Ntrials = currData$numChildren,
               rural = 1-currData$urban,
               idx = as.numeric(currData$admin1),
               idx2 = as.numeric(currData$admin1),
               idxEps = 1:length(currData$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
    dat$rural = c(rep(c(0,1), each = 47), dat$rural)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idx2 = c(rep(1:47, 2), dat$idx2)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
    
    # Run model
    result = inla(formula = formula,
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat,
                  control.compute = list(config = TRUE))
    
    if(includeCluster)
      sampClusterSigmaSRS = 1:Nsim
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
      
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
          muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
          sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
          sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
        }
      }
      sampCountySRSDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
      sampCountySRSDatMod[,,i] = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
      sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampCounty[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMat = cbind(sampCounty[, j], clusterSigma)
          sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat)
        }
      }
      sampCountySRSDat[,,i] = sampCounty
      sampCountySRSDatMod[,,i] = logit(sampCountyMod)
    }
  }
  
  # Go through datasets for overSampDat
  sampCountyOverSampDat = array(NA, dim = c(47, Nsim, maxDataSets))
  sampCountyOverSampDatMod = array(NA, dim = c(47, Nsim, maxDataSets))
  for(i in 1:maxDataSets) {
    # Extract data
    currData = overSampDat$clustDat[[i]]
    currData$admin1 = factor(currData$admin1)
    
    # INLA data
    dat = list(y = currData$died,
               Ntrials = currData$numChildren,
               rural = 1-currData$urban,
               idx = as.numeric(currData$admin1),
               idx2 = as.numeric(currData$admin1),
               idxEps = 1:length(currData$died))
    
    # Add unobserved data to make sampling easier
    dat$y = c(rep(NA, 47*2), dat$y)
    dat$Ntrials = c(rep(1, 47*2), dat$Ntrials)
    dat$rural = c(rep(c(0,1), each = 47), dat$rural)
    dat$idx = c(rep(1:47, 2), dat$idx)
    dat$idx2 = c(rep(1:47, 2), dat$idx2)
    dat$idxEps = c(rep(NA, 47*2), dat$idxEps)
    
    # Run model
    result = inla(formula = formula,
                  family="binomial",
                  Ntrials = Ntrials,
                  data=dat,
                  control.compute = list(config = TRUE))
    
    if(includeCluster) {
      sampClusterSigmaOverSampDat = 1:Nsim
    }
    
    # Simulate from posterior
    if(includeUrbanRural) {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampRural = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrban = matrix(NA, nrow = 47, ncol = Nsim)
      sampRuralMod = matrix(NA, nrow = 47, ncol = Nsim)
      sampUrbanMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampRural[, j] = samp[[j]]$latent[47 + (1:47)]
        sampUrban[, j] = samp[[j]]$latent[1:47]
        
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMatRural = cbind(sampRural[, j], clusterSigma)
          muSigmaMatUrban = cbind(sampUrban[, j], clusterSigma)
          sampRuralMod[, j] = logitNormMean(muSigmaMat = muSigmaMatRural)
          sampUrbanMod[, j] = logitNormMean(muSigmaMat = muSigmaMatUrban)
        }
      }
      sampCountyOverSampDat[,,i] = logit(expit(sampUrban)*urbRatio + expit(sampRural)*(1-urbRatio))
      sampCountyOverSampDatMod[,,i] = logit(sampUrbanMod*urbRatio + sampRuralMod*(1-urbRatio))
    } else {
      samp = inla.posterior.sample(n = Nsim, result = result)
      sampCounty = matrix(NA, nrow = 47, ncol = Nsim)
      sampCountyMod = matrix(NA, nrow = 47, ncol = Nsim)
      for(j in 1:Nsim){
        sampCounty[, j] = samp[[j]]$latent[1:47]
        if(includeCluster) {
          # if cluster effect is included, must debias predictions in each modeled strata
          clusterSigma = sqrt(1/samp[[j]]$hyperpar[3])
          muSigmaMat = cbind(sampCounty[, j], clusterSigma)
          sampCountyMod[, j] = logitNormMean(muSigmaMat = muSigmaMat)
        }
      }
      sampCountyOverSampDat[,,i] = sampCounty
      sampCountyOverSampDatMod[,,i] = logit(sampCountyMod)
    }
  }
  
  processSamples = function(samp){
    # 80% credible intervals
    CI = t(apply(X = samp, MARGIN = 1, FUN = quantile, probs = c(0.1, 0.5, 0.9)))
    mm = rowMeans(samp)
    ss = apply(X = samp, MARGIN = 1, FUN = sd)
    
    return(list(logit = list(CI = CI,
                             mean = mm,
                             stddev = ss),
                prob = list(CI = exp(CI))))
  }
  
  ## SRSdata
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDat)[3])
  for(i in 1:dim(sampCountySRSDat)[3]){
    tmp = processSamples(sampCountySRSDat[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resSRSdat = list(Q10 = Q10,
                   Q50 = Q50,
                   Q90 = Q90,
                   mean = mm,
                   stddev = ss)
  
  if(includeCluster) {
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampCountySRSDatMod)[3])
    for(i in 1:dim(sampCountySRSDatMod)[3]){
      tmp = processSamples(sampCountySRSDatMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resSRSdatMod = list(Q10 = Q10,
                        Q50 = Q50,
                        Q90 = Q90,
                        mean = mm,
                        stddev = ss)
  } else {
    resSRSdatMod = NULL
  }
  
  ## overSampDat
  Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  mm = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  ss = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDat)[3])
  for(i in 1:dim(sampCountyOverSampDat)[3]){
    tmp = processSamples(sampCountyOverSampDat[,,i])
    Q10[,i] = tmp$logit$CI[,1]
    Q50[,i] = tmp$logit$CI[,2]
    Q90[,i] = tmp$logit$CI[,3]
    mm[,i] = tmp$logit$mean
    ss[,i] = tmp$logit$stddev
  }
  resDatOverSamp = list(Q10 = Q10,
                        Q50 = Q50,
                        Q90 = Q90,
                        mean = mm,
                        stddev = ss)
  
  if(includeCluster) {
    Q10 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    Q50 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    Q90 = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    mm = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    ss = matrix(NA, nrow = 47, ncol = dim(sampCountyOverSampDatMod)[3])
    for(i in 1:dim(sampCountyOverSampDatMod)[3]){
      tmp = processSamples(sampCountyOverSampDatMod[,,i])
      Q10[,i] = tmp$logit$CI[,1]
      Q50[,i] = tmp$logit$CI[,2]
      Q90[,i] = tmp$logit$CI[,3]
      mm[,i] = tmp$logit$mean
      ss[,i] = tmp$logit$stddev
    }
    resDatOverSampMod = list(Q10 = Q10,
                             Q50 = Q50,
                             Q90 = Q90,
                             mean = mm,
                             stddev = ss)
  } else {
    resDatOverSampMod = NULL
  }
  
  # Full result
  designRes = list(SRSdat = resSRSdat,
                   overSampDat = resDatOverSamp)
  # save(file = 'kenyaSpatialDesignResultNew.RData', designRes = designRes)
  # save(file = paste0('kenyaSpatialDesignResultNewTausq0UrbRur', 
  #                      includeUrbanRural, '.RData'), designRes = designRes)
  testText = ifelse(test, "Test", "")
  save(file = paste0('kenyaSpatialDesignResultNewTausq', round(tausq, 4), 'UrbRur',
                     includeUrbanRural, 'Cluster', includeCluster, testText, '.RData'), 
       designRes = designRes)
  
  # include the debiased results if cluster effect is included
  if(includeCluster) {
    designRes = list(SRSdat = resSRSdatMod,
                     overSampDat = resDatOverSampMod)
    
    save(file = paste0('kenyaSpatialDesignResultNewTausq', round(tausq, 4), 'UrbRur',
                       includeUrbanRural, 'Cluster', includeCluster, 'debiased', testText, '.RData'), 
         designRes = designRes)
  }
  
  invisible(NULL)
}


# leave one county of data out at a time in order to validate the mercer model versus county level direct estimates
validateBYM2Dat = function(directLogitEsts, directLogitVars, directVars, dat=ed, includeUrbanRural=TRUE, includeCluster=TRUE, 
                           saveResults=TRUE, fileNameRoot="Ed", counties=sort(unique(poppc$County))) {
  
  # fit the full model once, calculating certain validation scores and generating in sample predictions
  print("Fitting full model integrating over hyperparameter uncertainty")
  modelFitFull = runBYM2Dat(dat=dat, includeUrbanRural=includeUrbanRural, includeCluster=includeCluster, saveResults=saveResults, fileNameRoot=fileNameRoot, 
                            doValidation=TRUE, doPredsAtPostMean=FALSE, getPosteriorDensity=FALSE, directLogitEsts=directLogitEsts)
  cpo = modelFitFull$mod$cpo$cpo
  cpoFailure = modelFitFull$mod$cpo$failure
  dic = modelFitFull$mod$dic$dic
  waic = modelFitFull$mod$waic$waic
  countyPredsUrbanInSample = modelFitFull$designRes$predictionsClusterUrban
  countyPredsRuralInSample = modelFitFull$designRes$predictionsClusterRural
  
  # now do leave one county out across validation
  for(i in 1:length(counties)) {
    if(i %% 10 == 1)
      print(paste0("Fitting model with data from county ", i, "/", length(counties), " left out"))
    thisCountyName = counties[i]
    
    # fit model, get all predictions for each areal level and each posterior sample
    fit = runBYM2Dat(dat=dat, includeUrbanRural=includeUrbanRural, includeCluster=includeCluster, 
                     fileNameRoot=fileNameRoot, saveResults=FALSE, 
                     previousResult=modelFitFull$mod, doValidation=FALSE, predCountyI=i)
    
    # get predictive distribution for the left out county
    resUrban = fit$predictionsClusterUrban
    resRural = fit$predictionsClusterRural
    
    if(i == 1) {
      countyPredsUrban = resUrban
      countyPredsRural = resRural
    } else {
      countyPredsUrban = rbind(countyPredsUrban, resUrban)
      countyPredsRural = rbind(countyPredsRural, resRural)
    }
  }
  
  ## since deciding to switch validation to the cluster level, we decided to not take weight averages of the estimates
  # calculate weights further predictions based on inverse direct estimate variances, calculate validation scores
  # weights = 1 / directVars
  # weights = weights / sum(weights)
  
  ## now calculate scoring rules on the cluster level
  countyI = match(dat$admin1, counties)
  isUrban = dat$urban
  
  # calculate scoring rules when leaving out individual counties
  thisLogitMean = countyPredsRural$mean[countyI]
  thisLogitMean[isUrban] = countyPredsUrban$mean[countyI][isUrban]
  thisLogitVar = countyPredsRural$stddev[countyI]^2
  thisLogitVar[isUrban] = countyPredsUrban$stddev[countyI][isUrban]^2
  theseScores = getValidationScores(dat$y / dat$n, 
                                    thisLogitMean, thisLogitVar, 
                                    usePearson=FALSE, n=dat$n, 
                                    urbanVec=dat$urban, filterType="leftOutCounty")
  
  # calculate scoring rules when leaving out clusters
  thisLogitMeanInSample = countyPredsRuralInSample$mean[countyI]
  thisLogitMeanInSample[isUrban] = countyPredsUrbanInSample$mean[countyI][isUrban]
  thisLogitVarInSample = countyPredsRuralInSample$stddev[countyI]^2
  thisLogitVarInSample[isUrban] = countyPredsUrbanInSample$stddev[countyI][isUrban]^2
  theseScoresInSample = getValidationScores(dat$y / dat$n, 
                                            thisLogitMeanInSample, thisLogitVarInSample, 
                                            usePearson=FALSE, n=dat$n, 
                                            urbanVec=dat$urban, filterType="inSample")
  theseScoresInSample$scores = cbind(WAIC=waic, DIC=dic, theseScoresInSample$scores)
  theseScoresInSample$allResults = cbind(WAIC=waic, DIC=dic, theseScoresInSample$allResults)
  
  # save (if necessary) and return results
  bym2ResultsInSample = theseScoresInSample
  bym2ResultsLeaveOutCounty = theseScores
  bym2ResultsLeaveOutCluster = data.frame(CPO=mean(cpo, na.rm=TRUE))
  cpoFailure = mean(cpoFailure[!is.na(cpoFailure)] != 0)
  
  bym2Results = list(countyPredsUrbanCluster=countyPredsUrban, countyPredsRuralCluster=countyPredsRural, 
                     countyPredsUrbanClusterInSample=countyPredsUrbanInSample, countyPredsRuralClusterInSample=countyPredsRuralInSample, 
                     modelFitFull=modelFitFull, 
                     bym2ResultsInSample=bym2ResultsInSample, 
                     bym2ResultsLeaveOutCounty=bym2ResultsLeaveOutCounty, 
                     bym2ResultsLeaveOutCluster=bym2ResultsLeaveOutCluster, 
                     cpoFailure=cpoFailure)
  
  if(saveResults)
    save(bym2Results, 
         file=paste0('bym2', fileNameRoot, 'ValidationAllUrbRur', includeUrbanRural, 'Cluster', includeCluster, '.RData'))
  
  bym2Results
}









