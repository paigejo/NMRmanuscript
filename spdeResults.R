# code for generating table for Jon's talk with model scores

# generate predictive results given the simulated datasets simDataMulti.RData
# includeClustEffect, urbanEffect: control whether or knot to include a 
# cluster random effect in urban fixed effect in the model
# genRegionLevel, keepPixelPreds, genEALevel: control whether to include 
# regional, pixel, and enumeration area level predictions
resultsSPDE = function(nPostSamples=1000, test=FALSE, nTest=2, verbose=TRUE, 
                       includeClustEffect=TRUE, int.strategy="eb", 
                       genRegionLevel=TRUE, keepPixelPreds=TRUE, kmres=5, 
                       genEALevel=TRUE, urbanEffect=TRUE, tausq=0, 
                       saveResults=!test && is.null(maxDataSets), margVar=.15^2, gamma=-1, 
                       beta0=-1.75, loadProgress=FALSE, continuousOnly=TRUE, strictPrior=FALSE, 
                       maxDataSets=NULL, integrateOutCluster=TRUE, seed=123, 
                       strictSpatialPrior=FALSE, effRange=150) {
  
  # Load data
  # load("simDataMulti.RData") # overSampDat, SRSDat
  # load a different 1 of these depending on whether a cluster effect should be included 
  # in the simulation of the data or not (tausq is the cluster effect variance)
  # load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  # load and relevant data
  rangeText = ifelse(effRange == 150, "", "Range50")
  if(!test)
    load(paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", rangeText, ".RData"))
  else
    load(paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0Test", rangeText, ".RData"))
  
  # modify the pixel indices if necessary
  if(kmres != 5) {
    overSampDat = addPixelIToEaDat2(overSampDat, kmres)
    SRSDat = addPixelIToEaDat2(SRSDat, kmres)
  }
  
  eaDat = overSampDat$eaDat
  clustSRS = SRSDat$clustDat
  clustOverSamp = overSampDat$clustDat
  
  if(test)
    maxDataSets = nTest
  if(!is.null(maxDataSets)) {
    clustSRS = lapply(1:maxDataSets, function(i) {clustSRS[[i]]})
    clustOverSamp = lapply(1:maxDataSets, function(i) {clustOverSamp[[i]]})
  }
  
  # set random number seed for generating random number seeds...
  if(!is.null(seed))
    set.seed(seed)
  
  # in the parallel case, we must generate random numbers for each data set, so do it even in the sequential case for consistency
  randomSeedsSRS = sample(1:2000000, length(clustSRS), replace=FALSE)
  randomSeedsOverSamp = sample(1:2000000, length(clustOverSamp), replace=FALSE)
  
  # SRS results are correct, so load those and recompute overSamp results
  # load(paste0("resultsSPDETausq", round(tausq, 4), 
  #             "urbanEffect", as.character(urbanEffect), ".RData"))
  # spdeSRS = resultsSPDEHelper(clustSRS, eaDat, nPostSamples = nPostSamples, verbose=verbose,
  #                             includeClustEffect=includeClustEffect, int.strategy=int.strategy,
  #                             genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds,
  #                             genEALevel=genEALevel, urbanEffect=urbanEffect, 
  #                             genCountLevel=genCountLevel, exactAggregation=exactAggregation)
  # spdeOverSamp = resultsSPDEHelper(clustOverSamp, eaDat, nPostSamples = nPostSamples, verbose=verbose, 
  #                                  includeClustEffect=includeClustEffect, int.strategy=int.strategy, 
  #                                  genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds, 
  #                                  genEALevel=genEALevel, urbanEffect=urbanEffect, 
  #                                  genCountLevel=genCountLevel, exactAggregation=exactAggregation)
  
  strictPriorText = ifelse(strictPrior, "strictPrior", "")
  testText = ifelse(test, "Test", "")
  # integrateClusterText = ifelse(integrateOutCluster, "IntClust", "")
  integrateClusterText = ""
  fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                    round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                    "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, 
                    testText, integrateClusterText, rangeText, ".RData")
  fileNameTemp = paste0("resultsSPDEBetaTemp", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                    round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                    "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, 
                    testText, integrateClusterText, rangeText, ".RData")
  
  if(!loadProgress) {
    print("Generating SRS results")
    spdeSRS = resultsSPDEHelper3(clustSRS, eaDat, nPostSamples = nPostSamples, verbose=verbose,
                                 includeClustEffect=includeClustEffect, int.strategy=int.strategy,
                                 genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds,
                                 genEALevel=genEALevel, urbanEffect=urbanEffect, kmres=kmres, 
                                 continuousOnly=continuousOnly, strictPrior=strictPrior, 
                                 integrateOutCluster=integrateOutCluster, randomSeeds=randomSeedsSRS, 
                                 strictSpatialPrior=strictSpatialPrior)
    
    # save our progress as we go
    if(saveResults)
      save(spdeSRS, file=fileNameTemp)
  }
  else {
    # load our previous progress if necessary
    print("Loading SRS results")
    load(fileNameTemp)
  }
  
  print("Generating urban oversampled results")
  spdeOverSamp = resultsSPDEHelper3(clustOverSamp, eaDat, nPostSamples = nPostSamples, verbose=verbose, 
                                    includeClustEffect=includeClustEffect, int.strategy=int.strategy, 
                                    genRegionLevel=genRegionLevel, keepPixelPreds=keepPixelPreds, 
                                    genEALevel=genEALevel, urbanEffect=urbanEffect, kmres=kmres, 
                                    continuousOnly=continuousOnly, strictPrior=strictPrior, 
                                    integrateOutCluster=integrateOutCluster, randomSeeds=randomSeedsOverSamp, 
                                    strictSpatialPrior=strictSpatialPrior)
  print(paste0("Saving final results under: ", fileName))
  if(saveResults)
    save(spdeSRS, spdeOverSamp, file=fileName)
  print(paste0("Finished saving final results under: ", fileName))
  
  list(spdeSRS=spdeSRS, spdeOverSamp=spdeOverSamp, randomSeedsSRS=randomSeedsSRS, randomSeedsOverSamp=randomSeedsOverSamp)
}

# generate predictive results given a simulated dataset simDataMulti.RData
resultsSPDEHelper = function(clustDatMulti, eaDat, nPostSamples=100, verbose=TRUE, 
                             includeClustEffect=FALSE, int.strategy="eb", 
                             genRegionLevel=TRUE, keepPixelPreds=TRUE, genEALevel=TRUE, 
                             urbanEffect=TRUE, exactAggregation=TRUE, 
                             predictionType=c("mean", "median"), parClust=cl, calcCrps=TRUE) {
  
  # match the requested prediction type with one of the possible options
  predictionType = match.arg(predictionType)
  
  if(calcCrps) {
    # compute truth based on superpopulation
    
    # region
    regions = sort(unique(eaDat$region))
    counties = sort(unique(eaDat$admin1))
    truthByRegion = getTruth("region", eaDat)
    truthByCounty = getTruth("county", eaDat)
    truthByPixel = getTruth("pixel", eaDat)
    truthByEa = getTruth("EA", eaDat)
  } else {
    truthByRegion = NULL
    truthByCounty = NULL
    truthByPixel = NULL
    truthByEa = NULL
  }
  
  nsim = length(clustDatMulti)
  
  # get prediction locations from population grid
  # popGrid = makeInterpPopGrid()
  load("popGrid.RData")
  predCoords = cbind(popGrid$east, popGrid$north)
  predUrban = popGrid$urban
  if(genEALevel) {
    # Must predict at enumeration areas as well. Include enumeration areas as 
    # first rows of prediction coordinates and prediction urban/rural
    predCoords = rbind(cbind(eaDat$east, eaDat$north), predCoords)
    predUrban = c(eaDat$urban, predUrban)
  }
  # we only care about the probability, not counts, so not used except for the purposes 
  # of calling inla:
  # predNs = rep(25, nrow(predCoords))
  predNs = rep(1, nrow(predCoords))
  
  # first make a function for combining the results that can work either in parallel or serial
  combineResults = function(...) {
    results = list(...)[[1]]
    countyResults = do.call("rbind", lapply(results, function(x) {x$countyResults}))
    regionResults = do.call("rbind", lapply(results, function(x) {x$regionResults}))
    pixelResults = do.call("rbind", lapply(results, function(x) {x$pixelResults}))
    eaResults = do.call("rbind", lapply(results, function(x) {x$eaResults}))
    list(countyResults=countyResults, regionResults=regionResults, 
         pixelResults=pixelResults, eaResults=eaResults)
  }
  
  # now make a function for generating the results that can work either in parallel or serial
  mainFunction = function(i, doSink=FALSE) {
    if(doSink)
      sink("log.txt", append=TRUE)
    
    print(paste0("iteration ", i, "/", nsim))
    
    # get the simulated sample
    clustDat = clustDatMulti[[i]]
    
    # get observations from dataset
    obsCoords = cbind(clustDat$east, clustDat$north)
    obsNs = clustDat$n
    obsCounts = clustDat$y
    obsUrban = clustDat$urban
    
    # fit model, get all predictions for each areal level and each posterior sample
    fit = fitSPDEModel(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                       predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=nPostSamples, 
                       verbose = verbose, clusterEffect=includeClustEffect, 
                       int.strategy=int.strategy, genRegionLevel=genRegionLevel, 
                       keepPixelPreds=keepPixelPreds, genEALevel=genEALevel, 
                       urbanEffect=urbanEffect, link=1, predictionType=predictionType, 
                       exactAggregation=exactAggregation, genCountLevel=genCountLevel, 
                       eaDat=eaDat, truthByCounty=truthByCounty, truthByRegion=truthByRegion, 
                       truthByPixel=truthByPixel)
    print(paste0("Fit completed: iteration ", i, "/", nsim))
    countyPredMat = fit$countyPredMat
    regionPredMat = fit$regionPredMat
    pixelPredMat = fit$pixelPredMat
    eaPredMat = fit$eaPredMat
    eaMarginals = fit$eaMarginals
    
    ## calculate model fit properties for each level of predictions
    
    # enumeration area level
    if(genEALevel) {
      # uncertainties and coverages that don't include binomial variation
      thislowerEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.1)}))
      thisupperEA = logit(apply(eaPredMat, 1, function(x) {quantile(x, probs=.9)}))
      thisvar.estEA = apply(logit(eaPredMat), 1, var)
      thisRightRejectEA = NA
      thisLeftRejectEA = NA
      
      # get the marginal "binomial" densities at each location
      # n = 25
      n = max(eaDat$n)
      binomProb = function(k, i) {
        inla.emarginal(function(logitP) {dbinom(k, eaDat$n[i], expit(logitP))}, eaMarginals[[i]])
      }
      
      ## make highest density coverage interval on count scale
      # generate the "binomial" pmfs for each marginal
      probs = matrix(mapply(binomProb, 
                            k = matrix(0:n, nrow=length(eaMarginals), ncol=n + 1, byrow = TRUE), 
                            i = matrix(1:length(eaMarginals), nrow=length(eaMarginals), ncol=n + 1)), 
                     nrow=length(eaMarginals), ncol=n + 1)
      
      # now generate the summary statistics about the marginals
      if(predictionType == "mean")
        thisu1mEABinom = rowSums(sweep(probs, 2, STATS = 0:n, "*"))
      else {
        thisu1mEABinom = rowSums(sweep(probs, 2, STATS = 0:n, "*"))
        warning("Using median predictions is not supported for discrete distributions. Used mean prediction instead")
      }
      
      # calculate PREDICTIVE variance on the probability scale
      # thisvar.estEA = rowSums(sweep(probs, 2, STATS = (0:n)^2, "*")) - thisu1mEA^2
      predVar = function(eaMarginal, power=2) {
        inla.emarginal(function(logitP) {expit(logitP)^power}, eaMarginal)
      }
      thisvar.estEA = sapply(eaMarginals, predVar) - sapply(eaMarginals, predVar, power=1)^2
      
      # get credible intervals
      intervals = apply(probs, 1, getQuantileInterval)
      thislowerEABinom = intervals[1,]
      thisupperEABinom = intervals[2,]
      thisRightRejectEABinom = intervals[3,]
      thisLeftRejectEABinom = intervals[4,]
    }
    else {
      thisu1mEA = NA
      thislowerEA = NA
      thisupperEA = NA
      thisvar.estEA = NA
      
      thislowerEABinom = NA
      thisupperEABinom = NA
      thisRightRejectEABinom = NA
      thisLeftRejectEABinom = NA
    }
    print(paste0("EA results generated: iteration ", i, "/", nsim))
    
    # pixel level
    if(keepPixelPreds) {
      # first get credible interval when not including binomial variation
      thislowerPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.1)}))
      thisupperPixel = logit(apply(pixelPredMat, 1, function(x) {quantile(x, probs=.9)}))
      
      # now get rest of the scoring rules when including binomial variation
      thisu1mPixel = pixelPredMat$preds
      thisvar.estPixel = pixelPredMat$vars
      thisupperPixel = pixelPredMat$upper
      thislowerPixel = pixelPredMat$lower
      thisRightRejectPixel = pixelPredMat$rightRejectProb
      thisLeftRejectPixel = pixelPredMat$leftRejectProb
    }
    else {
      thisu1mPixel = NA
      thislowerPixel = NA
      thisupperPixel = NA
      thisvar.estPixel = NA
      thisRightRejectPixel = NA
      thisLeftRejectPixel = NA
    }
    print(paste0("Pixel results generated: iteration ", i, "/", nsim))
    
    # County level (county level results are always included)
    if(!genCountLevel) {
      if(predictionType == "mean")
        thisu1mCounty = rowMeans(countyPredMat)
      else
        thisu1mCounty = apply(countyPredMat, 1, median)
      thislowerCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.1)}))
      thisupperCounty = logit(apply(countyPredMat, 1, function(x) {quantile(x, probs=.9)}))
      thisvar.estCounty = apply(logit(countyPredMat), 1, var)
      thisRightRejectCounty = NA
      thisLeftRejectCounty = NA
    } else {
      thisu1mCounty = countyPredMat$preds
      thisvar.estCounty = countyPredMat$vars
      thisupperCounty = countyPredMat$upper
      thislowerCounty = countyPredMat$lower
      thisRightRejectCounty = countyPredMat$rightRejectProb
      thisLeftRejectCounty = countyPredMat$leftRejectProb
    }
    print(paste0("County results generated: iteration ", i, "/", nsim))
    
    # region level
    if(genRegionLevel) {
      if(!genCountLevel) {
        if(predictionType == "mean")
          thisu1mRegion = rowMeans(regionPredMat)
        else
          thisu1mRegion = apply(regionPredMat, 1, median)
        thislowerRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.1)}))
        thisupperRegion = logit(apply(regionPredMat, 1, function(x) {quantile(x, probs=.9)}))
        thisvar.estRegion = apply(logit(regionPredMat), 1, var)
        thisRightRejectRegion = NA
        thisLeftRejectRegion = NA
      } else {
        thisu1mRegion = regionPredMat$preds
        thisvar.estRegion = regionPredMat$vars
        thisupperRegion = regionPredMat$upper
        thislowerRegion = regionPredMat$lower
        thisRightRejectRegion = regionPredMat$rightRejectProb
        thisLeftRejectRegion = regionPredMat$leftRejectProb
      }
    }
    else {
      thisu1mRegion = NA
      thislowerRegion = NA
      thisupperRegion = NA
      thisvar.estRegion = NA
      thisRightRejectRegion = NA
      thisLeftRejectRegion = NA
    }
    print(paste0("Region results generated: iteration ", i, "/", nsim))
    
    # collect results in a list, one element for each simulation
    thisCountyResults = data.frame(list(admin1=counties, 
                                        u1m.spde=thisu1mCounty, lower.spde=thislowerCounty, 
                                        upper.spde=thisupperCounty, logit.est.spde=logit(thisu1mCounty), 
                                        var.est.spde=thisvar.estCounty, 
                                        leftReject.spde=thisLeftRejectCounty, 
                                        rightReject.spde=thisRightRejectCounty))
    thisRegionResults = data.frame(list(region=regions, 
                                        u1m.spde=thisu1mRegion, lower.spde=thislowerRegion, 
                                        upper.spde=thisupperRegion, logit.est.spde=logit(thisu1mRegion), 
                                        var.est.spde=thisvar.estRegion, 
                                        leftReject.spde=thisLeftRejectRegion, 
                                        rightReject.spde=thisRightRejectRegion))
    thisPixelResults = data.frame(list(u1m.spde=thisu1mPixel, lower.spde=thislowerPixel, 
                                       upper.spde=thisupperPixel, logit.est.spde=logit(thisu1mPixel), 
                                       var.est.spde=thisvar.estPixel, 
                                       leftReject.spde=thisLeftRejectPixel, 
                                       rightReject.spde=thisRightRejectPixel))
    thisEaResults = data.frame(list(u1m.spde=thisu1mEA, lower.spde=thislowerEA, 
                                    upper.spde=thisupperEA, logit.est.spde=logit(thisu1mEA), 
                                    var.est.spde=thisvar.estEA, 
                                    leftReject.spde=thisLeftRejectEA, 
                                    rightReject.spde=thisRightRejectEA))
    list(countyResults=thisCountyResults, regionResults=thisRegionResults, 
         pixelResults=thisPixelResults, eaResults=thisEaResults)
  }
  
  # compute Bias & MSE & mean(Var) & 80\% coverage for each simulation
  if(is.null(parClust)) {
    # sequential version
    countyResults = as.list(1:nsim)
    regionResults = as.list(1:nsim)
    pixelResults = as.list(1:nsim)
    eaResults = as.list(1:nsim)
    for(i in 1:nsim) {
      results = mainFunction(i, FALSE)
      
      # collect results in a list, one element for each simulation
      countyResults[[i]] = results$CountyResults
      regionResults[[i]] = results$regionResults
      pixelResults[[i]] = results$pixelResults
      eaResults[[i]] = results$eaResults
    }
  } else {
    # parallel version
    
    results = foreach(i = 1:nsim, .combine=combineResults, .verbose=TRUE, .multicombine=TRUE) %dopar% {
      mainFunction(i, TRUE)
    }
    
    # separate results into the different aggregation levels
    countyResults = results$countyResults
    regionResults = results$regionResults
    pixelResults = results$pixelResults
    eaResults = results$eaResults
  }
  list(countyResults=countyResults, regionResults=regionResults, pixelResults=pixelResults, 
       eaResults=eaResults)
}

# Same as resultsSPDEHelper, but updated to use fitSPDEModel2 instead of fitSPDEModel, with 
# reorganized scoring rule calculations
resultsSPDEHelper2 = function(clustDatMulti, eaDat, nPostSamples=100, verbose=FALSE, 
                              includeClustEffect=FALSE, int.strategy="eb", 
                              genRegionLevel=TRUE, keepPixelPreds=TRUE, genEALevel=TRUE, 
                              urbanEffect=TRUE, kmres=5, nSamplePixel=10, 
                              predictionType=c("mean", "median"), parClust=cl, calcCrps=TRUE, 
                              significance=.8) {
  
  # match the requested prediction type with one of the possible options
  predictionType = match.arg(predictionType)
  
  if(calcCrps) {
    # compute truth based on superpopulation
    print("Calculating empirical rates at the desired aggregation levels")
    
    # region
    regions = sort(unique(eaDat$region))
    counties = sort(unique(eaDat$admin1))
    truthByCounty = getTruth("county", eaDat)
    if(genRegionLevel)
      truthByRegion = getTruth("region", eaDat)
    else
      truthByRegion = NULL
    if(keepPixelPreds)
      truthByPixel = getTruth("pixel", eaDat)
    else
      truthByPixel = NULL
    if(genEALevel)
      truthByEa = getTruth("EA", eaDat)
    else
      truthByEa = NULL
  } else {
    truthByRegion = NULL
    truthByCounty = NULL
    truthByPixel = NULL
    truthByEa = NULL
  }
  
  nsim = length(clustDatMulti)
  
  # get prediction locations from population grid
  if(kmres == 5)
    load("popGrid.RData")
  else
    popGrid = makeInterpPopGrid(kmres)
  
  predCoords = cbind(popGrid$east, popGrid$north)
  predUrban = popGrid$urban
  if(genEALevel) {
    # Must predict at enumeration areas as well. Include enumeration areas as 
    # first rows of prediction coordinates and prediction urban/rural
    predCoords = rbind(cbind(eaDat$east, eaDat$north), predCoords)
    predUrban = c(eaDat$urban, predUrban)
  }
  # we only care about the probability, not counts, so not used except for the purposes 
  # of calling inla:
  # predNs = rep(25, nrow(predCoords))
  predNs = rep(1, nrow(predCoords))
  
  # first make a function for combining the results that can work either in parallel or serial
  combineResults = function(...) {
    print("Combining results...")
    
    results = list(...)[[1]]
    
    # scoring rules
    scoresEaInexact = do.call("rbind", lapply(results, function(x) {x$scoresEaInexact}))
    scoresEaExact = do.call("rbind", lapply(results, function(x) {x$scoresEaExact}))
    scoresEaExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresEaExactBVar}))
    scoresPixelInexact = do.call("rbind", lapply(results, function(x) {x$scoresPixelInexact}))
    scoresPixelExact = do.call("rbind", lapply(results, function(x) {x$scoresPixelExact}))
    scoresPixelExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresPixelExactBVar}))
    scoresCountyInexact = do.call("rbind", lapply(results, function(x) {x$scoresCountyInexact}))
    scoresCountyExact = do.call("rbind", lapply(results, function(x) {x$scoresCountyExact}))
    scoresCountyExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresCountyExactBVar}))
    scoresRegionInexact = do.call("rbind", lapply(results, function(x) {x$scoresRegionInexact}))
    scoresRegionExact = do.call("rbind", lapply(results, function(x) {x$scoresRegionExact}))
    scoresRegionExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresRegionExactBVar}))
    
    # parameter estimates
    interceptSummary = do.call("rbind", lapply(results, function(x) {x$interceptSummary}))
    interceptSummary = colMeans(interceptSummary)
    if(urbanEffect) {
      urbanSummary = do.call("rbind", lapply(results, function(x) {x$urbanSummary}))
      urbanSummary = colMeans(urbanSummary)
    } else
      urbanSummary = NULL
    rangeSummary = do.call("rbind", lapply(results, function(x) {x$rangeSummary}))
    rangeSummary = colMeans(rangeSummary)
    sdSummary = do.call("rbind", lapply(results, function(x) {x$sdSummary}))
    sdSummary = colMeans(sdSummary)
    varSummary = do.call("rbind", lapply(results, function(x) {x$varSummary}))
    varSummary = colMeans(varSummary)
    if(includeClustEffect) {
      nuggetVarSummary = do.call("rbind", lapply(results, function(x) {x$nuggetVarSummary}))
      nuggetVarSummary = colMeans(nuggetVarSummary)
      nuggetSDSummary = do.call("rbind", lapply(results, function(x) {x$nuggetSDSummary}))
      nuggetSDSummary = colMeans(nuggetSDSummary)
    } else {
      nuggetVarSummary = NULL
      nuggetSDSummary = NULL
    }
    
    list(scoresEaInexact=scoresEaInexact, scoresEaExact=scoresEaExact, scoresEaExactBVar=scoresEaExactBVar, 
         scoresPixelInexact=scoresPixelInexact, scoresPixelExact=scoresPixelExact, scoresPixelExactBVar=scoresPixelExactBVar, 
         scoresCountyInexact=scoresCountyInexact, scoresCountyExact=scoresCountyExact, scoresCountyExactBVar=scoresCountyExactBVar, 
         scoresRegionInexact=scoresRegionInexact, scoresRegionExact=scoresRegionExact, scoresRegionExactBVar=scoresRegionExactBVar, 
         interceptSummary=interceptSummary, urbanSummary=urbanSummary, 
         rangeSummary=rangeSummary, varSummary=varSummary, sdSummary=sdSummary, 
         nuggetVarSummary=nuggetVarSummary, nuggetSDSummary=nuggetSDSummary)
  }
  
  # now make a function for generating the results that can work either in parallel or serial
  mainFunction = function(i, doSink=FALSE) {
    if(doSink)
      sink("log.txt", append=TRUE)
    
    print(paste0("iteration ", i, "/", nsim))
    
    # get the simulated sample
    clustDat = clustDatMulti[[i]]
    
    # get observations from dataset
    obsCoords = cbind(clustDat$east, clustDat$north)
    obsNs = clustDat$n
    obsCounts = clustDat$y
    obsUrban = clustDat$urban
    
    # fit model, get all predictions for each areal level and each posterior sample
    fit = fitSPDEModel2(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                        predUrban, clusterIndices=clustDat$eaIs, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=nPostSamples, 
                        verbose = verbose, clusterEffect=includeClustEffect, 
                        int.strategy=int.strategy, genRegionLevel=genRegionLevel, 
                        keepPixelPreds=keepPixelPreds, genEALevel=genEALevel, 
                        urbanEffect=urbanEffect, link=1, predictionType=predictionType, 
                        eaDat=eaDat, truthByCounty=truthByCounty, truthByRegion=truthByRegion, 
                        truthByPixel=truthByPixel, nSamplePixel=nSamplePixel, 
                        significance=significance)
    print(paste0("Fit completed: iteration ", i, "/", nsim))
    countyPreds = fit$countyPreds
    regionPreds = fit$regionPreds
    pixelPreds = fit$pixelPreds
    eaPreds = fit$eaPreds
    
    ## calculate model fit properties for each level of predictions
    
    # enumeration area level
    if(genEALevel) {
      print(paste0("EA results generating: iteration ", i, "/", nsim))
      
      ### generate enumeration area level predictions
      # calculate logit scale estimates and variances
      
      # scoring rules for all aggregation models
      scoresEaInexact = getScoresSPDE(truthByEa$truth, truthByEa$n, eaPreds$logitEst, 
                                      eaPreds$est, NULL, bVar=FALSE, probMat=eaPreds$eaPredMat, 
                                      nsim=nSamplePixel, logitVar=eaPreds$logitVar, 
                                      logitL=eaPreds$logitLower, logitU=eaPreds$logitUpper)
      cat(".")
      scoresEaExact = scoresEaInexact
      cat(".")
      scoresEaExactBVar = getScoresSPDE(truthByEa$truth, truthByEa$n, eaPreds$logitEst, 
                                        eaPreds$est, NULL, bVar=TRUE, probMat=eaPreds$eaPredMat, 
                                        nsim=nSamplePixel, logitVar=eaPreds$logitVar, 
                                        logitL=eaPreds$logitLower, logitU=eaPreds$logitUpper)
      cat(".")
    }
    else {
      scoresEaInexact = NA
      scoresEaExact = NA
      scoresEaExactBVar = NA
    }
    
    # pixel level
    if(keepPixelPreds) {
      print(paste0("Pixel results generating: iteration ", i, "/", nsim))
      
      # calculate estimates, but convert to the logit scale (estimates are the same with or without binomial variation)
      out = sapply(pixelPreds$pixelMarginalsWithData, function(m) {inla.emarginal(function(logitP) {c(logitP, logitP^2, expit(logitP), expit(logitP)^2)}, m)})
      logitEst = pixelPreds$logitEst
      logitVar = pixelPreds$logitVar
      est = pixelPreds$est
      logitLower = pixelPreds$logitLower
      logitUpper = pixelPreds$logitUpper
      
      # scoring rules for all aggregation models
      scoresPixelInexact = getScoresSPDE(truthByPixel$truth, truthByPixel$n, logitEst, 
                                         logitEst, est, bVar=FALSE, probMat=pixelPreds$pixelPredMatInexact, 
                                         nsim=nSamplePixel, logitVar=pixelPreds$logitVar, 
                                         logitL=pixelPreds$logitLower, logitU=pixelPreds$logitUpper)
      cat(".")
      scoresPixelExact = getScoresSPDE(truthByPixel$truth, truthByPixel$n, logitEst, 
                                       logitEst, est, bVar=FALSE, probMat=pixelPreds$pixelPredMatExact, 
                                       nsim=nSamplePixel)
      cat(".")
      scoresPixelExactBVar = getScoresSPDE(truthByPixel$truth, truthByPixel$n, logitEst, 
                                           logitEst, est, bVar=TRUE, pmfs=pixelPreds$pixelPredMatExactB, 
                                           nsim=nSamplePixel)
      cat(".")
    }
    else {
      scoresPixelInexact = NA
      scoresPixelExact = NA
      scoresPixelExactBVar = NA
    }
    
    # County level (county level results are always included)
    print(paste0("County results generating: iteration ", i, "/", nsim))
    
    # first generate the county estimates
    thisu1mCountyInexact = logit(rowMeans(countyPreds$countyPredMatInexact))
    thisu1mCountyExact = logit(rowMeans(countyPreds$countyPredMatExact))
    thisu1mCountyExactBVar = logit(sapply(countyPreds$countyPredMatExactB, 
                                          function(pmf) {sum((0:(length(pmf) - 1)) * pmf) / length(pmf)}))
    
    # scoring rules for all aggregation models
    scoresCountyInexact = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyInexact, 
                                       expit(thisu1mCountyInexact), NULL, bVar=FALSE, probMat=countyPreds$countyPredMatInexact)
    cat(".")
    scoresCountyExact = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyExact, 
                                     expit(thisu1mCountyExact), NULL, bVar=FALSE, probMat=countyPreds$countyPredMatExact)
    cat(".")
    scoresCountyExactBVar = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyExactBVar, 
                                         expit(thisu1mCountyExact), NULL, bVar=TRUE, pmfs=countyPreds$countyPredMatExactB)
    cat(".")
    
    # region level
    if(genRegionLevel) {
      print(paste0("Region results generating: iteration ", i, "/", nsim))
      
      # first generate the region estimates
      thisu1mRegionInexact = logit(rowMeans(regionPreds$regionPredMatInexact))
      thisu1mRegionExact = logit(rowMeans(regionPreds$regionPredMatExact))
      thisu1mRegionExactB = logit(sapply(regionPreds$regionPredMatExactB, 
                                         function(pmf) {sum((0:(length(pmf) - 1)) * pmf) / length(pmf)}))
      
      # scoring rules for all aggregation models
      scoresRegionInexact = getScoresSPDE(truthByRegion$truth, truthByRegion$n, thisu1mRegionInexact, 
                                          expit(thisu1mRegionInexact), NULL, bVar=FALSE, probMat=regionPreds$regionPredMatInexact)
      cat(".")
      scoresRegionExact = getScoresSPDE(truthByRegion$truth, truthByRegion$n, thisu1mRegionExact, 
                                        expit(thisu1mRegionExact), NULL, bVar=FALSE, probMat=regionPreds$regionPredMatExact)
      cat(".")
      scoresRegionExactBVar = getScoresSPDE(truthByRegion$truth, truthByRegion$n, thisu1mRegionExactB, 
                                            expit(thisu1mRegionExact), NULL, bVar=TRUE, pmfs=regionPreds$regionPredMatExactB)
      cat(".")
    }
    else {
      scoresRegionInexact = NA
      scoresRegionExact = NA
      scoresRegionExactBVar = NA
    }
    
    ## collect parameter means, sds, and quantiles
    # for fixed effects
    print(paste0("Parameter summaries generating: iteration ", i, "/", nsim))
    interceptQuants = inla.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.fixed[[1]])
    interceptMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.fixed[[1]])
    interceptSummary = c(interceptMoments[1], sqrt(interceptMoments[2] - interceptMoments[1]^2), interceptMoments[2] - interceptMoments[1]^2, interceptQuants, interceptQuants[3] - interceptQuants[1])
    
    urbanSummary = NULL
    if(urbanEffect) {
      urbanQuants = inla.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.fixed[[2]])
      urbanMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.fixed[[2]])
      urbanSummary = c(urbanMoments[1], sqrt(urbanMoments[2] - urbanMoments[1]^2), urbanMoments[2] - urbanMoments[1]^2, urbanQuants, urbanQuants[3] - urbanQuants[1])
    }
    
    # for hyperparameters
    rangeQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[1]])
    sdQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[2]])
    varQuants = sdQuants^2
    varMarg = my.tmarginal(function(x) {x^2}, fit$mod$marginals.hyperpar[[2]])
    rangeMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.hyperpar[[1]])
    rangeSummary = c(rangeMoments[1], sqrt(rangeMoments[2] - rangeMoments[1]^2), rangeMoments[2] - rangeMoments[1]^2, rangeQuants, rangeQuants[3] - rangeQuants[1])
    sdMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.hyperpar[[2]])
    sdSummary = c(sdMoments[1], sqrt(sdMoments[2] - sdMoments[1]^2), sdMoments[2] - sdMoments[1]^2, sdQuants, sdQuants[3] - sdQuants[1])
    varMoments = inla.emarginal(function(x) {c(x, x^2)}, varMarg)
    varSummary = c(varMoments[1], sqrt(varMoments[2] - varMoments[1]^2), varMoments[2] - varMoments[1]^2, varQuants, varQuants[3] - varQuants[1])
    nuggetVarSummary = NULL
    nuggetSDSummary = NULL
    if(includeClustEffect) {
      nuggetPrecQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[3]])
      nuggetVarQuants = 1/nuggetPrecQuants
      nuggetSDQuants = sqrt(nuggetVarQuants)
      nuggetVarMarg = my.tmarginal(function(x) {1/x}, fit$mod$marginals.hyperpar[[3]])
      nuggetSDMarg = my.tmarginal(function(x) {1/sqrt(x)}, fit$mod$marginals.hyperpar[[3]])
      nuggetVarMoments = inla.emarginal(function(x) {c(x, x^2)}, nuggetVarMarg)
      nuggetSDMoments = inla.emarginal(function(x) {c(x, x^2)}, nuggetSDMarg)
      
      nuggetVarSummary = c(nuggetVarMoments[1], sqrt(nuggetVarMoments[2] - nuggetVarMoments[1]^2), nuggetVarMoments[2] - nuggetVarMoments[1]^2, nuggetVarQuants, nuggetVarQuants[3] - nuggetVarQuants[1])
      nuggetSDSummary = c(nuggetSDMoments[1], sqrt(nuggetSDMoments[2] - nuggetSDMoments[1]^2), nuggetSDMoments[2] - nuggetSDMoments[1]^2, nuggetSDQuants, nuggetSDQuants[3] - nuggetSDQuants[1])
    }
    
    res = list(scoresEaInexact=scoresEaInexact, scoresEaExact=scoresEaExact, scoresEaExactBVar=scoresEaExactBVar, 
               scoresPixelInexact=scoresPixelInexact, scoresPixelExact=scoresPixelExact, scoresPixelExactBVar=scoresPixelExactBVar, 
               scoresCountyInexact=scoresCountyInexact, scoresCountyExact=scoresCountyExact, scoresCountyExactBVar=scoresCountyExactBVar, 
               scoresRegionInexact=scoresRegionInexact, scoresRegionExact=scoresRegionExact, scoresRegionExactBVar=scoresRegionExactBVar, 
               interceptSummary=interceptSummary, urbanSummary=urbanSummary, 
               rangeSummary=rangeSummary, varSummary=varSummary, sdSummary=sdSummary, 
               nuggetVarSummary=nuggetVarSummary, nuggetSDSummary=nuggetSDSummary)
    
    print(paste0("Completed iteration ", i, "/", nsim))
    res
  }
  
  # compute Bias & MSE & mean(Var) & 80\% coverage for each simulation
  if(is.null(parClust)) {
    # sequential version
    scoresEaInexact=as.list(1:nsim)
    scoresEaExact=as.list(1:nsim)
    scoresEaExactBVar=as.list(1:nsim)
    scoresPixelInexact=as.list(1:nsim)
    scoresPixelExact=as.list(1:nsim)
    scoresPixelExactBVar=as.list(1:nsim)
    scoresCountyInexact=as.list(1:nsim)
    scoresCountyExact=as.list(1:nsim)
    scoresCountyExactBVar=as.list(1:nsim)
    scoresRegionInexact=as.list(1:nsim)
    scoresRegionExact=as.list(1:nsim)
    scoresRegionExactBVar=as.list(1:nsim)
    
    interceptSummary=as.list(1:nsim)
    urbanSummary=as.list(1:nsim)
    rangeSummary=as.list(1:nsim)
    varSummary=as.list(1:nsim)
    sdSummary=as.list(1:nsim)
    nuggetVarSummary=as.list(1:nsim)
    nuggetSDSummary=as.list(1:nsim)
    for(i in 1:nsim) {
      results = mainFunction(i, FALSE)
      
      # collect results in a list, one element for each simulation
      scoresEaInexact[[i]]=results$scoresEaInexact
      scoresEaExact[[i]]=results$scoresEaExact
      scoresEaExactBVar[[i]]=results$scoresEaExactBVar
      scoresPixelInexact[[i]]=results$scoresPixelInexact
      scoresPixelExact[[i]]=results$scoresPixelExact
      scoresPixelExactBVar[[i]]=results$scoresPixelExactBVar
      scoresCountyInexact[[i]]=results$scoresCountyInexact
      scoresCountyExact[[i]]=results$scoresCountyExact
      scoresCountyExactBVar[[i]]=results$scoresCountyExactBVar
      scoresRegionInexact[[i]]=results$scoresRegionInexact
      scoresRegionExact[[i]]=results$scoresRegionExact
      scoresRegionExactBVar[[i]]=results$scoresRegionExactBVar
      
      interceptSummary[[i]]=results$interceptSummary
      urbanSummary[[i]]=results$urbanSummary
      rangeSummary[[i]]=results$rangeSummary
      varSummary[[i]]=results$varSummary
      sdSummary[[i]]=results$sdSummary
      nuggetVarSummary[[i]]=results$nuggetVarSummary
      nuggetSDSummary[[i]]=results$nuggetSDSummary
    }
    
    print("Combining results...")
    scoresEaInexact = do.call("rbind",  scoresEaInexact)
    scoresEaExact = do.call("rbind", scoresEaExact)
    scoresEaExactBVar = do.call("rbind", scoresEaExactBVar)
    scoresPixelInexact = do.call("rbind", scoresPixelInexact)
    scoresPixelExact = do.call("rbind", scoresPixelExact)
    scoresPixelExactBVar = do.call("rbind", scoresPixelExactBVar)
    scoresCountyInexact = do.call("rbind", scoresCountyInexact)
    scoresCountyExact = do.call("rbind", scoresCountyExact)
    scoresCountyExactBVar = do.call("rbind", scoresCountyExactBVar)
    scoresRegionInexact = do.call("rbind", scoresRegionInexact)
    scoresRegionExact = do.call("rbind", scoresRegionExact)
    scoresRegionExactBVar = do.call("rbind", scoresRegionExactBVar)
    
    interceptSummary = do.call("rbind", interceptSummary)
    interceptSummary = colMeans(interceptSummary)
    if(urbanEffect) {
      urbanSummary = do.call("rbind", urbanSummary)
      urbanSummary = colMeans(urbanSummary)
    } else {
      urbanSummary = NULL
    }
    rangeSummary = do.call("rbind", rangeSummary)
    rangeSummary = colMeans(rangeSummary)
    sdSummary = do.call("rbind", sdSummary)
    sdSummary = colMeans(sdSummary)
    varSummary = do.call("rbind", varSummary)
    varSummary = colMeans(varSummary)
    if(includeClustEffect) {
      nuggetVarSummary = do.call("rbind", nuggetVarSummary)
      nuggetVarSummary = colMeans(nuggetVarSummary)
      nuggetSDSummary = do.call("rbind", nuggetSDSummary)
      nuggetSDSummary = colMeans(nuggetSDSummary)
    } else {
      nuggetVarSummary = NULL
      nuggetSDSummary = NULL
    }
  } else {
    # parallel version
    
    results = foreach(i = 1:nsim, .combine=combineResults, .verbose=TRUE, .multicombine=TRUE) %dopar% {
      mainFunction(i, TRUE)
    }
    
    # separate results into the different aggregation levels
    scoresEaInexact=results$scoresEaInexact
    scoresEaExact=results$scoresEaExact
    scoresEaExactBVar=results$scoresEaExactBVar
    scoresPixelInexact=results$scoresPixelInexact
    scoresPixelExact=results$scoresPixelExact
    scoresPixelExactBVar=results$scoresPixelExactBVar
    scoresCountyInexact=results$scoresCountyInexact
    scoresCountyExact=results$scoresCountyExact
    scoresCountyExactBVar=results$scoresCountyExactBVar
    scoresRegionInexact=results$scoresRegionInexact
    scoresRegionExact=results$scoresRegionExact
    scoresRegionExactBVar=results$scoresRegionExactBVar
    
    interceptSummary=results$interceptSummary
    if(urbanEffect) {
      urbanSummary=results$urbanSummary
    } else {
      urbanSummary = NULL
    }
    rangeSummary=results$rangeSummary
    varSummary=results$varSummary
    sdSummary=results$sdSummary
    if(includeClustEffect) {
      nuggetVarSummary=results$nuggetVarSummary
      nuggetSDSummary=results$nuggetSDSummary
    } else {
      nuggetVarSummary=NULL
      nuggetSDSummary=NULL
    }
  }
  list(scoresEaInexact=scoresEaInexact, scoresEaExact=scoresEaExact, scoresEaExactBVar=scoresEaExactBVar, 
       scoresPixelInexact=scoresPixelInexact, scoresPixelExact=scoresPixelExact, scoresPixelExactBVar=scoresPixelExactBVar, 
       scoresCountyInexact=scoresCountyInexact, scoresCountyExact=scoresCountyExact, scoresCountyExactBVar=scoresCountyExactBVar, 
       scoresRegionInexact=scoresRegionInexact, scoresRegionExact=scoresRegionExact, scoresRegionExactBVar=scoresRegionExactBVar, 
       interceptSummary=interceptSummary, urbanSummary=urbanSummary, 
       rangeSummary=rangeSummary, varSummary=varSummary, sdSummary=sdSummary, 
       nuggetVarSummary=nuggetVarSummary, nuggetSDSummary=nuggetSDSummary)
}

# Same as resultsSPDEHelper, but updated to use fitSPDEModel2 instead of fitSPDEModel, with 
# reorganized scoring rule calculations
resultsSPDEHelper3 = function(clustDatMulti, eaDat, nPostSamples=100, verbose=FALSE, 
                              includeClustEffect=FALSE, int.strategy="eb", 
                              genRegionLevel=TRUE, keepPixelPreds=TRUE, genEALevel=TRUE, 
                              urbanEffect=TRUE, kmres=5, nSamplePixel=nPostSamples, 
                              predictionType=c("mean", "median"), parClust=cl, calcCrps=TRUE, 
                              significance=.8, continuousOnly=FALSE, strictPrior=TRUE, 
                              integrateOutCluster=TRUE, adjustPopSurface=TRUE, randomSeeds=NULL, 
                              strictSpatialPrior=FALSE) {
  
  # generate random seeds for each data set
  if(is.null(randomSeeds))
    randomSeeds = sample(1:2000000, length(clustDatMulti), replace=FALSE)
  
  # match the requested prediction type with one of the possible options
  predictionType = match.arg(predictionType)
  
  if(calcCrps) {
    # compute truth based on superpopulation
    print("Calculating empirical rates at the desired aggregation levels")
    
    # region
    regions = sort(unique(eaDat$region))
    counties = sort(unique(eaDat$admin1))
    truthByCounty = getTruth("county", eaDat)
    if(genRegionLevel)
      truthByRegion = getTruth("region", eaDat)
    else
      truthByRegion = NULL
    if(keepPixelPreds)
      truthByPixel = getTruth("pixel", eaDat)
    else
      truthByPixel = NULL
    if(genEALevel)
      truthByEa = getTruth("EA", eaDat)
    else
      truthByEa = NULL
  } else {
    truthByRegion = NULL
    truthByCounty = NULL
    truthByPixel = NULL
    truthByEa = NULL
  }
  
  nsim = length(clustDatMulti)
  
  # get prediction locations from population grid
  popGridAdjusted = NULL
  if(kmres == 5 && adjustPopSurface) {
    load("popGridAdjusted.RData")
    popGridAdjusted = popGrid
    load("popGrid.RData")
  }
  else if(kmres == 5 && !adjustPopSurface)
    load("popGrid.RData")
  else
    popGrid = makeInterpPopGrid(kmres, adjustPopSurface)
  
  predCoords = cbind(popGrid$east, popGrid$north)
  predUrban = popGrid$urban
  if(genEALevel) {
    # Must predict at enumeration areas as well. Include enumeration areas as 
    # first rows of prediction coordinates and prediction urban/rural
    predCoords = rbind(cbind(eaDat$east, eaDat$north), predCoords)
    predUrban = c(eaDat$urban, predUrban)
  }
  # we only care about the probability, not counts, so not used except for the purposes 
  # of calling inla:
  # predNs = rep(25, nrow(predCoords))
  predNs = rep(1, nrow(predCoords))
  
  # first make a function for combining the results that can work either in parallel or serial
  combineResults = function(...) {
    print("Combining results...")
    
    # not sure why the [[1]] is necessary in the non-parallel case
    if(!doParallel)
      results = list(...)[[1]]
    else
      results = list(...)
    
    # scoring rules
    scoresEaExact = do.call("rbind", lapply(results, function(x) {x$scoresEaExact}))
    scoresEaExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresEaExactBVar}))
    scoresPixelInexact = do.call("rbind", lapply(results, function(x) {x$scoresPixelInexact}))
    scoresPixelExact = do.call("rbind", lapply(results, function(x) {x$scoresPixelExact}))
    scoresPixelExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresPixelExactBVar}))
    scoresCountyInexact = do.call("rbind", lapply(results, function(x) {x$scoresCountyInexact}))
    scoresCountyInexactUnintegrated = do.call("rbind", lapply(results, function(x) {x$scoresCountyInexactUnintegrated}))
    scoresCountyInexactUnadjusted = do.call("rbind", lapply(results, function(x) {x$scoresCountyInexactUnadjusted}))
    scoresCountyExact = do.call("rbind", lapply(results, function(x) {x$scoresCountyExact}))
    scoresCountyExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresCountyExactBVar}))
    scoresRegionInexact = do.call("rbind", lapply(results, function(x) {x$scoresRegionInexact}))
    scoresRegionExact = do.call("rbind", lapply(results, function(x) {x$scoresRegionExact}))
    scoresRegionExactBVar = do.call("rbind", lapply(results, function(x) {x$scoresRegionExactBVar}))
    
    # parameter estimates
    interceptSummary = do.call("rbind", lapply(results, function(x) {x$interceptSummary}))
    interceptSummary = colMeans(interceptSummary)
    if(urbanEffect) {
      urbanSummary = do.call("rbind", lapply(results, function(x) {x$urbanSummary}))
      urbanSummary = colMeans(urbanSummary)
    } else
      urbanSummary = NULL
    rangeSummary = do.call("rbind", lapply(results, function(x) {x$rangeSummary}))
    rangeSummary = colMeans(rangeSummary)
    sdSummary = do.call("rbind", lapply(results, function(x) {x$sdSummary}))
    sdSummary = colMeans(sdSummary)
    varSummary = do.call("rbind", lapply(results, function(x) {x$varSummary}))
    varSummary = colMeans(varSummary)
    if(includeClustEffect) {
      nuggetVarSummary = do.call("rbind", lapply(results, function(x) {x$nuggetVarSummary}))
      nuggetVarSummary = colMeans(nuggetVarSummary)
      nuggetSDSummary = do.call("rbind", lapply(results, function(x) {x$nuggetSDSummary}))
      nuggetSDSummary = colMeans(nuggetSDSummary)
    } else {
      nuggetVarSummary = NULL
      nuggetSDSummary = NULL
    }
    
    list(scoresEaExact=scoresEaExact, scoresEaExactBVar=scoresEaExactBVar, 
         scoresPixelInexact=scoresPixelInexact, scoresPixelExact=scoresPixelExact, scoresPixelExactBVar=scoresPixelExactBVar, 
         scoresCountyInexact=scoresCountyInexact, scoresCountyInexactUnintegrated=scoresCountyInexactUnintegrated, scoresCountyInexactUnadjusted=scoresCountyInexactUnadjusted, scoresCountyExact=scoresCountyExact, scoresCountyExactBVar=scoresCountyExactBVar, 
         scoresRegionInexact=scoresRegionInexact, scoresRegionExact=scoresRegionExact, scoresRegionExactBVar=scoresRegionExactBVar, 
         interceptSummary=interceptSummary, urbanSummary=urbanSummary, 
         rangeSummary=rangeSummary, varSummary=varSummary, sdSummary=sdSummary, 
         nuggetVarSummary=nuggetVarSummary, nuggetSDSummary=nuggetSDSummary)
  }
  
  # now make a function for generating the results that can work either in parallel or serial
  mainFunction = function(i, doSink=FALSE) {
    if(doSink)
      sink("log.txt", append=TRUE)
    
    set.seed(randomSeeds[i])
    
    print(paste0("iteration ", i, "/", nsim))
    
    # get the simulated sample
    clustDat = clustDatMulti[[i]]
    
    # get observations from dataset
    obsCoords = cbind(clustDat$east, clustDat$north)
    obsNs = clustDat$numChildren
    obsCounts = clustDat$died
    obsUrban = clustDat$urban
    
    # fit model, get all predictions for each areal level and each posterior sample
    fit = fitSPDEModel3(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                        predUrban, clusterIndices=clustDat$eaIs, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=nPostSamples, 
                        verbose = verbose, clusterEffect=includeClustEffect, 
                        int.strategy=int.strategy, genRegionLevel=genRegionLevel, 
                        keepPixelPreds=keepPixelPreds, genEALevel=genEALevel, 
                        urbanEffect=urbanEffect, link=1, predictionType=predictionType, 
                        eaDat=eaDat, nSamplePixel=nSamplePixel, significance=significance, 
                        continuousOnly=continuousOnly, strictPrior=strictPrior, 
                        integrateOutCluster=integrateOutCluster, popGridAdjusted=popGridAdjusted, 
                        adjustPopSurface=adjustPopSurface, strictSpatialPrior=strictSpatialPrior)
    print(paste0("Fit completed: iteration ", i, "/", nsim))
    countyPreds = fit$countyPreds
    regionPreds = fit$regionPreds
    pixelPreds = fit$pixelPreds
    eaPreds = fit$eaPreds
    
    if(continuousOnly) {
      warning("Since continuousOnly is set to TRUE, setting genEALevel to FALSE")
      genEALevel = FALSE
    }
    
    ## calculate model fit properties for each level of predictions
    
    # enumeration area level
    if(genEALevel) {
      print(paste0("EA results generating: iteration ", i, "/", nsim))
      
      ### generate enumeration area level predictions
      # calculate estimates, but convert to the logit scale (estimates are the same with or without binomial variation)
      thisu1mEaExact = logit(rowMeans(eaPreds$eaPredMat))
      
      # scoring rules for all aggregation models
      cat(".")
      scoresEaExact = getScoresSPDE(truthByEa$truth, truthByEa$n, thisu1mEaExact, 
                                       expit(thisu1mEaExact), bVar=FALSE, probMat=eaPreds$eaPredMat)
      cat(".")
      scoresEaExactBVar = getScoresSPDE(truthByEa$truth, truthByEa$n, thisu1mEaExact, 
                                           expit(thisu1mEaExact), bVar=TRUE, probMat=eaPreds$eaPredMat)
      cat(".")
    }
    else {
      scoresEaExact = NA
      scoresEaExactBVar = NA
    }
    
    # pixel level
    if(keepPixelPreds) {
      print(paste0("Pixel results generating: iteration ", i, "/", nsim))
      # TODO: fix NA problem with CRPS
      # calculate estimates, but convert to the logit scale (estimates are the same with or without binomial variation)
      thisu1mPixelInexact = logit(rowMeans(pixelPreds$pixelPredMatInexact))
      if(!continuousOnly)
        thisu1mPixelExact = logit(rowMeans(pixelPreds$pixelPredMatExact))
      
      # scoring rules for all aggregation models
      scoresPixelInexact = getScoresSPDE(truthByPixel$truth, truthByPixel$n, thisu1mPixelInexact, 
                                         expit(thisu1mPixelInexact), bVar=FALSE, probMat=pixelPreds$pixelPredMatInexact)
      cat(".")
      if(!continuousOnly) {
        scoresPixelExact = getScoresSPDE(truthByPixel$truth, truthByPixel$n, thisu1mPixelExact, 
                                         expit(thisu1mPixelExact), bVar=FALSE, probMat=pixelPreds$pixelPredMatExact)
        cat(".")
        scoresPixelExactBVar = getScoresSPDE(truthByPixel$truth, truthByPixel$n, thisu1mPixelExact, 
                                             expit(thisu1mPixelExact), bVar=TRUE, pmfs=pixelPreds$pixelPredMatExactB)
        cat(".")
      }
      else {
        scoresPixelExact = NA
        scoresPixelExactBVar = NA
      }
    }
    else {
      scoresPixelInexact = NA
      scoresPixelExact = NA
      scoresPixelExactBVar = NA
    }
    
    # County level (county level results are always included)
    print(paste0("County results generating: iteration ", i, "/", nsim))
    
    # first generate the county estimates
    thisu1mCountyInexact = logit(rowMeans(countyPreds$countyPredMatInexact))
    thisu1mCountyInexactUnintegrated = logit(rowMeans(countyPreds$countyPredMatInexactUnintegrated))
    thisu1mCountyInexactUnadjusted = logit(rowMeans(countyPreds$countyPredMatInexactUnadjusted))
    if(!continuousOnly)
      thisu1mCountyExact = logit(rowMeans(countyPreds$countyPredMatExact))
    
    # scoring rules for all aggregation models
    scoresCountyInexact = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyInexact, 
                                        expit(thisu1mCountyInexact), NULL, bVar=FALSE, probMat=countyPreds$countyPredMatInexact)
    scoresCountyInexactUnintegrated = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyInexactUnintegrated, 
                                        expit(thisu1mCountyInexactUnintegrated), NULL, bVar=FALSE, probMat=countyPreds$countyPredMatInexactUnintegrated)
    scoresCountyInexactUnadjusted = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyInexactUnadjusted, 
                                                    expit(thisu1mCountyInexactUnadjusted), NULL, bVar=FALSE, probMat=countyPreds$countyPredMatInexactUnadjusted)
    cat(".")
    if(!continuousOnly) {
      scoresCountyExact = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyExact, 
                                        expit(thisu1mCountyExact), NULL, bVar=FALSE, probMat=countyPreds$countyPredMatExact)
      cat(".")
      scoresCountyExactBVar = getScoresSPDE(truthByCounty$truth, truthByCounty$n, thisu1mCountyExact, 
                                            expit(thisu1mCountyExact), NULL, bVar=TRUE, pmfs=countyPreds$countyPredMatExactB)
      cat(".")
    }
    else {
      scoresCountyExact = NA
      scoresCountyExactBVar = NA
    }
    
    # region level
    if(genRegionLevel) {
      print(paste0("Region results generating: iteration ", i, "/", nsim))
      
      # first generate the region estimates
      thisu1mRegionInexact = logit(rowMeans(regionPreds$regionPredMatInexact))
      if(!continuousOnly)
        thisu1mRegionExact = logit(rowMeans(regionPreds$regionPredMatExact))
      
      # scoring rules for all aggregation models
      scoresRegionInexact = getScoresSPDE(truthByRegion$truth, truthByRegion$n, thisu1mRegionInexact, 
                                          expit(thisu1mRegionInexact), NULL, bVar=FALSE, probMat=regionPreds$regionPredMatInexact)
      cat(".")
      if(!continuousOnly) {
        scoresRegionExact = getScoresSPDE(truthByRegion$truth, truthByRegion$n, thisu1mRegionExact, 
                                          expit(thisu1mRegionExact), NULL, bVar=FALSE, probMat=regionPreds$regionPredMatExact)
        cat(".")
        scoresRegionExactBVar = getScoresSPDE(truthByRegion$truth, truthByRegion$n, thisu1mRegionExact, 
                                              expit(thisu1mRegionExact), NULL, bVar=TRUE, pmfs=regionPreds$regionPredMatExactB)
        cat(".")
      }
      else {
        scoresRegionExact = NA
        scoresRegionExactBVar = NA
      }
    }
    else {
      scoresRegionInexact = NA
      scoresRegionExact = NA
      scoresRegionExactBVar = NA
    }
    
    ## collect parameter means, sds, and quantiles
    # for fixed effects
    print(paste0("Parameter summaries generating: iteration ", i, "/", nsim))
    interceptQuants = inla.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.fixed[[1]])
    interceptQuants = c(quant10=interceptQuants[1], quant50=interceptQuants[2], quant90=interceptQuants[3])
    interceptMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.fixed[[1]])
    interceptSummary = c(est=interceptMoments[1], sd=sqrt(interceptMoments[2] - interceptMoments[1]^2), var=interceptMoments[2] - interceptMoments[1]^2, interceptQuants, width=interceptQuants[3] - interceptQuants[1])
    
    urbanSummary = NULL
    if(urbanEffect) {
      urbanQuants = inla.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.fixed[[2]])
      urbanQuants = c(quant10=urbanQuants[1], quant50=urbanQuants[2], quant90=urbanQuants[3])
      urbanMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.fixed[[2]])
      urbanSummary = c(est=urbanMoments[1], sd=sqrt(urbanMoments[2] - urbanMoments[1]^2), var=urbanMoments[2] - urbanMoments[1]^2, urbanQuants, width=urbanQuants[3] - urbanQuants[1])
    }
    
    # for hyperparameters
    print(paste0("Hyperparameter summaries generating: iteration ", i, "/", nsim))
    rangeQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[1]])
    rangeQuants = c(quant10=rangeQuants[1], quant50=rangeQuants[2], quant90=rangeQuants[3])
    sdQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[2]])
    sdQuants = c(quant10=sdQuants[1], quant50=sdQuants[2], quant90=sdQuants[3])
    varQuants = sdQuants^2
    varMarg = my.tmarginal(function(x) {x^2}, fit$mod$marginals.hyperpar[[2]])
    rangeMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.hyperpar[[1]])
    rangeSummary = c(est=rangeMoments[1], sd=sqrt(rangeMoments[2] - rangeMoments[1]^2), var=rangeMoments[2] - rangeMoments[1]^2, rangeQuants, width=rangeQuants[3] - rangeQuants[1])
    sdMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.hyperpar[[2]])
    sdSummary = c(est=sdMoments[1], sd=sqrt(sdMoments[2] - sdMoments[1]^2), var=sdMoments[2] - sdMoments[1]^2, sdQuants, width=sdQuants[3] - sdQuants[1])
    varMoments = inla.emarginal(function(x) {c(x, x^2)}, varMarg)
    varSummary = c(est=varMoments[1], sd=sqrt(varMoments[2] - varMoments[1]^2), var=varMoments[2] - varMoments[1]^2, varQuants, width=varQuants[3] - varQuants[1])
    nuggetVarSummary = NULL
    nuggetSDSummary = NULL
    if(includeClustEffect) {
      print(paste0("Cluster hyperparameters summaries generating: iteration ", i, "/", nsim))
      nuggetPrecQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[3]])
      print(paste0("test1: iteration ", i, "/", nsim))
      nuggetVarQuants = 1/nuggetPrecQuants
      print(paste0("test2: iteration ", i, "/", nsim))
      nuggetVarQuants = c(quant10=nuggetVarQuants[3], quant50=nuggetVarQuants[2], quant90=nuggetVarQuants[1])
      print(paste0("test3: iteration ", i, "/", nsim))
      nuggetSDQuants = sqrt(nuggetVarQuants)
      print(paste0("test4: iteration ", i, "/", nsim))
      nuggetVarMarg = my.tmarginal(function(x) {1/x}, fit$mod$marginals.hyperpar[[3]])
      print(paste0("test5: iteration ", i, "/", nsim))
      nuggetSDMarg = my.tmarginal(function(x) {1/sqrt(x)}, fit$mod$marginals.hyperpar[[3]])
      print(paste0("test6: iteration ", i, "/", nsim))
      nuggetVarMoments = inla.emarginal(function(x) {c(x, x^2)}, nuggetVarMarg)
      print(paste0("test7: iteration ", i, "/", nsim))
      nuggetSDMoments = inla.emarginal(function(x) {c(x, x^2)}, nuggetSDMarg)
      print(paste0("test8: iteration ", i, "/", nsim))
      nuggetVarSummary = c(est=nuggetVarMoments[1], sd=sqrt(nuggetVarMoments[2] - nuggetVarMoments[1]^2), var=nuggetVarMoments[2] - nuggetVarMoments[1]^2, nuggetVarQuants, width=nuggetVarQuants[3] - nuggetVarQuants[1])
      print(paste0("test9: iteration ", i, "/", nsim))
      nuggetSDSummary = c(est=nuggetSDMoments[1], sd=sqrt(nuggetSDMoments[2] - nuggetSDMoments[1]^2), var=nuggetSDMoments[2] - nuggetSDMoments[1]^2, nuggetSDQuants, width=nuggetSDQuants[3] - nuggetSDQuants[1])
      print(paste0("test10: iteration ", i, "/", nsim))
    }
    
    print(paste0("Combining results: iteration ", i, "/", nsim))
    res = list(scoresEaExact=scoresEaExact, scoresEaExactBVar=scoresEaExactBVar, 
               scoresPixelInexact=scoresPixelInexact, scoresPixelExact=scoresPixelExact, scoresPixelExactBVar=scoresPixelExactBVar, 
               scoresCountyInexact=scoresCountyInexact, scoresCountyInexactUnintegrated=scoresCountyInexactUnintegrated, scoresCountyInexactUnadjusted=scoresCountyInexactUnadjusted, scoresCountyExact=scoresCountyExact, scoresCountyExactBVar=scoresCountyExactBVar, 
               scoresRegionInexact=scoresRegionInexact, scoresRegionExact=scoresRegionExact, scoresRegionExactBVar=scoresRegionExactBVar, 
               interceptSummary=interceptSummary, urbanSummary=urbanSummary, 
               rangeSummary=rangeSummary, varSummary=varSummary, sdSummary=sdSummary, 
               nuggetVarSummary=nuggetVarSummary, nuggetSDSummary=nuggetSDSummary)
    
    if(doSink)
      sink()
    print(paste0("Completed iteration ", i, "/", nsim))
    res
  }
  
  # compute Bias & MSE & mean(Var) & 80\% coverage for each simulation
  if(is.null(parClust)) {
    # # sequential version
    # scoresEaExact=as.list(1:nsim)
    # scoresEaExactBVar=as.list(1:nsim)
    # scoresPixelInexact=as.list(1:nsim)
    # scoresPixelExact=as.list(1:nsim)
    # scoresPixelExactBVar=as.list(1:nsim)
    # scoresCountyInexact=as.list(1:nsim)
    # scoresCountyExact=as.list(1:nsim)
    # scoresCountyExactBVar=as.list(1:nsim)
    # scoresRegionInexact=as.list(1:nsim)
    # scoresRegionExact=as.list(1:nsim)
    # scoresRegionExactBVar=as.list(1:nsim)
    # 
    # interceptSummary=as.list(1:nsim)
    # urbanSummary=as.list(1:nsim)
    # rangeSummary=as.list(1:nsim)
    # varSummary=as.list(1:nsim)
    # sdSummary=as.list(1:nsim)
    # nuggetVarSummary=as.list(1:nsim)
    # nuggetSDSummary=as.list(1:nsim)
    # for(i in 1:nsim) {
    #   results = mainFunction(i, FALSE)
    # 
    #   # collect results in a list, one element for each simulation
    #   scoresEaExact[[i]]=results$scoresEaExact
    #   scoresEaExactBVar[[i]]=results$scoresEaExactBVar
    #   scoresPixelInexact[[i]]=results$scoresPixelInexact
    #   scoresPixelExact[[i]]=results$scoresPixelExact
    #   scoresPixelExactBVar[[i]]=results$scoresPixelExactBVar
    #   scoresCountyInexact[[i]]=results$scoresCountyInexact
    #   scoresCountyExact[[i]]=results$scoresCountyExact
    #   scoresCountyExactBVar[[i]]=results$scoresCountyExactBVar
    #   scoresRegionInexact[[i]]=results$scoresRegionInexact
    #   scoresRegionExact[[i]]=results$scoresRegionExact
    #   scoresRegionExactBVar[[i]]=results$scoresRegionExactBVar
    # 
    #   interceptSummary[[i]]=results$interceptSummary
    #   urbanSummary[[i]]=results$urbanSummary
    #   rangeSummary[[i]]=results$rangeSummary
    #   varSummary[[i]]=results$varSummary
    #   sdSummary[[i]]=results$sdSummary
    #   nuggetVarSummary[[i]]=results$nuggetVarSummary
    #   nuggetSDSummary[[i]]=results$nuggetSDSummary
    # }
    # 
    # print("Combining results...")
    # scoresEaExact = do.call("rbind", scoresEaExact)
    # scoresEaExactBVar = do.call("rbind", scoresEaExactBVar)
    # scoresPixelInexact = do.call("rbind", scoresPixelInexact)
    # scoresPixelExact = do.call("rbind", scoresPixelExact)
    # scoresPixelExactBVar = do.call("rbind", scoresPixelExactBVar)
    # scoresCountyInexact = do.call("rbind", scoresCountyInexact)
    # scoresCountyExact = do.call("rbind", scoresCountyExact)
    # scoresCountyExactBVar = do.call("rbind", scoresCountyExactBVar)
    # scoresRegionInexact = do.call("rbind", scoresRegionInexact)
    # scoresRegionExact = do.call("rbind", scoresRegionExact)
    # scoresRegionExactBVar = do.call("rbind", scoresRegionExactBVar)
    # 
    # interceptSummary = do.call("rbind", interceptSummary)
    # interceptSummary = colMeans(interceptSummary)
    # if(urbanEffect) {
    #   urbanSummary = do.call("rbind", urbanSummary)
    #   urbanSummary = colMeans(urbanSummary)
    # } else {
    #   urbanSummary = NULL
    # }
    # rangeSummary = do.call("rbind", rangeSummary)
    # rangeSummary = colMeans(rangeSummary)
    # sdSummary = do.call("rbind", sdSummary)
    # sdSummary = colMeans(sdSummary)
    # varSummary = do.call("rbind", varSummary)
    # varSummary = colMeans(varSummary)
    # if(includeClustEffect) {
    #   nuggetVarSummary = do.call("rbind", nuggetVarSummary)
    #   nuggetVarSummary = colMeans(nuggetVarSummary)
    #   nuggetSDSummary = do.call("rbind", nuggetSDSummary)
    #   nuggetSDSummary = colMeans(nuggetSDSummary)
    # } else {
    #   nuggetVarSummary = NULL
    #   nuggetSDSummary = NULL
    # }
    
    # results = foreach(i = 1:nsim, .combine=combineResults, .verbose=TRUE, .multicombine=TRUE) %do% {
    #   mainFunction(i, FALSE)
    # }
    
    surveyResults = list()
    for(i in 1:nsim) {
      surveyResults = c(surveyResults, list(mainFunction(i, FALSE)))
    }
    results = combineResults(surveyResults)
  } else {
    # parallel version
    
    results = foreach(i = 1:nsim, .combine=combineResults, .verbose=TRUE, .multicombine=TRUE, .export=ls()) %dopar% {
      mainFunction(i, FALSE)
    }
    
    # # separate results into the different aggregation levels
    # scoresEaExact=results$scoresEaExact
    # scoresEaExactBVar=results$scoresEaExactBVar
    # scoresPixelInexact=results$scoresPixelInexact
    # scoresPixelExact=results$scoresPixelExact
    # scoresPixelExactBVar=results$scoresPixelExactBVar
    # scoresCountyInexact=results$scoresCountyInexact
    # scoresCountyExact=results$scoresCountyExact
    # scoresCountyExactBVar=results$scoresCountyExactBVar
    # scoresRegionInexact=results$scoresRegionInexact
    # scoresRegionExact=results$scoresRegionExact
    # scoresRegionExactBVar=results$scoresRegionExactBVar
    # 
    # interceptSummary=results$interceptSummary
    # if(urbanEffect) {
    #   urbanSummary=results$urbanSummary
    # } else {
    #   urbanSummary = NULL
    # }
    # rangeSummary=results$rangeSummary
    # varSummary=results$varSummary
    # sdSummary=results$sdSummary
    # if(includeClustEffect) {
    #   nuggetVarSummary=results$nuggetVarSummary
    #   nuggetSDSummary=results$nuggetSDSummary
    # } else {
    #   nuggetVarSummary=NULL
    #   nuggetSDSummary=NULL
    # }
  }
  
  # separate results into the different aggregation levels
  scoresEaExact=results$scoresEaExact
  scoresEaExactBVar=results$scoresEaExactBVar
  scoresPixelInexact=results$scoresPixelInexact
  scoresPixelExact=results$scoresPixelExact
  scoresPixelExactBVar=results$scoresPixelExactBVar
  scoresCountyInexact=results$scoresCountyInexact
  scoresCountyInexactUnintegrated=results$scoresCountyInexactUnintegrated
  scoresCountyInexactUnadjusted=results$scoresCountyInexactUnadjusted
  scoresCountyExact=results$scoresCountyExact
  scoresCountyExactBVar=results$scoresCountyExactBVar
  scoresRegionInexact=results$scoresRegionInexact
  scoresRegionExact=results$scoresRegionExact
  scoresRegionExactBVar=results$scoresRegionExactBVar
  
  interceptSummary=results$interceptSummary
  if(urbanEffect) {
    urbanSummary=results$urbanSummary
  } else {
    urbanSummary = NULL
  }
  rangeSummary=results$rangeSummary
  varSummary=results$varSummary
  sdSummary=results$sdSummary
  if(includeClustEffect) {
    nuggetVarSummary=results$nuggetVarSummary
    nuggetSDSummary=results$nuggetSDSummary
  } else {
    nuggetVarSummary=NULL
    nuggetSDSummary=NULL
  }
  
  list(scoresEaExact=scoresEaExact, scoresEaExactBVar=scoresEaExactBVar, 
       scoresPixelInexact=scoresPixelInexact, scoresPixelExact=scoresPixelExact, scoresPixelExactBVar=scoresPixelExactBVar, 
       scoresCountyInexact=scoresCountyInexact, scoresCountyInexactUnintegrated=scoresCountyInexactUnintegrated, scoresCountyInexactUnadjusted=scoresCountyInexactUnadjusted, scoresCountyExact=scoresCountyExact, scoresCountyExactBVar=scoresCountyExactBVar, 
       scoresRegionInexact=scoresRegionInexact, scoresRegionExact=scoresRegionExact, scoresRegionExactBVar=scoresRegionExactBVar, 
       interceptSummary=interceptSummary, urbanSummary=urbanSummary, 
       rangeSummary=rangeSummary, varSummary=varSummary, sdSummary=sdSummary, 
       nuggetVarSummary=nuggetVarSummary, nuggetSDSummary=nuggetSDSummary)
}

# Based on resultsSPDEHelper, but use for analyzing a single rate data set
resultsSPDEDat = function(clustDat=ed, nPostSamples=1000, verbose=FALSE, 
                          includeClustEffect=FALSE, int.strategy="eb", 
                          genRegionLevel=TRUE, keepPixelPreds=TRUE, 
                          urbanEffect=TRUE, kmres=5, nSamplePixel=nPostSamples, 
                          predictionType=c("mean", "median"), parClust=cl, 
                          significance=.8, previousResult=NULL, predCountyI=NULL, 
                          summarizeParameters=TRUE, targetPop=c("children", "women")) {
  
  # match the requested prediction type with one of the possible options
  predictionType = match.arg(predictionType)
  
  predCoords = cbind(popGrid$east, popGrid$north)
  predUrban = popGrid$urban
  
  # Must predict at clusters as well. Include clusters as 
  # first rows of prediction coordinates and prediction urban/rural
  predCoords = rbind(cbind(clustDat$east, clustDat$north), predCoords)
  predUrban = c(clustDat$urban, predUrban)
  
  # we only care about the probability, not counts, so not used except for the purposes 
  # of calling inla:
  # predNs = rep(25, nrow(predCoords))
  predNs = rep(1, nrow(predCoords))
  
  # get observations from dataset
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$n
  obsCounts = clustDat$y
  obsUrban = clustDat$urban
  
  # fit model, get all predictions for each areal level and each posterior sample 
  # (even though adjustPopSurface is FALSE by default, the predictions will be using the adjusted surface, 
  # passed in as the popGrid variable)
  fit = fitSPDEModel3(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                      predUrban, clusterIndices=1:nrow(clustDat), genCountyLevel=TRUE, nPostSamples=nPostSamples, 
                      verbose = verbose, clusterEffect=includeClustEffect, 
                      int.strategy=int.strategy, genRegionLevel=TRUE, counties=sort(unique(kenyaEAs$admin1)), 
                      keepPixelPreds=keepPixelPreds, genEALevel=TRUE, regions=sort(unique(kenyaEAs$region)), 
                      urbanEffect=urbanEffect, eaIndices=1:nrow(clustDat), 
                      eaDat=eaDat, nSamplePixel=nSamplePixel, 
                      significance=significance, onlyInexact=TRUE, allPixels=TRUE, 
                      newMesh=TRUE, previousResult=previousResult, predCountyI=predCountyI, integrateOutCluster=TRUE, 
                      targetPop=targetPop)
  
  countyPreds = fit$countyPreds
  regionPreds = fit$regionPreds
  pixelPreds = fit$pixelPreds
  clusterPreds = fit$eaPreds
  
  ## calculate model fit properties for each level of predictions
  generatePredictions = function(probMat) {
    if(is.null(probMat))
      return(NULL)
    
    pred = rowMeans(probMat)
    sds = apply(logit(probMat), 1, sd) # remember that the standard deviations are calculated on the logit scale
    lower = apply(probMat, 1, function(x) {quantile(probs=0.1, x)})
    upper = apply(probMat, 1, function(x) {quantile(probs=0.9, x)})
    
    results = list(pred=pred, sds=sds, lower=lower, upper=upper)
    results
  }
  
  # enumeration area level
  resultsCluster = generatePredictions(clusterPreds$eaPredMat)
  
  # pixel level
  resultsPixel = generatePredictions(pixelPreds$pixelPredMatInexact)
  
  # County level
  resultsCounty = generatePredictions(countyPreds$countyPredMatInexact)
  
  # region level
  resultsRegion = generatePredictions(regionPreds$regionPredMatInexact)
  
  if(summarizeParameters) {
    ## collect parameter means, sds, and quantiles
    # for fixed effects
    interceptQuants = inla.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.fixed[[1]])
    interceptQuants = c(quant10=interceptQuants[1], quant50=interceptQuants[2], quant90=interceptQuants[3])
    interceptMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.fixed[[1]])
    interceptSummary = c(est=interceptMoments[1], sd=sqrt(interceptMoments[2] - interceptMoments[1]^2), var=interceptMoments[2] - interceptMoments[1]^2, interceptQuants, width=interceptQuants[3] - interceptQuants[1])
    
    urbanSummary = NULL
    if(urbanEffect) {
      urbanQuants = inla.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.fixed[[2]])
      urbanQuants = c(quant10=urbanQuants[1], quant50=urbanQuants[2], quant90=urbanQuants[3])
      urbanMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.fixed[[2]])
      urbanSummary = c(est=urbanMoments[1], sd=sqrt(urbanMoments[2] - urbanMoments[1]^2), var=urbanMoments[2] - urbanMoments[1]^2, urbanQuants, width=urbanQuants[3] - urbanQuants[1])
    }
    
    # for hyperparameters
    rangeQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[1]])
    rangeQuants = c(quant10=rangeQuants[1], quant50=rangeQuants[2], quant90=rangeQuants[3])
    sdQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[2]])
    sdQuants = c(quant10=sdQuants[1], quant50=sdQuants[2], quant90=sdQuants[3])
    varQuants = sdQuants^2
    varMarg = my.tmarginal(function(x) {x^2}, fit$mod$marginals.hyperpar[[2]])
    rangeMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.hyperpar[[1]])
    rangeSummary = c(est=rangeMoments[1], sd=sqrt(rangeMoments[2] - rangeMoments[1]^2), var=rangeMoments[2] - rangeMoments[1]^2, rangeQuants, width=rangeQuants[3] - rangeQuants[1])
    sdMoments = inla.emarginal(function(x) {c(x, x^2)}, fit$mod$marginals.hyperpar[[2]])
    sdSummary = c(est=sdMoments[1], sd=sqrt(sdMoments[2] - sdMoments[1]^2), var=sdMoments[2] - sdMoments[1]^2, sdQuants, width=sdQuants[3] - sdQuants[1])
    varMoments = inla.emarginal(function(x) {c(x, x^2)}, varMarg)
    varSummary = c(est=varMoments[1], sd=sqrt(varMoments[2] - varMoments[1]^2), var=varMoments[2] - varMoments[1]^2, varQuants, width=varQuants[3] - varQuants[1])
    nuggetVarSummary = NULL
    nuggetSDSummary = NULL
    if(includeClustEffect) {
      nuggetPrecQuants = my.qmarginal(c(0.1, 0.5, 0.9), fit$mod$marginals.hyperpar[[3]])
      nuggetVarQuants = 1/nuggetPrecQuants
      nuggetVarQuants = c(quant10=nuggetVarQuants[3], quant50=nuggetVarQuants[2], quant90=nuggetVarQuants[1])
      nuggetSDQuants = sqrt(nuggetVarQuants)
      nuggetVarMarg = my.tmarginal(function(x) {1/x}, fit$mod$marginals.hyperpar[[3]])
      nuggetSDMarg = my.tmarginal(function(x) {1/sqrt(x)}, fit$mod$marginals.hyperpar[[3]])
      nuggetVarMoments = inla.emarginal(function(x) {c(x, x^2)}, nuggetVarMarg)
      nuggetSDMoments = inla.emarginal(function(x) {c(x, x^2)}, nuggetSDMarg)
      
      nuggetVarSummary = c(est=nuggetVarMoments[1], sd=sqrt(nuggetVarMoments[2] - nuggetVarMoments[1]^2), var=nuggetVarMoments[2] - nuggetVarMoments[1]^2, nuggetVarQuants, width=nuggetVarQuants[3] - nuggetVarQuants[1])
      nuggetSDSummary = c(est=nuggetSDMoments[1], sd=sqrt(nuggetSDMoments[2] - nuggetSDMoments[1]^2), var=nuggetSDMoments[2] - nuggetSDMoments[1]^2, nuggetSDQuants, width=nuggetSDQuants[3] - nuggetSDQuants[1])
    }
    
    list(resultsCluster=resultsCluster, resultsPixel=resultsPixel, resultsCounty=resultsCounty, resultsRegion=resultsRegion, 
         interceptSummary=interceptSummary, urbanSummary=urbanSummary, 
         rangeSummary=rangeSummary, varSummary=varSummary, sdSummary=sdSummary, 
         nuggetVarSummary=nuggetVarSummary, nuggetSDSummary=nuggetSDSummary, 
         pixelDraws=pixelPreds$pixelPredMatInexact, fit=fit)
  } else {
    list(resultsCluster=resultsCluster, resultsPixel=resultsPixel, resultsCounty=resultsCounty, resultsRegion=resultsRegion, fit=fit)
  }
}

# Based on resultsSPDEHelper, but use for analyzing a single data set given a set of direct estimates
validateSPDEDat = function(directLogitEsts, directLogitVars, directVars, 
                           clustDat=ed, nPostSamples=1000, verbose=FALSE, 
                           includeClustEffect=FALSE, int.strategy="eb", 
                           genRegionLevel=TRUE, keepPixelPreds=TRUE, 
                           urbanEffect=TRUE, kmres=5, nSamplePixel=nPostSamples, 
                           predictionType=c("mean", "median"), parClust=cl, 
                           significance=.8, saveResults=TRUE, fileNameRoot="Ed", 
                           loadPreviousFit=FALSE, targetPop=c("children", "women"), 
                           adjustPopSurface=TRUE) {
  # match the requested prediction type with one of the possible options
  predictionType = match.arg(predictionType)
  
  # if not supplied, get grid of population densities for pop-weighted integration
  if(kmres == 5) {
    if(adjustPopSurface) {
      if(targetPop == "children") {
        load("popGridAdjusted.RData")
        popGridAdjusted = popGrid
      }
      else {
        load("popGridAdjustedWomen.RData")
        popGridAdjusted = popGrid
      }
    }
    load("popGrid.RData")
  }
  else {
    popGrid = makeInterpPopGrid(kmres, FALSE, targetPop)
    if(adjustPopSurface)
      popGridAdjusted = makeInterpPopGrid(kmres, adjustPopSurface, targetPop)
  }
  
  predCoords = cbind(popGrid$east, popGrid$north)
  predUrban = popGrid$urban
  
  # Must predict at clusters as well. Include clusters as 
  # first rows of prediction coordinates and prediction urban/rural
  predCoords = rbind(cbind(clustDat$east, clustDat$north), predCoords)
  predUrban = c(clustDat$urban, predUrban)
  
  # we only care about the probability, not counts, so not used except for the purposes 
  # of calling inla:
  # predNs = rep(25, nrow(predCoords))
  predNs = rep(1, nrow(predCoords))
  
  # get observations from dataset
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$n
  obsCounts = clustDat$y
  obsUrban = clustDat$urban
  
  # first fit the full model (we will use this to initialize the model during the validation fits for each left out county)
  print("Fitting full model")
  fileName = paste0("resultsSPDE", fileNameRoot, "ValidationFull", "_includeClustEffect", includeClustEffect, 
                    "_urbanEffect", urbanEffect, ".RData")
  if(!loadPreviousFit) {
    out = fitSPDEModel3(obsCoords, obsNs=obsNs, obsCounts, obsUrban, predCoords, predNs=predNs, 
                        predUrban, clusterIndices=1:nrow(clustDat), genCountyLevel=FALSE, popGrid=popGrid, nPostSamples=nPostSamples, 
                        verbose = verbose, clusterEffect=includeClustEffect, 
                        int.strategy=int.strategy, genRegionLevel=FALSE, counties=sort(unique(kenyaEAs$admin1)), 
                        keepPixelPreds=keepPixelPreds, genEALevel=TRUE, regions=sort(unique(kenyaEAs$region)), 
                        urbanEffect=urbanEffect, eaIndices=1:nrow(clustDat), 
                        eaDat=eaDat, nSamplePixel=nSamplePixel, 
                        significance=significance, onlyInexact=TRUE, allPixels=TRUE, 
                        newMesh=TRUE, doValidation=TRUE, popGridAdjusted=popGridAdjusted)
    
    if(saveResults)
      save(out, file=fileName)
  }
  else {
    load(fileName)
  }
  cpo = out$mod$cpo$cpo
  cpoFailure = out$mod$cpo$failure
  dic = out$mod$dic$dic
  waic = out$mod$waic$waic
  modelFit = out$mod
  clusterPredsInSample = out$eaPreds$eaPredMat
  
  counties = sort(unique(poppc$County))
  allClusterI = c()
  for(i in 1:length(counties)) {
    if(i %% 10 == 1)
      print(paste0("Fitting model with data from county ", i, "/", length(counties), " left out"))
    thisCountyName = counties[i]
    
    thisCounty = as.character(clustDat$admin1) == thisCountyName
    thisObsCounts = obsCounts
    thisObsCounts[thisCounty] = NA
    
    # fit model, get all predictions for each areal level and each posterior sample
    fit = fitSPDEModel3(obsCoords, obsNs=obsNs, thisObsCounts, obsUrban, predCoords, predNs=predNs, 
                        predUrban, clusterIndices=1:nrow(clustDat), genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=nPostSamples, 
                        verbose = verbose, clusterEffect=includeClustEffect, 
                        int.strategy=int.strategy, genRegionLevel=FALSE, counties=sort(unique(kenyaEAs$admin1)), 
                        keepPixelPreds=keepPixelPreds, genEALevel=TRUE, regions=sort(unique(kenyaEAs$region)), 
                        urbanEffect=urbanEffect, eaIndices=1:nrow(clustDat), 
                        eaDat=eaDat, nSamplePixel=nSamplePixel, 
                        significance=significance, onlyInexact=TRUE, allPixels=TRUE, 
                        newMesh=TRUE, previousResult=modelFit, predCountyI=i, popGridAdjusted=popGridAdjusted)
    
    # thisCountyPreds = fit$countyPreds$countyPredMatInexact
    thisClusterPreds = fit$eaPreds$eaPredMat
    thisClusterI = which(clustDat$admin1 == counties[i])
    allClusterI = c(allClusterI, thisClusterI)
    
    if(i == 1) {
      clusterPreds = thisClusterPreds[thisClusterI,]
    } else {
      clusterPreds = rbind(clusterPreds, thisClusterPreds[thisClusterI,])
    }
  }
  
  # resort cluster predictions to correspond to the ordering of clustDat
  sortI = sort(allClusterI, index.return=TRUE)$ix
  clusterPreds = clusterPreds[sortI,]
  
  # # summary statistics for county level predictions
  # probMat = countyPreds
  # preds = rowMeans(probMat)
  # vars = apply(probMat, 1, var)
  # logitPreds = rowMeans(logit(probMat))
  # logitVars = apply(logit(probMat), 1, var)
  
  # summary statistics for cluster level predictions
  probMat = clusterPreds
  preds = rowMeans(probMat)
  vars = apply(probMat, 1, var)
  logitPreds = rowMeans(logit(probMat))
  logitVars = apply(logit(probMat), 1, var)
  
  probMatInSample = clusterPredsInSample
  predsInSample = rowMeans(probMatInSample)
  varsInSample = apply(probMatInSample, 1, var)
  logitPredsInSample = rowMeans(logit(probMatInSample))
  logitVarsInSample = apply(logit(probMatInSample), 1, var)
  
  # calculate validation scoring rules
  # weights = 1 / directVars
  # weights = weights / sum(weights)
  theseScores = getValidationScores(clustDat$y / clustDat$n, logitPreds, logitVars, 
                                    ests=preds, probMat=probMat, usePearson=FALSE, n=clustDat$n, 
                                    urbanVec=clustDat$urban, filterType="leftOutCounty")
  
  theseScoresInSample = getValidationScores(clustDat$y / clustDat$n, logitPredsInSample, logitVarsInSample, 
                                    ests=predsInSample, probMat=probMatInSample, usePearson=FALSE, n=clustDat$n, 
                                    urbanVec=clustDat$urban, filterType="inSample")
  theseScoresInSample$scores = cbind(WAIC=waic, DIC=dic, theseScoresInSample$scores)
  theseScoresInSample$allResults = cbind(WAIC=waic, DIC=dic, theseScoresInSample$allResults)
  
  # compile scores c("MSE", "CPO", "CRPS", "logScore")
  spdeResultsInSample = theseScoresInSample
  spdeResultsLeaveOutCounty = theseScores
  spdeResultsLeaveOutCluster = data.frame(CPO=mean(cpo, na.rm=TRUE))
  cpoFailure = mean(cpoFailure[!is.na(cpoFailure)] != 0)
  
  spdeResults = list(fullModelFit=modelFit, probMatInSample=probMatInSample, probMat=probMat, 
                     spdeResultsInSample=spdeResultsInSample, spdeResultsLeaveOutCounty=spdeResultsLeaveOutCounty, 
                     spdeResultsLeaveOutCluster=spdeResultsLeaveOutCluster, 
                     cpoFailure=cpoFailure)
  
  # save and return results
  fileName = paste0("resultsSPDE", fileNameRoot, "ValidationAll", "_includeClustEffect", includeClustEffect, 
                    "_urbanEffect", urbanEffect, ".RData")
  fileNameCompact = paste0("resultsSPDE", fileNameRoot, "ValidationAll", "_includeClustEffect", includeClustEffect, 
                           "_urbanEffect", urbanEffect, "compact.RData")
  if(saveResults) {
    save(spdeResults, file=fileName)
    
    # also save a compact version, without the full model fit
    temp = spdeResults
    spdeResults$fullModelFit = NULL
    save(spdeResults, file=fileNameCompact)
    
    spdeResults = temp
  }
  
  spdeResults
}

# compute true proportion of women that completed secondary education in each county in the order of the 
# input counties variable from the full dataset eaDat
getTruthByCounty = function(eaDat, counties=as.character(unique(ed$admin1))) {
  # get true rate for each county
  # diedByCounty = aggregate(eaDat$died, by=list(eaDat$admin1), FUN="sum")
  # theseCounties = diedByCounty$Group.1
  # diedByCounty = diedByCounty$V1
  countByCounty = tapply(eaDat$y, eaDat$admin1, sum)
  theseCounties = rownames(countByCounty)
  # nByCounty = aggregate(eaDat$n, by=list(eaDat$admin1), FUN="sum")$x
  nByCounty = tapply(eaDat$n, eaDat$admin1, sum)
  rate = countByCounty/nByCounty
  
  # sort to be in the order of the input countries variable
  sortI = match(counties, theseCounties)
  list(counties=counties, rate=rate[sortI], n=nByCounty)
}

# same as inla.tmarginal, except can fit splines to the marginals on log scale, which may 
# be necessary for highly peaked distributions near 0
my.tmarginal = function (fun, marginal, n = 2048L, h.diff = .Machine$double.eps^(1/3), 
                         method = c("quantile", "linear"), logX=TRUE) 
{
  ff = match.fun(fun)
  is.mat = is.matrix(marginal)
  m = inla.smarginal(marginal)
  r = range(m$x)
  method = match.arg(method)
  if (INLA:::inla.strcasecmp(method, "quantile")) {
    x = my.qmarginal((1:n)/(n + 1), marginal, logX=logX)
  }
  else if (INLA:::inla.strcasecmp(method, "linear")) {
    x = seq(r[1], r[2], length = n)
  }
  else {
    stop("unknown method")
  }
  xx = ff(x)
  fd = INLA:::inla.deriv.func(ff)
  log.dens = inla.dmarginal(x, marginal, log = FALSE)/abs(fd(x))
  if (xx[1] > xx[n]) {
    xx[1:n] = xx[n:1]
    log.dens[1:n] = log.dens[n:1]
  }
  if (is.mat) {
    ret = cbind(x = xx, y = log.dens)
  }
  else {
    ret = list(x = xx, y = log.dens)
  }
  if (FALSE) {
    class(ret) = "inla.marginal"
    attr(ret, "inla.tag") = paste(attr(marginal, "inla.tag"), 
                                  "transformed")
  }
  return(ret)
}

# same as inla.qmarginal, except can fit splines to the marginals on log scale, which may 
# be necessary for highly peaked distributions near 0
my.qmarginal = function(p, marginal, len = 2048L, logX=TRUE) 
{
  f = INLA:::inla.sfmarginal(inla.smarginal(marginal))
  xx = seq(f$range[1], f$range[2], length = len)
  
  d = cumsum(exp(f$fun(xx)))
  d = d/d[length(d)]
  eps = .Machine$double.eps * 1000
  for (val in c(0, 1)) {
    is.val = which(abs(d - val) <= eps)
    if (length(is.val) > 1) {
      is.val = is.val[-1]
      d = d[-is.val]
      xx = xx[-is.val]
    }
  }
  
  if(logX) {
    fq = splinefun(d, log(xx), method = "hyman")
    n = length(p)
    pp = pmin(pmax(p, rep(0, n)), rep(1, n))
    return(exp(fq(pp)))
  }
  else {
    fq = splinefun(d, xx, method = "hyman")
    n = length(p)
    pp = pmin(pmax(p, rep(0, n)), rep(1, n))
    return(fq(pp))
  }
}



