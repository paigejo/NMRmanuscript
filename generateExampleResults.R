# function for analyzing example datasets

# first name elements of ed to be the same as the corresponding elements of the simulated datasets

generateExampleResults = function(dat=ed, resultNameRoot="Ed") {
  resultNameRootLower = tolower(resultNameRoot)
  
  ##### first generate results for direct and naive models
  startIDat = 1
  for(i in 1:nrow(dat)) {
    if(i %% 100 == 1)
      print(paste0("i: ", i))
    
    tmpDat = extendDataDat(dat[i,], v001 = i)
    
    # initialize the data frames on the first iteration only
    if(i == 1) {
      resDat = as.data.frame(matrix(nrow=sum(dat$n), ncol=ncol(tmpDat) + 1))
    }
    
    # append to the data frames
    names(resDat) = c(names(tmpDat), "regionRural")
    endIDat = startIDat + nrow(tmpDat) - 1
    resDat[startIDat:endIDat, 1:ncol(tmpDat)] = tmpDat
    
    # update row index
    startIDat = endIDat + 1
  }
  
  # add in RegionRural interaction
  resDat$regionRural <- with(resDat, interaction(admin1, urban), drop=TRUE)
  
  # save the resulting data frame
  save(resDat, file=paste0("data4direct", resultNameRoot, ".RData"))
  
  # start analysis using naive approach and computing direct estimates
  directEstDat = naiveDat = list()
  
  # analyse the unstratified sampling scenario
  dat_obj2 = resDat
  res2 = defineSurveyDat(dat_obj2, 
                         stratVar=dat_obj2$regionRural,
                         useSamplingWeights = TRUE)
  directEstResults = res2
  defineSurvey
  resnA2 = run_naiveDat(dat_obj2)
  
  naiveResults = resnA2
  
  # save results
  save(directEstResults, naiveResults, file=paste0("resultsDirectNaive", resultNameRoot, ".RData"))
  
  ##### run Mercer et al. model
  source("mercer.R")
  tmpResults = mercer_u1m2(directEstResults$logit.est, directEstResults$var.est, 
                          graph.path = "Kenyaadm1.graph")
  
  res = data.frame(admin1=directEstResults$admin1,
                   est.mercer=expit(tmpResults$summary.linear.predictor$mean),
                   lower.mercer=tmpResults$summary.linear.predictor$"0.1quant",
                   upper.mercer=tmpResults$summary.linear.predictor$"0.9quant",
                   logit.est.mercer=tmpResults$summary.linear.predictor$mean, 
                   var.est.mercer=(tmpResults$summary.linear.predictor$sd)^2)
  mercerResults = res
  
  ## collect parameter estimates
  fixedParameters = tmpResults$summary.fixed
  
  # BYM2 hyperparameter phi
  hyperparameters = tmpResults$summary.hyperpar
  
  ## transformed hyperparameters
  # sample the hyperparameters, using the marginals to improve the sampling
  out = inla.hyperpar.sample(1000, tmpResults, improve.marginals=TRUE)
  transformFunction = function(x) {c(x[2], 1/x[1], 1/x[1]*x[2], 1/x[1]*(1-x[2]), sqrt(1/x[1]), sqrt(1/x[1]*x[2]), sqrt(1/x[1]*(1-x[2])))}
  transformedOut = apply(out, 1, transformFunction)
  
  # now calculate the summary statistics of the transformed BYM2 hyperparameters
  parNames = c("BYM2 Phi", "BYM2 Tot. Var", "BYM2 Spatial Var", "BYM2 iid Var", "BYM2 Tot. SD", "BYM2 Spatial SD", "BYM2 iid SD")
  rownames(transformedOut) = parNames
  mercerParEst = rowMeans(transformedOut)
  mercerParSD = apply(transformedOut, 1, sd)
  mercerPar10 = apply(transformedOut, 1, quantile, probs=.1)
  mercerPar50 = apply(transformedOut, 1, quantile, probs=.5)
  mercerPar90 = apply(transformedOut, 1, quantile, probs=.9)
  mercerParResults = cbind(mercerParEst, mercerParSD, mercerPar10, mercerPar50, mercerPar90)
  mercerParResults = rbind(Intercept=as.numeric(tmpResults$summary.fixed[1:5]), mercerParResults)
  colnames(mercerParResults) = c("Est", "SD", "Q10", "Q50", "Q90")
  
  save(mercerResults, mercerParResults, file=paste0("resultsMercer", resultNameRoot, ".RData"))
  
  ##### run BYM models
  source("designBased.R")
  runBYM2Dat(dat, includeUrbanRural = FALSE, includeCluster = FALSE, fileNameRoot=resultNameRoot)
  runBYM2Dat(dat, includeUrbanRural = FALSE, includeCluster = TRUE, fileNameRoot=resultNameRoot)
  runBYM2Dat(dat, includeUrbanRural = TRUE, includeCluster = FALSE, fileNameRoot=resultNameRoot)
  runBYM2Dat(dat, includeUrbanRural = TRUE, includeCluster = TRUE, fileNameRoot=resultNameRoot)
  
  ##### run SPDE 
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    spdeResults = do.call("resultsSPDEDat", args)
    includeClustEffect = args$includeClustEffect
    urbanEffect = args$urbanEffect
    fileName = paste0("resultsSPDE", resultNameRootLower, "_includeClustEffect", includeClustEffect, 
                      "_urbanEffect", urbanEffect, ".RData")
    save(spdeResults, file=fileName)
  }
  invisible(NULL)
}