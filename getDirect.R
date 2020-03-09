source("neonatalSimStudyWeighted.R")

# load("simDataMulti.RData")
# load a different 1 of these depending on whether a cluster effect should be included 
# in the simulation of the data or not (tausq is the cluster effect variance)
# load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
# load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")

getDirectNaive = function(tausq=0.1^2, test=FALSE, loadResults=FALSE, big=FALSE, margVar=0.15^2, gamma=-1, effRange=150) {
  bigText = ifelse(big, "Big", "")
  rangeText = ifelse(effRange == 150, "", "Range50")
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0", rangeText, bigText, ".RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0", rangeText, "Test", bigText, ".RData"))
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
  # load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
  
  if(!loadResults) {
    data4directSRS = list()
    data4directoverSamp = list()
    
    # extend the original simulated datasets to binary format
    # where each child gets one row in the dataset.
    # for(j in 1:length(SRSDat[[2]])){
    #   print(paste0("j: ", j))
    #   nChildrenSRS = sapply(SRSDat$clustDat, function(x) {x$numChildren})
    #   nChildrenOverSamp = sapply(overSampDat$clustDat, function(x) {x$numChildren})
    #   startISRS = 1
    #   startIOverSamp = 1
    #   for(i in 1:nrow(SRSDat[[2]][[j]])){
    #     if(i %% 100 == 1)
    #       print(paste0("i: ", i, " / ", nrow(SRSDat[[2]][[j]])))
    #     
    #     tempSRS = extendData(SRSDat[[2]][[j]][i,], v001 = i)
    #     tempOversamp = extendData(overSampDat[[2]][[j]][i,], v001 = i)
    #     
    #     if(i == 1) {
    #       resSRS = data.frame(matrix(nrow = nChildrenSRS[j], ncol=ncol(tempSRS)))
    #       resoverSamp = data.frame(matrix(nrow = nChildrenOverSamp[j], ncol=ncol(tempOversamp)))
    #     }
    #     resSRS[startISRS:(startISRS - 1 + nrow(tempSRS)),] = tempSRS
    #     resoverSamp[startIOverSamp:(startIOverSamp - 1 + nrow(tempOversamp)),] = tempOversamp
    #     
    #     startISRS = startISRS + nrow(tempSRS)
    #     startIOverSamp = startIOverSamp + nrow(tempOversamp)
    #   }
    #   data4directSRS[[j]] = resSRS
    #   data4directoverSamp[[j]] = resoverSamp
    # }
    
    
    for(j in 1:length(SRSDat[[2]])){
      print(paste0("j: ", j))
      startISRS = 1
      startIoverSamp = 1
      for(i in 1:nrow(SRSDat[[2]][[j]])) {
        if(i %% 100 == 1)
          print(paste0("i: ", i))

        tmpSRS = extendData(SRSDat[[2]][[j]][i,], v001 = i)
        tmpoverSamp = extendData(overSampDat[[2]][[j]][i,], v001 = i)

        # initialize the data frames on the first iteration only
        if(i == 1) {
          resSRS = as.data.frame(matrix(nrow=sum(SRSDat[[2]][[j]]$numChildren), ncol=ncol(tmpSRS) + 1))
          resoverSamp = as.data.frame(matrix(nrow=sum(overSampDat[[2]][[j]]$numChildren), ncol=ncol(tmpoverSamp) + 1))
        }

        # append to the data frames
        names(resSRS) = c(names(tmpSRS), "regionRural")
        names(resoverSamp) = c(names(tmpoverSamp), "regionRural")
        endISRS = startISRS + nrow(tmpSRS) - 1
        endIoverSamp = startIoverSamp + nrow(tmpoverSamp) - 1
        resSRS[startISRS:endISRS, 1:ncol(tmpSRS)] = tmpSRS
        resoverSamp[startIoverSamp:endIoverSamp, 1:ncol(tmpoverSamp)] = tmpoverSamp

        # update row index
        startISRS = endISRS + 1
        startIoverSamp = endIoverSamp + 1
      }
      
      # add in RegionRural interaction
      resSRS$regionRural <- with(resSRS, interaction(admin1, urbanRural), drop=TRUE)
      resoverSamp$regionRural <- with(resoverSamp, interaction(admin1, urbanRural), drop=TRUE)
      
      if(any(is.na(resSRS)) || any(is.na(resoverSamp)))
        stop()
      data4directSRS[[j]] = resSRS
      data4directoverSamp[[j]] = resoverSamp
    }
    
    tc = ifelse(tausq == 0, "", paste0("Tausq", round(tausq, 4)))
    if(!test)
      save(data4directSRS, data4directoverSamp, 
           file=paste0("data4direct", "simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                       "HHoldVar0urbanOverSamplefrac0", bigText, rangeText, ".RData"))
    else
      save(data4directSRS, data4directoverSamp, 
           file=paste0("data4direct", "simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                       "HHoldVar0urbanOverSamplefrac0", bigText, "Test", rangeText, ".RData"))
    # save(data4directSRS, data4directoverSamp, file="data4directTausq0.01.RData")
    # save(data4directSRS, data4directoverSamp, file="data4directTest.RData")
    # save(data4directSRS, data4directoverSamp, file="data4directTausq0.01Test.RData")
  } else {
    tc = ifelse(tausq == 0, "", paste0("Tausq", round(tausq, 4)))
    if(!test)
      load(paste0("data4direct", "simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                  "HHoldVar0urbanOverSamplefrac0", bigText, rangeText, ".RData"))
    else
      load(paste0("data4direct", "simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                  "HHoldVar0urbanOverSamplefrac0", bigText, "Test", rangeText, ".RData"))
  }
  
  # start analysis using naive approach and computing direct estimates
  directEstSRS = naiveSRS = list()
  directEstoverSamp = naiveoverSamp = list()
  
  for(i in 1:length(data4directoverSamp)){
    print(i)
    # analyse the oversampled scenario
    childBirths_obj = data4directoverSamp[[i]]
    res = defineSurvey(childBirths_obj, 
                       stratVar=childBirths_obj$regionRural,
                       useSamplingWeights = TRUE)
    directEstoverSamp[[i]] = res
    
    resnA = run_naive(childBirths_obj)
    
    naiveoverSamp[[i]] = resnA
    
    # analyse the unstratified sampling scenario
    childBirths_obj2 = data4directSRS[[i]]
    res2 = defineSurvey(childBirths_obj2, 
                        stratVar=childBirths_obj2$regionRural,
                        useSamplingWeights = TRUE)
    directEstSRS[[i]] = res2
    
    resnA2 = run_naive(childBirths_obj2)
    
    naiveSRS[[i]] = resnA2
  }
  
  if(!test)
    save(directEstSRS, directEstoverSamp, naiveSRS, naiveoverSamp,
         file=paste0("resultsDirectNaiveBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                     "HHoldVar0urbanOverSamplefrac0", bigText, rangeText, ".RData"))
  else
    save(directEstSRS, directEstoverSamp, naiveSRS, naiveoverSamp,
         file=paste0("resultsDirectNaiveBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                     "HHoldVar0urbanOverSamplefrac0Test", bigText, rangeText, ".RData"))
  # save(directEstSRS, directEstoverSamp, naiveSRS, naiveoverSamp,file="resultsDirectNaiveTausq0.01.RData")
  # save(directEstSRS, directEstoverSamp, naiveSRS, naiveoverSamp,file="resultsDirectNaiveTausq0Test.RData")
  # save(directEstSRS, directEstoverSamp, naiveSRS, naiveoverSamp,file="resultsDirectNaiveTausq0.01Test.RData")
  
  invisible(NULL)
}