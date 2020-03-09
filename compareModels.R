library(xtable)
source("scores.R")

# This model produces an xtable summarizing the scoring rules aggregated at the given level of interest. 
# Same as runCompareModels, except uses the updated scoring rule function getScores, so it produces 
# scoring rules based on the predictive distribution with and without binomial variation.
# test: whether or not to use the test data set based results or the true
#       simulation study results. The test data set results are based on fewer
#       observations, requiring fewer computations and should be used to test the code
# tausq: the spatial nugget in the simulated data, can only be the default or 0
# resultType: the level of aggregation over which scoring rules should be computed
# sampling: whether or not to use the simple random sample or urban over sampled results
# recomputeTruth: by default the should be set to TRUE, unless the user calls this function 
#                 twice in a row with the same set of inputs
# modelsI: which model results should be included. Models come in the order:
#          c("naive", "direct", "mercer", "bym", "bymMod", "bymNoUrb", "bymNoUrbMod", "bymNoClust", "bymNoUrbClust", "spde", "spdeNoUrb"), 
#          so the default value of this argument is to include results for the naive, 
#          direct, and mercer models. 1:11 is all models
# produceFigures: whether or not to produce figures based on the results
# printIEvery: how often to print progress
# maxDataSets: if not null, set this to a small integer to only score this many datasets 
#              for testing purposes
# nsim: the number of points with which to approximate the distribution of probability on the logit scale
# saveResults: whether or not to save a disk image of the results
# loadResults: if TRUE loads the saved disk image rather than recomputing results
# xtable.args: the arguments passed to xtable for printing results
# tableFormat: if "1", binomial scores are considered the same model as non-binomial scores. 
#             i "2", binomial scores are put on extra rows of the printed table
runCompareModels2 = function(test=FALSE, tausq=.1^2, margVar=.15^2, gamma=-1, 
                             beta0=-1.75, effRange=150, resultType=c("county", "pixel", "EA"), 
                             sampling=c("SRS", "oversamp"), recomputeTruth=TRUE, modelsI=1:21, 
                             produceFigures=FALSE, big=FALSE, printIEvery=50, 
                             maxDataSets=NULL, nsim=10, saveResults=TRUE, loadResults=TRUE, 
                             xtable.args=list(digits=c(0, 2, 2, 2, 2, 1, 2), display=rep("f", 7), auto=TRUE), 
                             tableFormat=c("2", "1"), colScale=c(10^4, 10^5, 10^5, 10^3, 100, 100), 
                             colUnits=c(" ($\\times 10^{-4}$)", " ($\\times 10^{-5}$)", " ($\\times 10^{-5}$)", 
                                        " ($\\times 10^{-3}$)", " ($\\times 10^{-2}$)", " ($\\times 10^{-2}$)"), 
                             colDigits=c(2, 2, 2, 2, 1, 2), counties=sort(unique(poppc$admin1)), 
                             loadTempProgress=FALSE, includeBVarResults=FALSE, continuousSPDEonly=TRUE, 
                             strictPriors=FALSE, doFancyTables=FALSE, printScoreTable=TRUE, printParTable=TRUE) {
  
  # match the arguments with their correct values
  resultType = match.arg(resultType)
  sampling = match.arg(sampling)
  tableFormat = match.arg(tableFormat)
  
  # test to make sure only the naive and direct models are included if big is set to TRUE
  if(!identical(modelsI, 1:2) && big)
    stop("if big is set to TRUE, only naive and direct results can be included")
  
  # load the true superpopulation
  testText = ifelse(test, "Test", "")
  bigText = ifelse(big, "Big", "")
  strictPriorText = ifelse(strictPriors, "strictPrior", "")
  rangeText = ifelse(effRange == 150, "", "Range50")
  if(!test)
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0", rangeText, bigText, ".RData"))
  else
    load(paste0("simDataMultiBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                "HHoldVar0urbanOverSamplefrac0Test", rangeText, bigText, ".RData"))
  eaDat = SRSDat$eaDat
  
  if(sampling == "SRS") {
    clustDat = SRSDat
  } else {
    clustDat = overSampDat
  }
  maxDataSets = ifelse(is.null(maxDataSets), length(clustDat$clustDat), maxDataSets)
  
  # allModels = c("naive", "direct", "mercer", "bym", "bymMod", "bymNoUrb", "bymNoUrbMod", "bymNoClust", "bymNoUrbClust", "spde", "spdeNoUrb")
  # allNames = c("Naive", "Direct ", "Smoothed Direct", "BYM (no urban/cluster)", "BYM (no urban)", "BYM (no cluster)", "BYM", "SPDE (no urban)", "SPDE")
  # allNamesBinomial = c("Naive Binom.", "Direct Binom.", "Mercer et al. Binom.", "BYM Binom. (no urb/clust)", "BYM Binom. (no urb)", "BYM Binom. (no clust)", "BYM Binom.", "SPDE Binom. (no urb)", "SPDE Binom.")
  # BYM models are in order of complexity: no urban/cluster, no urban, no cluster, full
  allNames = c("Naive", "Direct", "Smoothed Direct", "BYM2 ucA", "BYM2 uCA", "BYM2 uCA'", "BYM2 UcA", "BYM2 UCA", "BYM2 UCA'", 
               "BYM2 uca", "BYM2 uCa", "BYM2 uCa'", "BYM2 Uca", "BYM2 UCa", "BYM2 UCa'", 
               "SPDE uc", "SPDE uC", "SPDE Uc", "SPDE UC", "SPDE uC'", "SPDE UC'")
  allNamesBinomial = paste0(allNames, " Bin.")
  allModels = allNames
  models = allModels[modelsI]
  
  # this string carries all the information about the run
  runId = paste0("Beta-1.75margVar", round(margVar, 4), rangeText, "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                 "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                 "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
  
  # compute the results if necessary
  if(!loadResults && !loadTempProgress) {
    
    ## generate the true values at the county level for the 
    ## given settings if these are new settings
    counties = sort(unique(eaDat$admin1))
    if(recomputeTruth){
      truth = getTruth(resultType, eaDat)
      if(resultType == "county") {
        save(truth, file="truthbycounty.RData")
      } else if(resultType == "pixel") {
        save(truth, file="truthbyPixel.RData")
      } else if(resultType == "EA") {
        save(truth, file="truthbyEA.RData")
      }
    } else {
      # load the truth
      load(paste0("truthby", resultType, ".RData"))
    }
    
    ## Now pick which models to include in the results
    
    # load data
    tauText = ifelse(tausq == 0, "0", "0.01")
    testText = ifelse(test, "Test", "")
    if("Naive" %in% models || "Direct" %in% models) {
      out = load(paste0("resultsDirectNaiveBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                        "HHoldVar0urbanOverSamplefrac0", testText, bigText, rangeText, ".RData"))
      if(sampling == "SRS") {
        directEst = directEstSRS
        naive = naiveSRS
      }
      else {
        directEst = directEstoverSamp
        naive = naiveoverSamp
      }
    }
    if("Smoothed Direct" %in% models) {
      if(!test)
        load(paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                    "HHoldVar0urbanOverSamplefrac0", strictPriorText, rangeText, ".RData"))
      else
        load(paste0("resultsMercerBeta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                    "HHoldVar0urbanOverSamplefrac0", strictPriorText, rangeText, "Test.RData"))
      if(sampling == "SRS") {
        mercer = mercerSRS
        mercerPar = mercerSRSPar
      }
      else {
        mercer = merceroverSamp
        mercerPar = merceroverSampPar
      }
    }
    if("BYM2 ucA" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = FALSE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      # out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurFALSEClusterTRUE", strictPriorText, testText, ".RData"))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbClust = designRes
    }
    if("BYM2 uCA" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrb = designRes
    }
    if("BYM2 uCa'" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbMod = designRes
    }
    if("BYM2 UcA" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = FALSE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResNoClust = designRes
    }
    if("BYM2 UCA'" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResMod = designRes
    }
    if("BYM2 UCA" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = FALSE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResTemp = designRes
    }
    if("BYM2 uca" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = FALSE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      # out = load(paste0("kenyaSpatialDesignResultNewTausq", tauText, "UrbRurFALSEClusterTRUE", strictPriorText, testText, ".RData"))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbClustPopAgg = designRes
    }
    if("BYM2 uCa" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbPopAgg = designRes
    }
    if("BYM2 uCa'" %in% models) {
      includeUrbanRural = FALSE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$SRSdatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$SRSdatClusterUrban
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel")
          designRes[[1]] = designRes$overSampDatPixelUrban
        else if(resultType == "EA")
          designRes[[1]] = designRes$overSampDatClusterUrban
      }
      designResNoUrbModPopAgg = designRes
    }
    if("BYM2 Uca" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = FALSE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResNoClustPopAgg = designRes
    }
    if("BYM2 UCa'" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "debiasedMaxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResModPopAgg = designRes
    }
    if("BYM2 UCa" %in% models) {
      includeUrbanRural = TRUE
      includeCluster = TRUE
      aggregateByPopulation = TRUE
      load(paste0('bym2Beta-1.75margVar', round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 'UrbRur',
                  includeUrbanRural, 'Cluster', includeCluster, "aggByPop", aggregateByPopulation, "maxDataSets", 100, strictPriorText, testText, rangeText, '.RData'))
      if(sampling == "SRS") {
        designRes$overSampDat = NULL
        designRes$overSampDatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$SRSdatPixelUrban
          designRes[[1]] = designRes$SRSdatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$SRSdatClusterUrban
          designRes[[1]] = designRes$SRSdatClusterRural
        }
      }
      else {
        designRes$SRSdat = NULL
        designRes$SRSdatPar = NULL
        if(resultType == "pixel") {
          designRes[[2]] = designRes$overSampDatPixelUrban
          designRes[[1]] = designRes$overSampDatPixelRural
        }
        else if(resultType == "EA") {
          designRes[[2]] = designRes$overSampDatClusterUrban
          designRes[[1]] = designRes$overSampDatClusterRural
        }
      }
      designResPopAgg = designRes
    }
    if("BYM2 UCA" %in% models)
      designRes = designResTemp
    if("SPDE uc" %in% models) {
      urbanEffect = FALSE
      includeClustEffect = FALSE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, rangeText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spdeNoUrbClust = spdeSRS
      else
        spdeNoUrbClust = spdeOverSamp
    }
    if("SPDE uC" %in% models) {
      urbanEffect = FALSE
      includeClustEffect = TRUE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, rangeText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spdeNoUrb = spdeSRS
      else
        spdeNoUrb = spdeOverSamp
    }
    if("SPDE Uc" %in% models) {
      urbanEffect = TRUE
      includeClustEffect = FALSE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, rangeText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spdeNoClust = spdeSRS
      else
        spdeNoClust = spdeOverSamp
    }
    if("SPDE UC" %in% models) {
      urbanEffect = TRUE
      includeClustEffect = TRUE
      testText = ifelse(test, "Test", "")
      fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                        round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
                        "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, rangeText, ".RData")
      out = load(fileName)
      if(sampling == "SRS")
        spde = spdeSRS
      else
        spde = spdeOverSamp
    }
    # if("SPDE uC'" %in% models) {
    #   urbanEffect = FALSE
    #   includeClustEffect = TRUE
    #   testText = ifelse(test, "Test", "")
    #   fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
    #                     round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
    #                     "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, 
    #                     ".RData")
    #   out = load(fileName)
    #   if(sampling == "SRS")
    #     spdeNoUrbMod = spdeSRS
    #   else
    #     spdeNoUrbMod = spdeOverSamp
    # }
    # if("SPDE UC'" %in% models) {
    #   urbanEffect = TRUE
    #   includeClustEffect = TRUE
    #   testText = ifelse(test, "Test", "")
    #   fileName = paste0("resultsSPDEBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
    #                     round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0", 
    #                     "urbanEffect", urbanEffect, "clustEffect", includeClustEffect, strictPriorText, testText, 
    #                     ".RData")
    #   out = load(fileName)
    #   if(sampling == "SRS")
    #     spdeMod = spdeSRS
    #   else
    #     spdeMod = spdeOverSamp
    # }
    
    # relabel direct, naive, and mercer county names
    counties=sort(unique(poppc$admin1))
    for(i in 1:maxDataSets) {
      if("Smoothed Direct" %in% models)
        mercer[[i]]$admin1 = counties
      if("Direct" %in% models)
        directEst[[i]]$admin1 = counties
      if("Naive" %in% models)
        naive[[i]]$admin1 = counties
    }
    
    # commented out the below code, since it's not clear within strata predictions are necessary
    # # modify discrete estimates stratified by urban/rural if resultType is at finer scale than county
    # # (currently, further resolving discrete estimates is only supported for the direct estimates)
    # if(resultType == "EA") {
    #   filterByUrban = function(i) {
    #     urbanI = 1:nrow(directEstSRS[[i]])
    #   }
    #   directEstSRS = lapply(1:100, )
    # }
    
    # function for generating discrete model pixel and EA level predictions 
    # given the county level predictions. In the case that two result tables are given, 
    # the first is assumed to be for rural areas, and the second is for urban areas. 
    # Also, urbanVec is a vector of urban logical labels for the desired aggregation level
    getSubLevelResults = function(resultTable, resultTableUrban=NULL, urbanVec=NULL) {
      if(!is.null(resultTableUrban) && !is.null(urbanVec)) {
        tabRural = getSubLevelResults(resultTableRural)
        tabUrban = getSubLevelResults(resultTableUrban)
        tabRural[urbanVec,] = tabUrban[UrbanVec,]
        return(tabRural)
      }
      
      # get the order of the counties that resultTable is in
      counties = sort(unique(eaDat$admin1))
      
      # convert results to the desired aggregation level if necessary for discrete models:
      if(resultType == "pixel") {
        # compute pixel results
        pixelToAdmin = match(popGrid$admin1, counties)
        resultTable = resultTable[pixelToAdmin,]
      }
      else if(resultType == "EA") {
        # compute EA level results
        eaToAdmin = match(eaDat$admin1, counties)
        resultTable = resultTable[eaToAdmin,]
      }
      
      # modify result row and column table names according to aggregation level
      whichName = which(names(resultTable) == "admin1")
      names(resultTable)[whichName] = resultType
      if((resultType == "EA") && (is.data.frame(resultTable))) {
        # resultTable[[resultType]] = resultTable$eaI
        resultTable[[resultType]] = 1:nrow(resultTable)
      } else if((resultType == "pixel") && (is.data.frame(resultTable))) {
        # resultTable[[resultType]] = resultTable$eaI
        resultTable[[resultType]] = 1:nrow(resultTable)
      }
      
      if(resultType == "pixel")
        resultTable = resultTable[as.numeric(as.character(truth$pixel)),]
      
      if(is.data.frame(resultTable) && resultType == "county" && is.null(resultTable[[resultType]]))
        resultTable[[resultType]] = counties
      
      resultTable
    }
    
    # convert the truth to the desired aggregation level
    if(resultType == "county")
      truth = getSubLevelResults(truth)
    
    # calculate the binomial n (number of children) for each prediction unit
    numChildren = truth$numChildren
    
    # compute scores
    scoresDirect = scoresNaive = scoresMercer = scoresBYMNoUrb = scoresBYM = scoresBYMNoUrbMod = scoresBYMMod = scoresBYMNoUrbClust = scoresBYMNoClust = 
      scoresBYMNoUrbPopAgg = scoresBYMPopAgg = scoresBYMNoUrbModPopAgg = scoresBYMModPopAgg = scoresBYMNoUrbClustPopAgg = scoresBYMNoClustPopAgg = 
      scoresSPDENoUrbClust = scoresSPDENoUrb = scoresSPDENoClust = scoresSPDE = data.frame()
    
    # convert results to the desired aggregation level 
    # not including urban or cluster effect
    if("BYM2 ucA" %in% models) {
      designResNoUrbClust[[1]]$Q10 = getSubLevelResults(designResNoUrbClust[[1]]$Q10)
      designResNoUrbClust[[1]]$Q50 = getSubLevelResults(designResNoUrbClust[[1]]$Q50)
      designResNoUrbClust[[1]]$Q90 = getSubLevelResults(designResNoUrbClust[[1]]$Q90)
      designResNoUrbClust[[1]]$mean = getSubLevelResults(designResNoUrbClust[[1]]$mean)
      designResNoUrbClust[[1]]$stddev = getSubLevelResults(designResNoUrbClust[[1]]$stddev)
    }
    
    # not including urban effect
    if("BYM2 uCA" %in% models) {
      designResNoUrb[[1]]$Q10 = getSubLevelResults(designResNoUrb[[1]]$Q10)
      designResNoUrb[[1]]$Q50 = getSubLevelResults(designResNoUrb[[1]]$Q50)
      designResNoUrb[[1]]$Q90 = getSubLevelResults(designResNoUrb[[1]]$Q90)
      designResNoUrb[[1]]$mean = getSubLevelResults(designResNoUrb[[1]]$mean)
      designResNoUrb[[1]]$stddev = getSubLevelResults(designResNoUrb[[1]]$stddev)
    }
    
    # including urban effect
    if("BYM2 UcA" %in% models) {
      if(resultType != "county") {
        designResNoClust[[1]]$Q10 = getSubLevelResults(designResNoClust[[1]]$Q10, designResNoClust[[2]]$Q10, truth$urban)
        designResNoClust[[1]]$Q50 = getSubLevelResults(designResNoClust[[1]]$Q50, designResNoClust[[2]]$Q50, truth$urban)
        designResNoClust[[1]]$Q90 = getSubLevelResults(designResNoClust[[1]]$Q90, designResNoClust[[2]]$Q90, truth$urban)
        designResNoClust[[1]]$mean = getSubLevelResults(designResNoClust[[1]]$mean, designResNoClust[[2]]$mean, truth$urban)
        designResNoClust[[1]]$stddev = getSubLevelResults(designResNoClust[[1]]$stddev, designResNoClust[[2]]$stddev, truth$urban)
      }
      else {
        designResNoClust[[1]]$Q10 = getSubLevelResults(designResNoClust[[1]]$Q10)
        designResNoClust[[1]]$Q50 = getSubLevelResults(designResNoClust[[1]]$Q50)
        designResNoClust[[1]]$Q90 = getSubLevelResults(designResNoClust[[1]]$Q90)
        designResNoClust[[1]]$mean = getSubLevelResults(designResNoClust[[1]]$mean)
        designResNoClust[[1]]$stddev = getSubLevelResults(designResNoClust[[1]]$stddev)
      }
    }
    
    # including urban effect
    if("BYM2 UCA" %in% models) {
      if(resultType != "county") {
        designRes[[1]]$Q10 = getSubLevelResults(designRes[[1]]$Q10, designRes[[2]]$Q10, truth$urban)
        designRes[[1]]$Q50 = getSubLevelResults(designRes[[1]]$Q50, designRes[[2]]$Q50, truth$urban)
        designRes[[1]]$Q90 = getSubLevelResults(designRes[[1]]$Q90, designRes[[2]]$Q90, truth$urban)
        designRes[[1]]$mean = getSubLevelResults(designRes[[1]]$mean, designRes[[2]]$mean, truth$urban)
        designRes[[1]]$stddev = getSubLevelResults(designRes[[1]]$stddev, designRes[[2]]$stddev, truth$urban)
      }
      else {
        designRes[[1]]$Q10 = getSubLevelResults(designRes[[1]]$Q10)
        designRes[[1]]$Q50 = getSubLevelResults(designRes[[1]]$Q50)
        designRes[[1]]$Q90 = getSubLevelResults(designRes[[1]]$Q90)
        designRes[[1]]$mean = getSubLevelResults(designRes[[1]]$mean)
        designRes[[1]]$stddev = getSubLevelResults(designRes[[1]]$stddev)
      }
    }
    
    # not including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 uCA'" %in% models) {
      designResNoUrbMod[[1]]$Q10 = getSubLevelResults(designResNoUrbMod[[1]]$Q10)
      designResNoUrbMod[[1]]$Q50 = getSubLevelResults(designResNoUrbMod[[1]]$Q50)
      designResNoUrbMod[[1]]$Q90 = getSubLevelResults(designResNoUrbMod[[1]]$Q90)
      designResNoUrbMod[[1]]$mean = getSubLevelResults(designResNoUrbMod[[1]]$mean)
      designResNoUrbMod[[1]]$stddev = getSubLevelResults(designResNoUrbMod[[1]]$stddev)
    }
    
    # including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 UCA'" %in% models) {
      if(resultType != "county") {
        designResMod[[1]]$Q10 = getSubLevelResults(designResMod[[1]]$Q10, designResMod[[2]]$Q10, truth$urban)
        designResMod[[1]]$Q50 = getSubLevelResults(designResMod[[1]]$Q50, designResMod[[2]]$Q50, truth$urban)
        designResMod[[1]]$Q90 = getSubLevelResults(designResMod[[1]]$Q90, designResMod[[2]]$Q90, truth$urban)
        designResMod[[1]]$mean = getSubLevelResults(designResMod[[1]]$mean, designResMod[[2]]$mean, truth$urban)
        designResMod[[1]]$stddev = getSubLevelResults(designResMod[[1]]$stddev, designResMod[[2]]$stddev, truth$urban)
      }
      else {
        designResMod[[1]]$Q10 = getSubLevelResults(designResMod[[1]]$Q10)
        designResMod[[1]]$Q50 = getSubLevelResults(designResMod[[1]]$Q50)
        designResMod[[1]]$Q90 = getSubLevelResults(designResMod[[1]]$Q90)
        designResMod[[1]]$mean = getSubLevelResults(designResMod[[1]]$mean)
        designResMod[[1]]$stddev = getSubLevelResults(designResMod[[1]]$stddev)
      }
    }
    
    if("BYM2 uca" %in% models) {
      designResNoUrbClustPopAgg[[1]]$Q10 = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$Q10)
      designResNoUrbClustPopAgg[[1]]$Q50 = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$Q50)
      designResNoUrbClustPopAgg[[1]]$Q90 = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$Q90)
      designResNoUrbClustPopAgg[[1]]$mean = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$mean)
      designResNoUrbClustPopAgg[[1]]$stddev = getSubLevelResults(designResNoUrbClustPopAgg[[1]]$stddev)
    }
    
    # not including urban effect
    if("BYM2 uCa" %in% models) {
      designResNoUrbPopAgg[[1]]$Q10 = getSubLevelResults(designResNoUrbPopAgg[[1]]$Q10)
      designResNoUrbPopAgg[[1]]$Q50 = getSubLevelResults(designResNoUrbPopAgg[[1]]$Q50)
      designResNoUrbPopAgg[[1]]$Q90 = getSubLevelResults(designResNoUrbPopAgg[[1]]$Q90)
      designResNoUrbPopAgg[[1]]$mean = getSubLevelResults(designResNoUrbPopAgg[[1]]$mean)
      designResNoUrbPopAgg[[1]]$stddev = getSubLevelResults(designResNoUrbPopAgg[[1]]$stddev)
    }
    
    # including urban effect
    if("BYM2 Uca" %in% models) {
      if(resultType != "county") {
        designResNoClustPopAgg[[1]]$Q10 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q10, designResNoClustPopAgg[[2]]$Q10, truth$urban)
        designResNoClustPopAgg[[1]]$Q50 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q50, designResNoClustPopAgg[[2]]$Q50, truth$urban)
        designResNoClustPopAgg[[1]]$Q90 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q90, designResNoClustPopAgg[[2]]$Q90, truth$urban)
        designResNoClustPopAgg[[1]]$mean = getSubLevelResults(designResNoClustPopAgg[[1]]$mean, designResNoClustPopAgg[[2]]$mean, truth$urban)
        designResNoClustPopAgg[[1]]$stddev = getSubLevelResults(designResNoClustPopAgg[[1]]$stddev, designResNoClustPopAgg[[2]]$stddev, truth$urban)
      }
      else {
        designResNoClustPopAgg[[1]]$Q10 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q10)
        designResNoClustPopAgg[[1]]$Q50 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q50)
        designResNoClustPopAgg[[1]]$Q90 = getSubLevelResults(designResNoClustPopAgg[[1]]$Q90)
        designResNoClustPopAgg[[1]]$mean = getSubLevelResults(designResNoClustPopAgg[[1]]$mean)
        designResNoClustPopAgg[[1]]$stddev = getSubLevelResults(designResNoClustPopAgg[[1]]$stddev)
      }
    }
    
    # including urban effect
    if("BYM2 UCa" %in% models) {
      if(resultType != "county") {
        designResPopAgg[[1]]$Q10 = getSubLevelResults(designResPopAgg[[1]]$Q10, designResPopAgg[[2]]$Q10, truth$urban)
        designResPopAgg[[1]]$Q50 = getSubLevelResults(designResPopAgg[[1]]$Q50, designResPopAgg[[2]]$Q50, truth$urban)
        designResPopAgg[[1]]$Q90 = getSubLevelResults(designResPopAgg[[1]]$Q90, designResPopAgg[[2]]$Q90, truth$urban)
        designResPopAgg[[1]]$mean = getSubLevelResults(designResPopAgg[[1]]$mean, designResPopAgg[[2]]$mean, truth$urban)
        designResPopAgg[[1]]$stddev = getSubLevelResults(designResPopAgg[[1]]$stddev, designResPopAgg[[2]]$stddev, truth$urban)
      }
      else {
        designResPopAgg[[1]]$Q10 = getSubLevelResults(designResPopAgg[[1]]$Q10)
        designResPopAgg[[1]]$Q50 = getSubLevelResults(designResPopAgg[[1]]$Q50)
        designResPopAgg[[1]]$Q90 = getSubLevelResults(designResPopAgg[[1]]$Q90)
        designResPopAgg[[1]]$mean = getSubLevelResults(designResPopAgg[[1]]$mean)
        designResPopAgg[[1]]$stddev = getSubLevelResults(designResPopAgg[[1]]$stddev)
      }
    }
    
    # not including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 uCa'" %in% models) {
      designResNoUrbModPopAgg[[1]]$Q10 = getSubLevelResults(designResNoUrbModPopAgg[[1]]$Q10)
      designResNoUrbModPopAgg[[1]]$Q50 = getSubLevelResults(designResNoUrbModPopAgg[[1]]$Q50)
      designResNoUrbModPopAgg[[1]]$Q90 = getSubLevelResults(designResNoUrbModPopAgg[[1]]$Q90)
      designResNoUrbModPopAgg[[1]]$mean = getSubLevelResults(designResNoUrbModPopAgg[[1]]$mean)
      designResNoUrbModPopAgg[[1]]$stddev = getSubLevelResults(designResNoUrbModPopAgg[[1]]$stddev)
    }
    
    # including urban effect, modified to be debiased using marginal rather than conditional effect as prediction
    if("BYM2 UCa'" %in% models) {
      if(resultType != "county") {
        designResModPopAgg[[1]]$Q10 = getSubLevelResults(designResModPopAgg[[1]]$Q10, designResModPopAgg[[2]]$Q10, truth$urban)
        designResModPopAgg[[1]]$Q50 = getSubLevelResults(designResModPopAgg[[1]]$Q50, designResModPopAgg[[2]]$Q50, truth$urban)
        designResModPopAgg[[1]]$Q90 = getSubLevelResults(designResModPopAgg[[1]]$Q90, designResModPopAgg[[2]]$Q90, truth$urban)
        designResModPopAgg[[1]]$mean = getSubLevelResults(designResModPopAgg[[1]]$mean, designResModPopAgg[[2]]$mean, truth$urban)
        designResModPopAgg[[1]]$stddev = getSubLevelResults(designResModPopAgg[[1]]$stddev, designResModPopAgg[[2]]$stddev, truth$urban)
      }
      else {
        designResModPopAgg[[1]]$Q10 = getSubLevelResults(designResModPopAgg[[1]]$Q10)
        designResModPopAgg[[1]]$Q50 = getSubLevelResults(designResModPopAgg[[1]]$Q50)
        designResModPopAgg[[1]]$Q90 = getSubLevelResults(designResModPopAgg[[1]]$Q90)
        designResModPopAgg[[1]]$mean = getSubLevelResults(designResModPopAgg[[1]]$mean)
        designResModPopAgg[[1]]$stddev = getSubLevelResults(designResModPopAgg[[1]]$stddev)
      }
    }
    
    for(i in c(1:maxDataSets)) { # for problem fitting mercerSRS for SRS sampling, tausq=0
      # for(i in 1:100) {
      if((i %% printIEvery == 0) || (i == 1))
        print(i)
      resultName = paste0(resultType, "Results")
      if(resultType == "EA")
        resultName = "eaResults"
      
      # convert results to the desired aggregation level
      if("Direct" %in% models)
        directEsti = getSubLevelResults(directEst[[i]])
      if("Naive" %in% models)
        naivei = getSubLevelResults(naive[[i]])
      if("Smoothed Direct" %in% models)
        merceri = getSubLevelResults(mercer[[i]])
      if(resultType != "county") {
        # if("SPDE uc" %in% models)
        #   spdeNoUrbClusti = spdeNoUrbClust[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        # if("SPDE uC" %in% models)
        #   spdeNoUrbi = spdeNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        # if("SPDE Uc" %in% models)
        #   spdeNoClusti = spdeNoClust[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        # if("SPDE UC" %in% models)
        #   spdei = spde[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("Direct" %in% models)
          directEsti = getSubLevelResults(directEst[[i]])
      } else {
        # if("SPDE uc" %in% models)
        #   spdeNoUrbClusti = spdeNoUrbClust[[resultName]][[i]]
        # if("SPDE uC" %in% models)
        #   spdeNoUrbi = spdeNoUrb[[resultName]][[i]]
        # if("SPDE Uc" %in% models)
        #   spdeNoClusti = spdeNoClust[[resultName]][[i]]
        # if("SPDE UC" %in% models)
        #   spdei = spde[[resultName]][[i]]
      }
      
      if(resultType == "EA") {
        # set first row of spde results to be the EA index
        # if("SPDE uc" %in% models) {
        #   spdeNoUrbClusti[[resultType]] = 1:nrow(spdeNoUrbClusti)
        #   
        #   whichName = which(names(spdeNoUrbClusti) == "EA")
        #   spdeNoUrbClusti = cbind(spdeNoUrbClusti[,whichName], spdeNoUrbClusti[,-whichName])
        #   
        #   names(spdeNoUrbClusti)[1] = "EA"
        # }
        # if("SPDE uC" %in% models) {
        #   spdeNoUrbi[[resultType]] = 1:nrow(spdeNoUrbi)
        #   
        #   whichName = which(names(spdeNoUrbi) == "EA")
        #   spdeNoUrbi = cbind(spdeNoUrbi[,whichName], spdeNoUrbi[,-whichName])
        #   
        #   names(spdeNoUrbi)[1] = "EA"
        # }
        # if("SPDE Uc" %in% models) {
        #   spdeNoClusti[[resultType]] = 1:nrow(spdeNoClusti)
        #   
        #   whichName = which(names(spdeNoClusti) == "EA")
        #   spdeNoClusti = cbind(spdeNoClusti[,whichName], spdeNoClusti[,-whichName])
        #   
        #   names(spdeNoClusti)[1] = "EA"
        # }
        # if("SPDE UC" %in% models) {
        #   spdei[[resultType]] = 1:nrow(spdei)
        #   
        #   whichName = which(names(spdei) == "EA")
        #   spdei = cbind(spdei[,whichName], spdei[,-whichName])
        #   
        #   names(spdei)[1] = "EA"
        # }
      }
      
      # for spde results, modify the name of the results
      # modify result row and column table names according to aggregation level
      if(resultType == "county") {
        # if("SPDE uc" %in% models) {
        #   whichName = which(names(spdeNoUrbClusti) == "admin1")
        #   names(spdeNoUrbClusti)[whichName] = resultType
        # }
        # if("SPDE uC" %in% models) {
        #   whichName = which(names(spdeNoUrbi) == "admin1")
        #   names(spdeNoUrbi)[whichName] = resultType
        # }
        # if("SPDE Uc" %in% models) {
        #   whichName = which(names(spdeNoClusti) == "admin1")
        #   names(spdeNoClusti)[whichName] = resultType
        # }
        # if("SPDE UC" %in% models) {
        #   # with urban effect:
        #   whichName = which(names(spdei) == "admin1")
        #   names(spdei)[whichName] = resultType
        # }
      }
      
      if(resultType == "pixel") {
        # set first row of spde results to be the pixel index
        # if("SPDE uc" %in% models) {
        #   spdeNoUrbClusti[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdeNoUrbClusti) == "pixel")
        #   spdeNoUrbClusti = cbind(spdeNoUrbClusti[,whichName], spdeNoUrbClusti[,-whichName])
        #   
        #   names(spdeNoUrbClusti)[1] = "pixel"
        # }
        # if("SPDE uC" %in% models) {
        #   spdeNoUrbi[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdeNoUrbi) == "pixel")
        #   spdeNoUrbi = cbind(spdeNoUrbi[,whichName], spdeNoUrbi[,-whichName])
        #   
        #   names(spdeNoUrbi)[1] = "pixel"
        # }
        # if("SPDE Uc" %in% models) {
        #   spdeNoClusti[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdeNoClusti) == "pixel")
        #   spdeNoClusti = cbind(spdeNoClusti[,whichName], spdeNoClusti[,-whichName])
        #   
        #   names(spdeNoClusti)[1] = "pixel"
        # }
        # if("SPDE UC" %in% models) {
        #   spdei[[resultType]] = truth$pixel
        #   
        #   whichName = which(names(spdei) == "pixel")
        #   spdei = cbind(spdei[,whichName], spdei[,-whichName])
        #   
        #   names(spdei)[1] = "pixel"
        # }
      }
      
      # TODO: fix the below code block
      # change names of table variables in spde model with no urban effect to reflect that
      # if("SPDE uc" %in% models)
      #   names(spdeNoUrbClusti)[2:6] = paste0(names(spdeNoUrbClusti)[2:6], " uc")
      # if("SPDE uC" %in% models)
      #   names(spdeNoUrbi)[2:6] = paste0(names(spdeNoUrbi)[2:6], "NoUrb")
      # if("SPDE uC" %in% models)
      #   names(spdeNoUrbi)[2:6] = paste0(names(spdeNoUrbi)[2:6], "NoUrb")
      # if("SPDE uC" %in% models)
      #   names(spdeNoUrbi)[2:6] = paste0(names(spdeNoUrbi)[2:6], "NoUrb")
      
      if("Direct" %in% models) {
        allres = merge(truth, directEsti, by=resultType)
        colnames(allres) = c(resultType, "truth", paste(colnames(allres)[3:8], "Direct", sep=""))
      } else {
        stop("direct estimates must be included at this point in order to name the estimate table columns")
      }
      if("Naive" %in% models)
        allres = merge(allres, naivei, by=resultType)
      if("Smoothed Direct" %in% models)
        allres = merge(allres, merceri, by=resultType)
      # if("SPDE uc" %in% models)
      #   allres = merge(allres, spdeNoUrbClusti, by=resultType)
      # if("SPDE uC" %in% models)
      #   allres = merge(allres, spdeNoUrbi, by=resultType)
      # if("SPDE Uc" %in% models)
      #   allres = merge(allres, spdeNoClusti, by=resultType)
      # if("SPDE UC" %in% models)
      #   allres = merge(allres, spdei, by=resultType)
      
      # set whether or not to calculate scores on logit scale depending on result type
      thisTruth = allres$truth
      
      # if(is.null(useLogit)) {
      #   useLogit = FALSE
      #   if(resultType != "EA" && resultType != "pixel") {
      #     thisTruth = logit(thisTruth)
      #     useLogit=TRUE
      #   }
      # } else if(useLogit == TRUE) {
      #   thisTruth = logit(thisTruth)
      # }
      
      if("Direct" %in% models) {
        my.scoresdirect = getScores(thisTruth, numChildren, allres$logit.estDirect, allres$var.estDirect, nsim=nsim)
        scoresDirect <- rbind(scoresDirect,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresdirect))
      }
      if("Naive" %in% models) {
        my.scoresnaive = getScores(thisTruth, numChildren, allres$logit.est, allres$var.est, nsim=nsim)
        scoresNaive <- rbind(scoresNaive,
                             cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresnaive))
      }
      if("Smoothed Direct" %in% models) {
        my.scoresmercer = getScores(thisTruth, numChildren, allres$logit.est.mercer, allres$var.est.mercer, nsim=nsim)
        scoresMercer <- rbind(scoresMercer,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresmercer))
      }
      if("BYM2 ucA" %in% models) {
        my.scoresbymNoUrbClust = getScores(thisTruth, numChildren, designResNoUrbClust[[1]]$mean[,i], (designResNoUrbClust[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbClust <- rbind(scoresBYMNoUrbClust,
                                     cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbClust))
      }
      if("BYM2 uCA" %in% models) {
        my.scoresbymNoUrb = getScores(thisTruth, numChildren, designResNoUrb[[1]]$mean[,i], (designResNoUrb[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrb <- rbind(scoresBYMNoUrb,
                                cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrb))
      }
      if("BYM2 uCa'" %in% models) {
        my.scoresbymNoUrbMod = getScores(thisTruth, numChildren, designResNoUrbMod[[1]]$mean[,i], (designResNoUrbMod[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbMod <- rbind(scoresBYMNoUrbMod,
                                   cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbMod))
      }
      if("BYM2 UcA" %in% models) {
        my.scoresbymNoClust = getScores(thisTruth, numChildren, designResNoClust[[1]]$mean[,i], (designResNoClust[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoClust <- rbind(scoresBYMNoClust,
                                  cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoClust))
      }
      if("BYM2 UCA" %in% models) {
        my.scoresbym = getScores(thisTruth, numChildren, designRes[[1]]$mean[,i], (designRes[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYM <- rbind(scoresBYM,
                           cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbym))
      }
      if("BYM2 UCA'" %in% models) {
        my.scoresbymMod = getScores(thisTruth, numChildren, designResMod[[1]]$mean[,i], (designResMod[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMMod <- rbind(scoresBYMMod,
                              cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymMod))
      }
      if("BYM2 uca" %in% models) {
        my.scoresbymNoUrbClustPopAgg = getScores(thisTruth, numChildren, designResNoUrbClustPopAgg[[1]]$mean[,i], (designResNoUrbClustPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbClustPopAgg <- rbind(scoresBYMNoUrbClustPopAgg,
                                           cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbClustPopAgg))
      }
      if("BYM2 uCa" %in% models) {
        my.scoresbymNoUrbPopAgg = getScores(thisTruth, numChildren, designResNoUrbPopAgg[[1]]$mean[,i], (designResNoUrbPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbPopAgg <- rbind(scoresBYMNoUrbPopAgg,
                                      cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbPopAgg))
      }
      if("BYM2 uCa'" %in% models) {
        my.scoresbymNoUrbModPopAgg = getScores(thisTruth, numChildren, designResNoUrbModPopAgg[[1]]$mean[,i], (designResNoUrbModPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoUrbModPopAgg <- rbind(scoresBYMNoUrbModPopAgg,
                                         cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoUrbModPopAgg))
      }
      if("BYM2 Uca" %in% models) {
        my.scoresbymNoClustPopAgg = getScores(thisTruth, numChildren, designResNoClustPopAgg[[1]]$mean[,i], (designResNoClustPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMNoClustPopAgg <- rbind(scoresBYMNoClustPopAgg,
                                        cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymNoClustPopAgg))
      }
      if("BYM2 UCa" %in% models) {
        my.scoresbymPopAgg = getScores(thisTruth, numChildren, designResPopAgg[[1]]$mean[,i], (designResPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMPopAgg <- rbind(scoresBYMPopAgg,
                                 cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymPopAgg))
      }
      if("BYM2 UCa'" %in% models) {
        my.scoresbymModPopAgg = getScores(thisTruth, numChildren, designResModPopAgg[[1]]$mean[,i], (designResModPopAgg[[1]]$stddev[,i])^2, nsim=nsim)
        scoresBYMModPopAgg <- rbind(scoresBYMModPopAgg,
                                    cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresbymModPopAgg))
      }
      if("SPDE uC" %in% models) {
        # stop("determine if the spde code should compute all of these directly")
        # my.biasspdeNoUrb = bias(thisTruth, allres$logit.est.spdeNoUrb, logit=useLogit, my.var=allres$var.est.spdeNoUrb)
        # my.msespdeNoUrb = mse(thisTruth, allres$logit.est.spdeNoUrb, logit=useLogit, my.var=allres$var.est.spdeNoUrb)
        # my.dssspdeNoUrb = dss(thisTruth, allres$logit.est.spdeNoUrb, allres$var.est.spdeNoUrb)
        # # the below line used to be commented for some reason
        # my.crpsspdeNoUrb = crpsNormal(thisTruth, allres$logit.est.spdeNoUrb, allres$var.est.spdeNoUrb, logit=useLogit, n=numChildren)
        # my.crpsspdeNoUrb = allres$crps.spdeNoUrb
        # my.coveragespdeNoUrb = coverage(thisTruth, allres$lower.spdeNoUrb, allres$upper.spdeNoUrb, logit=useLogit)
        # my.lengthspdeNoUrb = intervalWidth(allres$lower.spdeNoUrb, allres$upper.spdeNoUrb, logit=useLogit)
        # scoresSPDENoUrb <- rbind(scoresSPDENoUrb,
        #                          cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresspdeNoUrb))
        # scoresSPDENoUrb = rbind(scoresSPDENoUrb # TODO: fix this
      }
      if("SPDE UC" %in% models) {
        # stop("determine if the spde code should compute all of these directly")
        # my.biasspde = bias(thisTruth, allres$logit.est.spde, logit=useLogit, my.var=allres$var.est.spde)
        # my.msespde = mse(thisTruth, allres$logit.est.spde, logit=useLogit, my.var=allres$var.est.spde)
        # my.dssspde = dss(thisTruth, allres$logit.est.spde, allres$var.est.spde)
        # # the below line needs to be commented for some reason
        # my.crpsspde = crpsNormal(thisTruth, allres$logit.est.spde, allres$var.est.spde, logit=useLogit, n=numChildren)
        # my.crpsspde = allres$crps.spde
        # my.coveragespde = coverage(thisTruth, allres$lower.spde, allres$upper.spde, logit=useLogit)
        # my.lengthspde = intervalWidth(allres$lower.spde, allres$upper.spde, logit=useLogit)
        # scoresSPDE <- rbind(scoresSPDE,
        #                     cbind(data.frame(dataset=i, region=allres[[resultType]]), my.scoresspde))
      }
    }
    
    # save progress
    runId = paste0("Beta-1.75margVar", round(margVar, 4), rangeText, "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                   "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                   "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
    
    # first collect all the results. Save everything except for the postprocessing arguments: 
    # produceFigures, digits
    objectNames = ls()
    objectNames = objectNames[-match(c("produceFigures", "xtable.args", "tableFormat", "colScale", 
                                       "colUnits", "colDigits"), objectNames)]
    save(list=objectNames, file=paste0("scoresTemp", runId, ".RData"))
  }
  else {
    # in this case, we have already computed the results so just load them into the environment
    print("Loading results...")
    temp = doFancyTables
    if(loadResults)
      load(paste0("scores", runId, ".RData"))
    else if(loadTempProgress) {
      temp = loadTempProgress
      load(paste0("scoresTemp", runId, ".RData"))
      loadTempProgress = temp
    }
    
    doFancyTables = temp
    allNames = c("Naive", "Direct", "Smoothed Direct", "BYM2 ucA", "BYM2 uCA", "BYM2 uCA'", "BYM2 UcA", "BYM2 UCA", "BYM2 UCA'", 
                 "BYM2 uca", "BYM2 uCa", "BYM2 uCa'", "BYM2 Uca", "BYM2 UCa", "BYM2 UCa'", 
                 "SPDE uc", "SPDE uC", "SPDE Uc", "SPDE UC", "SPDE uC'", "SPDE UC'")
    allNamesBinomial = paste0(allNames, " Bin.")
    allModels = allNames
    models = allModels[modelsI]
    
    # also, don't resave these results that we've already saved
    if(!loadTempProgress)
      saveResults = FALSE
    else
      saveResults = TRUE
  }
  
  # final table: 
  if("Naive" %in% models)
    naive = apply(scoresNaive[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("Direct" %in% models)
    direct = apply(scoresDirect[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("Smoothed Direct" %in% models)
    mercer = apply(scoresMercer[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 ucA" %in% models)
    bymNoUrbClust = apply(scoresBYMNoUrbClust[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCA" %in% models)
    bymNoUrb = apply(scoresBYMNoUrb[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCa'" %in% models)
    bymNoUrbMod = apply(scoresBYMNoUrbMod[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UcA" %in% models)
    bymNoClust = apply(scoresBYMNoClust[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCA" %in% models)
    bym = apply(scoresBYM[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCA'" %in% models)
    bymMod = apply(scoresBYMMod[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uca" %in% models)
    bymNoUrbClustPopAgg = apply(scoresBYMNoUrbClustPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCa" %in% models)
    bymNoUrbPopAgg = apply(scoresBYMNoUrbPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 uCa'" %in% models)
    bymNoUrbModPopAgg = apply(scoresBYMNoUrbModPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 Uca" %in% models)
    bymNoClustPopAgg = apply(scoresBYMNoClustPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCa" %in% models)
    bymPopAgg = apply(scoresBYMPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("BYM2 UCa'" %in% models)
    bymModPopAgg = apply(scoresBYMModPopAgg[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("SPDE uc" %in% models) {
    theseNames = names(spdeNoUrbClust)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeNoUrbClustScores = apply(spdeNoUrbClust[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeNoUrbClustScores = rbind(spdeNoUrbClustScores, apply(spdeNoUrbClust[[i]], 2, mean))
        }
      }
    }
  }
  if("SPDE uC" %in% models || "SPDE uC'" %in% models) {
    theseNames = names(spdeNoUrb)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeNoUrbScores = apply(spdeNoUrb[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeNoUrbScores = rbind(spdeNoUrbScores, apply(spdeNoUrb[[i]], 2, mean))
        }
      }
    }
  }
  if("SPDE Uc" %in% models) {
    theseNames = names(spdeNoClust)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeNoClustScores = apply(spdeNoClust[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeNoClustScores = rbind(spdeNoClustScores, apply(spdeNoClust[[i]], 2, mean))
        }
      }
    }
  }
  if("SPDE UC" %in% models || "SPDE UC'" %in% models) {
    theseNames = names(spde)
    namesI = grepl(tolower(resultType), tolower(theseNames))
    first = TRUE
    for(i in 1:length(theseNames)) {
      if(namesI[i]) {
        if(first) {
          spdeScores = apply(spde[[i]], 2, mean)
          first = FALSE
        }
        else {
          spdeScores = rbind(spdeScores, apply(spde[[i]], 2, mean))
        }
      }
    }
  }
  idx = 1:9
  tab = c()
  if("Naive" %in% models)
    tab = rbind(tab, c(naive[idx]))
  if("Direct" %in% models)
    tab = rbind(tab, c(direct[idx]))
  if("Smoothed Direct" %in% models)
    tab = rbind(tab, c(mercer[idx]))
  if("BYM2 ucA" %in% models)
    tab = rbind(tab, c(bymNoUrbClust[idx]))
  if("BYM2 uCA" %in% models)
    tab = rbind(tab, c(bymNoUrb[idx]))
  if("BYM2 uCa'" %in% models)
    tab = rbind(tab, c(bymNoUrbMod[idx]))
  if("BYM2 UcA" %in% models)
    tab = rbind(tab, c(bymNoClust[idx]))
  if("BYM2 UCA" %in% models)
    tab = rbind(tab, c(bym[idx]))
  if("BYM2 UCA'" %in% models)
    tab = rbind(tab, c(bymMod[idx]))
  if("BYM2 uca" %in% models)
    tab = rbind(tab, c(bymNoUrbClustPopAgg[idx]))
  if("BYM2 uCa" %in% models)
    tab = rbind(tab, c(bymNoUrbPopAgg[idx]))
  if("BYM2 uCa'" %in% models)
    tab = rbind(tab, c(bymNoUrbModPopAgg[idx]))
  if("BYM2 Uca" %in% models)
    tab = rbind(tab, c(bymNoClustPopAgg[idx]))
  if("BYM2 UCa" %in% models)
    tab = rbind(tab, c(bymPopAgg[idx]))
  if("BYM2 UCa'" %in% models)
    tab = rbind(tab, c(bymModPopAgg[idx]))
  # SPDE models are in format 2
  # if("SPDE uc" %in% models)
  #   tab = rbind(tab, c(spdeNoUrbClust[idx]))
  # if("SPDE uC" %in% models)
  #   tab = rbind(tab, c(spdeNoUrb[idx]))
  # if("SPDE Uc" %in% models)
  #   tab = rbind(tab, c(spdeNoClust[idx]))
  # if("SPDE UC" %in% models)
  #   tab = rbind(tab, c(spde[idx]))
  colnames(tab) = c("Bias", "Var", "MSE", "CRPS", "CRPS Bin.", "80\\% Cvg", "80\\% Cvg Bin.", "CI Width", "CI Width Bin.")
  finalNames = allNames[modelsI]
  finalNamesBinomial = allNamesBinomial[modelsI]
  
  if(tableFormat == "1")
    rownames(tab) = finalNames
  else {
    # separate out the binomial and non-binomial models into separate rows
    binomialColumns = c(1:3, 5, 7, 9)
    otherColumns = c(1:3, 4, 6, 8)
    binomialTable = tab[, binomialColumns]
    otherTable = tab[, otherColumns]
    tab = otherTable[numeric(0),]
    thisFinalNames = c()
    
    for(i in 1:nrow(otherTable)) {
      if(includeBVarResults) {
        tab = rbind(tab, 
                    otherTable[i,], 
                    binomialTable[i,])
        thisFinalNames = c(thisFinalNames, finalNames[i], finalNamesBinomial[i])
      }
      else {
        tab = rbind(tab, 
                    otherTable[i,])
        thisFinalNames = c(thisFinalNames, finalNames[i])
      }
    }
    
    # add in SPDE models if necessary
    if("SPDE uc" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE uc")
        tab = rbind(tab, spdeNoUrbClustScores[1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE uc", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeNoUrbClustScores)
      }
    }
    if("SPDE uC" %in% models || "SPDE uC'" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE uC", "SPDE uC'")
        tab = rbind(tab, spdeNoUrbScores[2:1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE uC", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeNoUrbScores)
      }
    }
    if("SPDE Uc" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE Uc")
        tab = rbind(tab, spdeNoClustScores[1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE Uc", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeNoClustScores)
      }
    }
    if("SPDE UC" %in% models || "SPDE UC'" %in% models) {
      if(continuousSPDEonly) {
        thisFinalNames = c(thisFinalNames, "SPDE UC", "SPDE UC'")
        tab = rbind(tab, spdeScores[2:1,])
      }
      else {
        thisFinalNames = c(thisFinalNames, paste0("SPDE UC", c(" Cts.", " Discrete", " Exact")))
        tab = rbind(tab, spdeScores)
      }
    }
    # if("SPDE uC'" %in% models) {
    #   if(continuousSPDEonly) {
    #     thisFinalNames = c(thisFinalNames, "SPDE uC'")
    #     tab = rbind(tab, spdeNoUrbModScores[1,])
    #   }
    #   else {
    #     thisFinalNames = c(thisFinalNames, paste0("SPDE uC'", c(" Cts.", " Discrete", " Exact")))
    #     tab = rbind(tab, spdeNoUrbModScores)
    #   }
    # }
    # if("SPDE UC" %in% models) {
    #   if(continuousSPDEonly) {
    #     thisFinalNames = c(thisFinalNames, "SPDE UC'")
    #     tab = rbind(tab, spdeModScores[1,])
    #   }
    #   else {
    #     thisFinalNames = c(thisFinalNames, paste0("SPDE UC'", c(" Cts.", " Discrete", " Exact")))
    #     tab = rbind(tab, spdeModScores)
    #   }
    # }
    
    rownames(tab) = thisFinalNames
  }
  
  # round the columns of tab and modify the column names to include the scale
  unroundedTab = tab
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(tab[,i] * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
    unroundedTab[,i] = as.numeric(unroundedTab[,i] * colScale[i])
    colnames(unroundedTab)[i] = paste0(colnames(unroundedTab)[i], colUnits[i])
  }
  
  # remove "Direct Bin." model, since direct estimates already account for Binomial variation
  # if("Direct" %and% models) {
  #   rowI = thisFinalNames == "Direct Bin."
  #   tab = tab[!rowI,]
  # }
  
  # remove redundant rows of BYM2 model
  require(stringr)
  modelTypes = word(models, 1)
  modelTypes[modelTypes == "Smoothed"] = "Smoothed Direct"
  uniqueModelTypes = unique(modelTypes)
  modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
  modelVariations = word(models, 2)
  modelVariations[is.na(modelVariations)] = ""
  modelVariations[modelVariations == "Direct"] = ""
  if("BYM2 ucA" %in% models)
    models[models == "BYM2 ucA"] = "BYM2 uc"
  if("BYM2 uCA" %in% models)
    models[models == "BYM2 uCA"] = "BYM2 uC"
  if("BYM2 uCA'" %in% models)
    models[models == "BYM2 uCA'"] = "BYM2 uC'"
  if("BYM2 uca" %in% models) {
    tab = tab[models != "BYM2 uca",]
    unroundedTab = unroundedTab[models != "BYM2 uca",]
    models = models[models != "BYM2 uca"]
  }
  if("BYM2 uCa" %in% models) {
    tab = tab[models != "BYM2 uCa",]
    unroundedTab = unroundedTab[models != "BYM2 uCa",]
    models = models[models != "BYM2 uCa"]
  }
  if("BYM2 uCa'" %in% models) {
    tab = tab[models != "BYM2 uCa'",]
    unroundedTab = unroundedTab[models != "BYM2 uCa'",]
    models = models[models != "BYM2 uCa'"]
  }
  if(identical(modelsI, 1:21)) {
    models = models[c(1:12, 13, 14, 17, 15, 16, 18)]
  }
  rownames(tab) = models
  rownames(unroundedTab) = models
  
  # recalculate model types and variations
  modelTypes = word(models, 1)
  modelTypes[modelTypes == "Smoothed"] = "Smoothed Direct"
  uniqueModelTypes = unique(modelTypes)
  modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
  modelVariations = word(models, 2)
  modelVariations[is.na(modelVariations)] = ""
  modelVariations[modelVariations == "Direct"] = ""
  
  if(!doFancyTables && printScoreTable)
    print(do.call("xtable", c(list(tab), xtable.args)), 
          include.colnames=TRUE,
          hline.after=0, 
          math.style.exponents=TRUE, 
          sanitize.text.function=function(x){x})
  
  if(doFancyTables && printScoreTable) {
    require(stringr)
    require(dplyr)
    require(kableExtra)
    
    options(knitr.table.format = "latex")
    
    # bold the best entries of each column, italicize worst entries of each column
    centers = c(rep(0, 4), 80, 0)
    columnBest = apply(cbind(abs(tab[,1]), tab[,2:4], abs(tab[,5]-80), tab[,6]), 2, min)
    columnWorst = apply(cbind(abs(tab[,1]), tab[,2:4], abs(tab[,5]-80), tab[,6]), 2, max)
    dat = data.table(tab)
    test = dat %>% mutate(Bias = cell_spec(tab[,1], "latex", bold=abs(tab[,1] - centers[1]) <= columnBest[1], italic = abs(tab[,1] - centers[1]) >= columnWorst[1], 
                                           monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          Var = cell_spec(tab[,2], "latex", bold=abs(tab[,2] - centers[2]) <= columnBest[2], italic = abs(tab[,2] - centers[2]) >= columnWorst[2], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          MSE = cell_spec(tab[,3], "latex", bold=abs(tab[,3] - centers[3]) <= columnBest[3], italic = abs(tab[,3] - centers[3]) >= columnWorst[3], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          CRPS = cell_spec(tab[,4], "latex", bold=abs(tab[,4] - centers[4]) <= columnBest[4], italic = abs(tab[,4] - centers[4]) >= columnWorst[4], 
                                           monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          CVG = cell_spec(tab[,5], "latex", bold=abs(tab[,5] - centers[5]) <= columnBest[5], italic = abs(tab[,5] - centers[5]) >= columnWorst[5], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          Width = cell_spec(tab[,6], "latex", bold=abs(tab[,6] - centers[6]) <= columnBest[6], italic = abs(tab[,6] - centers[6]) >= columnWorst[6], 
                                            monospace=FALSE, underline=FALSE, strikeout=FALSE)) %>%
      select(Bias, Var, MSE, CRPS, CVG, Width)
    
    # revert the column names to their true values, set the model variations to be the values in the first column
    colnames(test) = colnames(tab)
    test = cbind(" "=modelVariations, test)
    rownames(test)=NULL
    
    # group the rows by the type of model
    fullTab = test %>%
      kable("latex", escape = F, booktabs = T) %>% kable_styling()
    for(i in 1:length(uniqueModelTypes)) {
      startR = min(modelTypeGroups[[i]])
      endR = max(modelTypeGroups[[i]])
      fullTab = fullTab %>% pack_rows(uniqueModelTypes[i], startR, endR, latex_gap_space = "2em")
    }
    if(printScoreTable)
      print(fullTab)
  }
  
  ## append parameter tables from each smoothing model if necessary
  anySmoothingModels = as.numeric("Naive" %in% models) + as.numeric("Direct" %in% models)
  anySmoothingModels = length(models) > anySmoothingModels
  parTab = c()
  if(anySmoothingModels) {
    parRowNames = c()
    if("Smoothed Direct" %in% models) {
      parTab = rbind(parTab, mercerPar)
      parRowNames = c(parRowNames, rep("Smoothed Direct", nrow(mercerPar)))
    }
    if("BYM2 ucA" %in% models || "BYM2 uca" %in% models || "BYM2 uc" %in% models) {
      parTab = rbind(parTab, designResNoUrbClust[[length(designResNoUrbClust)]])
      parRowNames = c(parRowNames, rep("BYM2 uc", nrow(designResNoUrbClust[[length(designResNoUrbClust)]])))
    }
    if("BYM2 uCA" %in% models || "BYM2 uCA'" %in% models || "BYM2 uCa" %in% models || "BYM2 uCa'" %in% models  || "BYM2 uC" %in% models || "BYM2 uC'" %in% models) {
      parTab = rbind(parTab, designResNoUrb[[length(designResNoUrb)]])
      parRowNames = c(parRowNames, rep("BYM2 uC", nrow(designResNoUrb[[length(designResNoUrb)]])))
    }
    if("BYM2 UcA" %in% models || "BYM2 Uca" %in% models || "BYM2 Uca'" %in% models || "BYM2 Ucb'" %in% models) {
      parTab = rbind(parTab, designResNoClust[[length(designResNoClust)]])
      parRowNames = c(parRowNames, rep("BYM2 Uc", nrow(designResNoClust[[length(designResNoClust)]])))
    }
    if("BYM2 UCA" %in% models || "BYM2 UCa" %in% models || "BYM2 UCA'" %in% models || "BYM2 UCa'" %in% models) {
      parTab = rbind(parTab, designRes[[length(designRes)]])
      parRowNames = c(parRowNames, rep("BYM2 UC", nrow(designRes[[length(designRes)]])))
    }
    spdeParIndices = c(1:2, 4:6) # leave out variance and width
    if("SPDE uc" %in% models) {
      thisParTab = matrix(spdeNoUrbClust$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spdeNoUrbClust$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spdeNoUrbClust$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spdeNoUrbClust$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spdeNoUrbClust$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spdeNoUrbClust$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spdeNoUrbClust$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spdeNoUrbClust$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrbClust$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE uc", nrow(thisParTab)))
    }
    if("SPDE uC" %in% models || "SPDE uC'" %in% models) {
      thisParTab = matrix(spdeNoUrb$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spdeNoUrb$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spdeNoUrb$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spdeNoUrb$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spdeNoUrb$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spdeNoUrb$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spdeNoUrb$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spdeNoUrb$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoUrb$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE uC", nrow(thisParTab)))
    }
    if("SPDE Uc" %in% models) {
      thisParTab = matrix(spdeNoClust$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spdeNoClust$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spdeNoClust$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spdeNoClust$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spdeNoClust$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spdeNoClust$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spdeNoClust$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spdeNoClust$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spdeNoClust$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE Uc", nrow(thisParTab)))
    }
    if("SPDE UC" %in% models || "SPDE UC'" %in% models) {
      thisParTab = matrix(spde$interceptSummary[spdeParIndices], nrow=1)
      theseRowNames = "Intercept"
      if(!is.null(spde$urbanSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, Urban=spde$urbanSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Urban")
      }
      if(!is.null(spde$rangeSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$rangeSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Range")
      }
      if(!is.null(spde$varSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$varSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial Var")
      }
      if(!is.null(spde$sdSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$sdSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Spatial SD")
      }
      if(!is.null(spde$nuggetVarSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$nuggetVarSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster Var")
      }
      if(!is.null(spde$nuggetSDSummary[spdeParIndices])) {
        thisParTab = rbind(thisParTab, spde$nuggetSDSummary[spdeParIndices])
        theseRowNames = c(theseRowNames, "Cluster SD")
      }
      colnames(thisParTab) = names(parTab)
      rownames(thisParTab) = theseRowNames
      parTab = rbind(parTab, thisParTab)
      parRowNames = c(parRowNames, rep("SPDE UC", nrow(thisParTab)))
    }
    # if("SPDE uC'" %in% models) {
    #   thisParTab = matrix(spdeNoUrbMod$interceptSummary[spdeParIndices], nrow=1)
    #   theseRowNames = "Intercept"
    #   if(!is.null(spdeNoUrbMod$urbanSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, Urban=spdeNoUrbMod$urbanSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Urban")
    #   }
    #   if(!is.null(spdeNoUrbMod$rangeSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$rangeSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Range")
    #   }
    #   if(!is.null(spdeNoUrbMod$varSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$varSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial Var")
    #   }
    #   if(!is.null(spdeNoUrbMod$sdSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$sdSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial SD")
    #   }
    #   if(!is.null(spdeNoUrbMod$nuggetVarSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$nuggetVarSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster Var")
    #   }
    #   if(!is.null(spdeNoUrbMod$nuggetSDSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeNoUrbMod$nuggetSDSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster SD")
    #   }
    #   colnames(thisParTab) = names(parTab)
    #   rownames(thisParTab) = theseRowNames
    #   parTab = rbind(parTab, thisParTab)
    #   parRowNames = c(parRowNames, rep("SPDE uC'", nrow(thisParTab)))
    # }
    # if("SPDE UC" %in% models) {
    #   thisParTab = matrix(spdeMod$interceptSummary[spdeParIndices], nrow=1)
    #   theseRowNames = "Intercept"
    #   if(!is.null(spdeMod$urbanSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, Urban=spdeMod$urbanSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Urban")
    #   }
    #   if(!is.null(spdeMod$rangeSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$rangeSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Range")
    #   }
    #   if(!is.null(spdeMod$varSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$varSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial Var")
    #   }
    #   if(!is.null(spdeMod$sdSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$sdSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Spatial SD")
    #   }
    #   if(!is.null(spdeMod$nuggetVarSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$nuggetVarSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster Var")
    #   }
    #   if(!is.null(spdeMod$nuggetSDSummary[spdeParIndices])) {
    #     thisParTab = rbind(thisParTab, spdeMod$nuggetSDSummary[spdeParIndices])
    #     theseRowNames = c(theseRowNames, "Cluster SD")
    #   }
    #   colnames(thisParTab) = names(parTab)
    #   rownames(thisParTab) = theseRowNames
    #   parTab = rbind(parTab, thisParTab)
    #   parRowNames = c(parRowNames, rep("SPDE UC'", nrow(thisParTab)))
    # }
    
    # add the model to the row names, remove the numbers at the end of the duplicated row names, print out the aggregated parameter table
    rownames(parTab) = paste(parRowNames, rownames(parTab))
    for(i in 1:nrow(parTab)) {
      lastCharacter = substr(rownames(parTab)[i], nchar(rownames(parTab)[i]), nchar(rownames(parTab)[i]))
      if(grepl("\\d", lastCharacter))
        rownames(parTab)[i] = substr(rownames(parTab)[i], 1, nchar(rownames(parTab)[i])-1)
    }
    
    if(!doFancyTables && printParTable)
      print(xtable(parTab, digits=2, display=c("s", rep("fg", ncol(parTab)))))
    
    if(doFancyTables && printParTable) {
      # now make the fancy table using kableExtra by grouping the rows by models and model variations
      # fancyParTable = xtable2kable(xtable(parTab, digits=2, display=c("s", rep("fg", ncol(parTab)))))
      # fancyParTable = xtable2kable(xtable(parTab, digits=2, display=c("s", rep("e", ncol(parTab)))))
      
      # determine the model type and model variation groupings
      SmoothDirectI = (1:nrow(parTab))[grepl("Smoothed Direct", rownames(parTab))]
      BYM2I = (1:nrow(parTab))[grepl("BYM2", rownames(parTab))]
      BYM2I = setdiff(BYM2I, SmoothDirectI)
      SPDEI = (1:nrow(parTab))[grepl("SPDE", rownames(parTab))]
      modelTypes = rep("Smoothed Direct", nrow(parTab))
      modelTypes[BYM2I] = "BYM2"
      modelTypes[SPDEI] = "SPDE"
      uniqueModelTypes = unique(modelTypes)
      modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
      modelVariations = rep("", nrow(parTab))
      modelVariations[BYM2I] = word(rownames(parTab)[BYM2I], 2)
      modelVariations[SPDEI] = word(rownames(parTab)[SPDEI], 2)
      
      # fix BYM2 switch up between total variance and phi when no cluster effect is included:
      # switch data, but keep the perimeter names at the same rows
      noCluster = grepl("uc", rownames(parTab)) | grepl("Uc", rownames(parTab))
      badRows = which(noCluster & grepl("Phi", rownames(parTab)))
      currentNames = rownames(parTab)
      switchRows = parTab[badRows+1,]
      parTab[badRows+1,] = parTab[badRows,]
      parTab[badRows,] = switchRows
      rownames(parTab) = currentNames
      
      # determine the parameter names
      require("tm")
      parNames = trimws(removeWords(rownames(parTab), c(uniqueModelTypes, unique(modelVariations))))
      
      # round intercept, phi, and range parameters to 3, 3, and 0 digits respectively
      parTab[parNames == "Intercept",] = round(parTab[parNames == "Intercept",], digits=3)
      parTab[parNames == "Phi",] = round(parTab[parNames == "Phi",], digits=3)
      parTab[parNames == "Range",] = round(parTab[parNames == "Range",], digits=0)
      parTab[parNames == "Urban",] = round(parTab[parNames == "Urban",], digits=3)
      
      # round everything else to approximately three significant figures, paying special note 
      # to the non-SD quantities, which should be rounded out to the same decimal
      otherPar = (parNames != "Intercept") & (parNames != "Phi") & (parNames != "Range") & (parNames != "Urban")
      parTab[otherPar,2] = signif(parTab[otherPar,2], digits=3)
      parTab[otherPar,-2] = t(apply(parTab[otherPar,-2], 1, roundToFirstSigFigs, digits=3))
      
      # add the grouping variables and parameter names as new columns
      parTab = cbind(" "=modelTypes, " "=modelVariations, " "=parNames, parTab)
      rownames(parTab) = NULL
      
      # group the rows by model type and model variation
      row_group_label_fonts <-list(list(bold = T, italic = T), list(bold = F, italic = F))
      print(kable(parTab, "latex", booktabs = T, escape=FALSE, format.args=list(drop0trailing=TRUE, scientific=FALSE), 
                  longtable=TRUE, caption = "Longtable") %>%
              collapse_rows(1:2, row_group_label_position ='stack', latex_hline ='custom', custom_latex_hline = 1:2, 
                            row_group_label_fonts = row_group_label_fonts) %>%
              kable_styling(latex_options =c("repeat_header")))
    }
  }
  
  runId = paste0("Beta-1.75margVar", round(margVar, 4), rangeText, "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                 "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                 "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets)
  if(saveResults) {
    # first collect all the results. Save everything except for the postprocessing arguments: 
    # produceFigures, digits
    objectNames = ls()
    objectNames = objectNames[-match(c("produceFigures", "xtable.args", "tableFormat", "colScale", 
                                       "colUnits", "colDigits"), objectNames)]
    save(list=objectNames, file=paste0("scores", runId, ".RData"))
  }
  
  if(produceFigures) {
    # compare all five
    pdf(paste0("figures/biasbyregion", runId, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(bias~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    boxplot(bias~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    boxplot(bias~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    # axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    abline(h=0, lwd=2, col=2)
    dev.off()
    
    pdf(paste0("figures/crpsbyregion", runId, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(crps~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    if("Naive" %in% models)
      boxplot(crps~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    if("Smoothed Direct" %in% models)
      boxplot(crps~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    if("BYM" %in% models)
      boxplot(crps~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    if("SPDE" %in% models)
      boxplot(crps~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    dev.off()
  }
  
  list(tab=tab, parTab=parTab, unroundedTab=unroundedTab)
}

# This model produces an xtable summarizing the scoring rules aggregated at the given level of interest
# test: whether or not to use the test data set based results or the true
#       simulation study results. The test data set results are based on fewer
#       observations, requiring fewer computations and should be used to test the code
# tausq: the spatial nugget in the simulated data, can only be the default or 0
# resultType: the level of aggregation over which scoring rules should be computed
# sampling: whether or not to use the simple random sample or urban over sampled results
# recomputeTruth: by default the should be set to TRUE, unless the user calls this function 
#                 twice in a row with the same set of inputs
# modelsI: which model results should be included. Models come in the order:
#          c("naive", "direct", "mercer", "bym", "bymNoUrb", "spde", "spdeNoUrb"), 
#          so the default value of this argument is to include results for the naive, 
#          direct, and mercer models. 1:7 is all models
# produceFigures: whether or not to produce figures based on the results
# logitScale: whether or not to produce scoring rule results at the logit or count scale. 
#           By default this is to use loaded scale results unless the resultType is EA or pixel
# printIEvery: how often to print progress
runCompareModels = function(test=FALSE, tausq=.1^2, resultType=c("county", "pixel", "EA"), 
                            sampling=c("SRS", "oversamp"), recomputeTruth=TRUE, modelsI=1:3, 
                            produceFigures=FALSE, logitScale=NULL, big=FALSE, printIEvery=50) {
  
  # match the arguments with their correct values
  resultType = match.arg(resultType)
  sampling = match.arg(sampling)
  useLogit = logitScale
  
  # test to make sure only the naive and direct models are included if big is set to TRUE
  if(!identical(modelsI, 1:2) && big)
    stop("if big is set to TRUE, only naive and direct results can be included")
  
  testText = ifelse(test, "Test", "")
  bigText = ifelse(big, "Big", "")
  if(tausq == .1^2 || tausq == .01) {
    out = load(paste0("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0", testText, bigText, ".RData"))
  }
  else {
    if(tausq != 0)
      stop("tausq can only be equal to .1^2 or 0")
    
    out = load(paste0("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0", testText, bigText, ".RData"))
  }
  eaDat = SRSDat$eaDat
  
  if(sampling == "SRS") {
    clustDat = SRSDat
  } else {
    clustDat = overSampDat
  }
  
  ## generate the true values at the county level for the 
  ## given settings if these are new settings
  if(recomputeTruth){
    if(resultType == "county") {
      # compute truth based on superpopulation
      regions = sort(unique(eaDat$admin1))
      truthbycounty <- rep(NA, 47)
      childrenPerCounty = truthbycounty
      
      for(i in 1:47){
        super = eaDat[eaDat$admin1 == regions[i],]
        childrenPerCounty[i] = sum(super$numChildren)
        truthbycounty[i] <- sum(super$died)/childrenPerCounty[i]
      }
      truth = data.frame(admin1=regions, truth=truthbycounty, numChildren=childrenPerCounty)
      save(truth, file="truthbycounty.RData")
    } else if(resultType == "pixel") {
      counties = sort(unique(eaDat$admin1))
      eaToPixel = eaDat$pixelI
      childrenPerPixel = tapply(eaDat$numChildren, list(pixel=eaDat$pixelI), sum)
      urbanPixel = tapply(eaDat$urban, list(pixel=eaDat$pixelI), function(x) {mean(x[1])})
      deathsPerPixel = tapply(eaDat$died, list(pixel=eaDat$pixelI), sum)
      regions = names(childrenPerPixel) # these are the pixels with enumeration areas in them
      pixelToAdmin = match(popGrid$admin1[as.numeric(regions)], counties)
      
      truth = data.frame(pixel=regions, truth=deathsPerPixel / childrenPerPixel, countyI=pixelToAdmin, urban=urbanPixel, numChildren=childrenPerPixel)
      save(truth, file="truthbyPixel.RData")
    } else if(resultType == "EA") {
      truth = data.frame(EA = 1:nrow(eaDat), truth = eaDat$died/25, urban=eaDat$urban, numChildren=eaDat$numChildren)
      save(truth, file="truthbyEA.RData")
    }
  } else {
    # load the truth
    load(paste0("truthby", resultType, ".RData"))
  }
  
  ## Now pick which models to include in the results
  allModels = c("naive", "direct", "mercer", "bym", "bymNoUrb", "bymNoClust", "bymNoUrbClust", "spde", "spdeNoUrb")
  allNames = c("Naive", "Direct estimates", "Smoothed Direct", "Model-based BYM", "Model-based BYM (no urban effect)", "Model-based BYM (no cluster effect)", "Model-based BYM (no urban or cluster effect)", "SPDE", "SPDE (no urban effect)")
  models = allModels[modelsI]
  
  # load data
  if(tausq == .1^2 || tausq == 0.01) {
    if(test) {
      if("naive" %in% models || "direct" %in% models)
        out = load(paste0("resultsDirectNaiveTausq0.01Test", bigText, ".RData")) # this is the only case that uses the big dataset
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.01test.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterTRUETest.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterTRUETest.RData")
      if("bymNoUrbClust" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterFALSETest.RData")
        designResNoUrb = designRes
      }
      if("bymNoClust" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterFALSETest.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0.01urbanEffectFALSETest.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0.01urbanEffectTRUETest.RData")
    } else {
      if("naive" %in% models || "direct" %in% models)
        out = load(paste0("resultsDirectNaiveTausq0.01", bigText, ".RData"))
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.01.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterTRUE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterTRUE.RData")
      if("bymNoUrbClust" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurFALSEClusterFALSE.RData")
        designResNoUrb = designRes
      }
      if("bymNoClust" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0.01UrbRurTRUEClusterFALSE.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0.01urbanEffectFALSE.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0.01urbanEffectTRUE.RData")
    }
  } else if(tausq == 0) {
    if(test) {
      if("naive" %in% models || "direct" %in% models)
        out = load(paste0("resultsDirectNaiveTausq0Test", bigText, ".RData")) # this is the only case that uses the big dataset
      if("mercer" %in% models)
        out = load("resultsMercerTausq0test.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurFATestLSE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurTRUETest.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0urbanEffectFALSETest.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0urbanEffectTRUETest.RData")
    } else {
      if("naive" %in% models || "direct" %in% models)
        out = load(paste0("resultsDirectNaiveTausq0", bigText, ".RData")) # this is the only case that uses the big dataset
      if("mercer" %in% models)
        out = load("resultsMercerTausq0.RData")
      if("bymNoUrb" %in% models) {
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurFALSE.RData")
        designResNoUrb = designRes
      }
      if("bym" %in% models)
        out = load("kenyaSpatialDesignResultNewTausq0UrbRurTRUE.RData")
      if("spdeNoUrb" %in% models) {
        out = load("resultsSPDETausq0urbanEffectFALSE.RData")
        spdeSRSNoUrb = spdeSRS
        spdeOverSampNoUrb = spdeOverSamp
      }
      if("spde" %in% models)
        out = load("resultsSPDETausq0urbanEffectTRUE.RData")
    }
  }
  
  # commented out the below code, since it's not clear within strata predictions are necessary
  # # modify discrete estimates stratified by urban/rural if resultType is at finer scale than county
  # # (currently, further resolving discrete estimates is only supported for the direct estimates)
  # if(resultType == "EA") {
  #   filterByUrban = function(i) {
  #     urbanI = 1:nrow(directEstSRS[[i]])
  #   }
  #   directEstSRS = lapply(1:100, )
  # }
  
  # function for generating discrete model pixel and EA level predictions 
  # given the county level predictions
  getSubLevelResults = function(resultTable) {
    # get the order of the counties that resultTable is in
    counties = sort(unique(eaDat$admin1))
    
    # convert results to the desired aggregation level if necessary for discrete models:
    if(resultType == "pixel") {
      # compute pixel results
      pixelToAdmin = match(popGrid$admin1, counties)
      resultTable = resultTable[pixelToAdmin,]
    }
    else if(resultType == "EA") {
      # compute EA level results
      eaToAdmin = match(eaDat$admin1, counties)
      resultTable = resultTable[eaToAdmin,]
    }
    
    # modify result row and column table names according to aggregation level
    whichName = which(names(resultTable) == "admin1")
    names(resultTable)[whichName] = resultType
    if((resultType == "EA") && (is.data.frame(resultTable))) {
      # resultTable[[resultType]] = resultTable$eaI
      resultTable[[resultType]] = 1:nrow(resultTable)
    } else if((resultType == "pixel") && (is.data.frame(resultTable))) {
      # resultTable[[resultType]] = resultTable$eaI
      resultTable[[resultType]] = 1:nrow(resultTable)
    }
    
    if(resultType == "pixel")
      resultTable = resultTable[as.numeric(as.character(truth$pixel)),]
    
    resultTable
  }
  
  # convert the truth to the desired aggregation level
  if(resultType == "county")
    truth = getSubLevelResults(truth)
  
  # calculate the binomial n (number of children) for each prediction unit
  numChildren = truth$numChildren
  
  # compute scores
  scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMNoUrbSRS = scoresBYMSRS = scoresSPDENoUrbSRS = scoresSPDESRS = data.frame()
  scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMNoUrboverSamp = scoresBYMoverSamp = scoresSPDENoUrboverSamp = scoresSPDEoverSamp = data.frame()
  # scoresDirectSRS = scoresNaiveSRS = scoresMercerSRS = scoresBYMSRS = data.frame()
  # scoresDirectoverSamp = scoresNaiveoverSamp = scoresMerceroverSamp = scoresBYMoverSamp = data.frame()
  
  # convert results to the desired aggregation level
  # not including urban effect
  if("bymNoUrb" %in% models) {
    designResNoUrb[[1]]$Q10 = getSubLevelResults(designResNoUrb[[1]]$Q10)
    designResNoUrb[[1]]$Q50 = getSubLevelResults(designResNoUrb[[1]]$Q50)
    designResNoUrb[[1]]$Q90 = getSubLevelResults(designResNoUrb[[1]]$Q90)
    designResNoUrb[[1]]$mean = getSubLevelResults(designResNoUrb[[1]]$mean)
    designResNoUrb[[1]]$stddev = getSubLevelResults(designResNoUrb[[1]]$stddev)
    designResNoUrb$overSampDat$Q10 = getSubLevelResults(designResNoUrb$overSampDat$Q10)
    designResNoUrb$overSampDat$Q50 = getSubLevelResults(designResNoUrb$overSampDat$Q50)
    designResNoUrb$overSampDat$Q90 = getSubLevelResults(designResNoUrb$overSampDat$Q90)
    designResNoUrb$overSampDat$mean = getSubLevelResults(designResNoUrb$overSampDat$mean)
    designResNoUrb$overSampDat$stddev = getSubLevelResults(designResNoUrb$overSampDat$stddev)
  }
  
  # including urban effect
  if("bym" %in% models) {
    designRes[[1]]$Q10 = getSubLevelResults(designRes[[1]]$Q10)
    designRes[[1]]$Q50 = getSubLevelResults(designRes[[1]]$Q50)
    designRes[[1]]$Q90 = getSubLevelResults(designRes[[1]]$Q90)
    designRes[[1]]$mean = getSubLevelResults(designRes[[1]]$mean)
    designRes[[1]]$stddev = getSubLevelResults(designRes[[1]]$stddev)
    designRes$overSampDat$Q10 = getSubLevelResults(designRes$overSampDat$Q10)
    designRes$overSampDat$Q50 = getSubLevelResults(designRes$overSampDat$Q50)
    designRes$overSampDat$Q90 = getSubLevelResults(designRes$overSampDat$Q90)
    designRes$overSampDat$mean = getSubLevelResults(designRes$overSampDat$mean)
    designRes$overSampDat$stddev = getSubLevelResults(designRes$overSampDat$stddev)
  }
  
  for(i in c(1:length(clustDat$clustDat))) { # for problem fitting mercerSRS for SRS sampling, tausq=0
    # for(i in 1:100) {
    if((i %% printIEvery == 0) || (i == 1))
      print(i)
    resultName = paste0(resultType, "Results")
    if(resultType == "EA")
      resultName = "eaResults"
    
    # convert results to the desired aggregation level
    if(sampling == "SRS") {
      if("direct" %in% models)
        directEstSRSi = getSubLevelResults(directEstSRS[[i]])
      if("naive" %in% models)
        naiveSRSi = getSubLevelResults(naiveSRS[[i]])
      if("mercer" %in% models)
        mercerSRSi = getSubLevelResults(mercerSRS[[i]])
      if(resultType != "county") {
        if("spdeNoUrb" %in% models)
          spdeSRSNoUrbi = spdeSRSNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("spde" %in% models)
          spdeSRSi = spdeSRS[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("direct" %in% models)
          directEstSRSi = getSubLevelResults(directEstSRS[[i]])
      } else {
        if("spdeNoUrb" %in% models)
          spdeSRSNoUrbi = spdeSRSNoUrb[[resultName]][[i]]
        if("spde" %in% models)
          spdeSRSi = spdeSRS[[resultName]][[i]]
      }
    } else if(sampling == "oversamp") {
      if("direct" %in% models) {
        directEstoverSampi = getSubLevelResults(directEstoverSamp[[i]])
        # directEstoverSampUrbani = getSubLevelResults(directEstoverSampUrb[[i]])
        # directEstoverSampRurali = getSubLevelResults(directEstoverSampRural[[i]])
      }
      if("naive" %in% models)
        naiveoverSampi = getSubLevelResults(naiveoverSamp[[i]])
      if("mercer" %in% models)
        merceroverSampi = getSubLevelResults(merceroverSamp[[i]])
      if(resultType != "county") {
        if("spdeNoUrb" %in% models)
          spdeOverSampNoUrbi = spdeOverSampNoUrb[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
        if("spde" %in% models)
          spdeOverSampi = spdeOverSamp[[resultName]][[i]][as.numeric(as.character(truth[[resultType]])),]
      } else {
        if("spdeNoUrb" %in% models)
          spdeOverSampNoUrbi = spdeOverSampNoUrb[[resultName]][[i]]
        if("spde" %in% models)
          spdeOverSampi = spdeOverSamp[[resultName]][[i]]
      }
    }
    
    if(resultType == "EA") {
      # set first row of spde results to be the EA index
      if(sampling == "SRS") {
        if("spdeNoUrb" %in% models) {
          spdeSRSNoUrbi[[resultType]] = 1:nrow(spdeSRSNoUrbi)
          
          whichName = which(names(spdeSRSNoUrbi) == "EA")
          spdeSRSNoUrbi = cbind(spdeSRSNoUrbi[,whichName], spdeSRSNoUrbi[,-whichName])
          
          names(spdeSRSNoUrbi)[1] = "EA"
        }
        if("spde" %in% models) {
          spdeSRSi[[resultType]] = 1:nrow(spdeSRSi)
          
          whichName = which(names(spdeSRSi) == "EA")
          spdeSRSi = cbind(spdeSRSi[,whichName], spdeSRSi[,-whichName])
          
          names(spdeSRSi)[1] = "EA"
        }
      } else {
        if("spdeNoUrb" %in% models) {
          spdeOverSampNoUrbi[[resultType]] = 1:nrow(spdeOverSampNoUrbi)
          
          whichName = which(names(spdeOverSampNoUrbi) == "EA")
          spdeOverSampNoUrbi = cbind(spdeOverSampNoUrbi[,whichName], spdeOverSampNoUrbi[,-whichName])
          
          names(spdeOverSampNoUrbi)[1] = "EA"
        }
        if("spde" %in% models) {
          spdeOverSampi[[resultType]] = 1:nrow(spdeOverSampi)
          
          whichName = which(names(spdeOverSampi) == "EA")
          spdeOverSampi = cbind(spdeOverSampi[,whichName], spdeOverSampi[,-whichName])
          
          names(spdeOverSampi)[1] = "EA"
        }
      }
    }
    
    # for spde results, modify the name of the results
    # modify result row and column table names according to aggregation level
    if(resultType == "county") {
      # without urban effect:
      if("spdeNoUrb" %in% models) {
        if(sampling == "SRS") {
          whichName = which(names(spdeSRSNoUrbi) == "admin1")
          names(spdeSRSNoUrbi)[whichName] = resultType
        } else if (sampling == "oversamp") {
          whichName = which(names(spdeOverSampNoUrbi) == "admin1")
          names(spdeOverSampNoUrbi)[whichName] = resultType
        }
      }
      
      if("spde" %in% models) {
        # with urban effect:
        if(sampling == "SRS") {
          whichName = which(names(spdeSRSi) == "admin1")
          names(spdeSRSi)[whichName] = resultType
        } else if(sampling == "oversamp") {
          whichName = which(names(spdeOverSampi) == "admin1")
          names(spdeOverSampi)[whichName] = resultType
        }
      }
    }
    
    if(resultType == "pixel") {
      # set first row of spde results to be the pixel index
      if(sampling == "SRS") {
        if("spdeNoUrb" %in% models) {
          spdeSRSNoUrbi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeSRSNoUrbi) == "pixel")
          spdeSRSNoUrbi = cbind(spdeSRSNoUrbi[,whichName], spdeSRSNoUrbi[,-whichName])
          
          names(spdeSRSNoUrbi)[1] = "pixel"
        }
        if("spde" %in% models) {
          spdeSRSi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeSRSi) == "pixel")
          spdeSRSi = cbind(spdeSRSi[,whichName], spdeSRSi[,-whichName])
          
          names(spdeSRSi)[1] = "pixel"
        }
      } else {
        if("spdeNoUrb" %in% models) {
          spdeOverSampNoUrbi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeOverSampNoUrbi) == "pixel")
          spdeOverSampNoUrbi = cbind(spdeOverSampNoUrbi[,whichName], spdeOverSampNoUrbi[,-whichName])
          
          names(spdeOverSampNoUrbi)[1] = "pixel"
        }
        if("spde" %in% models) {
          spdeOverSampi[[resultType]] = truth$pixel
          
          whichName = which(names(spdeOverSampi) == "pixel")
          spdeOverSampi = cbind(spdeOverSampi[,whichName], spdeOverSampi[,-whichName])
          
          names(spdeOverSampi)[1] = "pixel"
        }
      }
    }
    
    # change names of table variables in spde model with no urban effect to reflect that
    if(sampling == "SRS") {
      if("spdeNoUrb" %in% models)
        names(spdeSRSNoUrbi)[2:6] = paste0(names(spdeSRSNoUrbi)[2:6], "NoUrb")
    } else {
      if("spdeNoUrb" %in% models)
        names(spdeOverSampNoUrbi)[2:6] = paste0(names(spdeOverSampNoUrbi)[2:6], "NoUrb")
    }
    
    if(sampling == "SRS") {
      if("direct" %in% models) {
        allresSRS = merge(truth, directEstSRSi, by=resultType)
        colnames(allresSRS) = c(resultType, "truth", paste(colnames(allresSRS)[3:8], "direct", sep=""))
      } else {
        stop("direct estimates must be included at this point in order to name the estimate table columns")
      }
      if("naive" %in% models)
        allresSRS = merge(allresSRS, naiveSRSi, by=resultType)
      if("mercer" %in% models)
        allresSRS = merge(allresSRS, mercerSRSi, by=resultType)
      if("spdeNoUrb" %in% models)
        allresSRS = merge(allresSRS, spdeSRSNoUrbi, by=resultType)
      if("spde" %in% models)
        allresSRS = merge(allresSRS, spdeSRSi, by=resultType)
    }
    
    if(sampling == "oversamp") {
      if("direct" %in% models) {
        allresoverSamp = merge(truth, directEstoverSampi, by=resultType)
        colnames(allresoverSamp) = c(resultType, "truth", paste(colnames(allresoverSamp)[3:8], "direct", sep=""))
      } else {
        stop("direct estimates must be included at this point in order to name the estimate table columns")
      }
      if("naive" %in% models)
        allresoverSamp = merge(allresoverSamp, naiveoverSampi, by=resultType)
      if("mercer" %in% models)
        allresoverSamp = merge(allresoverSamp, merceroverSampi, by=resultType)
      if("spdeNoUrb" %in% models)
        allresoverSamp = merge(allresoverSamp, spdeOverSampNoUrbi, by=resultType)
      if("spde" %in% models)
        allresoverSamp = merge(allresoverSamp, spdeOverSampi, by=resultType)
    }
    
    # set whether or not to calculate scores on logit scale depending on result type
    if(sampling == "SRS")
      thisTruth = allresSRS$truth
    else
      thisTruth = allresoverSamp$truth
    
    if(is.null(useLogit)) {
      useLogit = FALSE
      if(resultType != "EA" && resultType != "pixel") {
        thisTruth = logit(thisTruth)
        useLogit=TRUE
      }
    } else if(useLogit == TRUE) {
      thisTruth = logit(thisTruth)
    }
    
    # SRS setting 
    if(sampling == "SRS") {
      if("direct" %in% models) {
        my.biasSRSdirect = bias(thisTruth, allresSRS$logit.estDirect, logit=useLogit, my.var=allresSRS$var.estDirect)
        my.mseSRSdirect = mse(thisTruth, allresSRS$logit.estDirect, logit=useLogit, my.var=allresSRS$var.estDirect)
        my.dssSRSdirect = dss(thisTruth, allresSRS$logit.estDirect, allresSRS$var.estDirect)
        my.crpsSRSdirect = crpsNormal(thisTruth, allresSRS$logit.estDirect, allresSRS$var.estDirect, logit=useLogit, n=numChildren)
        my.coverageSRSdirect = coverage(thisTruth, allresSRS$lowerdirect, allresSRS$upperdirect, doLogit=useLogit)
        my.lengthSRSdirect = intervalWidth(allresSRS$lowerdirect, allresSRS$upperdirect, logit=useLogit)
      }
      
      if("naive" %in% models) {
        my.biasSRSnaive = bias(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
        my.mseSRSnaive = mse(thisTruth, allresSRS$logit.est, logit=useLogit, my.var=allresSRS$var.est)
        my.dssSRSnaive = dss(thisTruth, allresSRS$logit.est, allresSRS$var.est)
        my.crpsSRSnaive = crpsNormal(thisTruth, allresSRS$logit.est, allresSRS$var.est, logit=useLogit, n=numChildren)
        my.coverageSRSnaive = coverage(thisTruth, allresSRS$lower, allresSRS$upper, doLogit=useLogit)
        my.lengthSRSnaive = intervalWidth(allresSRS$lower, allresSRS$upper, logit=useLogit)
      }
      
      if("mercer" %in% models) {
        my.biasSRSmercer = bias(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
        my.mseSRSmercer = mse(thisTruth, allresSRS$logit.est.mercer, logit=useLogit, my.var=allresSRS$var.est.mercer)
        my.dssSRSmercer = dss(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer)
        my.crpsSRSmercer = crpsNormal(thisTruth, allresSRS$logit.est.mercer, allresSRS$var.est.mercer, logit=useLogit, n=numChildren)
        my.coverageSRSmercer = coverage(thisTruth, allresSRS$lower.mercer, allresSRS$upper.mercer, doLogit=useLogit)
        my.lengthSRSmercer = intervalWidth(allresSRS$lower.mercer, allresSRS$upper.mercer, logit=useLogit)
      }
      
      if("bymNoUrb" %in% models) {
        my.biasSRSbymNoUrb = bias(thisTruth, designResNoUrb[[1]]$mean[,i], logit=useLogit, my.var=(designResNoUrb[[1]]$stddev[,i])^2)
        my.mseSRSbymNoUrb = mse(thisTruth, designResNoUrb[[1]]$mean[,i], logit=useLogit, my.var=(designResNoUrb[[1]]$stddev[,i])^2)
        my.dssSRSbymNoUrb = dss(thisTruth, designResNoUrb[[1]]$mean[,i], (designResNoUrb[[1]]$stddev[,i])^2)
        my.crpsSRSbymNoUrb = crpsNormal(thisTruth, designResNoUrb[[1]]$mean[,i], (designResNoUrb[[1]]$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageSRSbymNoUrb = coverage(thisTruth, designResNoUrb[[1]]$Q10[,i],designResNoUrb[[1]]$Q90[,i], doLogit=useLogit)
        my.lengthSRSbymNoUrb = intervalWidth(designResNoUrb[[1]]$Q10[,i], designResNoUrb[[1]]$Q90[,i], logit=useLogit)
      }
      
      if("bym" %in% models) {
        my.biasSRSbym = bias(thisTruth, designRes[[1]]$mean[,i], logit=useLogit, my.var=(designRes[[1]]$stddev[,i])^2)
        my.mseSRSbym = mse(thisTruth, designRes[[1]]$mean[,i], logit=useLogit, my.var=(designRes[[1]]$stddev[,i])^2)
        my.dssSRSbym = dss(thisTruth, designRes[[1]]$mean[,i], (designRes[[1]]$stddev[,i])^2)
        my.crpsSRSbym = crpsNormal(thisTruth, designRes[[1]]$mean[,i], (designRes[[1]]$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageSRSbym = coverage(thisTruth, designRes[[1]]$Q10[,i],designRes[[1]]$Q90[,i], doLogit=useLogit)
        my.lengthSRSbym = intervalWidth(designRes[[1]]$Q10[,i], designRes[[1]]$Q90[,i], logit=useLogit)
      }
      
      if("spdeNoUrb" %in% models) {
        my.biasSRSspdeNoUrb = bias(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
        my.mseSRSspdeNoUrb = mse(thisTruth, allresSRS$logit.est.spdeNoUrb, logit=useLogit, my.var=allresSRS$var.est.spdeNoUrb)
        my.dssSRSspdeNoUrb = dss(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb)
        # the below line used to be commented for some reason
        my.crpsSRSspdeNoUrb = crpsNormal(thisTruth, allresSRS$logit.est.spdeNoUrb, allresSRS$var.est.spdeNoUrb, logit=useLogit, n=numChildren)
        my.crpsSRSspdeNoUrb = allresSRS$crps.spdeNoUrb
        my.coverageSRSspdeNoUrb = coverage(thisTruth, allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, doLogit=useLogit)
        my.lengthSRSspdeNoUrb = intervalWidth(allresSRS$lower.spdeNoUrb, allresSRS$upper.spdeNoUrb, logit=useLogit)
      }
      
      if("spde" %in% models) {
        my.biasSRSspde = bias(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
        my.mseSRSspde = mse(thisTruth, allresSRS$logit.est.spde, logit=useLogit, my.var=allresSRS$var.est.spde)
        my.dssSRSspde = dss(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde)
        # the below line needs to be commented for some reason
        my.crpsSRSspde = crpsNormal(thisTruth, allresSRS$logit.est.spde, allresSRS$var.est.spde, logit=useLogit, n=numChildren)
        my.crpsSRSspde = allresSRS$crps.spde
        my.coverageSRSspde = coverage(thisTruth, allresSRS$lower.spde, allresSRS$upper.spde, doLogit=useLogit)
        my.lengthSRSspde = intervalWidth(allresSRS$lower.spde, allresSRS$upper.spde, logit=useLogit)
      }
      
      if("direct" %in% models) {
        scoresDirectSRS <- rbind(scoresDirectSRS,
                                 data.frame(dataset=i, region=allresSRS[[resultType]],
                                            bias=my.biasSRSdirect,
                                            mse=my.mseSRSdirect,
                                            dss=my.dssSRSdirect,
                                            coverage=my.coverageSRSdirect,
                                            var=mean(allresSRS$var.estdirect),
                                            crps=my.crpsSRSdirect, 
                                            length=my.lengthSRSdirect))
      }
      
      if("naive" %in% models) {
        scoresNaiveSRS <- rbind(scoresNaiveSRS,
                                data.frame(dataset=i, region=allresSRS[[resultType]],
                                           bias=my.biasSRSnaive,
                                           mse=my.mseSRSnaive,
                                           dss=my.dssSRSnaive,
                                           coverage=my.coverageSRSnaive,
                                           var=mean(allresSRS$var.est),
                                           crps=my.crpsSRSnaive, 
                                           length=my.lengthSRSnaive))
      }
      
      if("mercer" %in% models) {
        scoresMercerSRS <- rbind(scoresMercerSRS,
                                 data.frame(dataset=i, region=allresSRS[[resultType]],
                                            bias=my.biasSRSmercer,
                                            mse=my.mseSRSmercer,
                                            dss=my.dssSRSmercer,
                                            coverage=my.coverageSRSmercer,
                                            var=mean(allresSRS$var.est.mercer),
                                            crps=my.crpsSRSmercer, 
                                            length=my.lengthSRSmercer))
      }
      
      if("bymNoUrb" %in% models) {
        scoresBYMNoUrbSRS <- rbind(scoresBYMNoUrbSRS,
                                   data.frame(dataset=i, region=allresSRS[[resultType]],
                                              bias=my.biasSRSbymNoUrb,
                                              mse=my.mseSRSbymNoUrb,
                                              dss=my.dssSRSbymNoUrb,
                                              coverage=my.coverageSRSbymNoUrb,
                                              var=mean((designResNoUrb[[1]]$stddev[,i])^2),
                                              crps=my.crpsSRSbymNoUrb, 
                                              length=my.lengthSRSbymNoUrb))
      }
      
      if("bym" %in% models) {
        scoresBYMSRS <- rbind(scoresBYMSRS,
                              data.frame(dataset=i, region=allresSRS[[resultType]],
                                         bias=my.biasSRSbym,
                                         mse=my.mseSRSbym,
                                         dss=my.dssSRSbym,
                                         coverage=my.coverageSRSbym,
                                         var=mean((designRes[[1]]$stddev[,i])^2),
                                         crps=my.crpsSRSbym, 
                                         length=my.lengthSRSbym))
      }
      
      if("spdeNoUrb" %in% models) {
        scoresSPDENoUrbSRS <- rbind(scoresSPDENoUrbSRS, 
                                    data.frame(dataset=i, region=allresSRS[[resultType]], 
                                               bias=my.biasSRSspdeNoUrb, 
                                               mse=my.mseSRSspdeNoUrb,
                                               dss=my.dssSRSspdeNoUrb,
                                               coverage=my.coverageSRSspdeNoUrb, 
                                               var=mean(allresSRS$var.est.spdeNoUrb),
                                               crps=my.crpsSRSspdeNoUrb, 
                                               length=my.lengthSRSspdeNoUrb))
      }
      
      if("spde" %in% models) {
        scoresSPDESRS <- rbind(scoresSPDESRS, 
                               data.frame(dataset=i, region=allresSRS[[resultType]], 
                                          bias=my.biasSRSspde, 
                                          mse=my.mseSRSspde,
                                          dss=my.dssSRSspde,
                                          coverage=my.coverageSRSspde, 
                                          var=mean(allresSRS$var.est.spde),
                                          crps=my.crpsSRSspde, 
                                          length=my.lengthSRSspde))
      }
      
    } else if(sampling == "oversamp") {
      # oversampling setting
      if("direct" %in% models) {
        my.biasoverSampdirect = bias(thisTruth, allresoverSamp$logit.estdirect, logit=useLogit, my.var=allresoverSamp$var.estdirect)
        my.mseoverSampdirect = mse(thisTruth, allresoverSamp$logit.estdirect, logit=useLogit, my.var=allresoverSamp$var.estdirect)
        my.dssoverSampdirect = dss(thisTruth, allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect)
        my.crpsoverSampdirect = crpsNormal(thisTruth, allresoverSamp$logit.estdirect, allresoverSamp$var.estdirect, logit=useLogit, n=numChildren)
        my.coverageoverSampdirect = coverage(thisTruth, allresoverSamp$upperdirect, allresoverSamp$lowerdirect, doLogit=useLogit)
        my.lengthoverSampdirect = intervalWidth(allresoverSamp$lowerdirect, allresoverSamp$upperdirect, logit=useLogit)
      }
      
      if("naive" %in% models) {
        my.biasoverSampnaive = bias(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
        my.mseoverSampnaive = mse(thisTruth, allresoverSamp$logit.est, logit=useLogit, my.var=allresoverSamp$var.est)
        my.dssoverSampnaive = dss(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est)
        my.crpsoverSampnaive = crpsNormal(thisTruth, allresoverSamp$logit.est, allresoverSamp$var.est, logit=useLogit, n=numChildren)
        my.coverageoverSampnaive = coverage(thisTruth, allresoverSamp$upper, allresoverSamp$lower, doLogit=useLogit)
        my.lengthoverSampnaive = intervalWidth(allresoverSamp$lower, allresoverSamp$upper, logit=useLogit)
      }
      
      if("mercer" %in% models) {
        my.biasoverSampmercer = bias(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
        my.mseoverSampmercer = mse(thisTruth, allresoverSamp$logit.est.mercer, logit=useLogit, my.var=allresoverSamp$var.est.mercer)
        my.dssoverSampmercer = dss(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer)
        my.crpsoverSampmercer = crpsNormal(thisTruth, allresoverSamp$logit.est.mercer, allresoverSamp$var.est.mercer, logit=useLogit, n=numChildren)
        my.coverageoverSampmercer = coverage(thisTruth, allresoverSamp$lower.mercer, allresoverSamp$upper.mercer, doLogit=useLogit)
        my.lengthoverSampmercer = intervalWidth(allresoverSamp$lower.mercer, allresoverSamp$upper.mercer, logit=useLogit)
      }
      
      if("bymNoUrb" %in% models) {
        my.biasoverSampbymNoUrb = bias(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
        my.mseoverSampbymNoUrb = mse(thisTruth, designResNoUrb$overSampDat$mean[,i], logit=useLogit, my.var=(designResNoUrb$overSampDat$stddev[,i])^2)
        my.dssoverSampbymNoUrb = dss(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2)
        my.crpsoverSampbymNoUrb = crpsNormal(thisTruth, designResNoUrb$overSampDat$mean[,i], (designResNoUrb$overSampDat$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageoverSampbymNoUrb = coverage(thisTruth, designResNoUrb$overSampDat$Q10[,i],designResNoUrb$overSampDat$Q90[,i], doLogit=useLogit)
        my.lengthoverSampbymNoUrb = intervalWidth(designResNoUrb$overSampDat$Q10[,i], designResNoUrb$overSampDat$Q90[,i], logit=useLogit)
      }
      
      if("bym" %in% models) {
        my.biasoverSampbym = bias(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
        my.mseoverSampbym = mse(thisTruth, designRes$overSampDat$mean[,i], logit=useLogit, my.var=(designRes$overSampDat$stddev[,i])^2)
        my.dssoverSampbym = dss(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2)
        my.crpsoverSampbym = crpsNormal(thisTruth, designRes$overSampDat$mean[,i], (designRes$overSampDat$stddev[,i])^2, logit=useLogit, n=numChildren)
        my.coverageoverSampbym = coverage(thisTruth, designRes$overSampDat$Q10[,i],designRes$overSampDat$Q90[,i], doLogit=useLogit)
        my.lengthoverSampbym = intervalWidth(designRes$overSampDat$Q10[,i], designRes$overSampDat$Q90[,i], logit=useLogit)
      }
      
      if("spdeNoUrb" %in% models) {
        my.biasoverSampspdeNoUrb = bias(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
        my.mseoverSampspdeNoUrb = mse(thisTruth, allresoverSamp$logit.est.spdeNoUrb, logit=useLogit, my.var=allresoverSamp$var.est.spdeNoUrb)
        my.dssoverSampspdeNoUrb = dss(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb)
        my.crpsoverSampspdeNoUrb = crpsNormal(thisTruth, allresoverSamp$logit.est.spdeNoUrb, allresoverSamp$var.est.spdeNoUrb, logit=useLogit, n=numChildren)
        my.coverageoverSampspdeNoUrb = coverage(thisTruth, allresoverSamp$lower.spdeNoUrb, allresoverSamp$upper.spdeNoUrb, doLogit=useLogit)
        my.lengthoverSampspdeNoUrb = intervalWidth(allresoverSamp$lower.spdeNoUrb, allresoverSamp$upper.spdeNoUrb, logit=useLogit)
      }
      
      if("spde" %in% models) {
        my.biasoverSampspde = bias(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
        my.mseoverSampspde = mse(thisTruth, allresoverSamp$logit.est.spde, logit=useLogit, my.var=allresoverSamp$var.est.spde)
        my.dssoverSampspde = dss(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde)
        my.crpsoverSampspde = crpsNormal(thisTruth, allresoverSamp$logit.est.spde, allresoverSamp$var.est.spde, logit=useLogit, n=numChildren)
        my.coverageoverSampspde = coverage(thisTruth, allresoverSamp$lower.spde, allresoverSamp$upper.spde, doLogit=useLogit)
        my.lengthoverSampspde = intervalWidth(allresoverSamp$lower.spde, allresoverSamp$upper.spde, logit=useLogit)
      }
      
      if("direct" %in% models) {
        scoresDirectoverSamp <- rbind(scoresDirectoverSamp, 
                                      data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                 bias=my.biasoverSampdirect, 
                                                 mse=my.mseoverSampdirect,
                                                 dss=my.dssoverSampdirect,
                                                 coverage=my.coverageoverSampdirect, 
                                                 var=mean(allresoverSamp$var.estdirect),
                                                 crps=my.crpsoverSampdirect, 
                                                 length=my.lengthoverSampdirect))
      }
      
      if("naive" %in% models) {
        scoresNaiveoverSamp<- rbind(scoresNaiveoverSamp, 
                                    data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                               bias=my.biasoverSampnaive, 
                                               mse=my.mseoverSampnaive,
                                               dss=my.dssoverSampnaive,
                                               coverage=my.coverageoverSampnaive, 
                                               var=mean(allresoverSamp$var.est),
                                               crps=my.crpsoverSampnaive, 
                                               length=my.lengthoverSampnaive))
      }
      
      if("mercer" %in% models) {
        scoresMerceroverSamp<- rbind(scoresMerceroverSamp, 
                                     data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                bias=my.biasoverSampmercer, 
                                                mse=my.mseoverSampmercer,
                                                dss=my.dssoverSampmercer,
                                                coverage=my.coverageoverSampmercer, 
                                                var=mean(allresoverSamp$var.est.mercer),
                                                crps=my.crpsoverSampmercer, 
                                                length=my.lengthoverSampmercer))
      }
      
      if("bymNoUrb" %in% models) {
        scoresBYMNoUrboverSamp <- rbind(scoresBYMNoUrboverSamp, 
                                        data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                   bias=my.biasoverSampbymNoUrb, 
                                                   mse=my.mseoverSampbymNoUrb,
                                                   dss=my.dssoverSampbymNoUrb,
                                                   coverage=my.coverageoverSampbymNoUrb, 
                                                   var=mean((designResNoUrb$overSampDat$stddev[,i])^2),
                                                   crps=my.crpsoverSampbymNoUrb, 
                                                   length=my.lengthoverSampbymNoUrb))
      }
      
      if("bym" %in% models) {
        scoresBYMoverSamp <- rbind(scoresBYMoverSamp, 
                                   data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                              bias=my.biasoverSampbym, 
                                              mse=my.mseoverSampbym,
                                              dss=my.dssoverSampbym,
                                              coverage=my.coverageoverSampbym, 
                                              var=mean((designRes$overSampDat$stddev[,i])^2),
                                              crps=my.crpsoverSampbym, 
                                              length=my.lengthoverSampbym))
      }
      
      if("spdeNoUrb" %in% models) {
        scoresSPDENoUrboverSamp <- rbind(scoresSPDENoUrboverSamp, 
                                         data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                                    bias=my.biasoverSampspdeNoUrb, 
                                                    mse=my.mseoverSampspdeNoUrb,
                                                    dss=my.dssoverSampspdeNoUrb,
                                                    coverage=my.coverageoverSampspdeNoUrb, 
                                                    var=mean(allresoverSamp$var.est.spdeNoUrb),
                                                    crps=my.crpsoverSampspdeNoUrb, 
                                                    length=my.lengthoverSampspdeNoUrb))
      }
      
      
      if("spde" %in% models) {
        scoresSPDEoverSamp <- rbind(scoresSPDEoverSamp, 
                                    data.frame(dataset=i, region=allresoverSamp[[resultType]], 
                                               bias=my.biasoverSampspde, 
                                               mse=my.mseoverSampspde,
                                               dss=my.dssoverSampspde,
                                               coverage=my.coverageoverSampspde, 
                                               var=mean(allresoverSamp$var.est.spde),
                                               crps=my.crpsoverSampspde, 
                                               length=my.lengthoverSampspde))
      }
      
    }
  }
  
  # the below code is commented out, since we want all scoring rules on the empirical proportion scale
  # if((resultType == "EA" || resultType == "pixel") && !useLogit) {
  #   # convert all variances to binomial scale
  #   N = 25
  #   if(sampling == "SRS") {
  #     scoresDirectSRS$var = scoresDirectSRS$var * N
  #     scoresNaiveSRS$var = scoresNaiveSRS$var * N
  #     scoresMercerSRS$var = scoresMercerSRS$var * N
  #     scoresBYMNoUrbSRS$var = scoresBYMNoUrbSRS$var * N
  #     scoresBYMSRS$var = scoresBYMSRS$var * N
  #     scoresSPDENoUrbSRS$var = scoresSPDENoUrbSRS$var * N
  #     scoresSPDESRS$var = scoresSPDESRS$var * N
  #   } else if(sampling == "oversamp") {
  #     scoresDirectoverSamp$var = scoresDirectoverSamp$var * N
  #     scoresNaiveoverSamp$var = scoresNaiveoverSamp$var * N
  #     scoresMerceroverSamp$var = scoresMerceroverSamp$var * N
  #     scoresBYMNoUrboverSamp$var = scoresBYMNoUrboverSamp$var * N
  #     scoresBYMoverSamp$var = scoresBYMoverSamp$var * N
  #     scoresSPDENoUrboverSamp$var = scoresSPDENoUrboverSamp$var * N
  #     scoresSPDEoverSamp$var = scoresSPDEoverSamp$var * N
  #   }
  # }
  
  # final table: 
  if("naive" %in% models)
    naive = apply(scoresNaiveSRS[, c("bias", "mse", "dss", "crps", "var", "coverage", "length")], 2, mean)
  if("direct" %in% models)
    direct = apply(scoresDirectSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("mercer" %in% models)
    mercer = apply(scoresMercerSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bymNoUrb" %in% models)
    bymNoUrb = apply(scoresBYMNoUrbSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("bym" %in% models)
    bym = apply(scoresBYMSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("spdeNoUrb" %in% models)
    spdeNoUrb = apply(scoresSPDENoUrbSRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  if("spde" %in% models)
    spde = apply(scoresSPDESRS[, c("bias", "var", "mse", "crps", "crpsB", "coverage", "coverageB", "length", "lengthB")], 2, mean)
  idx = 1:9
  tab = c()
  if("naive" %in% models)
    tab = rbind(tab, c(naive[idx]))
  if("direct" %in% models)
    tab = rbind(tab, c(direct[idx]))
  if("mercer" %in% models)
    tab = rbind(tab, c(mercer[idx]))
  if("bym" %in% models)
    tab = rbind(tab, c(bym[idx]))
  if("bymNoUrb" %in% models)
    tab = rbind(tab, c(bymNoUrb[idx]))
  if("bymNoClust" %in% models)
    tab = rbind(tab, c(bymNoClust[idx]))
  if("bymNoUrbClust" %in% models)
    tab = rbind(tab, c(bymNoUrbClust[idx]))
  if("spde" %in% models)
    tab = rbind(tab, c(spde[idx]))
  if("spdeNoUrb" %in% models)
    tab = rbind(tab, c(spdeNoUrb[idx]))
  
  rownames(tab) = allNames[modelsI]
  
  print(xtable(tab, digits=3), 
        only.contents=TRUE, 
        include.colnames=TRUE,
        hline.after=NULL)
  
  if(produceFigures) {
    # compare all five
    pdf(paste0("figures/biasbyregion", sampling, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(bias~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    boxplot(bias~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    boxplot(bias~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(bias~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    # axis(1, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    abline(h=0, lwd=2, col=2)
    dev.off()
    
    pdf(paste0("figures/crpsbyregion", sampling, ".pdf"), width=20, height=12)
    par(mar=c(10,4,2,1), las=2, cex.lab=2, cex.axis=1.4, cex.main=2)
    boxplot(crps~region, data=scoresDirect, at=seq(-1, 275, by=6), 
            col="yellow", xlim=c(0,279), names=FALSE, xaxt="n")
    boxplot(crps~region, data=scoresNaive, at=seq(0, 276, by=6), col="orange", xlim=c(0,230), add=TRUE)
    boxplot(crps~region, data=scoresMercer, at=seq(1, 277, by=6), col="green", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(crps~region, data=scoresBYM, at=seq(2, 278, by=6), col="lightblue", xlim=c(0,230), add=TRUE, xaxt="n")
    boxplot(crps~region, data=scoresSPDE, at=seq(3, 279, by=6), col="purple", xlim=c(0,230), add=TRUE, xaxt="n")
    # axis(2, at=seq(0.5, 276.5, by=6), labels=scoresDirect$region[1:47])
    legend("top", c("Direct estimates", "Naive", "Mercer", "BYM", "SPDE"),
           fill = c("yellow", "orange", "green", "lightblue", "purple"), ncol=4, cex=2)
    dev.off()
  }
  
  tab
}

getTruth = function(resultType = c("county", "region", "EA", "pixel"), eaDat) {
  resultType = match.arg(resultType)
  
  if(resultType == "county") {
    # compute truth based on superpopulation
    regions = sort(unique(eaDat$admin1))
    truthbycounty <- rep(NA, 47)
    childrenPerCounty = truthbycounty
    
    for(i in 1:47){
      super = eaDat[eaDat$admin1 == regions[i],]
      childrenPerCounty[i] = sum(super$numChildren)
      truthbycounty[i] <- sum(super$died)/childrenPerCounty[i]
    }
    truth = data.frame(admin1=regions, truth=truthbycounty, numChildren=childrenPerCounty)
  } else if(resultType == "pixel") {
    eaToPixel = eaDat$pixelI
    pixelsWithData = unique(eaToPixel)
    counties = sort(unique(eaDat$admin1))
    childrenPerPixel = aggregate(eaDat$numChildren, list(pixel=eaDat$pixelI), sum)
    urbanPixel = aggregate(eaDat$urban, list(pixel=eaDat$pixelI), function(x) {mean(x[1])})
    deathsPerPixel = aggregate(eaDat$died, list(pixel=eaDat$pixelI), sum)
    regions = childrenPerPixel$pixel # these are the pixels with enumeration areas in them
    # sort results by pixels with data
    sortI=match(pixelsWithData, regions)
    
    pixelToAdmin = match(popGrid$admin1[as.numeric(regions)], counties)
    
    truth = data.frame(pixel=regions, truth=deathsPerPixel$x / childrenPerPixel$x, countyI=pixelToAdmin, urban=urbanPixel$x, numChildren=childrenPerPixel$x)
    truth = truth[sortI,]
  } else if(resultType == "EA") {
    truth = data.frame(EA = 1:nrow(eaDat), truth = eaDat$died/eaDat$numChildren, urban=eaDat$urban, numChildren=eaDat$numChildren)
  } else if(resultType == "region") {
    regions = sort(unique(eaDat$region))
    truthbyregion <- rep(NA, 8)
    childrenPerRegion = truthbyregion
    
    for(i in 1:8){
      super = eaDat[eaDat$region == regions[i],]
      childrenPerRegion[i] = sum(super$numChildren)
      truthbyregion[i] <- sum(super$died)/childrenPerRegion[i]
    }
    truth = data.frame(admin1=regions, truth=truthbyregion, numChildren=childrenPerRegion)
  }
  
  truth
}

# takes all precomputed scoring rules from compareModels2 listed in compareModelCommandArgs.RData, and outputs results
# (use this for the parameter tables)
runCompareModelsAllLocal = function(indices=NULL, strictPriors=FALSE, doFancyTables=TRUE, printScoreTable=TRUE, 
                                    printParTable=TRUE, printBigResults=TRUE, spatialRange=150, spatialVar=0.15^2) {
  spatialRange = match.arg(as.character(spatialRange), choices=c(150, 50))
  spatialVar = match.arg(as.character(spatialVar), choices=c(0.15^2, 0.3^2))
  
  if(spatialRange == 150 && spatialVar == 0.15^2)
    load("compareModelCommandArgs.RData")
  else if(spatialRange == 50 && spatialVar == 0.3^2)
    stop("spatialRange == 50 && spatialVar == 0.3^2 not supported")
  else
    load("compareModelCommandArgsNew.RData")
  
  if(is.null(indices))
    indices = 1:length(compareModelCommandArgs)
  for(i in indices) {
    # get the arguments for the run, and specify that we want to load the precomputed results
    argList = compareModelCommandArgs[[i]]
    argList$loadResults = TRUE
    argList$strictPriors = strictPriors
    argList$doFancyTables = doFancyTables
    argList$printScoreTable = printScoreTable
    argList$printParTable = printParTable
    
    if(!printBigResults && argList$big == TRUE)
      next
    
    # get all elements from the list
    tausq = argList$tausq
    gamma = argList$gamma
    margVar = argList$margVar
    test = argList$test
    resultType = argList$resultType
    sampling = argList$sampling
    recomputeTruth = argList$recomputeTruth
    modelsI = argList$modelsI
    produceFigures = argList$produceFigures
    big = argList$big
    maxDataSets = argList$maxDataSets 
    nsim = argList$nsim
    range = argList$effRange
    if(is.null(range))
      range = 150
    
    # skip nonexistent populations and populations not in our scenario
    if(range != spatialRange && margVar != 0)
      next
    if(range == 50 && margVar == 0)
      argList$effRange = 150
    if(margVar != 0 && margVar != spatialVar)
      next
    
    # print out the population and design for this table
    if(margVar == 0 && gamma == 0 && tausq == 0)
      popText = "suc"
    else if(gamma == 0 && tausq == 0)
      popText = "Suc"
    else if(tausq == 0)
      popText = "SUc"
    else
      popText = "SUC"
    
    # skip nonexistent populations
    if(range == 50 && margVar == 0)
      next
    
    contextText = paste(popText, sampling)
    print(paste0("Printing table for ", contextText))
    
    # print out the table
    browser()
    do.call("runCompareModels2", argList)
  }
}

# (use this for the scoring rules)
runCompareModelsLocal2 = function(indices = NULL, strictPriors = FALSE, filterRows=c(1:3, 4, 6, 10, 12, 13:16), 
                                  incorrectlyAggregatedModels=TRUE, spatialRange=150, spatialVar=0.15^2) {
  spatialRange = match.arg(as.character(spatialRange), choices=c(150, 50))
  spatialVar = match.arg(as.character(spatialVar), choices=c(0.15^2, 0.3^2))
  
  if(spatialRange == 150 && spatialVar == 0.15^2)
    load("compareModelCommandArgs.RData")
  else if(spatialRange == 50 && spatialVar == 0.3^2)
    stop("spatialRange == 50 && spatialVar == 0.3^2 not supported")
  else
    load("compareModelCommandArgsNew.RData")
  rangeID = ifelse(spatialRange == 50, "Range50", "")
  spatialVarID = ifelse(spatialVar == 0.3^2, "margVar0.09", "")
  scenarioID = ifelse(spatialRange == 150 && spatialVar == 0.15^2, "", paste0(rangeID, "_", spatialVarID))
  
  if(is.null(indices))
    indices = 1:length(compareModelCommandArgs)
  
  fullTableSRS1 = c() # constant risk
  fullTableSRS2 = c() # constant plus spatial effect
  fullTableSRS3 = c() # all but cluster effect
  fullTableSRS4 = c() # all effects
  fullTableBigSRS1 = c() # constant risk
  fullTableBigSRS2 = c() # constant plus spatial effect
  fullTableBigSRS3 = c() # all but cluster effect
  fullTableBigSRS4 = c() # all effects
  fullTable1 = c() # constant risk
  fullTable2 = c() # constant plus spatial effect
  fullTable3 = c() # all but cluster effect
  fullTable4 = c() # all effects
  fullTableBig1 = c() # constant risk
  fullTableBig2 = c() # constant plus spatial effect
  fullTableBig3 = c() # all but cluster effect
  fullTableBig4 = c() # all effects
  for(i in indices) {
    # get the arguments for the run, and specify that we want to load the precomputed results
    argList = compareModelCommandArgs[[i]]
    argList$loadResults = TRUE
    argList$strictPriors = strictPriors
    argList$printScoreTable = FALSE
    argList$printParTable = FALSE
    
    # get all elements from the list
    tausq = argList$tausq
    gamma = argList$gamma
    margVar = argList$margVar
    test = argList$test
    resultType = argList$resultType
    sampling = argList$sampling
    recomputeTruth = argList$recomputeTruth
    modelsI = argList$modelsI
    produceFigures = argList$produceFigures
    big = argList$big
    maxDataSets = argList$maxDataSets 
    nsim = argList$nsim
    argList$xtable.args=list(digits=c(0, 2, 2, 3, 2, 1, 2), display=rep("f", 7), auto=TRUE)
    argList$colDigits = c(2, 2, 3, 2, 1, 2)
    range = argList$effRange
    if(is.null(range))
      range = 150
    
    # skip nonexistent populations and populations not in our scenario
    if(range != spatialRange && margVar != 0)
      next
    if(range == 50 && margVar == 0)
      argList$effRange = 150
    if(margVar != 0 && margVar != spatialVar)
      next
    
    # generate an informative id string to label the table we are about to print with
    testText = ifelse(test, "Test", "")
    bigText = ifelse(big, "Big", "")
    strictPriorText = ifelse(strictPriors, "strictPrior", "")
    runId = paste0("Beta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                   "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                   "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets, 
                   scenarioID)
    print(runId)
    
    # get the precomputed scoring rule results
    out = do.call("runCompareModels2", argList)
    theseScores = out$tab
    
    # write scores to the score table for the relevant population
    if(big && sampling == "oversamp") {
      if(margVar == 0 && gamma == 0 && tausq == 0)
        fullTableBig1 = theseScores
      else if(gamma == 0 && tausq == 0)
        fullTableBig2 = theseScores
      else if(tausq == 0)
        fullTableBig3 = theseScores
      else
        fullTableBig4 = theseScores
    } else if(!big && sampling == "oversamp") {
      if(margVar == 0 && gamma == 0 && tausq == 0)
        fullTable1 = theseScores
      else if(gamma == 0 && tausq == 0)
        fullTable2 = theseScores
      else if(tausq == 0)
        fullTable3 = theseScores
      else
        fullTable4 = theseScores
    }
    else if(big && sampling == "SRS") {
      if(margVar == 0 && gamma == 0 && tausq == 0)
        fullTableBigSRS1 = theseScores
      else if(gamma == 0 && tausq == 0)
        fullTableBigSRS2 = theseScores
      else if(tausq == 0)
        fullTableBigSRS3 = theseScores
      else
        fullTableBigSRS4 = theseScores
    } else if(!big && sampling == "SRS") {
      if(margVar == 0 && gamma == 0 && tausq == 0)
        fullTableSRS1 = theseScores
      else if(gamma == 0 && tausq == 0)
        fullTableSRS2 = theseScores
      else if(tausq == 0)
        fullTableSRS3 = theseScores
      else
        fullTableSRS4 = theseScores
    }
  }
  
  # overwrite the scoring rules for the non-smoothing models with the big scoring rules
  modelNames = rownames(fullTable1)
  nonSmoothingI = (modelNames == "Naive") | (modelNames == "Direct")
  fullTable1[nonSmoothingI,] = fullTableBig1
  fullTable2[nonSmoothingI,] = fullTableBig2
  fullTable3[nonSmoothingI,] = fullTableBig3
  fullTable4[nonSmoothingI,] = fullTableBig4
  fullTableSRS1[nonSmoothingI,] = fullTableBigSRS1
  fullTableSRS2[nonSmoothingI,] = fullTableBigSRS2
  fullTableSRS3[nonSmoothingI,] = fullTableBigSRS3
  fullTableSRS4[nonSmoothingI,] = fullTableBigSRS4
  
  # filter out only the desired models
  if(incorrectlyAggregatedModels)
    filterRows = c(1, 2, 3, 4, 6, 10, 12, 7, 9, 13, 15, 16, 18)
  fullTable1 = fullTable1[filterRows,]
  fullTable2 = fullTable2[filterRows,]
  fullTable3 = fullTable3[filterRows,]
  fullTable4 = fullTable4[filterRows,]
  fullTableSRS1 = fullTableSRS1[filterRows,]
  fullTableSRS2 = fullTableSRS2[filterRows,]
  fullTableSRS3 = fullTableSRS3[filterRows,]
  fullTableSRS4 = fullTableSRS4[filterRows,]
  modelNames = rownames(fullTable1)
  
  # remove unnecessary symbols from model names
  modelNames = gsub("'", "", modelNames)
  if(!incorrectlyAggregatedModels)
    modelNames[grepl("BYM2", modelNames)] = gsub("a", "", modelNames[grepl("BYM2", modelNames)])
  
  # given the table on the type of model, this function prints out the type of population and 
  # design and a fancy latex version of the table
  makeTable = function(fullTab, sampling=c("Unstratified", "Stratified"), popI=1, 
                       colDigits=c(1, 1, 1, 1, 0, 2)) {
    # first print out label for the table so we know which one it is
    sampling = match.arg(sampling)
    
    if(popI == 1)
      popText = "suc"
    if(popI == 2)
      popText = "Suc"
    if(popI == 3)
      popText = "SUc"
    if(popI == 4)
      popText = "SUC"
    
    contextText = paste(popText, sampling)
    print(paste0("Printing table for ", contextText))
    
     ## now generate the table
    # get the model variation and the model type
    require(stringr)
    modelTypes = word(modelNames, 1)
    modelTypes[modelTypes == "Smoothed"] = "Smoothed Direct"
    uniqueModelTypes = unique(modelTypes)
    modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
    modelVariations = word(modelNames, 2)
    modelVariations[is.na(modelVariations)] = ""
    modelVariations[modelVariations == "Direct"] = ""
    
    # separate the text that will be diagonal, vertical, and horizontal
    diagonalText = c(modelTypes[1:3], rep("", length(modelTypes) - 3))
    variationText = modelVariations
    categoryText = c("BYM2", "SPDE")
    
    # round each column of the table to the correct number of digits
    temp = data.frame(fullTab)
    # for constant risk population, add another dps to var for smoothed direct and bym2
    temp2 = as.matrix(do.call("cbind", lapply(1:ncol(temp), function(x) {round(temp[,x], digits=colDigits[x])})))
    if(popI == 1) {
      temp2[modelTypes %in% c("Smoothed Direct", "BYM2"),2] = round(temp[modelTypes %in% c("Smoothed Direct", "BYM2"),2], digits=colDigits[2]+1)
    }
    
    # put the rounded numbers in our data frame
    rownames(temp2) = rownames(fullTab)
    colnames(temp2) = colnames(fullTab)
    fullTab = temp2
    
    ## Generate the tables
    require(stringr)
    require(dplyr)
    require(kableExtra)
    
    options(knitr.table.format = "latex")
    
    # bold the best entries of each column, italicize worst entries of each column
    centers = c(rep(0, 4), 80, 0)
    columnBest = apply(cbind(abs(fullTab[,1]), fullTab[,2:4], abs(fullTab[,5]-80), fullTab[,6]), 2, min)
    columnWorst = apply(cbind(abs(fullTab[,1]), fullTab[,2:4], abs(fullTab[,5]-80), fullTab[,6]), 2, max)
    dat = data.table(fullTab)
    test = dat %>% mutate(Bias = cell_spec(fullTab[,1], "latex", bold=abs(fullTab[,1] - centers[1]) <= columnBest[1], italic = abs(fullTab[,1] - centers[1]) >= columnWorst[1], 
                                           monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          Var = cell_spec(fullTab[,2], "latex", bold=abs(fullTab[,2] - centers[2]) <= columnBest[2], italic = abs(fullTab[,2] - centers[2]) >= columnWorst[2], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          MSE = cell_spec(fullTab[,3], "latex", bold=abs(fullTab[,3] - centers[3]) <= columnBest[3], italic = abs(fullTab[,3] - centers[3]) >= columnWorst[3], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          CRPS = cell_spec(fullTab[,4], "latex", bold=abs(fullTab[,4] - centers[4]) <= columnBest[4], italic = abs(fullTab[,4] - centers[4]) >= columnWorst[4], 
                                           monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          CVG = cell_spec(fullTab[,5], "latex", bold=abs(fullTab[,5] - centers[5]) <= columnBest[5], italic = abs(fullTab[,5] - centers[5]) >= columnWorst[5], 
                                          monospace=FALSE, underline=FALSE, strikeout=FALSE), 
                          Width = cell_spec(fullTab[,6], "latex", bold=abs(fullTab[,6] - centers[6]) <= columnBest[6], italic = abs(fullTab[,6] - centers[6]) >= columnWorst[6], 
                                            monospace=FALSE, underline=FALSE, strikeout=FALSE)) %>%
      select(Bias, Var, MSE, CRPS, CVG, Width)
    
    # revert the column names to their true values, set the model variations to be the values in the first column
    colnames(test) = colnames(fullTab)
    test = cbind(" "=modelVariations, test)
    rownames(test)=NULL
    
    # move the column units to be in their own row, include the scoring rule names as a header above the table
    scoringRules = gsub(" \\(.*\\)","",colnames(fullTab))
    scoringRules = gsub("\\%", '\\\\%', scoringRules)
    columnUnits = str_extract(colnames(fullTab), "\\(.*\\)")
    colnames(test) = c(" ", columnUnits)
    numberColumns = rep(1, 7)
    names(numberColumns) = c(" ", scoringRules)
    
    # group the rows by the type of model
    fullTab = test %>%
      kable("latex", escape = F, booktabs = T) %>% kable_styling()
    for(i in 1:length(uniqueModelTypes)) {
      startR = min(modelTypeGroups[[i]])
      endR = max(modelTypeGroups[[i]])
      fullTab = fullTab %>% pack_rows(uniqueModelTypes[i], startR, endR, escape=FALSE, bold=TRUE, italic=TRUE)
    }
    
    print(add_header_above(fullTab, numberColumns, italic=FALSE, bold=TRUE, escape=FALSE, line=FALSE))
  }
  
  # print all of the tables to the console
  browser()
  makeTable(fullTableSRS1, sampling="Unstratified", popI=1, colDigits=c(1, 1, 2, 1, 0, 2)) # MSE is very small in this case so go out an extra digit
  makeTable(fullTableSRS2, sampling="Unstratified", popI=2, colDigits=c(1, 1, 2, 1, 0, 2))
  makeTable(fullTableSRS3, sampling="Unstratified", popI=3, colDigits=c(1, 1, 2, 1, 0, 2))
  makeTable(fullTableSRS4, sampling="Unstratified", popI=4, colDigits=c(1, 1, 2, 1, 0, 2))
  makeTable(fullTable1, sampling="Stratified", popI=1, colDigits=c(1, 1, 2, 1, 0, 2))
  makeTable(fullTable2, sampling="Stratified", popI=2, colDigits=c(1, 1, 2, 1, 0, 2))
  makeTable(fullTable3, sampling="Stratified", popI=3, colDigits=c(1, 1, 2, 1, 0, 2))
  makeTable(fullTable4, sampling="Stratified", popI=4, colDigits=c(1, 1, 2, 1, 0, 2))
}

# plot the scoring rules for each analysis and population model for a fixed type of survey design (the survey design being SRS or stratified)
plotCompareModelsAllLocal = function(strictPriors=FALSE, usePrecomputedResults=FALSE, saveResults=FALSE, 
                                     spatialRange=150, spatialVar=0.15^2) {
  spatialRange = match.arg(as.character(spatialRange), choices=c(150, 50))
  spatialVar = match.arg(as.character(spatialVar), choices=c(0.15^2, 0.3^2))
  
  # map the population type to a type of point plotted:
  ## constant risk: 1
  ## constant plus spatial: 2
  ## constant plus spatial plus urban: 0
  ## all effects: 5
  pch = c(1, 2, 0, 5)
  # cols = rainbow(4)
  cols = c("red1", "purple", "blue1", "green4")
  # cols = qualitative_hcl(4, h1=247, h2=54, c1=80, l1=61) # these colors are colorblind friendly
  
  if(spatialRange == 150 && spatialVar == 0.15^2)
    load("compareModelCommandArgs.RData")
  else if(spatialRange == 50 && spatialVar == 0.3^2)
    stop("spatialRange == 50 && spatialVar == 0.3^2 not supported")
  else
    load("compareModelCommandArgsNew.RData")
  rangeID = ifelse(spatialRange == 50, "Range50", "")
  spatialVarID = ifelse(spatialVar == 0.3^2, "margVar0.09", "")
  scenarioID = ifelse(spatialRange == 150 && spatialVar == 0.15^2, "", paste0(rangeID, "_", spatialVarID))
  
  indices = 1:length(compareModelCommandArgs)
  
  plotHelper = function(scoreI, goalVal=NULL, rangeIncludes=c(), scoreName="", filterRows=c(1:3, 4, 6, 10, 12, 13, 15, 16, 18), 
                        shareRange=FALSE, plotSRSLegend=FALSE, plotDHSLegend=TRUE, logScale=FALSE) {
    plotNameRoot = paste0(tolower(scoreName), "Plot")
    
    if(!usePrecomputedResults) {
      fullTableSRS1 = c() # constant risk
      fullTableSRS2 = c() # constant plus spatial effect
      fullTableSRS3 = c() # all but cluster effect
      fullTableSRS4 = c() # all effects
      fullTableBigSRS1 = c() # constant risk
      fullTableBigSRS2 = c() # constant plus spatial effect
      fullTableBigSRS3 = c() # all but cluster effect
      fullTableBigSRS4 = c() # all effects
      fullTable1 = c() # constant risk
      fullTable2 = c() # constant plus spatial effect
      fullTable3 = c() # all but cluster effect
      fullTable4 = c() # all effects
      fullTableBig1 = c() # constant risk
      fullTableBig2 = c() # constant plus spatial effect
      fullTableBig3 = c() # all but cluster effect
      fullTableBig4 = c() # all effects
      for(i in indices) {
        # get the arguments for the run, and specify that we want to load the precomputed results
        argList = compareModelCommandArgs[[i]]
        argList$loadResults = TRUE
        argList$strictPriors = strictPriors
        
        # get all elements from the list
        tausq = argList$tausq
        gamma = argList$gamma
        margVar = argList$margVar
        test = argList$test
        resultType = argList$resultType
        sampling = argList$sampling
        recomputeTruth = argList$recomputeTruth
        modelsI = argList$modelsI
        produceFigures = argList$produceFigures
        big = argList$big
        maxDataSets = argList$maxDataSets 
        nsim = argList$nsim
        range = argList$effRange
        if(is.null(range))
          range = 150
        
        # skip nonexistent populations and populations not in our scenario
        if(range != spatialRange && margVar != 0)
          next
        if(range == 50 && margVar == 0)
          argList$effRange = 150
        if(margVar != 0 && margVar != spatialVar)
          next
        
        # generate an informative id string to label the table we are about to print with
        testText = ifelse(test, "Test", "")
        bigText = ifelse(big, "Big", "")
        strictPriorText = ifelse(strictPriors, "strictPrior", "")
        runId = paste0("Beta-1.75margVar", round(margVar, 4), "tausq", round(tausq, 4), "gamma", round(gamma, 4), 
                       "HHoldVar0urbanOverSamplefrac0", strictPriorText, testText, bigText, sampling, 
                       "models", do.call("paste0", as.list(modelsI)), "nsim", nsim, "MaxDataSetI", maxDataSets, 
                       scenarioID)
        print(runId)
        
        # get the precomputed scoring rule results
        out = do.call("runCompareModels2", argList)
        theseScores = out$unroundedTab
        
        # append scores to the score table for the relevant population
        if(big && sampling == "oversamp") {
          if(margVar == 0 && gamma == 0 && tausq == 0)
            fullTableBig1 = cbind(fullTableBig1, theseScores[,scoreI])
          else if(gamma == 0 && tausq == 0)
            fullTableBig2 = cbind(fullTableBig2, theseScores[,scoreI])
          else if(tausq == 0)
            fullTableBig3 = cbind(fullTableBig3, theseScores[,scoreI])
          else
            fullTableBig4 = cbind(fullTableBig4, theseScores[,scoreI])
        } else if(!big && sampling == "oversamp") {
          if(margVar == 0 && gamma == 0 && tausq == 0)
            fullTable1 = cbind(fullTable1, theseScores[,scoreI])
          else if(gamma == 0 && tausq == 0)
            fullTable2 = cbind(fullTable2, theseScores[,scoreI])
          else if(tausq == 0)
            fullTable3 = cbind(fullTable3, theseScores[,scoreI])
          else
            fullTable4 = cbind(fullTable4, theseScores[,scoreI])
        }
        else if(big && sampling == "SRS") {
          if(margVar == 0 && gamma == 0 && tausq == 0)
            fullTableBigSRS1 = cbind(fullTableBigSRS1, theseScores[,scoreI])
          else if(gamma == 0 && tausq == 0)
            fullTableBigSRS2 = cbind(fullTableBigSRS2, theseScores[,scoreI])
          else if(tausq == 0)
            fullTableBigSRS3 = cbind(fullTableBigSRS3, theseScores[,scoreI])
          else
            fullTableBigSRS4 = cbind(fullTableBigSRS4, theseScores[,scoreI])
        } else if(!big && sampling == "SRS") {
          if(margVar == 0 && gamma == 0 && tausq == 0)
            fullTableSRS1 = cbind(fullTableSRS1, theseScores[,scoreI])
          else if(gamma == 0 && tausq == 0)
            fullTableSRS2 = cbind(fullTableSRS2, theseScores[,scoreI])
          else if(tausq == 0)
            fullTableSRS3 = cbind(fullTableSRS3, theseScores[,scoreI])
          else
            fullTableSRS4 = cbind(fullTableSRS4, theseScores[,scoreI])
        }
      }
      
      # overwrite the scoring rules for the non-smoothing models with the big scoring rules
      modelNames = rownames(fullTable1)
      nonSmoothingI = (modelNames == "Naive") | (modelNames == "Direct")
      fullTable1[nonSmoothingI,] = fullTableBig1
      fullTable2[nonSmoothingI,] = fullTableBig2
      fullTable3[nonSmoothingI,] = fullTableBig3
      fullTable4[nonSmoothingI,] = fullTableBig4
      fullTableSRS1[nonSmoothingI,] = fullTableBigSRS1
      fullTableSRS2[nonSmoothingI,] = fullTableBigSRS2
      fullTableSRS3[nonSmoothingI,] = fullTableBigSRS3
      fullTableSRS4[nonSmoothingI,] = fullTableBigSRS4
      
      # get the scoring rule name
      scoringRuleName = colnames(theseScores)[scoreI]
    }
    
    # save results if necessary
    if(usePrecomputedResults) {
      load(paste0("compareModelPlotVar", scoreI, scenarioID, ".RData"))
    }
    if(saveResults) {
      save(fullTable1, fullTable2, fullTable3, fullTable4, 
           fullTableSRS1, fullTableSRS2, fullTableSRS3, fullTableSRS4, scoringRuleName, 
           file=paste0("compareModelPlotVar", scoreI, scenarioID, ".RData"))
    }
    
    # filter out only the desired models
    fullTable1 = fullTable1[filterRows,]
    fullTable2 = fullTable2[filterRows,]
    fullTable3 = fullTable3[filterRows,]
    fullTable4 = fullTable4[filterRows,]
    fullTableSRS1 = fullTableSRS1[filterRows,]
    fullTableSRS2 = fullTableSRS2[filterRows,]
    fullTableSRS3 = fullTableSRS3[filterRows,]
    fullTableSRS4 = fullTableSRS4[filterRows,]
    modelNames = names(fullTable1)
    
    # remove unnecessary symbols from model names
    modelNames = gsub("'", "", modelNames)
    modelNames[grepl("BYM2", modelNames)] = gsub("a", "", modelNames[grepl("BYM2", modelNames)])
    
    # get the model variation and the model type
    require(stringr)
    modelTypes = word(modelNames, 1)
    modelTypes[modelTypes == "Smoothed"] = "Smoothed Direct"
    uniqueModelTypes = unique(modelTypes)
    modelTypeGroups = lapply(uniqueModelTypes, function(x) {(1:length(modelTypes))[modelTypes == x]})
    modelVariations = word(modelNames, 2)
    modelVariations[is.na(modelVariations)] = ""
    modelVariations[modelVariations == "Direct"] = ""
    
    # separate the text that will be diagonal, vertical, and horizontal
    diagonalText = c(modelTypes[1:3], rep("", length(modelTypes) - 3))
    variationText = modelVariations
    categoryText = c("BYM2", "SPDE")
    
    ## Generate the plots
    if(!shareRange)
      scoreRange = range(c(fullTable1, fullTable2, fullTable3, fullTable4, rangeIncludes))
    else
      scoreRange = range(c(fullTable1, fullTable2, fullTable3, fullTable4, fullTableSRS1, fullTableSRS2, fullTableSRS3, fullTableSRS4, rangeIncludes))
    print(paste0(scoreName, " range: (", scoreRange[1], ", ", scoreRange[2], ")"))
    
    thisLog = ""
    if(logScale)
      thisLog = "y"
    
    # plot the urban oversampled values
    tempModelNames = sort(factor(modelNames, labels=modelNames))
    # centers = seq(from=1, to=16, by=1)
    delta = .75
    # centers = rev(c(1:4, (5:13) + 1 * delta, 14 + 2 * delta, (15:16) + 3 * delta) * (16 / (16 + 3 * delta)))
    centers = rev(c(1:4 + .5 * delta, (5:8) + 1.5 * delta, 9 + 2.5 * delta, (10:11) + 3.5 * delta) * (11 / (11 + 4 * delta))) # would need a different one for a different filterRows
    centers = centers[1] + centers[length(centers)] - centers
    unabbreviatedTitle = gsub('\\\\%', "\\%", scoringRuleName)
    unabbreviatedTitle = gsub('Var', "Variance", unabbreviatedTitle)
    unabbreviatedTitle = gsub('Cvg', "Coverage", unabbreviatedTitle)
    unabbreviatedTitle = gsub("\\(", "(stratified design, ", unabbreviatedTitle)
    strictText = ifelse(strictPriors, "strictPrior", "")
    
    # browser()
    
    pdf(paste0("figures/", plotNameRoot, strictText, "Stratified", scenarioID, ".pdf"), width=6, height=5)
    # par(mar=c(4.1, 8.1, 5.1, 5.3), xpd=TRUE)
    par(mar=c(6.1, 4.1, 3.1, 6.3), xpd=TRUE)
    stripchart(fullTable1 ~ tempModelNames, cex=0, las=2, ylim=scoreRange, main="", 
               at=rev(centers), ylab="", axes=FALSE, vertical=TRUE, log=thisLog)
    
    title(TeX(unabbreviatedTitle))
    box()
    # axis(3, las=2)
    axis(2, las=2)
    axis(1, srt=60, at=centers, labels=rep("", length(centers)))
    # text(par("usr")[1] - 1, centers, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yShift = diff(par("usr")[3:4]) / 15
    # text(centers +  delta, par("usr")[3] - yShift, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yloc = par("usr")[3] - yShift
    if(logScale)
      yloc = 10^yloc
    text(centers + max(centers)/30, yloc, labels = diagonalText, srt = 45, pos = 2, xpd = TRUE)
    yloc = par("usr")[3] - 1.25 * yShift
    if(logScale)
      yloc = 10^yloc
    text(centers + 1.5 * max(centers)/30, yloc, labels = variationText, pos = 2, xpd = TRUE)
    spdeCenter = mean(centers[(length(centers) - 3):length(centers)])
    bym2Center = mean(centers[(length(centers) - 7):(length(centers) - 4)])
    yloc = par("usr")[3] - yShift * 2.75
    if(logScale)
      yloc = 10^yloc
    text(c(bym2Center, spdeCenter) + 2.75 * max(centers)/35, yloc, labels = categoryText, pos = 2, xpd = TRUE)
    
    if(!is.null(goalVal)) {
      segments(y0=goalVal, x0=par("usr")[1], y1=goalVal, x1=par("usr")[2], lty=2, col="black")
    }
    
    stripchart(fullTable1 ~ tempModelNames, col=cols[1], pch=pch[1], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTable2 ~ tempModelNames, col=cols[2], pch=pch[2], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTable3 ~ tempModelNames, col=cols[3], pch=pch[3], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTable4 ~ tempModelNames, col=cols[4], pch=pch[4], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    
    if(plotDHSLegend) {
      # pos = legend("right", c("suc", "Suc", "SUc", "SUC"), pch=pch, col=cols, horiz=FALSE, inset=c(-0.225,0), 
      #              title="Population\nmodel", bty="n")
      pos = legend("right", c(expression("Pop"[suc]), expression("Pop"[Suc]), expression("Pop"[SUc]), expression("Pop"[SUC])), pch=pch, col=cols, horiz=FALSE, inset=c(-0.315,0), 
                   title="Population\nmodel", bty="n", lwd=2, lty=NA)
      xleft <- pos$rect[["left"]]
      ytop <- pos$rect[["top"]]
      ybottom <- ytop - pos$rect[["h"]]
      xright <- xleft + pos$rect[["w"]]
      rect(xleft, ybottom + 0.25 * yShift, xright, ytop-1.35 * yShift)
    }
    
    dev.off()
    
    # plot the SRS values
    pdf(paste0("figures/", plotNameRoot, strictText, "Unstratified", scenarioID, ".pdf"), width=6, height=5)
    # par(mar=c(4.1, 8.1, 5.1, 5.3), xpd=TRUE)
    par(mar=c(6.1, 4.1, 3.1, 6.3), xpd=TRUE)
    
    if(!shareRange)
      scoreRange = range(c(fullTableSRS1, fullTableSRS2, fullTableSRS3, fullTableSRS4, rangeIncludes))
    else
      scoreRange = range(c(fullTable1, fullTable2, fullTable3, fullTable4, fullTableSRS1, fullTableSRS2, fullTableSRS3, fullTableSRS4, rangeIncludes))
    stripchart(fullTable1 ~ tempModelNames, cex=0, las=2, ylim=scoreRange, main="", 
               at=rev(centers), ylab="", axes=FALSE, vertical=TRUE, log=thisLog)
    
    unabbreviatedTitle = gsub("\\(stratified design, ", "(unstratified design, ", unabbreviatedTitle)
    # title(TeX(paste0(unabbreviatedTitle)), line=4)
    title(TeX(unabbreviatedTitle))
    box()
    # axis(3, las=2)
    axis(2, las=2)
    axis(1, srt=60, at=centers, labels=rep("", length(centers)))
    # text(par("usr")[1] - 1, centers, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yShift = diff(par("usr")[3:4]) / 15
    # text(centers +  delta, par("usr")[3] - yShift, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yloc = par("usr")[3] - yShift
    if(logScale)
      yloc = 10^yloc
    text(centers + max(centers)/30, yloc, labels = diagonalText, srt = 45, pos = 2, xpd = TRUE)
    yloc = par("usr")[3] - 1.25 * yShift
    if(logScale)
      yloc = 10^yloc
    text(centers + 1.5 * max(centers)/30, yloc, labels = variationText, pos = 2, xpd = TRUE)
    spdeCenter = mean(centers[(length(centers) - 3):length(centers)])
    bym2Center = mean(centers[(length(centers) - 7):(length(centers) - 4)])
    yloc = par("usr")[3] - yShift * 2.75
    if(logScale)
      yloc = 10^yloc
    text(c(bym2Center, spdeCenter) + 2.75 * max(centers)/35, yloc, labels = categoryText, pos = 2, xpd = TRUE)
    
    if(!is.null(goalVal)) {
      segments(y0=goalVal, x0=par("usr")[1], y1=goalVal, x1=par("usr")[2], lty=2, col="black")
    }
    
    # stripchart(fullTable1 ~ tempModelNames, col=cols[1], pch=pch[1], add=TRUE, 
    #            at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTableSRS1 ~ tempModelNames, col=cols[1], pch=pch[1], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTableSRS2 ~ tempModelNames, col=cols[2], pch=pch[2], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTableSRS3 ~ tempModelNames, col=cols[3], pch=pch[3], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTableSRS4 ~ tempModelNames, col=cols[4], pch=pch[4], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    
    if(plotSRSLegend && scoreName!="Bias") {
      pos = legend("right", c(expression("Pop"[suc]), expression("Pop"[Suc]), expression("Pop"[SUc]), expression("Pop"[SUC])), pch=pch, col=cols, horiz=FALSE, inset=c(-0.315,0), 
                   title="Population\nmodel", bty="n", lwd=2, lty=NA)
      xleft <- pos$rect[["left"]]
      ytop <- pos$rect[["top"]]
      ybottom <- ytop - pos$rect[["h"]]
      xright <- xleft + pos$rect[["w"]]
      rect(xleft, ybottom + 0.25 * yShift, xright, ytop-1.35 * yShift) 
    }
    
    dev.off()
    
    ##### Now constructe simplified versions of the above plots
    simpleI = c(2, 3, 7, 11)
    
    # separate the text that will be diagonal, vertical, and horizontal
    # TO DO: fixed model type indices
    diagonalText = modelTypes[c(2:3, 4, 8)]
    variationText = modelVariations[c(2:3, 7, 11)]
    diagonalText[3] = paste(diagonalText[3], variationText[3])
    diagonalText[4] = paste(diagonalText[4], variationText[4])
    
    ## Generate the plots
    if(!shareRange)
      scoreRange = range(c(fullTable1[simpleI], fullTable2[simpleI], fullTable3[simpleI], fullTable4[simpleI], rangeIncludes))
    else
      scoreRange = range(c(fullTable1[simpleI], fullTable2[simpleI], fullTable3[simpleI], fullTable4[simpleI], fullTableSRS1[simpleI], fullTableSRS2[simpleI], fullTableSRS3[simpleI], fullTableSRS4[simpleI], rangeIncludes))
    
    thisLog = ""
    if(logScale)
      thisLog = "y"
    
    # plot the urban oversampled values
    tempModelNames = sort(factor(diagonalText, labels=diagonalText))
    diagonalText[3] = expression("BYM2"[UC])
    diagonalText[4] = expression("SPDE"[UC])
    # centers = seq(from=1, to=16, by=1)
    delta = .5
    # centers = rev(c(1:4, (5:13) + 1 * delta, 14 + 2 * delta, (15:16) + 3 * delta) * (16 / (16 + 3 * delta)))
    # centers = rev(c(1:4 + .5 * delta, (5:8) + 1.5 * delta, 9 + 2.5 * delta, (10:11) + 3.5 * delta) * (11 / (11 + 4 * delta))) # would need a different one for a different filterRows
    centers = rev(seq(1 + delta, 4 - delta, l=4)) # would need a different one for a different filterRows
    centers = centers[1] + centers[length(centers)] - centers
    unabbreviatedTitle = gsub('\\\\%', "\\%", scoringRuleName)
    unabbreviatedTitle = gsub('Var', "Variance", unabbreviatedTitle)
    unabbreviatedTitle = gsub('Cvg', "Coverage", unabbreviatedTitle)
    unabbreviatedTitle = gsub("\\(", "(stratified design, ", unabbreviatedTitle)
    strictText = ifelse(strictPriors, "strictPrior", "")
    
    # browser()
    
    pdf(paste0("figures/", plotNameRoot, strictText, "StratifiedSimple", scenarioID, ".pdf"), width=6, height=5)
    # par(mar=c(4.1, 8.1, 5.1, 5.3), xpd=TRUE)
    par(mar=c(6.1, 4.1, 3.1, 6.3), xpd=TRUE)
    stripchart(fullTable1[simpleI] ~ tempModelNames, cex=0, las=2, ylim=scoreRange, main="", 
               at=rev(centers), ylab="", axes=FALSE, vertical=TRUE, log=thisLog)
    
    title(TeX(unabbreviatedTitle))
    box()
    # axis(3, las=2)
    axis(2, las=2)
    axis(1, srt=60, at=centers, labels=rep("", length(centers)))
    # text(par("usr")[1] - 1, centers, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yShift = diff(par("usr")[3:4]) / 15
    # text(centers +  delta, par("usr")[3] - yShift, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yloc = par("usr")[3] - yShift
    if(logScale)
      yloc = 10^yloc
    text(centers + max(centers)/30, yloc, labels = diagonalText, srt = 45, pos = 2, xpd = TRUE)
    
    if(!is.null(goalVal)) {
      segments(y0=goalVal, x0=par("usr")[1], y1=goalVal, x1=par("usr")[2], lty=2, col="black")
    }
    
    stripchart(fullTable1[simpleI] ~ tempModelNames, col=cols[1], pch=pch[1], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTable2[simpleI] ~ tempModelNames, col=cols[2], pch=pch[2], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTable3[simpleI] ~ tempModelNames, col=cols[3], pch=pch[3], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTable4[simpleI] ~ tempModelNames, col=cols[4], pch=pch[4], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    
    if(plotDHSLegend) {
      # pos = legend("right", c("suc", "Suc", "SUc", "SUC"), pch=pch, col=cols, horiz=FALSE, inset=c(-0.225,0), 
      #              title="Population\nmodel", bty="n")
      pos = legend("right", c(expression("Pop"[suc]), expression("Pop"[Suc]), expression("Pop"[SUc]), expression("Pop"[SUC])), pch=pch, lwd=2, col=cols, horiz=FALSE, inset=c(-0.315,0), 
                   title="Population\nmodel", bty="n", lty=NA)
      xleft <- pos$rect[["left"]]
      ytop <- pos$rect[["top"]]
      ybottom <- ytop - pos$rect[["h"]]
      xright <- xleft + pos$rect[["w"]]
      rect(xleft, ybottom + 0.25 * yShift, xright, ytop-1.35 * yShift)
    }
    
    dev.off()
    
    # plot the SRS values
    pdf(paste0("figures/", plotNameRoot, strictText, "UnstratifiedSimple", scenarioID, ".pdf"), width=6, height=5)
    # par(mar=c(4.1, 8.1, 5.1, 5.3), xpd=TRUE)
    par(mar=c(6.1, 4.1, 3.1, 6.3), xpd=TRUE)
    
    if(!shareRange)
      scoreRange = range(c(fullTableSRS1[simpleI], fullTableSRS2[simpleI], fullTableSRS3[simpleI], fullTableSRS4[simpleI], rangeIncludes))
    else
      scoreRange = range(c(fullTable1[simpleI], fullTable2[simpleI], fullTable3[simpleI], fullTable4[simpleI], fullTableSRS1[simpleI], fullTableSRS2[simpleI], fullTableSRS3[simpleI], fullTableSRS4[simpleI], rangeIncludes))
    stripchart(fullTable1[simpleI] ~ tempModelNames, cex=0, las=2, ylim=scoreRange, main="", 
               at=rev(centers), ylab="", axes=FALSE, vertical=TRUE, log=thisLog)
    
    unabbreviatedTitle = gsub("\\(stratified design, ", "(unstratified design, ", unabbreviatedTitle)
    # title(TeX(paste0(unabbreviatedTitle)), line=4)
    title(TeX(unabbreviatedTitle))
    box()
    # axis(3, las=2)
    axis(2, las=2)
    axis(1, srt=60, at=centers, labels=rep("", length(centers)))
    # text(par("usr")[1] - 1, centers, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yShift = diff(par("usr")[3:4]) / 15
    # text(centers +  delta, par("usr")[3] - yShift, labels = tempModelNames, srt = 45, pos = 2, xpd = TRUE)
    yloc = par("usr")[3] - yShift
    if(logScale)
      yloc = 10^yloc
    text(centers + max(centers)/30, yloc, labels = diagonalText, srt = 45, pos = 2, xpd = TRUE)
    
    if(!is.null(goalVal)) {
      segments(y0=goalVal, x0=par("usr")[1], y1=goalVal, x1=par("usr")[2], lty=2, col="black")
    }
    
    stripchart(fullTable1[simpleI] ~ tempModelNames, col=cols[1], pch=pch[1], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTableSRS2[simpleI] ~ tempModelNames, col=cols[2], pch=pch[2], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTableSRS3[simpleI] ~ tempModelNames, col=cols[3], pch=pch[3], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    stripchart(fullTableSRS4[simpleI] ~ tempModelNames, col=cols[4], pch=pch[4], add=TRUE, 
               at=jitter(centers, amount = .15), vertical=TRUE, lwd=2)
    
    if(plotSRSLegend) {
      pos = legend("right", c(expression("Pop"[suc]), expression("Pop"[Suc]), expression("Pop"[SUc]), expression("Pop"[SUC])), pch=pch, col=cols, horiz=FALSE, inset=c(-0.315,0), 
                   title="Population\nmodel", bty="n", lwd=2, lty=NA)
      xleft <- pos$rect[["left"]]
      ytop <- pos$rect[["top"]]
      ybottom <- ytop - pos$rect[["h"]]
      xright <- xleft + pos$rect[["w"]]
      rect(xleft, ybottom + 0.25 * yShift, xright, ytop-1.35 * yShift) 
    }
    
    dev.off()
  }
  
  # generate plots for each scoring rule (1-6: bias, variance, mse, crps, coverage, width)
  plotHelper(3, goalVal=0, rangeIncludes=c(0.0488351441837245, 1.96136848911688), scoreName="MSE", plotDHSLegend=FALSE, shareRange=TRUE, logScale=TRUE)
  plotHelper(1, goalVal=0, rangeIncludes=c(0, -61.4285011791729, 36.9516649422178), scoreName="Bias", plotDHSLegend=TRUE, shareRange=TRUE)
  plotHelper(2, goalVal=0, rangeIncludes=c(0, 17.2193571578031), scoreName="Var", plotDHSLegend=FALSE, shareRange=TRUE)
  plotHelper(4, goalVal=0, rangeIncludes=c(0, 7.64137139664491), scoreName="CRPS", shareRange=TRUE, plotDHSLegend=FALSE)
  plotHelper(5, goalVal=80, rangeIncludes=c(55.1914893617021, 100), scoreName="Cvg", shareRange=TRUE, plotDHSLegend=TRUE)
  plotHelper(6, goalVal=0, rangeIncludes=c(0, 2.73626752901516), scoreName="Width", shareRange=TRUE, plotDHSLegend=FALSE)
}













