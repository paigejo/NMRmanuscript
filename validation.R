library(kableExtra)

# validate the smoothing models by leaving out data from one county at a time
validateExample = function(dat=ed, resultNameRoot="Ed", directEstResults=NULL, counties=sort(unique(poppc$County)), 
                           startFrom=0, loadPreviousSPDEFits=FALSE) {
  
  if(is.null(directEstResults)) {
    print("Loading Direct estimates results...")
    out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData"))
    # directEstResults should now be loaded into the environment
  }
  
  # no need to run the direct and naive examples, since they cannot make predictions for a county that does not have data
  
  # calculate the variance of the direct estimates assuming Gaussianity on the logit scale to get the validation weights
  directMoments = logitNormMoments(cbind(directEstResults$logit.est, directEstResults$var.est))
  directEstResults$varProbScale = directMoments[,2]
  
  ##### run Mercer et al. model
  if(startFrom <= 0) {
    print("Validating Mercer model...")
    mercerResults = validateMercer(dat=dat, directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale)

    save(mercerResults, file=paste0("resultsMercer", resultNameRoot, "ValidationAll.RData"))
  }
  else {
    print("Loading Smoothed Direct results...")
    load(file=paste0("resultsMercer", resultNameRoot, "ValidationAll.RData"))
  }
  
  ##### run BYM models
  source("designBased.R")
  if(startFrom <= 1) {
    print("Validating BYM2 I model...")
    bym2I = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                            dat, saveResults=TRUE, includeUrbanRural = FALSE, includeCluster = FALSE, 
                            fileNameRoot=resultNameRoot)$bym2Results
  }
  else {
    print("Loading BYM2 I results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', FALSE, 'Cluster', FALSE, '.RData'))
    bym2I = bym2Results
  }
  
  if(startFrom <= 2) {
    print("Validating BYM2 II model...")
    bym2II = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                             dat, saveResults=TRUE, includeUrbanRural = FALSE, includeCluster = TRUE, 
                             fileNameRoot=resultNameRoot)
  }
  else {
    print("Loading BYM2 II results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', FALSE, 'Cluster', TRUE, '.RData'))
    bym2II = bym2Results
  }
  
  if(startFrom <= 3) {
    print("Validating BYM2 III model...")
    bym2III = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                              dat, saveResults=TRUE, includeUrbanRural = TRUE, includeCluster = FALSE, 
                              fileNameRoot=resultNameRoot)
  }
  else {
    print("Loading BYM2 III results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', TRUE, 'Cluster', FALSE, '.RData'))
    bym2III = bym2Results
  }
  
  if(startFrom <= 4) {
    print("Validating BYM2 IV model...")
    bym2IV = validateBYM2Dat(directEstResults$logit.est, directEstResults$var.est, directEstResults$varProbScale, 
                             dat, saveResults=TRUE, includeUrbanRural = TRUE, includeCluster = TRUE, 
                             fileNameRoot=resultNameRoot)
  }
  else {
    print("Loading BYM2 IV results...")
    load(paste0('bym2', resultNameRoot, 'ValidationAllUrbRur', TRUE, 'Cluster', TRUE, '.RData'))
    bym2IV = bym2Results
  }
  
  ##### run SPDE models
  if(resultNameRoot == "Ed")
    targetPop = "women"
  else
    targetPop = "children"
  if(startFrom <= 5) {
    print("Validating SPDE I model...")
    spdeI = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                            directVars=directEstResults$varProbScale, includeClustEffect = FALSE, urbanEffect = FALSE, saveResults=TRUE, 
                            loadPreviousFit=loadPreviousSPDEFits, fileNameRoot=resultNameRoot, targetPop=targetPop)
  }
  else {
    print("Loading SPDE I results...")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", FALSE, 
                      "_urbanEffect", FALSE, ".RData")
    load(fileName)
    spdeI = spdeResults
  }
  
  if(startFrom <= 6) {
    print("Validating SPDE II model...")
    spdeII = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                             directVars=directEstResults$varProbScale, includeClustEffect = TRUE, urbanEffect = FALSE, saveResults = TRUE, 
                             loadPreviousFit=loadPreviousSPDEFits, fileNameRoot=resultNameRoot, targetPop=targetPop)
  }
  else {
    print("Loading SPDE II results...")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", TRUE, 
                      "_urbanEffect", FALSE, ".RData")
    load(fileName)
    spdeII = spdeResults
  }
  
  if(startFrom <= 7) {
    print("Validating SPDE III model...")
    spdeIII = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                              directVars=directEstResults$varProbScale, includeClustEffect = FALSE, urbanEffect = TRUE, saveResults = TRUE, 
                              loadPreviousFit=loadPreviousSPDEFits, fileNameRoot=resultNameRoot, targetPop=targetPop)
  }
  else {
    print("Loading SPDE III results...")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", FALSE, 
                      "_urbanEffect", TRUE, ".RData")
    load(fileName)
    spdeIII = spdeResults
  }
  
  if(startFrom <= 8) {
    print("Validating SPDE IV model...")
    spdeIV = validateSPDEDat(clustDat = dat, directLogitEsts=directEstResults$logit.est, directLogitVars=directEstResults$var.est, 
                             directVars=directEstResults$varProbScale, includeClustEffect = TRUE, urbanEffect = TRUE, saveResults = TRUE, 
                             loadPreviousFit=loadPreviousSPDEFits, fileNameRoot=resultNameRoot, targetPop=targetPop)
  }
  else {
    print("Loading SPDE IV results...")
    # fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", TRUE, 
    #                   "_urbanEffect", TRUE, ".RData")
    fileName = paste0("resultsSPDE", resultNameRoot, "ValidationAll", "_includeClustEffect", TRUE, 
                      "_urbanEffect", TRUE, "compact.RData")
    load(fileName)
    spdeIV = spdeResults
  }
  
  ##### aggregate results
  print("Aggregating and saving results...")
  filterI = c(-10, -14)
  tabInSample = rbind(mercerResults$mercerResultsInSample$scores, 
                      bym2I$bym2ResultsInSample$scores, 
                      bym2II$bym2ResultsInSample$scores, 
                      bym2III$bym2ResultsInSample$scores, 
                      bym2IV$bym2ResultsInSample$scores, 
                      spdeI$spdeResultsInSample$scores, 
                      spdeII$spdeResultsInSample$scores, 
                      spdeIII$spdeResultsInSample$scores, 
                      spdeIV$spdeResultsInSample$scores)
  tabInSample=tabInSample[,filterI]
  sortI = c(1:3, 7, 10, 4, 8, 11, 5, 9, 12, 6)
  tabInSample=tabInSample[,sortI]
  rownames(tabInSample) = c("Smoothed Direct", "BYM2 I", "BYM2 II", "BYM2 III", "BYM2 IV", "SPDE I", "SPDE II", "SPDE III", "SPDE IV")
  PITsInSample = list(mercer=mercerResults$mercerResultsInSample$pit, 
                      bym2I=bym2I$bym2ResultsInSample$pit, 
                      bym2II=bym2II$bym2ResultsInSample$pit, 
                      bym2III=bym2III$bym2ResultsInSample$pit, 
                      bym2IV=bym2IV$bym2ResultsInSample$pit, 
                      spdeI=spdeI$spdeResultsInSample$pit, 
                      spdeII=spdeII$spdeResultsInSample$pit, 
                      spdeIII=spdeIII$spdeResultsInSample$pit, 
                      spdeIV=spdeIV$spdeResultsInSample$pit)
  
  # filterI = c(1:3, 4:7, 9:11)
  # tabLeaveOutCounty = rbind(mercerResults$mercerResultsLeaveOutCounty$scores, 
  #                           bym2I$bym2ResultsLeaveOutCounty$scores, 
  #                           bym2II$bym2ResultsLeaveOutCounty$scores, 
  #                           bym2II$bym2ResultsModLeaveOutCounty$scores, 
  #                           bym2III$bym2ResultsLeaveOutCounty$scores, 
  #                           bym2IV$bym2ResultsLeaveOutCounty$scores, 
  #                           bym2IV$bym2ResultsModLeaveOutCounty$scores, 
  #                           spdeI$spdeResultsLeaveOutCounty$scores, 
  #                           spdeII$spdeResultsLeaveOutCounty$scores, 
  #                           spdeIII$spdeResultsLeaveOutCounty$scores, 
  #                           spdeIV$spdeResultsLeaveOutCounty$scores)
  # rownames(tabLeaveOutCounty) = c("Smoothed Direct", "BYM2 I", "BYM2 II", "BYM2 III", "BYM2 IV", "SPDE I", "SPDE II", "SPDE III", "SPDE IV")
  PITsLeaveOutCounty = list(mercer=mercerResults$mercerResultsLeaveOutCounty$pit, 
                            bym2I=bym2I$bym2ResultsLeaveOutCounty$pit, 
                            bym2II=bym2II$bym2ResultsLeaveOutCounty$pit, 
                            bym2III=bym2III$bym2ResultsLeaveOutCounty$pit, 
                            bym2IV=bym2IV$bym2ResultsLeaveOutCounty$pit, 
                            spdeI=spdeI$spdeResultsLeaveOutCounty$pit, 
                            spdeII=spdeII$spdeResultsLeaveOutCounty$pit, 
                            spdeIII=spdeIII$spdeResultsLeaveOutCounty$pit, 
                            spdeIV=spdeIV$spdeResultsLeaveOutCounty$pit)
  
  tabLeaveOutCluster = c(mean(bym2I$bym2ResultsLeaveOutCluster$CPO, na.rm = TRUE),  # the mercer model isn't fit to clusters, so it cannot be included here
                         mean(bym2II$bym2ResultsLeaveOutCluster$CPO, na.rm = TRUE), 
                         mean(bym2III$bym2ResultsLeaveOutCluster$CPO, na.rm = TRUE), 
                         mean(bym2IV$bym2ResultsLeaveOutCluster$CPO, na.rm = TRUE), 
                         mean(spdeI$spdeResultsLeaveOutCluster$CPO, na.rm = TRUE), 
                         mean(spdeII$spdeResultsLeaveOutCluster$CPO, na.rm = TRUE), 
                         mean(spdeIII$spdeResultsLeaveOutCluster$CPO, na.rm = TRUE), 
                         mean(spdeIV$spdeResultsLeaveOutCluster$CPO, na.rm = TRUE))
  tabLeaveOutCluster = data.frame(CPO=tabLeaveOutCluster)
  rownames(tabLeaveOutCluster) = c("BYM2 I", "BYM2 II", "BYM2 III", "BYM2 IV", "SPDE I", "SPDE II", "SPDE III", "SPDE IV")
  PITsLeaveOutCluster = list(bym2I=bym2I$bym2ResultsLeaveOutCluster$pit, 
                             bym2II=bym2II$bym2ResultsLeaveOutCluster$pit, 
                             bym2III=bym2III$bym2ResultsLeaveOutCluster$pit, 
                             bym2IV=bym2IV$bym2ResultsLeaveOutCluster$pit, 
                             spdeI=spdeI$spdeResultsLeaveOutCluster$pit, 
                             spdeII=spdeII$spdeResultsLeaveOutCluster$pit, 
                             spdeIII=spdeIII$spdeResultsLeaveOutCluster$pit, 
                             spdeIV=spdeIV$spdeResultsLeaveOutCluster$pit)
  
  # also generate the full, unfiltered result tables with results broken down by urban/rural
  tabInSampleAll = rbind(mercerResults$mercerResultsInSample$allResults, 
                         bym2I$bym2ResultsInSample$allResults, 
                         bym2II$bym2ResultsInSample$allResults, 
                         bym2III$bym2ResultsInSample$allResults, 
                         bym2IV$bym2ResultsInSample$allResults, 
                         spdeI$spdeResultsInSample$allResults, 
                         spdeII$spdeResultsInSample$allResults, 
                         spdeIII$spdeResultsInSample$allResults, 
                         spdeIV$spdeResultsInSample$allResults)
  rownames(tabInSampleAll) = c("Smoothed Direct", "BYM2 I", "BYM2 II", "BYM2 III", "BYM2 IV", "SPDE I", "SPDE II", "SPDE III", "SPDE IV")
  tabLeaveOutCountyAll = rbind(mercerResults$mercerResultsLeaveOutCounty$allResults, 
                               bym2I$bym2ResultsLeaveOutCounty$allResults, 
                               bym2II$bym2ResultsLeaveOutCounty$allResults, 
                               bym2III$bym2ResultsLeaveOutCounty$allResults, 
                               bym2IV$bym2ResultsLeaveOutCounty$allResults, 
                               spdeI$spdeResultsLeaveOutCounty$allResults, 
                               spdeII$spdeResultsLeaveOutCounty$allResults, 
                               spdeIII$spdeResultsLeaveOutCounty$allResults, 
                               spdeIV$spdeResultsLeaveOutCounty$allResults)
  rownames(tabLeaveOutCountyAll) = c("Smoothed Direct", "BYM2 I", "BYM2 II", "BYM2 III", "BYM2 IV", "SPDE I", "SPDE II", "SPDE III", "SPDE IV")
  filterI = c(1:5, 7:9, 13:15)
  tabLeaveOutCounty = tabLeaveOutCountyAll[,filterI]
  sortI = c(1, 6, 9, 2, 7, 10, 3, 8, 11, 4:5)
  tabLeaveOutCounty = tabLeaveOutCounty[,sortI]
  tabLeaveOutClusterAll = tabLeaveOutCluster
  
  # save and return results
  validationResults = list(scoresInSample=tabInSample, pitInSample=PITsInSample, 
                           scoresLeaveOutCounty=tabLeaveOutCounty, pitLeaveOutCounty=PITsLeaveOutCounty, 
                           scoresLeaveOutCluster=tabLeaveOutCluster, pitLeaveOutCluster=PITsLeaveOutCluster, 
                           allResultsInSample=tabInSampleAll, 
                           allResultsLeaveOutCounty=tabLeaveOutCountyAll, 
                           allResultsLeaveOutCluster=tabLeaveOutClusterAll)
  save(validationResults, file=paste0("resultsValidationAll", resultNameRoot, ".RData"))
  invisible(validationResults)
}

# print out validation results and plot the PITs
printValidationResults = function(resultNameRoot="Ed") {
  require(stringr)
  require(dplyr)
  require(kableExtra)
  
  # first load the validation results
  out = load(paste0("resultsValidationAll", resultNameRoot, ".RData"))
  modelNames = rownames(validationResults$scoresInSample)
  
  # change the modelNames to the adjusted ones:
  badNames = c("IV", "III", "II", "I")
  goodNames = c("UC", "Uc", "uC", "uc")
  for(i in 1:length(badNames)) {
    badName = badNames[i]
    goodName = goodNames[i]
    whichI = grepl(badName, modelNames)
    modelNames[whichI] = gsub(badName, goodName, modelNames[whichI])
  }
  rownames(validationResults$scoresInSample) = modelNames
  rownames(validationResults$scoresLeaveOutCounty) = modelNames
  rownames(validationResults$scoresLeaveOutCluster) = modelNames[-1]
  
  modelCodes = gsub(" ", "", modelNames, fixed = TRUE)
  modelCodes[3] = "BYM2uClust"
  modelCodes[4] = "BYM2uClust'"
  modelCodes[5] = "BYM2Urbc"
  modelCodes[6] = "BYM2UrbClust"
  modelCodes[7] = "BYM2UrbClust'"
  modelCodes[9] = "SPDEuClust"
  modelCodes[10] = "SPDEUrbc"
  modelCodes[11] = "SPDEUrbClust"
  ## plot the PITs
  for(i in 1:length(modelNames)) {
    thisModelCode = modelCodes[i]
    pdf(paste0("figures/validationPITHistogram", resultNameRoot, "_", thisModelCode, ".pdf"), width=9, height=5)
    par(mfrow=c(1,2))
    hist(unlist(validationResults$pitInSample[[i]]), main=paste0(modelNames[i], " In Sample PITs"), 
         xlab="PIT", freq=FALSE, breaks=30, col="skyblue")
    
    hist(unlist(validationResults$pitLeaveOutCounty[[i]]), main=paste0(modelNames[i], " Leave Out County PITs"), 
         xlab="PIT", freq=FALSE, breaks=30, col="skyblue")
    dev.off()
    
    # only CPO is calculated when just leaving out individual clusters
    # if(i != 1) {
    #   hist(unlist(validationResults$pitLeaveOutCluster[[i]]), main=paste0("Histogram of ", modelNames[i], " Left Out Cluster PITs"), 
    #        xlab="PIT", freq=FALSE, breaks=30, col="skyblue")
    # }
  }
  
  ## modify the BYM2 model names 
  tab = validationResults$scoresInSample
  
  ## print out the scoring rule tables
  print("In sample results:")
  displayRow = "s"
  displayIC = rep("f", 2)
  displayMSE = rep("f", 3)
  displayVar = rep("f", 3)
  displayBias = rep("f", 3)
  displayCRPS = "f"
  display = c(displayRow, displayIC, displayMSE, displayVar, displayBias, displayCRPS)
  scalingPower = c(MSE=2, Var=3, Bias=3)
  if(resultNameRoot == "Mort") {
    scalingPower = c(MSE=4, Var=4, Bias=4)
  }
  colScaleIC = rep(1, 2)
  colScaleMSE = rep(10^scalingPower[1], 3)
  colScaleVar = rep(10^scalingPower[2], 3)
  colScaleBias = rep(10^scalingPower[3], 3)
  colScaleCRPS = 1
  colScale = c(colScaleIC, colScaleMSE, colScaleVar, colScaleBias, colScaleCRPS)
  colUnitsIC = rep("", 2)
  colUnitsMSE = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[1]), "}$)"), 3)
  colUnitsVar = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[2]), "}$)"), 3)
  colUnitsBias = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[3]), "}$)"), 3)
  colUnitsCRPS = ""
  colUnits = c(colUnitsIC, colUnitsMSE, colUnitsVar, colUnitsBias, colUnitsCRPS)
  digitsIC = rep(0, 2)
  digitsMSE = rep(1, 3)
  digitsVar = rep(1, 3)
  digitsBias = rep(1, 3)
  digitsCRPS = 2
  colDigits = c(digitsIC, digitsMSE, digitsVar, digitsBias, digitsCRPS)
  
  tab = validationResults$scoresInSample
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(tab[,i] * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
  }
  colnames(tab) = gsub(".", " ", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Urban", "Urb", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Rural", "Rur", colnames(tab), fixed=TRUE)
  tab = tab[,-1] # filter out WAIC since it doesn't work for spatially correlated data
  colDigits=colDigits[-1]
  display=display[-1]
  colUnits = colUnits[-1]
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  tab = t(tab)
  colnames(tab) = gsub("SPDE ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("BYM2 ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Smoothed Direct", " ", colnames(tab), fixed=TRUE)
  tab = tab[c(2:11, 1), ]
  colUnits = colUnits[c(2:11, 1)]
  # temp = xtable2kable(xtable(tab), 
  #                     include.colnames=TRUE,
  #                     hline.after=0, 
  #                     math.style.exponents=TRUE, 
  #                     sanitize.text.function=function(x){x})
  # print(add_header_above(kable_styling(temp), c(" " = 1, "Smoothed Direct" = 1, "BYM2"=4, "SPDE"=4)))
  
  # remove the smoothed direct results
  tab = tab[,-1]
  
  # bold the best entries of each row, italicize worst entries of each row
  centers = c(rep(0, nrow(tab)))
  rowBest = apply(abs(sweep(tab, 1, centers, "-")), 1, min, na.rm=TRUE)
  rowWorst = apply(abs(sweep(tab, 1, centers, "-")), 1, max, na.rm=TRUE)
  boldFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) <= rowBest[i]
    out
  }
  italicFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) >= rowWorst[i]
    out
  }
  test = t(sapply(1:nrow(tab), function(i) {cell_spec(tab[i,], "latex", bold=boldFun(i), italic=italicFun(i), 
                                             monospace=FALSE, underline=FALSE, strikeout=FALSE)}))
  
  # revert the column names to their true values
  colnames(test) = colnames(tab)
  scoreVariations = c("Avg", "Urban", "Rural")
  scoreVariations = c(rep(scoreVariations, 3), "Avg", "Avg")
  test = cbind(" "=scoreVariations, test)
  rownames(test)=NULL
  
  # group the rows by urbanicity
  fullTab = test %>%
    kable("latex", escape = F, booktabs = T, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
          align=c("l", rep("r", ncol(test) - 1))) %>% kable_styling()
  uniqueScoreTypes = c("MSE", "Var", "Bias", "CRPS", "DIC")
  uniqueColUnits = colUnits[c(1, 4, 7, 10:11)]
  for(i in 1:length(uniqueScoreTypes)) {
    uniqueScoreTypes[i] = paste0(uniqueScoreTypes[i], uniqueColUnits[i])
  }
  scoreTypeGroups = c(rep(uniqueScoreTypes[1:3], each=3), uniqueScoreTypes[4:5])
  
  for(i in 1:length(uniqueScoreTypes)) {
    startR = match(TRUE, scoreTypeGroups == uniqueScoreTypes[i])
    endR = nrow(tab) - match(TRUE, rev(scoreTypeGroups == uniqueScoreTypes[i])) + 1
    fullTab = fullTab %>% pack_rows(uniqueScoreTypes[i], startR, endR, latex_gap_space = "0.3em", escape=FALSE)
  }
  print(add_header_above(fullTab, c(" " = 1, "BYM2"=4, "SPDE"=4), italic=TRUE, bold=TRUE, escape=FALSE))
  
  # print(xtable(tab, digits=digits, display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  
  print("")
  print("Leave out county results:")
  
  scalingPower = c(MSE=2, Var=3, Bias=3)
  if(resultNameRoot == "Mort") {
    scalingPower = c(MSE=4, Var=4, Bias=4)
  }
  displayRow = "s"
  displayMSE = rep("f", 3)
  displayVar = rep("f", 3)
  displayBias = rep("f", 3)
  displayCPO = "f"
  displayCRPS = "f"
  display = c(displayRow, displayMSE, displayVar, displayBias, displayCPO, displayCRPS)
  colScaleMSE = rep(10^scalingPower[1], 3)
  colScaleVar = rep(10^scalingPower[2], 3)
  colScaleBias = rep(10^scalingPower[3], 3)
  colScaleCPO = 1
  colScaleCRPS = 1
  colScale = c(colScaleMSE, colScaleVar, colScaleBias, colScaleCPO, colScaleCRPS)
  colUnitsMSE = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[1]), "}$)"), 3)
  colUnitsVar = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[2]), "}$)"), 3)
  colUnitsBias = rep(paste0(" ($\\times 10^{-", as.character(scalingPower[3]), "}$)"), 3)
  colUnitsCPO = ""
  colUnitsCRPS = ""
  colUnits = c(colUnitsMSE, colUnitsVar, colUnitsBias, colUnitsCPO, colUnitsCRPS)
  digitsMSE = rep(1, 3)
  digitsVar = rep(1, 3)
  digitsBias = rep(1, 3)
  digitsCPO = 2
  digitsCRPS = 2
  if(resultNameRoot == "Mort") {
    digitsCPO = 3
    digitsCRPS = 3
  }
  colDigits = c(digitsMSE, digitsVar, digitsBias, digitsCPO, digitsCRPS)
  
  tab = data.frame(validationResults$scoresLeaveOutCounty)
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(unlist(tab[,i]) * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
  }
  colnames(tab) = gsub(".", " ", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Urban", "Urb", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Rural", "Rur", colnames(tab), fixed=TRUE)
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  tab = t(tab)
  colnames(tab) = gsub("SPDE ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("BYM2 ", "", colnames(tab), fixed=TRUE)
  colnames(tab) = gsub("Smoothed Direct", " ", colnames(tab), fixed=TRUE)
  # temp = xtable2kable(xtable(tab), 
  #                     include.colnames=TRUE,
  #                     hline.after=0, 
  #                     math.style.exponents=TRUE, 
  #                     sanitize.text.function=function(x){x})
  # print(add_header_above(kable_styling(temp), c(" " = 1, "Smoothed Direct" = 1, "BYM2"=4, "SPDE"=4)))
  
  # remove the smoothed direct results
  tab = tab[,-1]
  
  # bold the best entries of each row, italicize worst entries of each row
  centers = c(rep(0, nrow(tab)))
  rowWorst = apply(abs(sweep(tab, 1, centers, "-")), 1, max, na.rm=TRUE)
  rowBest = apply(abs(sweep(tab, 1, centers, "-")), 1, min, na.rm=TRUE)
  temp = rowWorst[10]
  rowWorst[10] = rowBest[10] # highest CPO is best, not the lowest
  rowBest[10] = temp
  boldFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowBest[i]
    out
  }
  italicFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowWorst[i]
    out
  }
  if(resultNameRoot == "Ed") {
    # bold the best model for each scoring rule, italicize the worst
    test = t(sapply(1:nrow(tab), function(i) {cell_spec(tab[i,], "latex", bold=boldFun(i), italic=italicFun(i), 
                                                        monospace=FALSE, underline=FALSE, strikeout=FALSE)}))
  }
  else {
    # don't bother bolding or italicizing, since all models perform essentially equally
    test = tab
  }
  
  # revert the column names to their true values
  colnames(test) = colnames(tab)
  scoreVariations = c("Avg", "Urban", "Rural")
  scoreVariations = c(rep(scoreVariations, 3), "Avg", "Avg")
  test = cbind(" "=scoreVariations, test)
  rownames(test)=NULL
  
  # group the rows by urbanicity
  fullTab = test %>%
    kable("latex", escape = F, booktabs = T, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
          align=c("l", rep("r", ncol(test) - 1))) %>% kable_styling()
  uniqueScoreTypes = c("MSE", "Var", "Bias", "CPO", "CRPS")
  uniqueColUnits = colUnits[c(1, 4, 7, 10:11)]
  for(i in 1:length(uniqueScoreTypes)) {
    uniqueScoreTypes[i] = paste0(uniqueScoreTypes[i], uniqueColUnits[i])
  }
  scoreTypeGroups = c(rep(uniqueScoreTypes[1:3], each=3), uniqueScoreTypes[4:5])
  
  for(i in 1:length(uniqueScoreTypes)) {
    startR = match(TRUE, scoreTypeGroups == uniqueScoreTypes[i])
    endR = nrow(tab) - match(TRUE, rev(scoreTypeGroups == uniqueScoreTypes[i])) + 1
    fullTab = fullTab %>% pack_rows(uniqueScoreTypes[i], startR, endR, latex_gap_space = "0.3em", escape=FALSE)
  }
  print(add_header_above(fullTab, c(" " = 1, "BYM2"=4, "SPDE"=4), italic=TRUE, bold=TRUE, escape=FALSE))
  
  print("")
  print("Leave out cluster results:")
  colScale = colScaleCPO
  colDigits = digitsCPO
  colUnits = colUnitsCPO
  display = c(displayRow, displayCPO)
  tab = data.frame(validationResults$scoresLeaveOutCluster)
  for(i in 1:ncol(tab)) {
    tab[,i] = as.numeric(round(unlist(tab[,i]) * colScale[i], digits=colDigits[i]))
    colnames(tab)[i] = paste0(colnames(tab)[i], colUnits[i])
  }
  # print(xtable(tab, digits=c(1,colDigits), display=display), 
  #       include.colnames=TRUE,
  #       hline.after=0, 
  #       math.style.exponents=TRUE, 
  #       sanitize.text.function=function(x){x}, 
  #       scalebox=0.5)
  
  tab = t(tab)
  
  # bold the best entries of each row, italicize worst entries of each row
  centers = c(rep(0, nrow(tab)))
  rowWorst = apply(abs(sweep(tab, 1, centers, "-")), 1, min, na.rm=TRUE)
  rowBest = apply(abs(sweep(tab, 1, centers, "-")), 1, max, na.rm=TRUE)
  boldFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowBest[i]
    out
  }
  italicFun = function(i) {
    vals = tab[i,]
    isNA = is.na(vals)
    out = rep(FALSE, length(vals))
    out[!isNA] = abs(tab[i,!isNA] - centers[i]) == rowWorst[i]
    out
  }
  test = t(sapply(1:nrow(tab), function(i) {cell_spec(tab[i,], "latex", bold=boldFun(i), italic=italicFun(i), 
                                                      monospace=FALSE, underline=FALSE, strikeout=FALSE)}))
  
  # revert the column names to their true values
  colnames(test) = word(colnames(tab), 2)
  rownames(test)= "CPO"
  
  # construct the kable version of the table
  fullTab = test %>%
    kable("latex", escape = F, booktabs = T, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
          align=c("l", rep("r", ncol(test) - 1))) %>% kable_styling()
  
  # group the rows by urbanicity
  print(add_header_above(fullTab, c(" "=1, "BYM2"=4, "SPDE"=4), italic=TRUE, bold=TRUE, escape=FALSE))
}





