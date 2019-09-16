library(colorspace)
library(mgcv)

# makes plots as well as parameter estimate tables for the example applications in the manuscript
makeAllPlots = function(dat=ed, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                        meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                        sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                        meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                        meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                        varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                        sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                        ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                        plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                        makeScreenSplitPlot=FALSE, sharedPredictionScale=TRUE, widthRange, widthTicks, 
                        widthTickLabels, widthRangeSPDE, widthTicksSPDE, widthTickLabelsSPDE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  plotDataVisualizations(dat, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                         meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                         sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                         meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                         meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                         varName, plotNameRoot, resultNameRoot, meanCols, sdCols, popCols, ncols, 
                         relativeCols, urbCols, plotUrbanMap, kenyaLatRange, kenyaLonRange)
  
  plotModelPredictions(dat, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                       meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                       sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                       meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                       meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                       varName, plotNameRoot, resultNameRoot, meanCols, sdCols, popCols, ncols, 
                       relativeCols, urbCols, plotUrbanMap, kenyaLatRange, kenyaLonRange, widthRange, widthTicks, 
                       widthTickLabels, widthRangeSPDE, widthTicksSPDE, widthTickLabelsSPDE)
  
  plotModelPredictionsPart2(dat, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                            meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                            sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                            meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                            meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                            varName, plotNameRoot, resultNameRoot, meanCols, sdCols, popCols, ncols, 
                            relativeCols, urbCols, plotUrbanMap, kenyaLatRange, kenyaLonRange, 
                            sharedPredictionScale=sharedPredictionScale, widthRange, widthTicks, 
                            widthTickLabels, widthRangeSPDE, widthTicksSPDE, widthTickLabelsSPDE)
  
  makePairPlots(dat, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                varName, plotNameRoot, resultNameRoot, meanCols, sdCols, popCols, ncols, 
                relativeCols, urbCols, plotUrbanMap, kenyaLatRange, kenyaLonRange, 
                makeScreenSplitPlot=makeScreenSplitPlot)
  
  ##### now print the parameter estimates:
  print("printing parameter estimates...")
  
  # BYM2
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    # add on interval width
    designRes$parameters$width = designRes$parameters$Q90 - designRes$parameters$Q10
    
    # reorder the columns TODO: fix the BYM2 code so this is necessary, so that variance and medians are included, and so that urban effect CI and points are switched
    designRes$parameters = designRes$parameters[,c(3, 4, 1, 2, 5)]
    
    # rename the columns
    colnames(designRes$parameters) = c("Est", "SD", "Q10", "Q90", "80% CI Width")
    
    print(paste0("Parameter summary table for BYM2 ", typeText, " model:"))
    print(xtable(designRes$parameters, digits=3))
  }
  
  # SPDE
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    parameters = spdeResults
    parameters$resultsPixel = NULL
    parameters$resultsCounty = NULL
    parameters$resultsRegion = NULL
    parameters$resultsCluster = NULL
    parameters$pixelDraws = NULL
    parameters$fit = NULL
    parameters = do.call("rbind", parameters)
    parameters = parameters[,-3]
    
    # modify the row names do not include the word "Summary"
    allNames = rownames(parameters)
    rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
    
    # rename the columns
    colnames(parameters) = c("Est", "SD", "Q10", "Q50", "Q90", "80% CI Width")
    
    print(paste0("Parameter summary table for SPDE ", typeText, " model:"))
    print(xtable(parameters, digit=3))
  }
  
  invisible(NULL)
}

plotDataVisualizations = function(dat=ed, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                                  meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                                  sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                                  meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                                  meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                                  varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                                  sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                                  ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                                  plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  print("generating data visualizations...")
  
  # plot the actual data
  png(file=paste0("figures/", resultNameRoot, "/clustersUrban", plotNameRoot, ".png"), width=500, height=500)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
  urban = dat$urban
  plot(dat$lon[!urban], dat$lat[!urban], pch=19, col="green", main=paste0("Urban vs. rural clusters"), xlim=kenyaLonRange, 
       ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.2, asp=1)
  points(dat$lon[urban], dat$lat[urban], pch=19, col="blue", cex=.2)
  # world(add=TRUE)
  plotMapDat(adm1)
  dev.off()
  
  # plot a map of urbanicity if requested (this can take ~10 minutes)
  if(plotUrbanMap) {
    # inside this if statement since it takes around ten minutes to run
    makeUrbanMap(kmres=1, savePlot=TRUE, lonLim=kenyaLonRange, latLim=kenyaLatRange)
  }
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", plotNameRoot, ".png"), width=500, height=500)
  par(oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
  # plot(cbind(dat$lon, dat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
  #      xlab="Longitude", ylab="Latitude", main=paste0("Empirical ", varName), asp=1)
  plot(cbind(dat$lon, dat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
       xlab="Longitude", ylab="Latitude", main="", asp=1)
  quilt.plot(dat$lon, dat$lat, dat$y / dat$n, nx=150, ny=150, col=meanCols, add=TRUE)
  # world(add=TRUE)
  plotMapDat(adm1)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", plotNameRoot, "Logit.png"), width=500, height=500)
  par(oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6))
  ticks = pretty(seq(0, max(dat$y / dat$n), l=10), n=10)
  ticks = logit(ticks[-c(1, 11)])
  varRange = expit(range(ticks))
  # par( oma=c( 0,0,0,5)) # save some room for the legend
  plot(cbind(dat$lon, dat$lat), type="n", main=paste0("Kenya empirical ", varName), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  quilt.plot(cbind(dat$lon, dat$lat), logit(dat$y / dat$n), col=meanCols, 
             nx=100, ny=100, add.legend=FALSE, add=TRUE, zlim=logit(varRange))
  plotMapDat(adm1)
  # world(add=TRUE)
  # par( oma=c(0,0,0,2))
  image.plot(zlim=logit(varRange), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=meanCols, add = TRUE, axis.args=list(at=ticks, labels=expit(ticks)))
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", plotNameRoot, "Discrete.png"), width=500, height=500)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
  totals = aggregate(dat$n, list(dat$admin1), FUN=sum)
  counts = aggregate(dat$y, list(dat$admin1), FUN=sum)
  # plotMapDat(adm1, plotVar=counts$x / totals$x, new = TRUE, main=paste0("Empirical ", varName), cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=counts$x / totals$x, new = TRUE, main="", cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/denominatorsDiscrete", plotNameRoot, ".png"), width=500, height=500)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
  plotMapDat(adm1, plotVar=totals$x, new = TRUE, main=paste0("Sample size"), cols=popCols, zlim=c(0, max(totals$x)), xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/empirical", varName, "DiscreteLogit.png"), width=500, height=500)
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
  totals = aggregate(dat$n, list(dat$admin1), FUN=sum)
  counts = aggregate(dat$y, list(dat$admin1), FUN=sum)
  plotMapDat(adm1, plotVar=counts$x / totals$x, new = TRUE, main=paste0("Empirical ", varName), col=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  dev.off()
  
  png(file="figures/populationDiscrete.png", width=500, height=500)
  countyPops = poppc$popTotal
  sortI = sort(poppc$County, index.return=TRUE)$ix
  countyPops = countyPops[sortI]
  par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
  zlim = range(countyPops)
  plotMapDat(adm1, plotVar=countyPops, new = TRUE, main=paste0("Total population"), col=popCols, zlim=log(zlim), scaleFun=log, scaleFunInverse=exp, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  # entirely for testing the extendData function
  out = load(paste0("data4direct", resultNameRoot, ".RData"))
  png(file=paste0("figures/", resultNameRoot, "/extended", varName, "Discrete.png"), width=500, height=500)
  testY = aggregate(resDat$y, list(resDat$admin1), FUN=sum)
  testN = aggregate(rep(1, nrow(resDat)), list(resDat$admin1), FUN=sum)
  plotMapDat(adm1, plotVar=testY$x / testN$x, new = TRUE, main=paste0("Extended ", varName), cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  popGrid = makeInterpPopGrid(kmRes=5, adjustPopSurface = TRUE)
  thisPop = popGrid$popOrig / 5^2
  png(file=paste0("figures/populationDensityAdjusted.png"), width=550, height=500)
  par(oma=c( 0,0,0,5), mar=c(5.1, 5.1, 4.1, 4), lab.cex=2.0, cex.lab=1.5, cex.axis=1.5, cex.main=2)
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=TeX("Adjusted population density (people/$km^2$)"), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  popRange = log(c(.1, max(thisPop)))
  popTicks = log(pretty(exp(popRange)))
  popTicks = c(log(.1), log(1), log(10), log(100), log(1000), popTicks[-c(1, 2, 4, 5, 6)])
  tickLabels = as.character(exp(popTicks))
  # tickLabels[7] = ""
  quilt.plot(cbind(popGrid$lon, popGrid$lat, log(thisPop)), nx=150, ny=150, add.legend=FALSE, add=TRUE, 
             zlim=range(popTicks), col=popCols)
  plotMapDat(adm1, lwd=.5)
  image.plot(zlim=popRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add = TRUE, axis.args=list(at=popTicks, labels=tickLabels), legend.mar = 0)
  dev.off()
  
  popGrid = makeInterpPopGrid(kmRes=5)
  thisPop = popGrid$popOrig
  totalPop = 43*10^6 # from DHS 2014 survey final report page 2
  thisPop = totalPop * thisPop / sum(thisPop) / 5^2 # population density per km^2
  png(file=paste0("figures/populationDensity.png"), width=550, height=500)
  par(oma=c( 0,0,0,5), mar=c(5.1, 5.1, 4.1, 4), lab.cex=2.0, cex.lab=1.5, cex.axis=1.5, cex.main=2)
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=TeX("Population density (people/$km^2$)"), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  popRange = log(c(1, max(thisPop)))
  popTicks = log(pretty(exp(popRange)))
  popTicks = c(log(1), log(10), log(100), log(1000), popTicks[-c(1, 2, 4, 6)])
  tickLabels = as.character(exp(popTicks))
  tickLabels[7] = ""
  quilt.plot(cbind(popGrid$lon, popGrid$lat, log(thisPop)), nx=150, ny=150, add.legend=FALSE, add=TRUE, 
             zlim=range(popTicks), col=popCols)
  plotMapDat(adm1, lwd=.5)
  image.plot(zlim=popRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
             col=popCols, add = TRUE, axis.args=list(at=popTicks, labels=tickLabels), legend.mar = 0)
  dev.off()
  
  png(file=paste0("figures/populationDensityDiscrete.png"), width=500, height=500)
  test = aggregate(thisPop*5^2, list(popGrid$admin1), FUN=sum)
  zlim = range(test$x)
  ticks = c(100000, 500000, 1000000, 2000000, 3000000, 4000000)
  par(mar=c(5.1, 5.1, 4.1, 4), lab.cex=2.0, cex.lab=1.5, cex.axis=1.5, cex.main=2)
  plotMapDat(adm1, plotVar=test$x, new = TRUE, main=TeX("Integrated population density"), 
             col=popCols, zlim=log(zlim), scaleFun=log, scaleFunInverse=exp, ticks = ticks, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  # plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " SD of logit predictions"), ylim=kenyaLatRange, 
  #      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  # quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
  #            nx=150, ny=150, add.legend=FALSE, add=TRUE, zlim=range(log(sdRange2)))
  # plotMapDat(adm1, lwd=.5)
  # points(dat$lon, dat$lat, pch=".")
  # image.plot(zlim=range(log(sdRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
  #            col=cols, add = TRUE, axis.args=list(at=log(sdTicks2), labels=sdTickLabels2), legend.mar = 0)
  
  tab = totals
  tab$x = counts$x / totals$x
  tab
}

plotModelPredictions = function(dat=ed, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                                meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                                sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                                meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                                meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                                varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                                sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                                ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                                plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), widthRange, widthTicks, 
                                widthTickLabels, widthRangeSPDE, widthTicksSPDE, widthTickLabelsSPDE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  ##### Naive and direct estimates
  print("plotting direct and naive estimates...")
  
  # save(directEstEd, naiveEd, file="resultsDirectNaiveEd.RData")
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData"))
  plotName = paste0(resultNameRoot, "Preds")
  png(file=paste0("figures/", resultNameRoot, "/naive", plotName, ".png"), width=1000, height=1200)
  par(mfrow=c(2,2))
  zlim = range(c(expit(naiveResults$upper),expit(naiveResults$lower) ))
  plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive SD of logit predictions"), zlim=sdRange, cols=sdCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  png(file=paste0("figures/", resultNameRoot, "/naive", plotName, "Logit.png"), width=1000, height=1200)
  par(mfrow=c(2,2))
  zlim = range(c(expit(naiveResults$upper),expit(naiveResults$lower) ))
  plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels)
  plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive SD of logit predictions"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  dev.off()
  png(file=paste0("figures/", resultNameRoot, "/direct", plotName, ".png"), width=800, height=1000)
  par(mfrow=c(2,2))
  zlim = range(c(expit(directEstResults$upper),expit(directEstResults$lower) ))
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), 
             zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct SD of logit predictions"), 
             zlim=sdRange, cols=sdCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), 
             zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), 
             zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  png(file=paste0("figures/", resultNameRoot, "/direct", plotName, "Logit.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct SD of logit predictions"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  dev.off()
  
  ##### Mercer et al. estimates
  print("plotting mercer estimates...")
  
  # save(mercerEd, file=paste0("resultsMercerEd.RData"))
  plotName = paste0("mercer", resultNameRoot, "Preds")
  out = load(paste0("resultsMercer", resultNameRoot, ".RData" ))
  zlim = range(c(expit(mercerResults$lower.mercer),expit(mercerResults$upper.mercer)))
  png(file=paste0("figures/", resultNameRoot, "/", plotName, ".png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Smoothed direct ", varName, " estimates"), 
             zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Smoothed direct SD of logit predictions"), 
             zlim=sdRange, cols=sdCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct 10th percentile"), 
             zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Smoothed direct 90th percentile"), 
             zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/", plotName, "Logit.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Smoothed direct ", varName, " estimates"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Smoothed direct SD of logit predictions"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct 10th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Smoothed direct 90th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
  dev.off()
  
  zlim = range(c(expit(mercerResults$lower.mercer),expit(mercerResults$upper.mercer)))
  png(file=paste0("figures/", resultNameRoot, "/", plotName, "2.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Smoothed direct ", varName, " estimates"), 
             zlim=meanRange2, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Smoothed direct SD of logit predictions"), 
             zlim=sdRange2, cols=sdCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct 10th percentile"), 
             zlim=meanRange2, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Smoothed direct 90th percentile"), 
             zlim=meanRange2, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  png(file=paste0("figures/", resultNameRoot, "/", plotName, "Logit2.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Smoothed direct ", varName, " estimates"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, cols=meanCols)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Smoothed direct SD of logit predictions"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, cols=sdCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct 10th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, cols=meanCols)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Smoothed direct 90th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
             zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, cols=meanCols)
  dev.off()
  
  ##### BYM2 estimates
  print("plotting BYM2 estimates...")
  
  # save(file = paste0('bym2EdUrbRur',includeUrbanRural, 'Cluster', includeCluster, '.RData'), 
  #      designRes = designRes)
  # if(includeCluster) {
  #   save(file = paste0('bym2EdUrbRur',includeUrbanRural, 'Cluster', includeCluster, 'debiased.RData'), 
  #        designRes = designRes)
  # }
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    zlim = range(c(expit(designRes$predictions$Q10),expit(designRes$predictions$Q90)))
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, ".png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), 
               zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), 
               zlim=sdRange, cols=sdCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), 
               zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), 
               zlim=meanRange, cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
               zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
               zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, cols=sdCols)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
               zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
               zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, cols=meanCols)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "2.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), cols=sdCols, zlim=sdRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit2.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      # typeText = paste0(typeText, "'")
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, ".png"), width=800, height=1000)
      par(mfrow=c(2,2))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), cols=sdCols, zlim=sdRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
      dev.off()
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit.png"), width=800, height=1000)
      par(mfrow=c(2,2))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
      dev.off()
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "2.png"), width=800, height=1000)
      par(mfrow=c(2,2))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), cols=sdCols, zlim=sdRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      dev.off()
      
      png(file=paste0("figures/", resultNameRoot, "/", nameRoot, "Logit2.png"), width=800, height=1000)
      par(mfrow=c(2,2))
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
      dev.off()
    }
  }
  
  ##### SPDE estimates
  print("plotting SPDE estimates...")
  
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    zlim = range(c(spdeResults$resultsCounty$lower,spdeResults$resultsCounty$upper))
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, ".png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeText, " SD of logit predictions"), cols=sdCols, zlim=sdRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeText, " 10th percentile"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeText, " 90th percentile"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "Logit.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeText, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeText, " 10th percentile"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeText, " 90th percentile"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "2.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeText, " SD of logit predictions"), cols=sdCols, zlim=sdRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeText, " 10th percentile"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeText, " 90th percentile"), cols=meanCols, zlim=meanRange2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "Logit2.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeText, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeText, " 10th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeText, " 90th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    ## plot continuous prediction surface
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, ".png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeText, " SD of logit predictions"), cols=sdCols, zlim=sdRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeText, " 10th percentile"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeText, " 90th percentile"), cols=meanCols, zlim=meanRange, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "Logit.png"), width=800, height=1000)
    par(mfrow=c(2,2))
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeText, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeText, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeText, " 10th percentile"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeText, " 90th percentile"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousLogit.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " estimates"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " SD of logit predictions"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(log(sdRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks), labels=sdTickLabels), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 10th percentile"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 90th percentile"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks), labels=meanTickLabels), legend.mar = 5)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousLogit2.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " estimates"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " SD of logit predictions"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(log(sdRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicks2), labels=sdTickLabels2), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 10th percentile"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 90th percentile"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange2)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 5)
    dev.off()
    
    png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousLogitSelf.png"), width=800, height=1000)
    par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " estimates"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " SD of logit predictions"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=sdCols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 10th percentile"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 90th percentile"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 5)
    dev.off()
  }
  
  ##### now put the predictions from each together in the same plot when relevant
  # 1: put naive and direct together with standard deviations and credible intervals
  ## make two by four plot with these two models and lower, mean, and upper bounds, and SDs
  # 2: put values from the smooth models including all effects together:
  ## mercer, BYM2 with urban and cluster effects, SPDE with urban and cluster effects
  ## make three by four plot with these three models and lower, mean, and upper bounds, and SDs
  # 3: plot test spde models together
  ## 4 x 4 plot showing estimates, sds, 10% and 90% CIs
  # 4: plot relative differences between the models, relative to the BYM2 model including urban and cluster effects that has been debiased
  ## 2 x 2 plot
  # 5: plot the BYM2 models together
  ## 4 x 6 plot
  # 6: plot relative differences between the BYM2 models
  ## 3 x 2 plot (one of the subplots will be empty)
  # 7: plot relative differences between the SPDE models on a continuous scale
  ## 2 x 2 plot, all relative to the model with both urban and cluster effects
  # 8: plot relative differences between the SPDE models on a discrete scale
  ## 2 x 2 plot
  # 9: make pairs plot of the full models with respect to each other
  
  ## Plot 1: direct and naive models
  print("plotting direct and naive estimates together...")
  
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData" ))
  
  png(file=paste0("figures/", resultNameRoot, "/fullDirectNaive", plotNameRoot, ".png"), width=800, height=1200)
  par(mfrow=c(4,2))
  plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive SD of logit predictions"), cols=sdCols, zlim=log(sdRangeND), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksND, tickLabels=sdTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct SD of logit predictions"), cols=sdCols, zlim=log(sdRangeND), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksND, tickLabels=sdTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  ## Plot 2: mercer, BYM2 with urban and cluster effects, SPDE with urban and cluster effects
  print("plotting smoothed estimates together...")
  
  png(file=paste0("figures/", resultNameRoot, "/fullSmoothed", plotNameRoot, ".png"), width=1000, height=1200)
  par(mfrow=c(4,3))
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  debiasedText = "Debiased"
  # typeTextBYM = "UC"
  typeTextBYM = "UC"
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  # typeTextSPDE = "UC"
  typeTextSPDE = "UC"
  
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Smoothed direct ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Smoothed direct SD of logit predictions"), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeTextBYM, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks2, tickLabels=sdTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct 10th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " 10th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " 10th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Smoothed direct 90th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " 90th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " 90th percentile"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  ## Plot 3: plot all of the SPDE plots together (4 x 4 plot)
  print("plotting SPDE models together...")
  
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousTogetherSelf.png"), width=1400, height=1400)
  # par(mfrow=c(4,4), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousTogetherSelf.png"), width=1400, height=700)
  par(mfrow=c(2,4), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " estimates"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    if(i == 4) {
      image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
                 col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
    }
  }
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    # don't plot the SDs.  Instead plut 80% interval widths
    # plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " SD of logit predictions"), ylim=kenyaLatRange, 
    #      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    # quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
    #            nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRangeSPDE)))
    # plotMapDat(adm1, lwd=.5)
    # points(dat$lon, dat$lat, pch=".")
    # image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
    #            col=sdCols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 80% credible interval width"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$upper-spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(widthRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    if(i == 4) {
      image.plot(zlim=range(log(widthRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
                 col=sdCols, add = TRUE, axis.args=list(at=log(widthTicksSPDE), labels=widthTickLabelsSPDE), legend.mar = 0)
    }
  }
  
  # don't plot the actual credible interval:
  # for(i in 1:length(argList)) {
  #   args = argList[[i]]
  #   includeUrban = args$urbanEffect
  #   includeCluster = args$includeClustEffect
  #   clusterText = ifelse(includeCluster, "", "NoClust")
  #   
  #   nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
  #                     "_urbanEffect", includeUrban)
  #   out = load(paste0("results", nameRoot, '.RData'))
  #   
  #   urbanText = ifelse(includeUrban, "", "noUrb")
  #   clusterText = ifelse(includeCluster, "", "NoClust")
  #   both = includeUrban && includeCluster
  #   notBothText = ifelse(both, "", " ")
  #   # typeText = romanText(i)
  #   if(i == 1)
  #     typeText = "uc"
  #   else if(i == 2)
  #     typeText = "uC"
  #   else if(i == 3)
  #     typeText = "Uc"
  #   else if(i == 4)
  #     typeText = "UC"
  #   
  #   plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 10th percentile"), ylim=kenyaLatRange, 
  #        xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  #   quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$lower), 
  #              nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
  #   plotMapDat(adm1, lwd=.5)
  #   points(dat$lon, dat$lat, pch=".")
  #   image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
  #              col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
  # }
  # 
  # for(i in 1:length(argList)) {
  #   args = argList[[i]]
  #   includeUrban = args$urbanEffect
  #   includeCluster = args$includeClustEffect
  #   clusterText = ifelse(includeCluster, "", "NoClust")
  #   
  #   nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
  #                     "_urbanEffect", includeUrban)
  #   out = load(paste0("results", nameRoot, '.RData'))
  #   
  #   urbanText = ifelse(includeUrban, "", "noUrb")
  #   clusterText = ifelse(includeCluster, "", "NoClust")
  #   both = includeUrban && includeCluster
  #   notBothText = ifelse(both, "", " ")
  #   # typeText = romanText(i)
  #   if(i == 1)
  #     typeText = "uc"
  #   else if(i == 2)
  #     typeText = "uC"
  #   else if(i == 3)
  #     typeText = "Uc"
  #   else if(i == 4)
  #     typeText = "UC"
  #   
  #   plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 90th percentile"), ylim=kenyaLatRange, 
  #        xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  #   quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$upper), 
  #              nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
  #   plotMapDat(adm1, lwd=.5)
  #   points(dat$lon, dat$lat, pch=".")
  #   image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
  #              col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 5)
  # }
  dev.off()
  
  ##### Plot 3.5: plot continuous cluster models together
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousTogetherSelf.png"), width=1400, height=1400)
  # par(mfrow=c(4,4), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousTogetherSelf.png"), width=1400, height=700)
  par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uC"
    else if(i == 2)
      typeText = "UC"
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " estimates"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    if(i == 2) {
      image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
                 col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
    }
  }
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uC"
    else if(i == 2)
      typeText = "UC"
    
    # don't plot the SDs.  Instead plut 80% interval widths
    # plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " SD of logit predictions"), ylim=kenyaLatRange, 
    #      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    # quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
    #            nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRangeSPDE)))
    # plotMapDat(adm1, lwd=.5)
    # points(dat$lon, dat$lat, pch=".")
    # image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
    #            col=sdCols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 80% credible interval width"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$upper-spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(widthRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    if(i == 2) {
      image.plot(zlim=range(log(widthRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
                 col=sdCols, add = TRUE, axis.args=list(at=log(widthTicksSPDE), labels=widthTickLabelsSPDE), legend.mar = 0)
    }
  }
  dev.off()
  
  ##### Plot 4: plot relative differences between the models
  print("plotting relative differences between models...")
  
  ## all together relative to the bym2
  # estimates
  # library("colorspace")
  # pal <-choose_palette()
  # cols=diverging_hcl(101, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
  # cols=
  
  nPerSide = (ncols - 1) / 2
  # cols=diverging_hcl(ncols, h1=265, h2=17, c1=100, l1=50, l2=92, p1=0.6)
  cols=relativeCols
  png(file=paste0("figures/", resultNameRoot, "/fullRelative", plotNameRoot, ".png"), width=800, height=1000)
  par(mfrow=c(2, 2))
  
  zlim = range(c((plotVar=naiveResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
                 (directEstResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
                 (mercerResults$est.mercer - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
                 (spdeResults$resultsCounty$pred - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean)))
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  plotMapDat(adm1, plotVar=(plotVar=naiveResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("Naive relative to BYM2 ", varName, " estimates"), zlim=zlim, col=cols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=(directEstResults$est - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("Direct relative to BYM2 ", varName, " estimates"), zlim=zlim, col=cols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=(mercerResults$est.mercer - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("Smoothed direct relative to BYM2 ", varName, " estimates"), zlim=zlim, col=cols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=(spdeResults$resultsCounty$pred - expit(designRes$predictions$mean)) / expit(designRes$predictions$mean), 
             new = TRUE, main=paste0("SPDE relative to BYM2 ", varName, " estimates"), zlim=zlim, col=cols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  dev.off()
  
  ## Plot 5: plot the BYM2 models together (four by six plot)
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  png(file=paste0("figures/", resultNameRoot, "/predsBYM2togetherSelf", resultNameRoot, "Logit.png"), width=1500, height=1000)
  par(mfrow=c(4,6))
  
  # plot the estimates
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), 
               zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2, 
               cols=meanCols, xlim=kenyaLonRange, ylim=kenyaLatRange)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeCluster
      if(!both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      # typeText = paste0(typeText, "'")
      
      
      plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " estimates"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
                 cols=meanCols, zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2)
    }
  }
  
  # plot the BYM2 standard deviations together
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
               zlim=log(sdRangeBYM2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksBYM2, tickLabels=sdTickLabelsBYM2, cols=sdCols)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeCluster
      if(!both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      # typeText = paste0(typeText, "'")
      
      plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeText, " SD of logit predictions"), 
                 cols=sdCols, zlim=log(sdRangeBYM2), scaleFun=log, scaleFunInverse=exp, ticks=sdTicksBYM2, 
                 tickLabels=sdTickLabelsBYM2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    }
  }
  
  # plot the BYM2 tenth percentiles together
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
               zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2, cols=meanCols)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeCluster
      if(!both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      # typeText = paste0(typeText, "'")
      
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeText, " 10th percentile"), 
                 cols=meanCols, zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, 
                 tickLabels=meanTickLabelsBYM2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    }
  }
  
  # plot the BYM2 ninetieth percentiles together
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), xlim=kenyaLonRange, ylim=kenyaLatRange, 
               zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, tickLabels=meanTickLabelsBYM2, cols=meanCols)
    
    if(includeCluster) {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeCluster
      if(!both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      # typeText = paste0(typeText, "'")
      
      plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeText, " 90th percentile"), 
                 cols=meanCols, zlim=logit(meanRangeBYM2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksBYM2, 
                 tickLabels=meanTickLabelsBYM2, xlim=kenyaLonRange, ylim=kenyaLatRange)
    }
  }
  dev.off()
  
  ## Plot 6: plot the BYM2 model relative differences versus the BYM2 UC model
  # 2 x 3 plot. Columns: debiased cluster, cluster, no cluster. Rows: urban, no urban
  argList = list(list(includeUrbanRural = TRUE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = FALSE, includeCluster = FALSE))
  png(file=paste0("figures/", resultNameRoot, "/fullRelativeBYM2", resultNameRoot, ".png"), width=1200, height=800)
  par(mfrow=c(2,3))
  
  # before the plotting begins, determine the range of relative differences
  zlim=c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    
    if(includeCluster) {
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      if(i == 1) {
        # get the values relative to which other predictions will be compared
        expected = expit(designRes$predictions$mean)
      }
      else {
        vals = (expit(designRes$predictions$mean) - expected) / expected
        zlim = range(c(zlim, vals))
      }
    }
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    vals = (expit(designRes$predictions$mean) - expected) / expected
    zlim = range(c(zlim, vals))
  }
  
  # now shift the color scales so that it center occurs at 0
  cols = relativeCols
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  
  # finally, plot the relative differences
  expectedTypeText = " full debiased"
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    
    if(includeCluster) {
      # gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeCluster
      if(!both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      # typeText = paste0(typeText, "'")
      
      if(i == 1) {
        expected = expit(designRes$predictions$mean)
      }
      vals = expit(designRes$predictions$mean)
      plotMapDat(adm1, plotVar=(vals - expected) / expected, new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " relative to full debiased"), 
                 cols=cols, zlim=zlim, xlim=kenyaLonRange, ylim=kenyaLatRange)
    }
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = expit(designRes$predictions$mean)
    plotMapDat(adm1, plotVar=(vals - expected) / expected, new = TRUE, main=paste0("BYM2 ", typeText, " ", varName, " relative to full debiased"), 
               zlim=zlim, cols=cols, xlim=kenyaLonRange, ylim=kenyaLatRange)
  }
  dev.off()
  
  ## Plot 7: plot the SPDE model relative differences versus the SPDE UC model on continuous scale
  # 2 x 2 plot. Columns: cluster, no cluster. Rows: urban, no urban
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE))
  png(file=paste0("figures/", resultNameRoot, "/fullRelativeSPDE", resultNameRoot, "Continuous.png"), width=800, height=1000)
  par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  
  # before the plotting begins, determine the range of relative differences
  zlim=c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    if(i == 1)
      expected = spdeResults$resultsPixel$pred
    
    vals = (spdeResults$resultsPixel$pred - expected) / expected
    zlim = range(c(zlim, vals))
  }
  
  # now shift the color scales so that it center occurs at 0
  cols = relativeCols
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  
  # finally, plot the relative differences
  expectedTypeText = " full"
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsPixel$pred
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " relative to SPDE UC"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), (vals-expected)/expected, 
               nx=150, ny=150, add.legend=TRUE, add=TRUE, col=cols, zlim=zlim)
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    # image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
    #            col=cols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
  }
  dev.off()
  
  ## Plot 8: plot the SPDE model relative differences versus the SPDE UC model on discrete scale
  # 2 x 2 plot. Columns: cluster, no cluster. Rows: urban, no urban
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE))
  png(file=paste0("figures/", resultNameRoot, "/fullRelativeSPDE", resultNameRoot, "Discrete.png"), width=800, height=1000)
  par(mfrow=c(2,2))
  
  # before the plotting begins, determine the range of relative differences
  zlim=c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    if(i == 1)
      expected = spdeResults$resultsCounty$pred
    
    vals = (spdeResults$resultsCounty$pred - expected) / expected
    zlim = range(c(zlim, vals))
  }
  
  # now shift the color scales so that it center occurs at 0
  cols = relativeCols
  lengthUp = zlim[2] - 0
  lengthDown = 0 - zlim[1]
  if(lengthUp >= lengthDown) {
    # shorten the blue part of the color scale
    propDown = lengthDown / lengthUp
    propUp = 1
  } else {
    propDown = 1
    propUp = lengthUp / lengthDown
  }
  # cols = pal(101) # keep an extra color in the middle for 0
  numDown = round(nPerSide * propDown)
  numUp = round(nPerSide * propUp)
  cols = c(cols[(nPerSide + 1-numDown):nPerSide], cols[nPerSide + 1], cols[(nPerSide + 2):(nPerSide + 1 + numUp)])
  
  # finally, plot the relative differences
  expectedTypeText = " full"
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsCounty$pred
    plotMapDat(adm1, plotVar=(vals-expected)/expected, new = TRUE, 
               main=paste0("SPDE ", typeText, " ", varName, " relative to SPDE UC"), 
               cols=cols, zlim=zlim, xlim=kenyaLonRange, ylim=kenyaLatRange)
  }
  dev.off()
  
  ##### Plot 8.5: all estimates and all standard deviations
  print("plotting all estimates together...")
  
  png(file=paste0("figures/", resultNameRoot, "/predsAll", plotNameRoot, ".png"), width=2300, height=1200)
  par(mfrow=c(2,4), lab.cex=2.0, cex.lab=2.5, cex.axis=2.0, cex.main=2.5, oma=c( 0,1,0,5), mar=c(5.1, 6.1, 5.1, 5))
  
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData" ))
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  debiasedText = "Debiased"
  # typeTextBYM = "UC"
  typeTextBYM = "UC"
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  # typeTextSPDE = "UC"
  typeTextSPDE = "UC"
  
  # plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Smoothed direct ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  # plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Smoothed direct SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeTextBYM, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  dev.off()
  
  ##### Plot 8.5.5: all 10th and 90th quantiles
  print("plotting all quantiles together...")
  
  png(file=paste0("figures/", resultNameRoot, "/predsAllQuantiles", plotNameRoot, ".png"), width=2300, height=1200)
  par(mfrow=c(2,4), lab.cex=2.0, cex.lab=2.5, cex.axis=2.0, cex.main=2.5, oma=c( 0,1,0,5), mar=c(5.1, 6.1, 5.1, 5))
  
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData" ))
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  debiasedText = "Debiased"
  # typeTextBYM = "UC"
  typeTextBYM = "UC"
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  # typeTextSPDE = "UC"
  typeTextSPDE = "UC"
  
  # plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Smoothed direct 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  # plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  dev.off()
}

plotModelPredictionsPart2 = function(dat=ed, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                                     meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                                     sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                                     meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                                     meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                                     varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                                     sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                                     ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                                     plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), 
                                     sharedPredictionScale=TRUE, widthRange, widthTicks, widthTickLabels, widthRangeSPDE, 
                                     widthTicksSPDE, widthTickLabelsSPDE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  ##### Plot 3.5: plot continuous cluster models together
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # png(file=paste0("figures/", resultNameRoot, "/preds", nameRoot, "ContinuousTogetherSelf.png"), width=1400, height=1400)
  # par(mfrow=c(4,4), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  
  print("plotting continuous SPDE 'C' model predictions together")
  png(file=paste0("figures/", resultNameRoot, "/preds", resultNameRoot, "ContinuousTogetherSelf.png"), width=700, height=700)
  par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uC"
    else if(i == 2)
      typeText = "UC"
    
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " estimates"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    if(i == 2) {
      image.plot(zlim=range(logit(meanRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
                 col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicksSPDE), labels=meanTickLabelsSPDE), legend.mar = 0)
    }
  }
  
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uC"
    else if(i == 2)
      typeText = "UC"
    
    # don't plot the SDs.  Instead plut 80% interval widths
    # plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " SD of logit predictions"), ylim=kenyaLatRange, 
    #      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    # quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$sds), 
    #            nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(sdRangeSPDE)))
    # plotMapDat(adm1, lwd=.5)
    # points(dat$lon, dat$lat, pch=".")
    # image.plot(zlim=range(log(sdRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
    #            col=sdCols, add = TRUE, axis.args=list(at=log(sdTicksSPDE), labels=sdTickLabelsSPDE), legend.mar = 0)
    plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " 80% credible interval width"), ylim=kenyaLatRange, 
         xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
    quilt.plot(cbind(popGrid$lon, popGrid$lat), log(spdeResults$resultsPixel$upper-spdeResults$resultsPixel$lower), 
               nx=150, ny=150, add.legend=FALSE, add=TRUE, col=sdCols, zlim=range(log(widthRangeSPDE)))
    plotMapDat(adm1, lwd=.5)
    points(dat$lon, dat$lat, pch=".")
    if(i == 2) {
      image.plot(zlim=range(log(widthRangeSPDE)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
                 col=sdCols, add = TRUE, axis.args=list(at=log(widthTicksSPDE), labels=widthTickLabelsSPDE), legend.mar = 0)
    }
  }
  dev.off()
  
  ##### Plot 8.5: all estimates and all standard deviations
  print("plotting all estimates together...")
  
  png(file=paste0("figures/", resultNameRoot, "/predsAll", plotNameRoot, ".png"), width=2300, height=1200)
  par(mfrow=c(2,4), lab.cex=2.0, cex.lab=2.5, cex.axis=2.0, cex.main=2.5, oma=c( 0,1,0,5), mar=c(5.1, 6.1, 5.1, 5))
  
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData" ))
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  debiasedText = "Debiased"
  # typeTextBYM = "UC"
  typeTextBYM = "UC"
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  # typeTextSPDE = "UC"
  typeTextSPDE = "UC"
  
  # plotMapDat(adm1, plotVar=naiveResults$est, new = TRUE, main=paste0("Naive ", varName, " estimates"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=directEstResults$est, new = TRUE, main=paste0("Direct ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=mercerResults$est.mercer, new = TRUE, main=paste0("Smoothed direct ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$pred, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks, tickLabels=meanTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  # plotMapDat(adm1, plotVar=sqrt(naiveResults$var.est), new = TRUE, main=paste0("Naive SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  # plotMapDat(adm1, plotVar=sqrt(directEstResults$var.est), new = TRUE, main=paste0("Direct SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  # plotMapDat(adm1, plotVar=sqrt(mercerResults$var.est.mercer), new = TRUE, main=paste0("Smoothed direct SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  # plotMapDat(adm1, plotVar=designRes$predictions$stddev, new = TRUE, main=paste0("BYM2 ", typeTextBYM, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  # plotMapDat(adm1, plotVar=spdeResults$resultsCounty$sds, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " SD of logit predictions"), cols=sdCols, zlim=log(sdRange), scaleFun=log, scaleFunInverse=exp, ticks=sdTicks, tickLabels=sdTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower)-expit(directEstResults$upper), new = TRUE, main=paste0("Direct 80% CI width"), cols=sdCols, zlim=log(widthRange), scaleFun=log, scaleFunInverse=exp, ticks=widthTicks, tickLabels=widthTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer)-expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct  80% CI width"), cols=sdCols, zlim=log(widthRange), scaleFun=log, scaleFunInverse=exp, ticks=widthTicks, tickLabels=widthTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90)-expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeTextBYM, "  80% CI width"), cols=sdCols, zlim=log(widthRange), scaleFun=log, scaleFunInverse=exp, ticks=widthTicks, tickLabels=widthTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper-spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeTextSPDE, "  80% CI width"), cols=sdCols, zlim=log(widthRange), scaleFun=log, scaleFunInverse=exp, ticks=widthTicks, tickLabels=widthTickLabels, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  dev.off()
  
  ##### Plot 8.5.5: all 10th and 90th quantiles
  print("plotting all quantiles together...")
  
  png(file=paste0("figures/", resultNameRoot, "/predsAllQuantiles", plotNameRoot, ".png"), width=2300, height=1200)
  par(mfrow=c(2,4), lab.cex=2.0, cex.lab=2.5, cex.axis=2.0, cex.main=2.5, oma=c( 0,1,0,5), mar=c(5.1, 6.1, 5.1, 5))
  
  out = load(paste0("resultsDirectNaive", resultNameRoot, ".RData" ))
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  debiasedText = "Debiased"
  # typeTextBYM = "UC"
  typeTextBYM = "UC"
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  # typeTextSPDE = "UC"
  typeTextSPDE = "UC"
  
  # plotMapDat(adm1, plotVar=expit(naiveResults$lower), new = TRUE, main=paste0("Naive 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(directEstResults$lower), new = TRUE, main=paste0("Direct 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(mercerResults$upper.mercer), new = TRUE, main=paste0("Smoothed direct 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q90), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$upper, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " 90th percentile"), cols=meanCols, zlim=log(meanRangeND), scaleFun=log, scaleFunInverse=exp, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  # plotMapDat(adm1, plotVar=expit(naiveResults$upper), new = TRUE, main=paste0("Naive 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(directEstResults$upper), new = TRUE, main=paste0("Direct 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(mercerResults$lower.mercer), new = TRUE, main=paste0("Smoothed direct 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$Q10), new = TRUE, main=paste0("BYM2 ", typeTextBYM, " 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  plotMapDat(adm1, plotVar=spdeResults$resultsCounty$lower, new = TRUE, main=paste0("SPDE ", typeTextSPDE, " 10th percentile"), cols=meanCols, zlim=logit(meanRangeND), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicksND, tickLabels=meanTickLabelsND, xlim=kenyaLonRange, ylim=kenyaLatRange, legend.width=5)
  
  dev.off()
}

printModelPredictionTables = function(dat=ed, resultNameRoot="Ed", counties=sort(unique(poppc$County)), 
                                      nDigitsPredictions=3, nDigitsSmoothedDirect=3, nDigitsBYM2=1, nDigitsSPDE=2, byRow=FALSE) {
  resultNameRootLower = tolower(resultNameRoot)
  
  ##### Mercer et al. estimates
  print("getting mercer estimates...")
  
  # save(mercerEd, file=paste0("resultsMercerEd.RData"))
  plotName = paste0("mercer", resultNameRoot, "Preds")
  out = load(paste0("resultsMercer", resultNameRoot, ".RData" ))
  
  tab = cbind(mercerResults$est.mercer, expit(mercerResults$lower.mercer), expit(mercerResults$upper.mercer))
  
  ##### BYM2 estimates
  print("getting BYM2 estimates...")
  
  includeUrban = TRUE
  includeCluster = TRUE
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  tab = cbind(tab, expit(designRes$predictions$mean), expit(designRes$predictions$Q10), expit(designRes$predictions$Q90))
  
  ##### SPDE estimates
  print("getting SPDE estimates...")
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  tab = cbind(tab, spdeResults$resultsCounty$pred, spdeResults$resultsCounty$lower, spdeResults$resultsCounty$upper)
  
  ##### now we construct the latex table
  tab = format(tab, digits=nDigitsPredictions) # round to the nearest 1 child per 1000
  tab = cbind(counties, tab)
  colnames(tab) = c("County", rep(c("Est", "Q10", "Q90"), 3))
  rownames(tab) = NULL
  
  require(kableExtra)
  fullTab = tab %>%
    kable("latex", escape = F, booktabs = TRUE, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
          align=c("l", rep("r", ncol(tab) - 1)), longtable=TRUE, caption = "Longtable") %>% 
    kable_styling(latex_options =c("repeat_header", "scale_down"))
  numberColumns = c(" "=1, "Smoothed Direct"=3, "BYM2 UC"=3, "SPDE UC"=3)
  print(add_header_above(fullTab, numberColumns, italic=TRUE, bold=TRUE, escape=FALSE, line=TRUE))
  
  ##### now construct the parameter tables
  ##### now print the parameter estimates:
  print("printing parameter estimates...")
  
  # smoothed direct
  out = load(paste0("resultsMercer", resultNameRoot, ".RData" ))
  parTable = round(mercerParResults, digits=nDigitsSmoothedDirect)
  
  # BYM2
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  # typeText = romanText(i)
  typeText = "UC"
  # add on interval width
  designRes$parameters$width = designRes$parameters$Q90 - designRes$parameters$Q10
  
  # rename the columns
  colnames(designRes$parameters) = c("Est", "SD", "Q10", "Q50", "Q90", "80% CI Width")
  
  # print(paste0("Parameter summary table for BYM2 ", typeText, " model:"))
  # print(xtable(designRes$parameters, digits=3))
  colnames(designRes$parameters[,1:5]) = colnames(parTable)
  parTable = rbind(parTable, 
                   format(as.matrix(designRes$parameters[,1:5]), digits=nDigitsBYM2, scientific=FALSE))
  
  # SPDE
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  # typeText = romanText(i)
  typeText = "UC"
  
  parameters = spdeResults
  parameters$resultsPixel = NULL
  parameters$resultsCounty = NULL
  parameters$resultsRegion = NULL
  parameters$resultsCluster = NULL
  parameters$pixelDraws = NULL
  parameters$fit = NULL
  parameters = do.call("rbind", parameters)
  parameters = parameters[,c(-3, -7)]
  
  # modify the row names do not include the word "Summary"
  allNames = rownames(parameters)
  # rownames(parameters) = unlist(sapply(allNames, strsplit, split="Summary"))
  
  if(byRow == FALSE)
    nonRangePar = format(parameters[-3,], digits=nDigitsSPDE, scientific=FALSE)
  else
    nonRangePar = t(apply(parameters[-3,], 1, format, digits=2, scientific=FALSE))
  formattedParameters = rbind(nonRangePar[1:2,], format(parameters[3,], digits=0), nonRangePar[3:nrow(nonRangePar),])
  rownames(formattedParameters) = c("Intercept", "Urban", "Range", "Spatial Var", "Spatial SD", "Nugget Var", "Nugget SD")
  
  # rename the columns
  colnames(parameters) = c("Est", "SD", "Q10", "Q50", "Q90")
  
  # print(paste0("Parameter summary table for SPDE ", typeText, " model:"))
  # print(xtable(parameters, digit=3))
  colnames(parameters) = colnames(parTable)
  parTable = rbind(parTable, formattedParameters)
  parTable = cbind(rownames(parTable), parTable)
  rownames(parTable) = NULL
  
  fullTab = kable(parTable, "latex", booktabs = T, escape=FALSE, format.args=list(drop0trailing=FALSE, scientific=FALSE), 
                  align=c("l", rep("r", ncol(parTable) - 1))) %>% 
    kable_styling()
  fullTab = fullTab %>% 
    pack_rows("Smoothed Direct", 1, 8, escape=FALSE, bold=TRUE, italic=TRUE) %>% 
    pack_rows("BYM2 UC", 9, 19, escape=FALSE, bold=TRUE, italic=TRUE) %>% 
    pack_rows("SPDE UC", 20, nrow(parTable), escape=FALSE, bold=TRUE, italic=TRUE)
  print(fullTab)
}

makeRedBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=10, h2=-115, c1=100, c2=100, l1=44, l2=59, p1=0, p2=2.3)
}

makeGreenBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=128, h2=250, c1=117, cmax=74, c2=107, l1=71, l2=55, p1=2, p2=2)
}

makeRedBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange)  || is.null(center)) {
    diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange
    propUp = (valRange[2] - center) / diff(valRange)
    propDown = 1 - propUp
    totalColors = ceiling(2 * max(propUp, propDown) * n)
    tempColors = makeRedBlueDivergingColors(totalColors, rev=rev)
    totalMissingColors = totalColors - n
    
    if(propUp >= propDown)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
  }
}

makeRedGrayBlueDivergingColors = function(n, valRange=NULL, center=NULL, rev=FALSE) {
  # library("colorspace")
  # pal <-choose_palette()
  if(is.null(valRange)  || is.null(center)) {
    diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=90, p1=0.9, rev=rev)
    # diverging_hcl(n, h1=10, h2=-115, c1=90, l1=40, l2=100, p1=0.9, p2=0.6)
  }
  else {
    # in this case we want white to be at the center of valRange
    propUp = (valRange[2] - center) / diff(valRange)
    propDown = 1 - propUp
    totalColors = ceiling(2 * max(propUp, propDown) * n)
    tempColors = makeRedGrayBlueDivergingColors(totalColors, rev=rev)
    totalMissingColors = totalColors - n
    
    if(propUp >= propDown)
      tempColors[-(1:totalMissingColors)]
    else
      tempColors[1:n]
  }
}

makeBlueSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  # sequential_hcl(n, h1=260, c1=80, l1=30, l2=90, p1=1.5, rev=TRUE)
  sequential_hcl(n, h1=245, c1=50, cmax=75, l1=20, l2=98, p1=0.8, rev=TRUE)
}

makeBlueYellowSequentialColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  sequential_hcl(n, h1=300, h2=75, c1=40, c2=95, l1=15, l2=90, p1=1.0, p2=1.1)
}

makeRedGreenDivergingColors = function(n) {
  # library("colorspace")
  # pal <-choose_palette()
  diverging_hcl(n, h1=265, h2=101, c1=100, l1=50, l2=92, p1=0.6, p2=1.5)
}

myPairs = function(x, labels, panel = points, ..., horInd = 1:nc, verInd = 1:nc, 
          lower.panel = panel, upper.panel = panel, diag.panel = NULL, 
          text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3, 
          cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1, 
          log = "", lims=NULL) 
{
  if (doText <- missing(text.panel) || is.function(text.panel)) 
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
                                                                 y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, i, j, xpd, bg, col = NULL, main, 
                        oma, ...) {
    if(!is.null(lims)) {
      x = lims[[i]]
      y = lims[[j]]
    }
    xpd <- NA
    if (side%%2L == 1L && xl[j]) 
      xpd <- FALSE
    if (side%%2L == 0L && yl[i]) 
      xpd <- FALSE
    if (side%%2L == 1L) 
      Axis(x, side = side, xpd = xpd, ...)
    else Axis(y, side = side, xpd = xpd, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])) 
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]]))) 
        stop("non-numeric argument to 'pairs'")
    }
  }
  else if (!is.numeric(x)) 
    stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2L) 
    stop("only one column in the argument to 'pairs'")
  if (!all(horInd >= 1L && horInd <= nc)) 
    stop("invalid argument 'horInd'")
  if (!all(verInd >= 1L && verInd <= nc)) 
    stop("invalid argument 'verInd'")
  if (doText) {
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels)) 
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels)) 
      doText <- FALSE
  }
  oma <- if ("oma" %in% nmdots) 
    dots$oma
  main <- if ("main" %in% nmdots) 
    dots$main
  if (is.null(oma)) 
    oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
  opar <- par(mfcol = c(length(horInd), length(verInd)), mar = rep.int(gap/2, 
                                                                       4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  xl <- yl <- logical(nc)
  if (is.numeric(log)) 
    xl[log] <- yl[log] <- TRUE
  else {
    xl[] <- grepl("x", log)
    yl[] <- grepl("y", log)
  }
  ni <- length(iSet <- if (row1attop) horInd else rev(horInd))
  nj <- length(jSet <- verInd)
  for (j in jSet) for (i in iSet) {
    l <- paste0(if (xl[j]) 
      "x"
      else "", if (yl[i]) 
        "y"
      else "")
    if(is.null(lims)) {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", ..., log = l)
    } else {
      localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, 
                type = "n", xlim=lims[[j]], ylim=lims[[i]], ..., log = l)
    }
    if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
      box()
      j.odd <- (match(j, jSet) + !row1attop)%%2L
      i.odd <- (match(i, iSet) + !row1attop)%%2L
      if (i == iSet[1L] && (!j.odd || !has.upper || !has.lower)) 
        localAxis(3L, x[, j], x[, i], j, i, ...)
      if (i == iSet[ni] && (j.odd || !has.upper || !has.lower)) 
        localAxis(1L, x[, j], x[, i], j, i, ...)
      if (j == jSet[1L] && (!i.odd || !has.upper || !has.lower)) 
        localAxis(2L, x[, j], x[, i], j, i, ...)
      if (j == jSet[nj] && (i.odd || !has.upper || !has.lower)) 
        localAxis(4L, x[, j], x[, i], j, i, ...)
      mfg <- par("mfg")
      if (i == j) {
        if (has.diag) 
          localDiagPanel(as.vector(x[, i]), ...)
        if (doText) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          xlp <- if (xl[i]) 
            10^0.5
          else 0.5
          ylp <- if (yl[j]) 
            10^label.pos
          else label.pos
          text.panel(xlp, ylp, labels[i], cex = cex.labels, 
                     font = font.labels)
        }
      }
      else if (i < j) 
        localLowerPanel(as.vector(x[, j]), as.vector(x[, 
                                                       i]), ...)
      else localUpperPanel(as.vector(x[, j]), as.vector(x[, 
                                                          i]), ...)
      if (any(par("mfg") != mfg)) 
        stop("the 'panel' function made a new plot")
    }
    else par(new = FALSE)
  }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots) 
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots) 
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main, 
          font = font.main)
  }
  invisible(NULL)
}

makePairPlots = function(dat=ed, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                         meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                         sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                         meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                         meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                         varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                         sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                         ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                         plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), makeScreenSplitPlot=FALSE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  print("printing pair plots...")
  
  ## calculate the colors of each of the counties based on urbanicity
  # first get the proportion of population of each county that
  popI = match(names(spdeResults$resultsCounty$pred), poppc$County)
  propUrban = poppc$pctUrb[popI] / 100
  
  ## Plot 9: make pair plots of the models
  # now get the color index of each county
  colI = cut(propUrban, breaks=seq(0 - .0001, 1, l=ncols+1), labels=FALSE)
  countyCols = urbCols[colI]
  
  ## do the same for regions
  popI = match(names(spdeResults$resultsRegion$pred), poppr$Region)
  propUrban = poppr$pctUrb[popI]
  
  ## Plot 9: make pair plots of the models
  # now get the color index of each county
  colI = cut(propUrban, breaks=seq(0 - .0001, 1, l=ncols+1), labels=FALSE)
  regionCols = urbCols[colI]
  
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlot", resultNameRoot, ".pdf"), width=6, height=6)
  # first load BYM2 UC and SPDE models
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  debiasedText = "Debiased"
  typeTextBYM = "UC"
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  typeTextSPDE = "UC"
  
  valMat = cbind(naiveResults$est, directEstResults$est, 
                 mercerResults$est.mercer, 
                 expit(designRes$predictions$mean), spdeResults$resultsCounty$pred)
  zlim = range(valMat)
  zlim2 = range(valMat[,3:5])
  
  # valMat = rbind(1:5, valMat)
  my_line <- function(x,y,...){
    if(diff(range(x)) >= .04)
      xlim = zlim
    else
      xlim = zlim2
    if(diff(range(y)) >= .04)
      ylim = zlim
    else
      ylim = zlim2
    # if(diff(range(c(x, y))) > 0.04)
    #   par(usr = c(zlim, zlim))
    # else
    #   par(usr = c(zlim2, zlim2))
    # par(usr = c(xlim, ylim))
    abline(a = 0,b = 1,...)
    points(x,y,..., col="blue")
    # points(x,y,..., col=countyCols)
  }
  
  # pairs(valMat, 
  #       c("Naive", "Direct", "Mercer et al.", "BYM2 UC", "SPDE UC"), 
  #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
  #       main=paste0("County ", varName, " estimate comparisons"))
  lims = c(list(zlim), list(zlim), list(zlim2), list(zlim2), list(zlim2))
  myPairs(valMat, 
          c("Naive", "Direct", "Smoothed direct", "BYM2 UC", "SPDE UC"), 
          pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
          main=paste0("County ", varName, " estimate comparisons"), 
          lims=lims)
  dev.off()
  
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlot", resultNameRoot, "UrbCol.pdf"), width=6, height=6)
  # first load BYM2 UC and SPDE models
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  debiasedText = "Debiased"
  typeTextBYM = "UC"
  
  includeUrban = TRUE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  typeTextSPDE = "UC"
  
  valMat = cbind(naiveResults$est, directEstResults$est, 
                 mercerResults$est.mercer, 
                 expit(designRes$predictions$mean), spdeResults$resultsCounty$pred)
  zlim = range(valMat)
  zlim2 = range(valMat[,3:5])
  
  # valMat = rbind(1:5, valMat)
  my_line <- function(x,y,...){
    if(diff(range(x)) >= .04)
      xlim = zlim
    else
      xlim = zlim2
    if(diff(range(y)) >= .04)
      ylim = zlim
    else
      ylim = zlim2
    # if(diff(range(c(x, y))) > 0.04)
    #   par(usr = c(zlim, zlim))
    # else
    #   par(usr = c(zlim2, zlim2))
    # par(usr = c(xlim, ylim))
    # points(x,y,..., col="blue")
    abline(a = 0,b = 1,...)
    points(x,y,..., col=countyCols)
  }
  
  # pairs(valMat, 
  #       c("Naive", "Direct", "Mercer et al.", "BYM2 UC", "SPDE UC"), 
  #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
  #       main=paste0("County ", varName, " estimate comparisons"))
  lims = c(list(zlim), list(zlim), list(zlim2), list(zlim2), list(zlim2))
  myPairs(valMat, 
          c("Naive", "Direct", "Smoothed direct", "BYM2 UC", "SPDE UC"), 
          pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
          main=paste0("County ", varName, " estimate comparisons"), 
          lims=lims, oma=c(3,3,6,7))
  image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=29, legend.mar=3.3, col=urbCols, add=TRUE, 
             legend.lab = "Urbanicity", legend.line=1.2, legend.width=.5, legend.shrink=.8, 
             legend.cex=.8, axis.args=list(cex.axis=.5, tck=-1, hadj=.8))
  dev.off()
  
  ## Plot 10: make pair plots of the BYM2 models
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotBYM2", resultNameRoot, ".pdf"), width=7, height=7)
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    # don't include BYM2 results without including cluster overdispersion
    if(!includeCluster) {
      valMat = cbind(valMat, expit(designRes$predictions$mean))
      labels = c(labels, paste0("BYM2 ", typeText))
    }
    else {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeCluster
      if(!both)
        debiasedText = "Debiased"
      else
        debiasedText = "debiased"
      # typeText = paste0(typeText, "'")
      
      valMat = cbind(valMat, expit(designRes$predictions$mean))
      labels = c(labels, paste0("BYM2 ", typeText))
    }
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col="blue")
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("BYM2 ", varName, " estimate comparisons"), 
        ylim=zlim, xlim=zlim)
  dev.off()
  
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotBYM2", resultNameRoot, "UrbCol.pdf"), width=7, height=7)
  argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
                 list(includeUrbanRural = FALSE, includeCluster = TRUE), 
                 list(includeUrbanRural = TRUE, includeCluster = FALSE), 
                 list(includeUrbanRural = TRUE, includeCluster = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$includeUrbanRural
    includeCluster = args$includeCluster
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster)
    out = load(paste0(nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    if(!includeCluster) {
      valMat = cbind(valMat, expit(designRes$predictions$mean))
      labels = c(labels, paste0("BYM2 ", typeText))
    }
    else {
      # also gather and plot the debiased results
      nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
      out = load(paste0(nameRoot, '.RData'))
      
      # typeText = paste0(typeText, "'")
      
      valMat = cbind(valMat, expit(designRes$predictions$mean))
      labels = c(labels, paste0("BYM2 ", typeText))
    }
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col=countyCols)
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.4, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("BYM2 ", varName, " estimate comparisons"), 
        ylim=zlim, xlim=zlim, oma=c(3,3,6,7), asp=1)
  image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=29, legend.mar=3.3, col=urbCols, add=TRUE, 
             legend.lab = "Urbanicity", legend.line=1.2, legend.width=.5, legend.shrink=.8, 
             legend.cex=.8, axis.args=list(cex.axis=.5, tck=-1, hadj=.8))
  dev.off()
  
  ## Plot 11: make pair plots of the SPDE models region estimates
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "Region.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsRegion$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col="blue")
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.8, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("SPDE region ", varName, " estimate comparisons"), 
        ylim=zlim, xlim=zlim)
  dev.off()
  
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "RegionUrbCol.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsRegion$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col=regionCols)
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.8, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("SPDE region ", varName, " estimate comparisons"), 
        ylim=zlim, xlim=zlim, oma=c(3,3,6,7), asp=1)
  image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=29, legend.mar=3.3, col=urbCols, add=TRUE, 
             legend.lab = "Urbanicity", legend.line=1.2, legend.width=.5, legend.shrink=.8, 
             legend.cex=.8, axis.args=list(cex.axis=.5, tck=-1, hadj=.8))
  dev.off()
  
  ## Plot 12: make pair plots of the SPDE models county estimates
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "County.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsCounty$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col="blue")
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("SPDE county ", varName, " estimate comparisons"), 
        ylim=zlim, xlim=zlim)
  dev.off()
  
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "CountyUrbCol.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsCounty$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col=countyCols)
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("SPDE county ", varName, " estimate comparisons"), 
        ylim=zlim, xlim=zlim, oma=c(3,3,6,7), asp=1)
  image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=29, legend.mar=3.3, col=urbCols, add=TRUE, 
             legend.lab = "Urbanicity", legend.line=1.2, legend.width=.5, legend.shrink=.8, 
             legend.cex=.8, axis.args=list(cex.axis=.5, tck=-1, hadj=.8))
  dev.off()
  
  ## Plot 13: make pair plots of the SPDE models pixel estimates
  png(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "Pixel.png"), width=1000, height=1000)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsPixel$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  urban = popGrid$urban
  my_line <- function(x,y,...){
    points(x[!urban],y[!urban],..., col="green")
    points(x[urban],y[urban],..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.1, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("Urban (blue) and rural (green) SPDE pixel ", varName, " estimate comparisons"), 
        ylim=zlim, xlim=zlim)
  dev.off()
  
  ## Plot 14: make pair plots of the SPDE models cluster estimates
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDE", resultNameRoot, "Cluster.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = FALSE, urbanEffect = TRUE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i)
    if(i == 1)
      typeText = "uc"
    else if(i == 2)
      typeText = "uC"
    else if(i == 3)
      typeText = "Uc"
    else if(i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsCluster$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  urban = dat$urban
  my_line <- function(x,y,...){
    points(x[!urban],y[!urban],..., col="green")
    points(x[urban],y[urban],..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  zlim = range(valMat)
  pairs(valMat, labels, 
        pch=19, cex=.15, lower.panel=my_line, upper.panel = my_line, 
        main=paste0("Urban/rural (blue/green) SPDE cluster ", varName, " estimates"), 
        ylim=zlim, xlim=zlim)
  dev.off()
  
  ## Plot 15: make pair plots of the SPDE models with and without urban effects at different aggregation levels
  png(file=paste0("figures/", resultNameRoot, "/pairPlotSPDEUrb", resultNameRoot, "All.png"), width=800, height=800)
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valList = list()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i*2)
    if(2*i == 1)
      typeText = "uc"
    else if(2*i == 2)
      typeText = "uC"
    else if(2*i == 3)
      typeText = "Uc"
    else if(2*i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsCluster$pred
    
    valList = c(valList, list(spdeResults$resultsCluster$pred))
    valList = c(valList, list(spdeResults$resultsPixel$pred))
    valList = c(valList, list(spdeResults$resultsCounty$pred))
    valList = c(valList, list(spdeResults$resultsRegion$pred))
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  urban = dat$urban
  my_line <- function(x,y,...){
    points(x[!urban],y[!urban],..., col="green")
    points(x[urban],y[urban],..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  zlim = range(sapply(valList, range))
  par(mfrow=c(2,2))
  clusterUrban = dat$urban
  pixelUrban = popGrid$urban
  plot(valList[[5]][!clusterUrban], valList[[1]][!clusterUrban], main=paste0("SPDE cluster ", varName, " estimates"), 
       ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], col="green", pch=19, cex=.15)
  points(valList[[5]][clusterUrban], valList[[1]][clusterUrban], col="blue", pch=19, cex=.15)
  abline(0, 1)
  legend("topleft", c("urban", "rural"), pch=19, col=c("blue", "green"))
  plot(valList[[6]][!pixelUrban], valList[[2]][!pixelUrban], main=paste0("SPDE pixel ", varName, " estimates"), 
       ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], col="green", pch=19, cex=.1)
  points(valList[[6]][pixelUrban], valList[[2]][pixelUrban], col="blue", pch=19, cex=.2)
  abline(0, 1)
  legend("topleft", c("urban", "rural"), pch=19, col=c("blue", "green"))
  plot(valList[[7]], valList[[3]], main=paste0("SPDE county ", varName, " estimates"), 
       ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], pch=19, cex=.5, col="blue")
  abline(0, 1)
  plot(valList[[8]], valList[[4]], main=paste0("SPDE region ", varName, " estimates"), 
       ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], pch=19, cex=.8, col="blue")
  abline(0, 1)
  dev.off()
  
  if(makeScreenSplitPlot) {
    png(file=paste0("figures/", resultNameRoot, "/pairPlotSPDEUrb", resultNameRoot, "AllUrbCol.png"), width=1000, height=1000)
    argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                   list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
    
    # collect the estimates and model labels
    valList = list()
    labels = c()
    for(i in 1:length(argList)) {
      args = argList[[i]]
      includeUrban = args$urbanEffect
      includeCluster = args$includeClustEffect
      clusterText = ifelse(includeCluster, "", "NoClust")
      
      nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                        "_urbanEffect", includeUrban)
      out = load(paste0("results", nameRoot, '.RData'))
      
      urbanText = ifelse(includeUrban, "", "noUrb")
      clusterText = ifelse(includeCluster, "", "NoClust")
      both = includeUrban && includeCluster
      notBothText = ifelse(both, "", " ")
      # typeText = romanText(i * 2)
      if(2*i == 1)
        typeText = "uc"
      else if(2*i == 2)
        typeText = "uC"
      else if(2*i == 3)
        typeText = "Uc"
      else if(2*i == 4)
        typeText = "UC"
      
      vals = spdeResults$resultsCluster$pred
      
      valList = c(valList, list(spdeResults$resultsCluster$pred))
      valList = c(valList, list(spdeResults$resultsPixel$pred))
      valList = c(valList, list(spdeResults$resultsCounty$pred))
      valList = c(valList, list(spdeResults$resultsRegion$pred))
      labels = c(labels, paste0("SPDE ", typeText))
    }
    
    # now construct the pair plot
    urban = dat$urban
    my_line <- function(x,y,...){
      points(x[!urban],y[!urban],..., col="green")
      points(x[urban],y[urban],..., col="blue")
      abline(a = 0,b = 1,...)
    }
    
    zlim = range(sapply(valList, range))
    allScreens = split.screen( rbind(c(0, .9,0,1), c(.9,1,0,1)))
    
    # now subdivide up the figure region into two parts
    split.screen(c(2,2), screen=allScreens[1])-> ind
    
    # first image
    screen( ind[1])
    clusterUrban = dat$urban
    pixelUrban = popGrid$urban
    plot(valList[[5]][!clusterUrban], valList[[1]][!clusterUrban], main=paste0("SPDE cluster ", varName, " estimates"), 
         ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], col=urbCols[1], pch=19, cex=.15)
    points(valList[[5]][clusterUrban], valList[[1]][clusterUrban], col=urbCols[ncols], pch=19, cex=.15)
    abline(0, 1)
    # legend("topleft", c("urban", "rural"), pch=19, col=c("blue", "green"))
    
    # second image
    screen( ind[2])
    plot(valList[[6]][!pixelUrban], valList[[2]][!pixelUrban], main=paste0("SPDE pixel ", varName, " estimates"), 
         ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], col=urbCols[1], pch=19, cex=.1)
    points(valList[[6]][pixelUrban], valList[[2]][pixelUrban], col=urbCols[ncols], pch=19, cex=.25)
    abline(0, 1)
    # legend("topleft", c("urban", "rural"), pch=19, col=c("blue", "green"))
    
    # third image
    screen( ind[3])
    plot(valList[[7]], valList[[3]], main=paste0("SPDE county ", varName, " estimates"), 
         ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], pch=19, cex=1, col=countyCols)
    abline(0, 1)
    
    # fourth image
    screen( ind[4])
    plot(valList[[8]], valList[[4]], main=paste0("SPDE region ", varName, " estimates"), 
         ylim=zlim, xlim=zlim, xlab=labels[2], ylab=labels[1], pch=19, cex=1.2, col=regionCols)
    abline(0, 1)
    
    # move to skinny region on right and draw the legend strip 
    screen( allScreens[2])
    
    # image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=29, legend.mar=3.3, col=urbCols, add=TRUE, 
    #            legend.lab = "Urbanicity", legend.line=1.2, legend.width=.5, legend.shrink=.8, 
    #            legend.cex=.8, axis.args=list(cex.axis=.5, tck=-1, hadj=.8))
    image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=29, col=urbCols, smallplot=c(0.05,.1, .1,.9), 
               legend.lab = "Urbanicity", legend.cex=1.2, graphics.reset=FALSE)
    dev.off()
    
    close.screen( all=TRUE)
    plot.new()
  }
  
  ## Plot 16: make pair plots of the SPDE models county estimates with and without urban effects versus direct estimates
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDEUrbDirect", resultNameRoot, "County.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c(directEstResults$est)
  labels = c("Direct")
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i * 2)
    if(2*i == 1)
      typeText = "uc"
    else if(2*i == 2)
      typeText = "uC"
    else if(2*i == 3)
      typeText = "Uc"
    else if(2*i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsCounty$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    points(x,y,..., col="blue")
    abline(a = 0,b = 1,...)
  }
  
  zlim = range(valMat)
  lims = c(list(range(valMat[,1])), list(range(valMat[,2])), list(range(valMat[,3])))
  # pairs(valMat, labels, 
  #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
  #       main=paste0("SPDE and direct county ", varName, " estimate comparisons"), 
  #       ylim=zlim, xlim=zlim)
  myPairs(valMat, labels, 
          pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
          main=paste0("SPDE and direct county ", varName, " estimate comparisons"), 
          lims=lims)
  dev.off()
  
  pdf(file=paste0("figures/", resultNameRoot, "/pairPlotSPDEUrbDirect", resultNameRoot, "CountyUrbCol.pdf"), width=6, height=6)
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE), 
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the estimates and model labels
  valMat = c(directEstResults$est)
  labels = c("Direct")
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i * 2)
    if(2*i == 1)
      typeText = "uc"
    else if(2*i == 2)
      typeText = "uC"
    else if(2*i == 3)
      typeText = "Uc"
    else if(2*i == 4)
      typeText = "UC"
    
    vals = spdeResults$resultsCounty$pred
    
    valMat = cbind(valMat, vals)
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  # now construct the pair plot
  my_line <- function(x,y,...){
    abline(a = 0,b = 1,...)
    points(x,y,..., col=countyCols)
  }
  
  zlim = range(valMat)
  lims = c(list(range(valMat[,1])), list(range(valMat[,2])), list(range(valMat[,3])))
  # pairs(valMat, labels, 
  #       pch=19, cex=.3, lower.panel=my_line, upper.panel = my_line, 
  #       main=paste0("SPDE and direct county ", varName, " estimate comparisons"), 
  #       ylim=zlim, xlim=zlim)
  myPairs(valMat, labels, 
          pch=19, cex=.5, lower.panel=my_line, upper.panel = my_line, 
          main=paste0("SPDE and direct county ", varName, " estimate comparisons"), 
          lims=lims, oma=c(3,3,6,7))
  image.plot(legend.only = TRUE, zlim=c(0,1), nlevel=29, legend.mar=3.3, col=urbCols, add=TRUE, 
             legend.lab = "Urbanicity", legend.line=1.2, legend.width=.5, legend.shrink=.8, 
             legend.cex=.8, axis.args=list(cex.axis=.5, tck=-1, hadj=.8))
  dev.off()
  
  ## Plot 17: make plot of SPDE differences in predictions at each pixel between full and noUrb models
  png(file=paste0("figures/", resultNameRoot, "/pairPlotSPDEUrb", resultNameRoot, "AllDiff.png"), width=600, height=600)
  argList = list(list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = FALSE),
                 list(clustDat = dat, includeClustEffect = TRUE, urbanEffect = TRUE))
  
  # collect the predictive draws and model labels
  valList = list()
  labels = c()
  for(i in 1:length(argList)) {
    args = argList[[i]]
    includeUrban = args$urbanEffect
    includeCluster = args$includeClustEffect
    clusterText = ifelse(includeCluster, "", "NoClust")
    
    nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster,
                      "_urbanEffect", includeUrban)
    out = load(paste0("results", nameRoot, '.RData'))
    
    if(i == 1)
      spdeResultsNoUrb = spdeResults
    
    urbanText = ifelse(includeUrban, "", "noUrb")
    clusterText = ifelse(includeCluster, "", "NoClust")
    both = includeUrban && includeCluster
    notBothText = ifelse(both, "", " ")
    # typeText = romanText(i * 2)
    if(2*i == 1)
      typeText = "uc"
    else if(2*i == 2)
      typeText = "uC"
    else if(2*i == 3)
      typeText = "Uc"
    else if(2*i == 4)
      typeText = "UC"
    
    
    vals = spdeResults$pixelDraws
    
    valList = c(valList, list(vals))
    labels = c(labels, paste0("SPDE ", typeText))
  }
  
  fullDraws = valList[[2]]
  noUrbDraws = valList[[1]]
  diffs = fullDraws - noUrbDraws
  y = spdeResults$resultsPixel$pred - spdeResultsNoUrb$resultsPixel$pred
  CI95s = apply(diffs, 1, quantile, probs=c(0.025, 0.975))
  
  # determine how transparent each credible interval should be
  alphaFun = function(x) {
    # make inversely proportional to the density of the points
    d = density(spdeResults$resultsPixel$pred)
    s = splinefun(d$x, d$y)
    testAlphas = 1 / s(x)
    minAlpha = min(testAlphas)
    testAlphas = testAlphas * (0.005 / minAlpha)
    # testAlphas[testAlphas >= .1] = .1
    testAlphas
  }
  alphaFun2 = function(x) {
    # make proportional to the inversely density of the points * the standard deviation of differences
    d = density(spdeResults$resultsPixel$pred)
    s1 = splinefun(d$x, d$y)
    
    # estimate the standard deviation
    sqDiffs = diff(sort(y[!popGrid$urban]))^2
    logSqRes = log(sqDiffs)
    
    # fit g from Wakefield Eq (11.51) using GAM
    xRural = spdeResults$resultsPixel$pred[!popGrid$urban]
    xRural = xRural[2:length(xRural)]
    mod = gam(logSqRes ~ s(xRural, bs="cr"))
    
    # get SD estimates for data locations
    gHats = predict(mod, list(xRural=x), type="response")
    sigmaHats = sqrt(exp(gHats))
    
    testAlphas = 1 / s1(x) * sigmaHats^0.5
    minAlpha = min(testAlphas)
    testAlphas = testAlphas * (0.005 / minAlpha)
    # testAlphas[testAlphas >= .1] = .1
    testAlphas
  }
  alphas = alphaFun2(spdeResults$resultsPixel$pred)
  
  # randomly sample pixels depending on population density
  sampleI = sample(1:ncol(CI95s), size=1000, prob = alphas / sum(alphas), replace = FALSE)
  y = y[sampleI]
  CI95s = CI95s[,sampleI]
  xsAll = spdeResults$resultsPixel$pred[sampleI]
  pixelUrban = popGrid$urban[sampleI]
  zlim = range(CI95s)
  
  # make the plot
  plot(xsAll, CI95s[1,], main=paste0("Difference in pixel level ", varName, " predictions (SPDE UC - SPDE uC)"),
       ylim=zlim, xlab=labels[1], ylab="Prediction differences (SPDE UC - SPDE uC)", type="n")
  
  # plot every pixel CI as a translucent line
  for(i in 1:ncol(CI95s)) {
    xs = c(xsAll[i], xsAll[i])
    ys = CI95s[,i]
    if(pixelUrban[i]) {
      thisColAlpha = rgb(0, 0, 1, 0.1)
      thislwd = 10
    }
    else {
      # if(xs[1] <= 0.15) {
      #   thisColAlpha = rgb(0, 1, 0, 0.005)
      # } else {
      #   thisColAlpha = rgb(0, 1, 0, 0.1)
      # }
      # thisColAlpha = rgb(0, 1, 0, alphas[i])
      thisColAlpha = rgb(0, 1, 0, 0.1)
      thislwd = 10
    }
    
    # plot the CI
    lines(xs, ys, col=thisColAlpha, lwd=thislwd)
  }
  
  # plot every pixel central estimate as a point
  for(i in 1:ncol(CI95s)) {
    xs = c(xsAll[i], xsAll[i])
    ys = CI95s[,i]
    if(pixelUrban[i]) {
      thisCol = rgb(0, 0, .7)
    }
    else {
      thisCol = rgb(0, .7, 0)
    }
    
    # plot the central estimates
    points(xs[1], y[i], col=thisCol, pch=19, cex=.1)
  }
  abline(0, 0, lty=2)
  legend("topleft", c("urban", "rural"), pch=19, col=c("blue", "green"))
  dev.off()
}

# produce plot requested by jon of continuous and discrete level predictions
makeJonPlot = function(dat=ed, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2, 
                       meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2, 
                       sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE, 
                       meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND, 
                       meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2, 
                       varName="SEP", plotNameRoot="Education", resultNameRoot="Ed", meanCols=makeRedBlueDivergingColors(64), 
                       sdCols=makeBlueYellowSequentialColors(64), popCols=makeBlueSequentialColors(64), 
                       ncols=29, relativeCols=makeRedGreenDivergingColors(ncols), urbCols=makeGreenBlueSequentialColors(ncols), 
                       plotUrbanMap=FALSE, kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0), makeScreenSplitPlot=FALSE) {
  plotNameRootLower = tolower(plotNameRoot)
  resultNameRootLower = tolower(resultNameRoot)
  
  png(file=paste0("figures/", resultNameRoot, "/JonPlot", plotNameRoot, ".png"), width=500, height=900)
  par(mfrow=c(2,1))
  includeUrban = FALSE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2', resultNameRoot, 'UrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
  out = load(paste0(nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  typeTextBYM = "Discrete model"
  
  includeUrban = FALSE
  includeCluster = TRUE
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDE", resultNameRootLower, "_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  urbanText = ifelse(includeUrban, "", "noUrb")
  clusterText = ifelse(includeCluster, "", "NoClust")
  both = includeUrban && includeCluster
  notBothText = ifelse(both, "", " ")
  typeTextSPDE = "Continuous model"
  
  zlim = range(c(expit(designRes$predictions$mean), spdeResults$resultsPixel$pred))
  # plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("Discrete ", varName, " estimates"), cols=meanCols, zlim=logit(meanRange2), scaleFun=logit, scaleFunInverse=expit, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  plotMapDat(adm1, plotVar=expit(designRes$predictions$mean), new = TRUE, main=paste0("Discrete Model"), cols=meanCols, zlim=zlim, ticks=meanTicks2, tickLabels=meanTickLabels2, xlim=kenyaLonRange, ylim=kenyaLatRange)
  
  # par(mfrow=c(2,2), oma=c( 0,0,0,1.5), mar=c(5.1, 4.1, 4.1, 6))
  # plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("SPDE ", typeText, " ", varName, " estimates"), ylim=kenyaLatRange, 
  #      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  # quilt.plot(cbind(popGrid$lon, popGrid$lat), logit(spdeResults$resultsPixel$pred), 
  #            nx=150, ny=150, add.legend=FALSE, add=TRUE, col=meanCols, zlim=range(logit(meanRange2)))
  # plotMapDat(adm1, lwd=.5)
  # points(dat$lon, dat$lat, pch=".")
  # image.plot(zlim=range(logit(meanRange2)), nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
  #            col=meanCols, add = TRUE, axis.args=list(at=logit(meanTicks2), labels=meanTickLabels2), legend.mar = 0)
  plot(cbind(popGrid$lon, popGrid$lat), type="n", main=paste0("Continuous Model"), ylim=kenyaLatRange, 
       xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", asp=1)
  # quilt.plot(cbind(popGrid$lon, popGrid$lat), spdeResults$resultsPixel$pred, 
  #            nx=150, ny=150, add=TRUE, col=meanCols, zlim=meanRange2, 
  #            axis.args=list(at=meanTicks2, labels=meanTickLabels2))
  quilt.plot(cbind(popGrid$lon, popGrid$lat), spdeResults$resultsPixel$pred, 
             nx=150, ny=150, add=TRUE, col=meanCols, zlim=zlim, legend.mar=7, 
             axis.args=list(at=meanTicks2, labels=meanTickLabels2))
  plotMapDat(adm1, lwd=.5)
  dev.off()
}

romanText = function(number) {
  if(number == 1) {
    "I"
  } else if(number == 2) {
    "II"
  } else if(number == 3) {
    "III"
  } else if(number == 4){
    "IV"
  } else {
    stop(paste0("unsupported roman-numeral number: ", number))
  }
}