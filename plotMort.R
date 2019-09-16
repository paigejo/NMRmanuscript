# script for plotting predictions for neonatal mortality in Kenya
source("plotGenerator.R")

##### before we make any plots, put all of them on the same scale
##### make multiple scales, two for estimates and quantiles, and two
##### for standard deviations. For each type of scale, make one that 
##### includes the naive and direct estimates, on anther that does not
# naive and direct
out = load("resultsDirectNaiveMort.RData")
meanRange = range(c(expit(naiveResults$upper),expit(naiveResults$lower) ))
meanRange = range(c(meanRange, expit(naiveResults$upper),expit(naiveResults$lower)))
zlim = range(c(expit(directEstResults$upper),expit(directEstResults$lower) ))
meanRange = range(c(meanRange, zlim))
sdRange = range(sqrt(naiveResults$var.est))
zlim = range(sqrt(directEstResults$var.est))
sdRange = range(c(sdRange, zlim))
zlim = range(c(expit(directEstResults$lower)-expit(directEstResults$upper) ))
widthRange = zlim

# make scales specialized for the direct and naive models
meanRangeND = meanRange
sdRangeND = sdRange
sdTicksND = pretty(sdRangeND)
sdTickLabelsND = as.character(sdTicksND)
meanTicksND = c(0.005, 0.01, pretty(c(.01, meanRangeND[2]), n=5))
meanTicksND = meanTicksND[-which(meanTicksND == 0)]
meanTickLabelsND = as.character(meanTicksND)

# mercer
out = load("resultsMercerMort.RData")
zlim = range(c(expit(mercerResults$lower.mercer),expit(mercerResults$upper.mercer)))
meanRange = range(c(meanRange, zlim))
sdRange = range(c(sdRange, sqrt(mercerResults$var.est.mercer)))
meanRange2 = zlim
sdRange2 = range(sqrt(mercerResults$var.est.mercer))
widthRange = range(c(widthRange, range(expit(mercerResults$upper.mercer)-expit(mercerResults$lower.mercer))))

# bym2
argList = list(list(includeUrbanRural = FALSE, includeCluster = FALSE), 
               list(includeUrbanRural = FALSE, includeCluster = TRUE), 
               list(includeUrbanRural = TRUE, includeCluster = FALSE), 
               list(includeUrbanRural = TRUE, includeCluster = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$includeUrbanRural
  includeCluster = args$includeCluster
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster)
  out = load(paste0(nameRoot, '.RData'))
  zlim = range(c(expit(designRes$predictions$Q10),expit(designRes$predictions$Q90)))
  meanRange = range(c(meanRange, zlim))
  sdRange = range(c(sdRange, designRes$predictions$stddev))
  meanRange2 = range(c(meanRange2, zlim))
  sdRange2 = range(c(sdRange2, designRes$predictions$stddev))
  widthRange = range(c(widthRange, expit(designRes$predictions$Q90)-expit(designRes$predictions$Q10)))
  
  if(i==1) {
    meanRangeBYM2 = zlim
    sdRangeBYM2 = range(designRes$predictions$stddev)
  } else {
    meanRangeBYM2 = range(c(meanRangeBYM2, zlim))
    sdRangeBYM2 = range(c(sdRangeBYM2, designRes$predictions$stddev))
  }
  
  if(includeCluster) {
    # also gather and plot the debiased results
    nameRoot = paste0('bym2MortUrbRur',includeUrban, 'Cluster', includeCluster, "debiased")
    out = load(paste0(nameRoot, '.RData'))
    
    zlim = range(c(expit(designRes$predictions$Q10),expit(designRes$predictions$Q90)))
    meanRange = range(c(meanRange, zlim))
    sdRange = range(c(sdRange, designRes$predictions$stddev))
    meanRange2 = range(c(meanRange2, zlim))
    sdRange2 = range(c(sdRange2, designRes$predictions$stddev))
    widthRange = range(c(widthRange, expit(designRes$predictions$Q90)-expit(designRes$predictions$Q10)))
    
    meanRangeBYM2 = range(c(meanRangeBYM2, zlim))
    sdRangeBYM2 = range(c(sdRangeBYM2, designRes$predictions$stddev))
  }
}

# spde
argList = list(list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = FALSE, urbanEffect = TRUE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = FALSE), 
               list(clustDat = mort, includeClustEffect = TRUE, urbanEffect = TRUE))
for(i in 1:length(argList)) {
  args = argList[[i]]
  includeUrban = args$urbanEffect
  includeCluster = args$includeClustEffect
  clusterText = ifelse(includeCluster, "", "NoClust")
  
  nameRoot = paste0("SPDEmort_includeClustEffect", includeCluster, 
                    "_urbanEffect", includeUrban)
  out = load(paste0("results", nameRoot, '.RData'))
  
  zlim = range(c(spdeResults$resultsCounty$lower,spdeResults$resultsCounty$upper))
  meanRange = range(c(meanRange, zlim))
  sdRange = range(c(sdRange, spdeResults$resultsCounty$sds))
  meanRange2 = range(c(meanRange2, zlim))
  sdRange2 = range(c(sdRange2, spdeResults$resultsCounty$sds))
  widthRange = range(c(widthRange, spdeResults$resultsCounty$upper-spdeResults$resultsCounty$lower))
  
  # get a range just for the SPDE continuous
  zlim = range(c(spdeResults$resultsPixel$lower,spdeResults$resultsPixel$upper))
  if(i==1) {
    meanRangeSPDE = zlim
    sdRangeSPDE = range(spdeResults$resultsPixel$sds)
    widthRangeSPDE = range(spdeResults$resultsPixel$upper-spdeResults$resultsPixel$lower)
  } else {
    meanRangeSPDE = range(c(meanRangeSPDE, zlim))
    sdRangeSPDE = range(c(sdRangeSPDE, spdeResults$resultsPixel$sds))
    widthRangeSPDE = range(c(widthRangeSPDE, spdeResults$resultsPixel$upper-spdeResults$resultsPixel$lower))
  }
}

# set plot tick marks to be reasonable on logit and log scales
sdTicks = pretty(c(.1, sdRange[2]))
sdTickLabels = as.character(sdTicks)
# sdTickLabels[c(5, 7)] = ""
meanTicks = pretty(c(.01, meanRange[2]), n=10)
meanTicks = meanTicks[-c(5, 7, 9, 11, 13)]
meanTickLabels = as.character(meanTicks)
widthTicks = pretty(widthRange, n=10)
widthTickLabels = as.character(widthTicks)
# meanTickLabels[c(5, 7, 9, 11, 13)] = ""
sdTicks2 = pretty(sdRange2)
sdTickLabels2 = as.character(sdTicks2)
meanTicks2 = pretty(meanRange2, n=10)
meanTicks2 = c(0.01, meanTicks2[-1])
meanTickLabels2 = as.character(meanTicks2)
meanTicksSPDE = pretty(meanRangeSPDE, n=11)
meanTicksSPDE = meanTicksSPDE[-1]
sdTicksSPDE = pretty(sdRangeSPDE, n=10)
sdTicksSPDE = c(0.05, sdTicksSPDE[-1])
meanTickLabelsSPDE = as.character(meanTicksSPDE)
sdTickLabelsSPDE = as.character(sdTicksSPDE)
meanTicksBYM2 = pretty(meanRangeBYM2, n=5)
sdTicksBYM2 = pretty(sdRangeBYM2, n=5)
meanTickLabelsBYM2 = as.character(meanTicksBYM2)
sdTickLabelsBYM2 = as.character(sdTicksBYM2)
widthTicksSPDE = pretty(widthRangeSPDE, n=10)
widthTickLabelsSPDE = as.character(widthTicksSPDE)

# add in a few extra tick marks
meanTicks = c(.01, meanTicks)
meanTickLabels = c("0.01", meanTickLabels)
sdTicks = c(0.05, sdTicks)
sdTickLabels = c("0.05", sdTickLabels)
sdTicks2 = c(0.005, 0.01, 0.05, sdTicks2)
sdTickLabels2 = c("0.005", "0.01", "0.05", sdTickLabels2)
meanTicksSPDE = c(0.005, 0.01, 0.05, meanTicksSPDE)
meanTickLabelsSPDE = as.character(meanTicksSPDE)

makeAllPlots(mort, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2,
             meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2,
             sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE,
             meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND,
             meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2,
             varName="NMR", plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)),
             relativeCols=rev(makeRedGreenDivergingColors(29)), plotUrbanMap=FALSE, makeScreenSplitPlot=TRUE, 
             sharedPredictionScale=TRUE, widthRange=widthRange, 
             widthTicks=widthTicks, widthTickLabels=widthTickLabels, widthRangeSPDE=widthRangeSPDE, 
             widthTicksSPDE=widthTicksSPDE, widthTickLabelsSPDE=widthTickLabelsSPDE)

printModelPredictionTables(mort, resultNameRoot="Mort", nDigitsPredictions=2, nDigitsSmoothedDirect=3)
# Browse[2]> diff(range(tab[,7]))
# [1] 0.008698917
# Browse[2]> median(tab[,9] - tab[,8])
# [1] 0.006962701
# plotModelPredictions(mort, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2,
#                      meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2,
#                      sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE,
#                      meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND,
#                      meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2,
#                      varName="NMR", plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)),
#                      relativeCols=rev(makeRedGreenDivergingColors(29)), plotUrbanMap=FALSE, widthRange=widthRange, 
#                      widthTicks=widthTicks, widthTickLabels=widthTickLabels, widthRangeSPDE=widthRangeSPDE, 
#                      widthTicksSPDE=widthTicksSPDE, widthTickLabelsSPDE=widthTickLabelsSPDE)

# makePairPlots(mort, meanRange, meanRange2, meanTicks, meanTicks2, meanTickLabels, meanTickLabels2,
#               meanRangeSPDE, meanTicksSPDE, meanTickLabelsSPDE, sdRange, sdRange2,
#               sdTicks, sdTicks2, sdTicksSPDE, sdTickLabels, sdTickLabels2, sdTickLabelsSPDE,
#               meanRangeND, meanTicksND, meanTickLabelsND, sdRangeND, sdTicksND, sdTickLabelsND,
#               meanRangeBYM2, meanTicksBYM2, meanTickLabelsBYM2, sdTicksBYM2, sdTickLabelsBYM2,
#               varName="NMR", plotNameRoot="Mort", resultNameRoot="Mort", meanCols=rev(makeRedBlueDivergingColors(64)),
#               relativeCols=rev(makeRedGreenDivergingColors(29)), plotUrbanMap=FALSE, makeScreenSplitPlot=FALSE)

