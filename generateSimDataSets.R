# load("simDat.RData")
# library(profvis)
# library(logitnorm)
source("setup.R")
# setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")

##### the code below does not use the same enumeration areas for each simulation, 
##### which was why was commented out


# 
# # number of datasets to be generated
# numData <- 100
# 
# set.seed(580252)
# my.seeds <- round(runif(numData)*100000)
# 
# dataSets = list()
# for(i in 1:numData){
#   print(i)
#   beta0 = -2
#   margVar = .15^2
#   tausq = .1^2
#   gamma = -.5
#   HHoldVar = 0
#   # HHoldVar = .3^2
#   tmp <- simDat(kenyaDat, beta0=beta0, margVar=margVar, tausq=tausq, gamma=gamma, seed=my.seeds[i])
#   # add the vector of sampling weights to the dataset 
#   #samplingWeight = 1/(table(tmp$clustDat$admin1)/table(tmp$eaDat$admin1))
#   #tmp$clustDat$samplingWeight = samplingWeight
#   dataSets[[i]] = tmp
# }
# 
# save(dataSets, file=paste0("simData4analysisBeta", round(beta0, 2), "margVar", round(margVar, 2), "tausq", 
#                            round(tausq, 2), "gamma", round(gamma, 2), "HHoldVar", HHoldVar, ".RData"))
# 

## CAUTION!!!!
#Warning messages:
# 1: In doTryCatch(return(expr), name, parentenv, handler) :
#  restarting interrupted promise evaluation
#2: In doTryCatch(return(expr), name, parentenv, handler) :
#  restarting interrupted promise evaluation



# unlike the above script, this script holds the EA data for each simulation fixed, 
# each simulation instead varying which clusters were sampled
# urbanOverSample: within any county, any individual urban EA is urbanOverSample times 
#                  as likely to be sampled than any individual rural EA to be a cluster
# nsim=100

# generate the empirical distributions and save them
# wd = getwd()
# setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
# empiricalDistributions = getSurveyEmpiricalDistributions2()
# save(empiricalDistributions, file="empiricalDistributions.RData")
# save(empiricalDistributions, file="~/git/U5MR/empiricalDistributions.RData")
# setwd(wd)


# simulate and save datasets used for the simulation study with the given model parameters
# NOTE: paired with the dataset using the passed parameters will be another dataset from the 
#       same model without a nugget/cluster effect
# nsim: number of surveys taken from the true latent population in the standard size survey collections
# nsimBig: number of surveys taken from the true latent population in the large size survey collections
# seeds: random number seeds used for making the latent population and generating surveys respectively
# beta0: latent gaussian model intercept
# margVar: marginal variance of the spatial field
# tausq: the nugget/cluster effect variance
# gamma: latent gaussian model urban effect
# HHoldVar: household effect variance
# effRange: spatial range
# urbanOverSamplefrac: the proportion with which to inflate the amount of urban samples in the surveys
generateSimDataSets = function(nsim=100, nsimBig = 250, seeds=c(580252, 1234), margVar=.15^2,
                               gamma = -1, HHoldVar = 0, effRange = 150, 
                               urbanOverSamplefrac = 0) {
  beta0 = -1.75
  tausq = .1^2
  set.seed(seeds[1])
  wd = getwd()
  setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
  
  rangeText = ""
  if(effRange != 150)
    rangeText = paste0("Range", effRange)
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                  round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                  "urbanOverSamplefrac", round(urbanOverSamplefrac, 4), rangeText)
  dataID0 = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                   round(0, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                   "urbanOverSamplefrac", round(urbanOverSamplefrac, 4), rangeText)
  
  # there should be 1 true data set, but many simulated cluster samples
  load("empiricalDistributions.RData")
  simulatedEAs = simDatEmpirical(empiricalDistributions, kenyaEAs, clustDat=NULL, nsim=1, 
                                 beta0=beta0, margVar=margVar, urbanOverSamplefrac=urbanOverSamplefrac, 
                                 tausq=tausq, gamma=gamma, HHoldVar=HHoldVar, effRange=effRange)
  kenyaEAs = simulatedEAs$eaDat
  kenyaEAs$eaIs = 1:nrow(kenyaEAs)
  kenyaEAsLong = kenyaEAs[rep(1:nrow(kenyaEAs), kenyaEAs$nHH),]
  
  set.seed(seeds[2])
  # simulate the cluster sampling and add to the data sets
  overSampClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, verbose=FALSE)
  clustList = genAndreaFormatFromEAIs(simulatedEAs$eaDat, overSampClustDat$eaIs, overSampClustDat$sampleWeights)
  overSampDat = list(eaDat=kenyaEAs, clustDat=clustList)
  
  overSampClustDatTest = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, urbanOverSamplefrac, fixedPerStrata=TRUE, nPerStrata=3, verbose=FALSE)
  clustListTest = genAndreaFormatFromEAIs(kenyaEAs, overSampClustDatTest$eaIs, overSampClustDatTest$sampleWeights)
  overSampDatTest = list(eaDat=kenyaEAs, clustDat=clustListTest)
  
  SRSClustDat = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, SRS=TRUE, verbose=FALSE)
  clustList = genAndreaFormatFromEAIs(kenyaEAs, SRSClustDat$eaIs, SRSClustDat$sampleWeights)
  SRSDat = list(eaDat=kenyaEAs, clustDat=clustList) # the only thing different is the sampling of the clusters
  
  SRSClustDatTest = simClustersEmpirical(kenyaEAs, kenyaEAsLong, nsimBig, NULL, fixedPerStrata=TRUE, nPerStrata=3, SRS=TRUE, verbose=FALSE)
  clustListTest = genAndreaFormatFromEAIs(kenyaEAs, SRSClustDatTest$eaIs, SRSClustDatTest$sampleWeights)
  SRSDatTest = list(eaDat=kenyaEAs, clustDat=clustListTest) # the only thing different is the sampling of the clusters
  
  # plot the first simulation of the over sampled and simple random sample data sets
  clustDat = SRSDat$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  eaDat = overSampDat$eaDat
  pdf(paste0("figures/exampleSRSSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "Big.RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "TestBig.RData"))
  out = load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # Now take only the first nsim simulations from the "big" dataset
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, ".RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "Test.RData"))
  
  # reload the data
  out = load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # Now simulate the data without a cluster effect but with the same underlying probability surface otherwise
  tausq = 0
  overSampDat$eaDat$trueProbDeath = overSampDat$eaDat$trueProbDeathNoNug
  SRSDat$eaDat$trueProbDeath = SRSDat$eaDat$trueProbDeathNoNug
  overSampDat$eaDat$died = rbinom(nrow(overSampDat$eaDat), overSampDat$eaDat$numChildren, overSampDat$eaDat$trueProbDeathNoNug)
  SRSDat$eaDat$died = overSampDat$eaDat$died
  
  overSampDatTest$eaDat$trueProbDeath = overSampDatTest$eaDat$trueProbDeathNoNug
  SRSDatTest$eaDat$trueProbDeath = SRSDatTest$eaDat$trueProbDeathNoNug
  overSampDatTest$eaDat$died = rbinom(nrow(overSampDatTest$eaDat), overSampDatTest$eaDat$numChildren, overSampDatTest$eaDat$trueProbDeathNoNug)
  SRSDatTest$eaDat$died = overSampDatTest$eaDat$died
  for(i in 1:nsimBig) {
    overSampDat$clustDat[[i]]$trueProbDeath = overSampDat$clustDat[[i]]$trueProbDeathNoNug
    SRSDat$clustDat[[i]]$trueProbDeath = SRSDat$clustDat[[i]]$trueProbDeathNoNug
    overSampDatTest$clustDat[[i]]$trueProbDeath = overSampDatTest$clustDat[[i]]$trueProbDeathNoNug
    SRSDatTest$clustDat[[i]]$trueProbDeath = SRSDatTest$clustDat[[i]]$trueProbDeathNoNug
    
    overSampDat$clustDat[[i]]$died = overSampDat$eaDat$died[overSampClustDat$eaIs[,i]]
    SRSDat$clustDat[[i]]$died = SRSDat$eaDat$died[SRSClustDat$eaIs[,i]]
    overSampDatTest$clustDat[[i]]$died = overSampDatTest$eaDat$died[overSampClustDatTest$eaIs[,i]]
    SRSDatTest$clustDat[[i]]$died = SRSDatTest$eaDat$died[SRSClustDatTest$eaIs[,i]]
  }
  
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, "Big.RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, "TestBig.RData"))
  load(paste0("simDataMulti", dataID0, "Big.RData"))
  
  # Again, take only the first nsim simulations from the "big"" dataset
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, ".RData"))
  overSampDat = overSampDatTest
  SRSDat = SRSDatTest
  overSampDat$clustDat = overSampDat$clustDat[1:nsim]
  SRSDat$clustDat = SRSDat$clustDat[1:nsim]
  save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID0, "Test.RData"))
  
  # clustDat = SRSDat$clustDat[[1]]
  # clustDat = SRSDatTest$clustDat[[1]]
  # clustDat = overSampDat$clustDat[[1]]
  clustDat = overSampDatTest$clustDat[[1]]
  eaDat = overSampDat$eaDat
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  
  pdf(paste0("figures/exampleOverSampTestSimulationNoNug", dataID0, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  setwd(wd)
  
  invisible(NULL)
}

generateAllDataSets = function() {
  generateSimDataSets(gamma=-1, margVar=.15^2)
  generateSimDataSets(gamma=-1, margVar=0)
  generateSimDataSets(gamma=0, margVar=.15^2)
  generateSimDataSets(gamma=0, margVar=0)
}

generateAllNewDataSets = function() {
  generateSimDataSets(gamma=-1, margVar=.15^2, effRange=50)
  generateSimDataSets(gamma=0, margVar=.15^2, effRange=50)
  generateSimDataSets(gamma=0, margVar=0, effRange=50)
  
  generateSimDataSets(gamma=-1, margVar=.3^2, effRange=150)
  generateSimDataSets(gamma=0, margVar=.3^2, effRange=150)
}

# simulate and save datasets used for the simulation study with the given model parameters
# NOTE: paired with the dataset using the passed parameters will be another dataset from the 
#       same model without a nugget/cluster effect
# nsim: number of surveys taken from the true latent population in the standard size survey collections
# nsimBig: number of surveys taken from the true latent population in the large size survey collections
# seeds: random number seeds used for making the latent population and generating surveys respectively
# beta0: latent gaussian model intercept
# margVar: marginal variance of the spatial field
# tausq: the nugget/cluster effect variance
# gamma: latent gaussian model urban effect
# HHoldVar: household effect variance
# effRange: spatial range
# urbanOverSamplefrac: the proportion with which to inflate the amount of urban samples in the surveys
plotSimDataSets = function(nsim=100, nsimBig = 250, seeds=c(580252, 1234), beta0 = -1.75, margVar = .15^2, 
                           tausq = .1^2, gamma = -1, HHoldVar = 0, effRange = 150, 
                           urbanOverSamplefrac = 0, colorScale=makeRedBlueDivergingColors(64, rev = TRUE), 
                           kenyaLatRange=c(-4.6, 5), kenyaLonRange=c(33.5, 42.0)) {
  set.seed(seeds[1])
  wd = getwd()
  setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
  
  rangeText = ""
  if(effRange != 150)
    rangeText = paste0("Range", effRange)
  
  # make strings representing the simulation with and without cluster effects
  dataID = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                  round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                  "urbanOverSamplefrac", round(urbanOverSamplefrac, 4), rangeText)
  dataID0 = paste0("Beta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
                   round(0, 4), "gamma", round(gamma, 4), "HHoldVar", HHoldVar, 
                   "urbanOverSamplefrac", round(urbanOverSamplefrac, 4), rangeText)
  
  # there should be 1 true data set, but many simulated cluster samples
  load("empiricalDistributions.RData")
  
  # save(overSampDat, SRSDat, file=paste0("simDataMulti", dataID, "Big.RData"))
  load(paste0("simDataMulti", dataID, "Big.RData"))
  
  # plot the first simulation of the over sampled and simple random sample data sets
  clustDat = SRSDat$clustDat[[1]]
  clustDat = overSampDat$clustDat[[1]]
  eaDat = overSampDat$eaDat
  pdf(paste0("figures/exampleSRSSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  png(file=paste0("figures/exampleOverSampSimulation", dataID, ".png"), width=900, height=500)
  par(oma=c( 0,0,0,2), mar=c(5.1, 4.1, 4.1, 6), mfrow =c(1, 2))
  
  plot(cbind(eaDat$lon, eaDat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
       xlab="Longitude", ylab="Latitude", main=paste0("Simulated NMRs"), asp=1)
  quilt.plot(eaDat$lon, eaDat$lat, eaDat$died/eaDat$numChildren, nx=150, ny=150, col=colorScale, zlim=zlim, add=TRUE)
  # world(add=TRUE)
  plotMapDat(adm1)
  
  plot(cbind(clustDat$lon, clustDat$lat), type="n", ylim=kenyaLatRange, xlim=kenyaLonRange, 
       xlab="Longitude", ylab="Latitude", main=paste0("Example Survey NMRs"), asp=1)
  quilt.plot(clustDat$lon, clustDat$lat, clustDat$died/clustDat$numChildren, nx=150, ny=150, col=colorScale, zlim=zlim, add=TRUE)
  # world(add=TRUE)
  plotMapDat(adm1)
  dev.off()
  
  clustDat = overSampDat$clustDat[[1]]
  pdf(paste0("figures/exampleUrbanOverSampSimulation", dataID, ".pdf"), width=8, height=8)
  par(mfrow =c(2, 2))
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = clustDat$numChildren
  obsCounts = clustDat$died
  zlim = c(0, quantile(c(eaDat$died/eaDat$numChildren, clustDat$died/clustDat$numChildren, 
                         eaDat$trueProbDeath), probs=.975))
  quilt.plot(eaDat$east, eaDat$north, eaDat$died/eaDat$numChildren, main="All Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, obsCounts/obsNs, main="Sample Empirical Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(eaDat$east, eaDat$north, eaDat$trueProbDeath, main="All True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, clustDat$trueProbDeath, main="Sample True Mortality Rates", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  dev.off()
  
  
  setwd(wd)
  
  invisible(NULL)
}

# function for calculating the average number of urban and rural clusters per simulated dataset
getAverageUrbanRuralClusters = function() {
  out = system("ls ~/Google\\ Drive/UW/Wakefield/WakefieldShared/U5MR/simDataMultiBeta*0Big.RData", intern=TRUE)
  numUrban = c()
  numRural = c()
  for(i in 1:length(out)) {
    print(paste0("Loading ", out[i], "... (", i, "/", length(out), ")"))
    load(out[i])
    numUrban = c(numUrban, sapply(SRSDat$clustDat, function(x){sum(x$urban)}))
    numRural = c(numRural, sapply(SRSDat$clustDat, function(x){sum(!x$urban)}))
  }
  
  print(paste0("Average number of urban clusters: ", mean(numUrban)))
  print(paste0("Average number of rural clusters: ", mean(numRural)))
  
  invisible(list(numUrban=numUrban, numRural=numRural))
}


