# simulation study scripts for UM5R study

# simulate enumeration areas from population data
simEAs = function(kenyaPop, numEAs=96251, totalKenyaPop=43.0 * 10^6, seed=123) {
  set.seed(seed)
  
  # sample without replacement proportionally to population
  EAIs = sample(1:nrow(kenyaPop), size=numEAs, replace=F, prob=kenyaPop$pop/sum(kenyaPop$pop))
  kenyaEAs = kenyaPop[EAIs,]
  decreaseFac = totalKenyaPop/sum(kenyaEAs$pop)
  kenyaEAs$pop = kenyaEAs$pop*decreaseFac
  # eaReg = getRegion(cbind(kenyaEAs$lon, kenyaEAs$lat), adm0)
  eaCountyDat = getRegion(cbind(kenyaEAs$lon, kenyaEAs$lat), adm1)
  eaCounty = eaCountyDat$regionID
  eaCountyNames = eaCountyDat$regionNames
  eaCountyNames[5] = "Elgeyo Marakwet"
  eaCountyNames[42] = "Trans-Nzoia"
  mapping = match(eaCountyNames, as.character(unique(mort$admin1)))
  eaCounty = unique(mort$admin1)[mapping[eaCounty]]
  kenyaEAs$admin1 = eaCounty
  kenyaEAs$region = countyToRegion(eaCounty)
  kenyaEAs$eaIs = EAIs
  
  kenyaEAs
}

# for kmres=2.5: nx=350, ny=430
makeUrbanMap = function(popGrid=NULL, kmres=2.5, savePlot=FALSE, fileName=ifelse(whiteRural, "figures/UrbanMapWhiteRural.png", "figures/urbanMap.png"), 
                        nx=850, ny=1050, width=500, height=500, lonLim=kenyaLonRange, latLim=kenyaLatRange, whiteRural=TRUE) {
  # get prediction locations from population grid
  if(is.null(popGrid)) {
    if(kmres == 5)
      load("popGrid.RData")
    else
      popGrid = makeInterpPopGrid(kmres)
  }
  
  # determine which points in Kenya are urban
  threshes = setThresholds2()
  popThreshes = sapply(1:nrow(popGrid), function(i) {threshes$threshes[threshes$counties == popGrid$admin1[i]]})
  urban = popGrid$popOrig > popThreshes
  
  if(savePlot) {
    png(file=fileName, width=width, height=height)
    par(oma=c( 0,0,0,0), mar=c(5.1, 4.1, 4.1, 4.1))
  }
  plot(popGrid$lon, popGrid$lat, xlab="Longitude", ylab="Latitude", main=TeX("Urbanicity"), xlim=lonLim, ylim=latLim, asp=1, type="n")
  # quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, 
  #            xlab="Longitude", ylab="Latitude", main=TeX("Urbanicity"), xlim=lonLim, ylim=latLim, asp=1)
  if(whiteRural)
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c(rgb(0, 0, 0, 0), "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE)
  else {
    quilt.plot(popGrid$lon, popGrid$lat, urban, col=c("green", "blue"), nx=850, ny=1050, add.legend = FALSE, add=TRUE)
    legend("bottomleft", c("urban", "rural"), col=c("blue", "green"), pch=19)
  }
  # world(add=TRUE)
  plotMapDat(adm1)
  if(savePlot) {
    dev.off()
  }
}

# simulate enumeration areas from population data.  stratified by urban 
# and rural and county
simEAs2 = function(kenyaPop, numEAs=96251, totalKenyaPop=43.0 * 10^6, seed=123, 
                   sampleByPop=TRUE, fixNumUrbanAtTruth=TRUE, thisEasPC=easpc) {
  set.seed(seed)
  
  # determine which points in Kenya are urban
  threshes = setThresholds2()
  popThreshes = sapply(1:nrow(kenyaPop), function(i) {threshes$threshes[threshes$counties == kenyaPop$admin1[i]]})
  urban = kenyaPop$popOrig > popThreshes
  
  counties = unique(kenyaPop$admin1)
  allSamples = c()
  allUrban = c()
  for(i in 1:length(counties)) {
    countyName = counties[i]
    countyI = kenyaPop$admin1 == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # urban samples
    thisCountyI = countyI & urban
    if(sampleByPop)
      probs = thisCountyI * kenyaPop$popOrig
    else
      probs = thisCountyI
    probs = probs/sum(probs)
    
    # number of rural and urban clusters to sample: either fix it or choose it 
    # randomly based on sample probabilities
    if(fixNumUrbanAtTruth) {
      numUrban = easpc$EAUrb[easpc$County == countyName]
      numRural = easpc$EARur[easpc$County == countyName]
    }
    else {
      numSamplesTotal = easpc$EATotal[easpc$County == countyName]
      countyProbs = countyI
      if(sampleByPop)
        countyProbs = countyProbs * kenyaPop$popOrig
      countyProbs = countyProbs/sum(countyProbs)
      sampleIs = sample(1:nrow(kenyaPop), numSamplesTotal, TRUE, prob=countyProbs)
      numUrban = sum(urban[sampleIs])
      numRural = numSamplesTotal - numUrban
    }
    allSamples = c(allSamples, sample(1:nrow(kenyaPop), numUrban, TRUE, prob=probs))
    allUrban = c(allUrban, rep(TRUE, numUrban))
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      if(sampleByPop)
        probs = thisCountyI * kenyaPop$popOrig
      else
        probs = thisCountyI
      probs = probs/sum(probs)
      allSamples = c(allSamples, sample(1:nrow(kenyaPop), numRural, TRUE, prob=probs))
      allUrban = c(allUrban, rep(FALSE, numRural))
    }
  }
  
  EAIs = allSamples
  kenyaEAs = kenyaPop[EAIs,]
  coordsX = runif(numEAs, min=kenyaEAs$lon - lonRes/2, max=kenyaEAs$lon + lonRes/2)
  coordsY = runif(numEAs, min=kenyaEAs$lat - latRes/2, max=kenyaEAs$lat + latRes/2)
  kenyaEAs$lon = coordsX
  kenyaEAs$lat = coordsY
  decreaseFac = totalKenyaPop/sum(kenyaEAs$pop)
  kenyaEAs$pop = kenyaEAs$pop*decreaseFac
  kenyaEAs$region = countyToRegion(kenyaEAs$admin1)
  kenyaEAs$urban = allUrban
  
  # project coordinates to northing/easting
  projCoords = projKenya(coordsX, coordsY)
  kenyaEAs$east = projCoords[,1]
  kenyaEAs$north = projCoords[,2]
  kenyaEAs$pixelIs = EAIs
  
  kenyaEAs
}

# get empirical mothers and children distribution per cluster
getSurveyEmpiricalDistributions = function(data=NULL, yearOfBirth = 2010) {
  # read in DHS data
  if(is.null(data))
    data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))
  
  # initialize lists of numbers for the empirical distributions
  clusters = unique(data$v001)
  nHouseholds = rep(0, length(clusters))
  nMothers = list()
  nChildren = list()
  
  for(i in 1:length(clusters)) {
    if(i %% 25 == 0)
      print(paste0("household ", i, "/", length(nHouseholds)))
    
    thisCluster = clusters[i]
    clusterData = data[data$v001 == thisCluster, ]
    households = unique(data$v002[data$v001 == thisCluster])
    nHouseholds[i] = max(households)
    
    for(j in 1:length(households)) {
      thisHousehold = households[j]
      thisData = clusterData[clusterData$v002 == thisHousehold, ]
      nMothers = c(nMothers, list(thisData$v138[1]))
      nChildren = c(nChildren, list(sum(thisData$b2 == yearOfBirth)))
    }
  }
  
  nMothers = unlist(nMothers)
  nChildren = unlist(nChildren)
  householdDistribution = ecdf(nHouseholds)
  motherDistribution = ecdf(nMothers)
  childrenDistribution = ecdf(nChildren)
  list(households=householdDistribution, mothers=motherDistribution, children=childrenDistribution)
}

# get empirical mothers and children distribution per cluster
# maxAge: the maximum age of the children being included in the distribution rounded down
getSurveyEmpiricalDistributions2 = function(data=NULL, dataDHS=NULL, maxAge=4) {
  
  # read in DHS data for the estimate of the number of households within each cluster
  print("reading in DHS data")
  if(is.null(dataDHS))
    dataDHS <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))

  # initialize lists of numbers for the empirical distributions
  print("calculating households per cluster empirical distribution")
  urbanDat = dataDHS[dataDHS$v025 == 1, ]
  ruralDat = dataDHS[dataDHS$v025 == 2, ]
  nHouseholds = aggregate(dataDHS$v002, data.frame(v001=dataDHS$v001), function(x) {max(x)})$x
  nHouseholdsUrban = aggregate(urbanDat$v002, data.frame(v001=urbanDat$v001), function(x) {max(x)})$x
  nHouseholdsRural = aggregate(ruralDat$v002, data.frame(v001=ruralDat$v001), function(x) {max(x)})$x
  
  nHouseholds=nHouseholds[(nHouseholds >= 50) & (nHouseholds <= 150)]
  nHouseholdsUrban=nHouseholdsUrban[(nHouseholdsUrban >= 50) & (nHouseholdsUrban <= 150)]
  nHouseholdsRural=nHouseholdsRural[(nHouseholdsRural >= 50) & (nHouseholdsRural <= 150)]

  # read in census data for the other empirical distributions
  print("reading and census data if necessary")
  if(is.null(data))
    data <- data.frame(read_dta("popSurvey2009/Population_2009KPHC_10PCT_STATA.dta"))
  data = data[, c("PROVINCE", "DISTRICT", "DIVISION", "LOCATION", "SUBLOC", "COUNTY", 
                  "RecreatedEANO", "HHNO", "EATYPE", "BARCODE", "P12", "P13")]
  
  print("beginning data table operations")
  dat = as.data.table(data)
  print(system.time(outMothers <- dat[, {
    .("nMothers"=uniqueN(P13[P13 != 0 & P12 <= maxAge]), "urban"=any(EATYPE != 1))
  }
  , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  
  dat0 = dat[P13 != 0 & P12 <= maxAge,]
  print(system.time(outChildren <- dat0[, {
    if(any(unlist(tapply(P13, P13, length, simplify = FALSE))) >= 15)
      print(.SD)
    .("nChildren"=unlist(tapply(P13, P13, length, simplify = FALSE)), "urban"=rep(any(EATYPE != 1), uniqueN(P13)))
  }
  , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  
  nMothers = outMothers[["nMothers"]]
  nMothersUrban = outMothers[urban == TRUE, nMothers]
  nMothersRural = outMothers[urban == FALSE, nMothers]
  nChildren = outChildren[["nChildren"]]
  nChildrenUrban = outChildren[urban == TRUE, nChildren]
  nChildrenRural = outChildren[urban == FALSE, nChildren]
  
  householdDistribution = ecdf(nHouseholds)
  motherDistribution = ecdf(nMothers)
  childrenDistribution = ecdf(nChildren)
  householdDistributionUrban = ecdf(nHouseholdsUrban)
  motherDistributionUrban = ecdf(nMothersUrban)
  childrenDistributionUrban = ecdf(nChildrenUrban)
  householdDistributionRural = ecdf(nHouseholdsRural)
  motherDistributionRural = ecdf(nMothersRural)
  childrenDistributionRural = ecdf(nChildrenRural)
  list(households=householdDistribution, mothers=motherDistribution, children=childrenDistribution,
       householdsUrban=householdDistributionUrban, mothersUrban=motherDistributionUrban, childrenUrban=childrenDistributionUrban,
       householdsRural=householdDistributionRural, mothersRural=motherDistributionRural, childrenRural=childrenDistributionRural)
}

# get empirical women distribution per household
# minAge: the minimum age of the women being included in the distribution rounded down
# maxAge: the maximum age of the women being included in the distribution rounded down
getSurveyEmpiricalDistributionsWomen = function(data=NULL, dataDHS=NULL, minAge=20, maxAge=29, 
                                                datawd="~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/") {
  require(haven)
  thiswd = getwd()
  setwd(datawd)
  
  # read in census data for the other empirical distributions
  print("reading in census data if necessary")
  if(is.null(data))
    data <- data.frame(read_dta("popSurvey2009/Population_2009KPHC_10PCT_STATA.dta"))
  # data = data[, c("PROVINCE", "DISTRICT", "DIVISION", "LOCATION", "SUBLOC", "COUNTY", 
  #                 "RecreatedEANO", "HHNO", "EATYPE", "BARCODE", "P12", "P13")]
  data = data[, c("PROVINCE", "DISTRICT", "DIVISION", "LOCATION", "SUBLOC", "COUNTY", 
                  "RecreatedEANO", "HHNO", "EATYPE", "BARCODE", "P11", "P12")]
  
  # P11: sex (1: male, 2: female)
  # P12: age in years completed
  
  print("beginning data table operations")
  dat = as.data.table(data)
  # print(system.time(outMothers <- dat[, {
  #   .("nMothers"=uniqueN(P13[P13 != 0 & P12 <= maxAge]), "urban"=any(EATYPE != 1))
  # }
  # , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  print(system.time(outWomen <- dat[, {
    .("nWomen"=length(P11[P11 != 1 & P12 <= maxAge & P12 >= minAge]), "urban"=any(EATYPE != 1))
  }
  , by = .(PROVINCE, DISTRICT, DIVISION, LOCATION, SUBLOC, COUNTY, RecreatedEANO, HHNO, BARCODE)]))
  
  nWomen = outWomen[["nWomen"]]
  nWomenUrban = outWomen[urban == TRUE, nWomen]
  nWomenRural = outWomen[urban == FALSE, nWomen]
  
  womenDistribution = ecdf(nWomen)
  womenDistributionUrban = ecdf(nWomenUrban)
  womenDistributionRural = ecdf(nWomenRural)
  
  setwd(thiswd)
  list(women=womenDistribution, womenUrban=womenDistributionUrban, womenRural=womenDistributionRural)
}

# simulate EAs household, mother, and child values in the clusters based on empirical distributions 
# from getSurveyEmpiricalDistributions. First simulate the "long form", where each row is 
# a household, then simulate the "short form", were each row is a EA
simEAEmpirical = function(empiricalDistributions, nEAs=423) {
  # simulate the long form
  householdValues = knots(empiricalDistributions$households)
  householdProbs = empiricalDistributions$households(householdValues) - empiricalDistributions$households(householdValues - 1)
  motherValues = knots(empiricalDistributions$mothers)
  motherProbs = empiricalDistributions$mothers(motherValues) - empiricalDistributions$mothers(motherValues - 1)
  childrenValues = knots(empiricalDistributions$children)
  childrenProbs = empiricalDistributions$children(childrenValues) - empiricalDistributions$children(childrenValues - 1)
  
  # simulate the number of households per numeration area, and mothers and children per household
  nHouseholds = sample(householdValues, nEAs, TRUE, householdProbs)
  nMothers = sample(motherValues, sum(nHouseholds), TRUE, motherProbs)
  nChildren = sample(childrenValues, sum(nHouseholds), TRUE, childrenProbs)
  
  # generate full data sets
  longform = data.frame(eaI=rep(1:nEAs, nHouseholds), mothers=nMothers, children=nChildren)
  totalMothers = aggregate(longform$mothers, data.frame(list(house=longform$eaI)), sum)$x
  totalChildren = aggregate(longform$children, data.frame(list(house=longform$eaI)), sum)$x
  shortform = data.frame(eaI=1:nEAs, mothers=totalMothers, children=totalChildren)
  
  list(longform=longform, shortform=shortform)
}

# same as simEAEmpirical, except use empirical distributions from getSurveyEmpiricalDistributions2, 
# which stratify by urban/rural
simEAEmpirical2 = function(empiricalDistributions, eaDat, eaDatLong) {
  urban = eaDat$urban
  nEAs = length(urban)
  nUrban = sum(urban)
  nRural = nEAs - nUrban
  
  # simulate the number of households per numeration area
  nHouseholdsUrban = recdf(nUrban, empiricalDistributions$householdsUrban)
  nHouseholdsRural = recdf(nRural, empiricalDistributions$householdsRural)
  nHouseholds = urban
  nHouseholds[urban] = nHouseholdsUrban
  nHouseholds[!urban] = nHouseholdsRural
  nMothersUrban = recdf(nUrban, empiricalDistributions$mothersUrban)
  nMothersRural = recdf(nRural, empiricalDistributions$mothersRural)
  nMothers = urban
  nMothers[urban] = nMothersUrban
  nMothers[!urban] = nMothersRural
  nChildrenUrban = recdf(nUrban, empiricalDistributions$childrenUrban)
  nChildrenRural = recdf(nRural, empiricalDistributions$childrenRural)
  nChildren = urban
  nChildren[urban] = nChildrenUrban
  nChildren[!urban] = nChildrenRural
  
  # generate full data sets
  longform = data.frame(eaI=rep(1:nEAs, nHouseholds), mothers=nMothers, children=nChildren)
  totalMothers = aggregate(longform$mothers, data.frame(list(house=longform$eaI)), sum)$x
  totalChildren = aggregate(longform$children, data.frame(list(house=longform$eaI)), sum)$x
  shortform = data.frame(eaI=1:nEAs, mothers=totalMothers, children=totalChildren)
  
  list(longform=longform, shortform=shortform)
}

saveSimDat = function(seed=123, HHoldVar=0) {
  # out = simDat(kenyaEAs, beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5, effRange=300, seed=seed)
  out = simDat(kenyaEAs, beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5, effRange=300, seed=seed, HHoldVar=HHoldVar)
  eaDat = out$eaDat
  clustDat = out$clustDat
  eaID = out$eaI
  save(eaDat, clustDat, eaID, file=paste0("simDatHH", round(HHoldVar, 3), ".RData"))
}

# set thresholds within each region based on number of urban and rural EAs
setThresholds = function() {
  warning("setThresholds is DEPRACATED.  Use setThresholds2 instead.")
  
  # get population data, take log
  pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values=TRUE)
  clusterPop = extract(pop, SpatialPoints(cbind(mort$lon, mort$lat)),method="bilinear")
  clusterLog10Pop = log10(clusterPop)
  
  # only bother with predictions over range of uncertainty
  doPred = (clusterPop != 0) & (clusterPop < 40000) & (clusterPop > 250) & (as.character(mort$admin1) != "Nairobi") & (as.character(mort$admin1) != "Mombasa")
  
  getCountyThresh = function(countyName) {
    # if Nairobi or Mombasa, always urban
    if((countyName == "Nairobi") || (countyName == "Mombasa"))
      return(-Inf)
    
    # get indices for this county
    thisDoPred = doPred & (as.character(mort$admin1) == countyName)
    
    # get data for county model
    modPop = clusterLog10Pop[thisDoPred]
    modUrban = mort$urban[thisDoPred]
    modAd = mort$admin1[thisDoPred]
    
    # fit model and compute threshold
    mod = glm(modUrban ~ modPop, family=binomial)
    thetas = coef(mod)
    popThresh = 10^(-thetas[1]/thetas[2])
    
    popThresh
  }
  
  # compute threshold for each county
  counties = as.character(unique(mort$admin1))
  threshes = sapply(counties, getCountyThresh)
  
  list(counties=counties, threshes=threshes)
}

# set thresholds within each county based on percent population urban
setThresholds2 = function(popRaster=NULL) {
  if(is.null(popRaster)) {
    popRaster
  }
  
  getCountyThresh = function(countyName) {
    # if Nairobi or Mombasa, always urban
    if((countyName == "Nairobi") || (countyName == "Mombasa"))
      return(-Inf)
    
    # do the setup
    thisCounty = as.character(kenyaPop$admin1) == countyName
    thisPop = kenyaPop$popOrig[thisCounty]
    thisTot = sum(thisPop)
    pctUrb = poppc$pctUrb[poppc$County == countyName]/100
    pctRural = 1 - pctUrb
    
    # objective function to minimize
    # objFun = function(thresh) {
    #   curPctUrb = sum(thisPop[thisPop > thresh])/thisTot
    #   (curPctUrb - pctUrb)^2
    # }
    
    # do optimization
    # out = optim(10, objFun)
    # thresh = out$par
    # out = optimize(objFun, c(.01, 50000))
    # thresh = out$par
    
    # calculate threshold by integrating ecdf via sorted value cumulative sum
    sortedPop = sort(thisPop)
    cumsumPop = cumsum(sortedPop)
    threshI = match(1, cumsumPop >= thisTot*pctRural)
    thresh = sortedPop[threshI]
    
    # print(paste0("pctUrb: ", pctUrb, "; resPctUrb: ", sum(thisPop[thisPop > thresh])/thisTot, "; thresh: ", thresh, "; obj: ", out$objective))
    thresh
  }
  
  # compute threshold for each county
  counties = as.character(unique(mort$admin1))
  threshes = sapply(counties, getCountyThresh)
  
  list(counties=counties, threshes=threshes)
}

# given a vector of populations with associated counties, and a threshes 
# object such as an output to a setThresholds function, this computes 
# a vector of length equal to the pops vector whether the point is urban.
setUrbanByThreshes = function(pops, counties, threshes) {
  allThreshes = sapply(1:length(pops), function(i) {threshes$threshes[threshes$counties == counties[i]]})
  urban = pops > allThreshes
  urban
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample
simClusters = function(eaDat, numClusters = 423, seed=123) {
  set.seed(seed)
  
  # set number of clusters empirically if not otherwise set, scale by a factor depending on how many clusters we sample
  thisclustpc = clustpc
  thisclustpc[,2:ncol(thisclustpc)] = round(thisclustpc[,2:ncol(thisclustpc)] * (numClusters / sum(thisclustpc[,4])))
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  counties = as.character(thisclustpc$County)
  urban = eaDat$urban
  allSamples = c()
  for(i in 1:nrow(thisclustpc)) {
    countyName = counties[i]
    countyI = as.character(eaDat$admin1) == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # number of rural and urban clusters to sample
    numUrban = thisclustpc$clustUrb[thisclustpc$County == countyName]
    numRural = thisclustpc$clustRur[thisclustpc$County == countyName]
    
    # urban samples
    thisCountyI = countyI & urban
    # probs = thisCountyI * eaDat$popOrig
    probs = thisCountyI
    probs = probs/sum(probs)
    allSamples = c(allSamples, sample(1:nrow(eaDat), numUrban, TRUE, prob=probs))
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      # probs = thisCountyI * eaDat$popOrig
      probs = thisCountyI
      probs = probs/sum(probs)
      allSamples = c(allSamples, sample(1:nrow(eaDat), numRural, TRUE, prob=probs))
    }
  }
  
  clustDat = eaDat[allSamples,]
  list(sampleI = allSamples, clustDat = clustDat)
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample.  Simulate same number of clusters per county (9)
# urbanProps: either a single proportion or a sequence of proportions of length 47 in order corresponding to counties
simClusters2 = function(eaDat, numClusters=423, urbanProps=NULL, counties=as.character(clustpc$County), seed=NULL, doJitter=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  
  # compute number of clusters to sample from each county
  numPerCounty = round(numClusters/47)
  if(round(numClusters/47) != numClusters/47)
    warning("number of clusters doens't divide evenly into number of counties")
  
  # set the sampling frame
  if(is.null(urbanProps)) {
    # set proportion of urban clusters empirically if not otherwise set, scale by a factor depending on how many clusters we sample
    thisclustpc = clustpc
    thisclustpc[,2:ncol(thisclustpc)] = round(thisclustpc[,2:ncol(thisclustpc)] * (numClusters / sum(thisclustpc$clustTotal)))
  }
  else {
    # set sampling frame based on user proportions that are urban
    
    if(length(urbanProps) == 1)
      urbanProps = rep(urbanProps, 47)
    
    thisclustpc = data.frame(list(County = counties, clustUrb=numPerCounty*urbanProps, clustRur = numPerCounty*(1-urbanProps), 
                                  clustTotal = numPerCounty))
  }
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  urban = eaDat$urban
  allSamples = c()
  sampleWeights = c()
  for(i in 1:nrow(thisclustpc)) {
    countyName = counties[i]
    countyI = as.character(eaDat$admin1) == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # number of rural and urban clusters to sample (sample numPerCounty total per county, so must rescale each)
    numUrban = thisclustpc$clustUrb[thisclustpc$County == countyName] * (numPerCounty/thisclustpc$clustTotal[thisclustpc$County == countyName])
    numRural = thisclustpc$clustRur[thisclustpc$County == countyName] * (numPerCounty/thisclustpc$clustTotal[thisclustpc$County == countyName])
    
    # make sure there's a whole number of rural and urban summing to numPerCounty (9 by default)
    if(round(numUrban) + round(numRural) != numPerCounty) {
      if(runif(1) > .5) {
        numUrban = ceiling(numUrban)
        numRural = floor(numRural)
      }
      else {
        numUrban = floor(numUrban)
        numRural = ceiling(numRural)
      }
    }
    else {
      numRural = round(numRural)
      numUrban = round(numUrban)
    }
    
    # urban samples
    thisCountyI = countyI & urban
    # probs = thisCountyI * eaDat$popOrig
    probs = thisCountyI
    probs = probs/sum(probs)
    newSamples = sample(1:nrow(eaDat), numUrban, FALSE, prob=probs)
    newWeights = rep(1/(numUrban/sum(eaDat$urban[thisCountyI])), numUrban)
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      # probs = thisCountyI * eaDat$popOrig
      probs = thisCountyI
      probs = probs/sum(probs)
      newSamples = c(newSamples, sample(1:nrow(eaDat), numRural, FALSE, prob=probs))
      newWeights = c(newWeights, rep(1/(numRural/sum(!eaDat$urban[thisCountyI])), numRural))
    }
    
    # update clusters and sampling weights for this county
    allSamples = c(allSamples, newSamples)
    sampleWeights = c(sampleWeights, newWeights)
  }
  
  clustDat = eaDat[allSamples,]
  clustDat$samplingWeight = sampleWeights
  
  # computer cluster locations if necessary
  if(doJitter) {
    stop("how do we want to jitter the data?")
  }
  
  list(sampleI = allSamples, clustDat = clustDat)
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample.  Simulate same number of clusters per county (9)
# urbanOverSample: oversample in urban areas by this ratio.  Any urban EA is urbanProp times more likely 
#            to be sampled than a rural EA. If urbanOverSample is 1 than no oversampling occurs
simClusters3 = function(eaDat, numClusters=423, urbanOverSample=1, nsim=1, seed=NULL) {
  if(!is.null(seed))
    set.seed(seed)
  
  # compute number of clusters to sample from each county
  numPerCounty = round(numClusters/47)
  if(round(numClusters/47) != numClusters/47)
    warning("number of clusters doens't divide evenly into number of counties")
  
  # set some basic variables
  thisclustpc = clustpc
  counties = clustpc$County
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  urban = eaDat$urban
  eaIs = matrix(nrow=numClusters, ncol=nsim)
  sampleWeights = matrix(nrow=numClusters, ncol=nsim)
  curRow = 1
  for(i in 1:nrow(thisclustpc)) {
    countyName = counties[i]
    countyI = as.character(eaDat$admin1) == countyName
    print(paste0("County ", i, "/", length(counties), ": ", countyName))
    
    # sampling probabilities depend on urban and rural strata
    thisUrban = countyI & urban
    thisRural = countyI & !urban
    probs = thisRural + thisUrban * urbanOverSample
    probs = probs/sum(probs)
    endRow = curRow - 1 + numPerCounty
    getSamples = function(i) {
      sample(1:nrow(eaDat), numPerCounty, FALSE, prob=probs)
    }
    thisEAIs = sapply(1:nsim, getSamples)
    
    eaIs[curRow:endRow, ] = thisEAIs
    sampleWeights[curRow:endRow, ] <- matrix(1/probs[thisEAIs], nrow=numPerCounty, ncol=nsim)
    curRow = endRow + 1
  }
  
  list(eaIs=eaIs, sampleWeights=sampleWeights)
}

# convert from output of simClusters3 to the result of generateSimDataSets.R
genAndreaFormatFromEAIs = function(eaDat, eaIs, sampleWeights) {
  lapply(1:ncol(eaIs), genClustDatFromEAIs, eaDat=eaDat, eaIs=eaIs, sampleWeights=sampleWeights)
}

# given the above function (simClusters2Mult), generates a dataframe following the format 
# of simClusters2 given a simulation index, i
genClustDatFromEAIs = function(eaDat, eaIs, sampleWeights, i) {
  clustDat = eaDat[eaIs[,i],]
  clustDat$samplingWeight = sampleWeights[,i]
  clustDat$eaIs = eaIs[,i]
  clustDat
}

# convert from output of simClusters3 to the result of generateSimDataSets.R
genAndreaFormatFromEAIsLong = function(eaDat, eaIs, eaDatLong, HHIs, sampleWeights) {
  lapply(1:ncol(eaIs), genClustDatFromEAIsLong, eaDat=eaDat, eaIs=eaIs, eaDatLong=eaDatLong, 
         HHIs=HHIs, sampleWeights=sampleWeights)
}

# given the above function (simClusters2Mult), generates a dataframe following the format 
# of simClusters2 given a simulation index, i
genClustDatFromEAIsLong = function(eaDat, eaIs, eaDatLong, HHIs, sampleWeights, i) {
  clustDatLong = as.data.table(eaDatLong[HHIs[,i],])
  
  # aggregate the long format to the short format
  # TODO: figure out what columns are in the cluster data and how to aggregate
  clustDat = clustDatLong[, .(lon=lon[1], lat=lat[1], pop=pop[1], popOrig=popOrig[1], 
                              admin1=admin1[1], region=region[1], 
                              east=east[1], north=north[1], urban=urban[1], nHH=nHH[1], 
                              numWomen=sum(numWomen), numChildren=sum(numChildren), 
                              died=sum(died), 
                              probsHH=probsHH[1], trueProbDeathNoNug=trueProbDeathNoNug[1], 
                              trueProbDeath=trueProbDeath[1]), 
                          by=eaIs]
  
  clustDat$samplingWeight = sampleWeights[,i]
  clustDat$eaIs = eaIs[,i]
  clustDat
}

# simulate clusters from enumeration areas.  Returns indices of rows in eaDat that are included 
# in the cluster sample.  Simulate same number of clusters per county (9)
# eaDatLong: long format enumeration area table, where each row is a household
# nPerStrata: if fixedPerStrata == TRUE, this is the number of clusters sampled per strata
# SRS: either SRS or PPS sampling
simClustersEmpirical = function(eaDat, eaDatLong, nsim=1, seed=NULL, urbanOverSamplefrac=0, 
                                nHHSampled=25, fixedPerStrata=FALSE, nPerStrata=3, SRS=FALSE, verbose=TRUE) {
  if(!is.null(seed))
    set.seed(seed)
  
  # get the number of clusters for each strata
  nUrbanClusters = clustpc$clustUrb
  nRuralClusters = clustpc$clustRur
  
  # make sure the ea tables are a data.tables
  if(!("data.table" %in% class(eaDatLong)))
    eaDatLong = as.data.table(eaDatLong)
  if(!("data.table" %in% class(eaDat)))
    eaDat = as.data.table(eaDat)
  eaDat = cbind(eaIs=1:nrow(eaDat), eaDat)
  
  if(!fixedPerStrata) {
    # set number of clusters empirically if not otherwise set, scale by a factor depending on how much to oversample urban areas
    thisclustpc = clustpc
    thisclustpc$clustUrb = round(thisclustpc$clustUrb * (1 + urbanOverSamplefrac))
    thisclustpc$clustUrb[thisclustpc$clustUrb >= thisclustpc$clustTotal] = thisclustpc$clustTotal[thisclustpc$clustUrb >= thisclustpc$clustTotal]
    thisclustpc$clustRur = thisclustpc$clustTotal - thisclustpc$clustUrb
    
    if(SRS) {
      # in this case, we want the naive model to be unbiased so we need the proportion of urban clusters and rural clusters to 
      # be representative of the population within each county
      thisclustpc$clustUrb = (poppc$pctUrb / 100) * thisclustpc$clustTotal
      thisclustpc$clustRur = (1 - poppc$pctUrb / 100) * thisclustpc$clustTotal
      
      # NOTE: the number of urban and rural clusters will be fractions, so we need to pick the number randomly for each county 
      #       so that on average the estimate within each county will be unbiased
    }
  } else {
    thisclustpc = clustpc
    thisclustpc$clustUrb = nPerStrata
    thisclustpc$clustRur = nPerStrata
    thisclustpc$clustRur[thisclustpc$County == "Mombasa"] = 0
    thisclustpc$clustRur[thisclustpc$County == "Nairobi"] = 0
  }
  
  # collect samples stratified by county and rural/urban based on empirical distribution
  counties = as.character(thisclustpc$County)
  urban = eaDat[["urban"]]
  
  # check if we need to randomly allocate clusters
  roundedUrban = round(thisclustpc$clustUrb, 0)
  roundedRural = round(thisclustpc$clustRur, 0)
  randomClusters = any(roundedUrban != thisclustpc$clustUrb) | any(roundedRural != thisclustpc$clustRur)
  
  # if we have random clusters, precalculate certain probabilities
  if(randomClusters) {
    extraUrban = thisclustpc$clustUrb - floor(thisclustpc$clustUrb)
    extraRural = thisclustpc$clustRur - floor(thisclustpc$clustRur)
    probUrbanTotal = sum(extraUrban) / sum(extraUrban + extraRural)
    probUrban = extraUrban / sum(extraUrban)
    probRural = extraRural / sum(extraRural)
  }
  
  generateSample = function() {
    thisclustpcI = thisclustpc
    
    if(randomClusters) {
      # we need to randomly allot clusters to strata so the naive case remains unbiased
      
      # calculate the number of clusters we need to randomly allot
      nExtra = round(sum(extraUrban + extraRural), 0)
      
      # calculate the number of extra urban and rural clusters
      nExtraUrbanTotal = rbinom(1, nExtra, probUrbanTotal)
      nExtraRuralTotal = nExtra - nExtraUrbanTotal
      
      # calculate the number of extra clusters in each strata conditional on the number of urban and rural
      nExtraUrban = rmultinom(1, nExtraUrbanTotal, probUrban)
      nExtraRural = rmultinom(1, nExtraRuralTotal, probRural)
      
      # add extra clusters to their corresponding strata
      thisclustpcI$clustUrb = floor(thisclustpcI$clustUrb) + nExtraUrban
      thisclustpcI$clustRur = floor(thisclustpcI$clustRur) + nExtraRural
      thisclustpcI$clustTotal = thisclustpcI$clustUrb + thisclustpcI$clustRur
    }
    
    allSamples = c()
    # allSamples = list()
    for(i in 1:nrow(thisclustpcI)) {
      countyName = counties[i]
      countyI = eaDat[,as.character(admin1) == countyName]
      if(verbose)
        print(paste0("County ", i, "/", length(counties), ": ", countyName))
      
      # number of rural and urban clusters to sample
      numUrban = thisclustpcI$clustUrb[thisclustpcI$County == countyName]
      numRural = thisclustpcI$clustRur[thisclustpcI$County == countyName]
      
      # urban samples
      thisCountyI = countyI & urban
      
      # probs = thisCountyI * eaDat$popOrig
      if(SRS)
        probs = thisCountyI
      else
        probs = thisCountyI * eaDat$nHH
      probs = probs/sum(probs)
      # allSamples = c(allSamples, sample(1:nrow(eaDat), numUrban, FALSE, prob=probs))
      zs = probs[thisCountyI] * numUrban # the inclusion probabilities
      if(any(zs >= 1-1e-6)) {
        warning("inclusion probabilities greater than one, modifying to be .999")
        zs[zs >= 1-1e-6] = .999
      }
      # newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPtille(zs))]
      newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPmidzuno(zs))]
      allSamples = c(allSamples, newSamples)
      
      # allSamples = c(allSamples, pps.sampling(eaDat$nHH[thisCountyI], numUrban, (1:nrow(eaDat))[thisCountyI], FALSE, prob=probs))
      
      # rural samples (if the county has any)
      if((countyName != "Mombasa") && (countyName != "Nairobi")) {
        thisCountyI = countyI & !urban
        # probs = thisCountyI * eaDat$popOrig
        if(SRS)
          probs = thisCountyI
        else
          probs = thisCountyI * eaDat$nHH
        probs = probs/sum(probs)
        
        zs = probs[thisCountyI] * numRural # the inclusion probabilities
        if(any(zs >= 1-1e-6)) {
          warning("inclusion probabilities greater than one, modifying to be .999")
          zs[zs >= 1-1e-6] = .999
        }
        # newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPtille(zs))]
        newSamples = (1:length(thisCountyI))[thisCountyI][as.logical(UPmidzuno(zs))]
        
        # allSamples = c(allSamples, sample(1:nrow(eaDat), numRural, FALSE, prob=probs))
        allSamples = c(allSamples, newSamples)
      }
    }
    sort(allSamples)
  }
  # the comment out line below works, accept it does not enable you to print out progress
  # eaIs = replicate(nsim, generateSample())
  eaIs = matrix(ncol=nsim, nrow=sum(thisclustpc$clustUrb) + sum(thisclustpc$clustRur))
  for(i in 1:nsim) {
    print(paste0("simulating data set ", i, "/", nsim))
    eaIs[, i] = generateSample()
  }
  
  ## decide which households within each cluster to sample
  eaDatLong = data.table(HHI=1:nrow(eaDatLong), eaDatLong)
  sampleHouseholds = function(sim) {
    assign("thisEAIs", eaIs[,sim], globalenv())
    clusterDatLong = eaDatLong[eaIs %in% get("thisEAIs", globalenv()), ]
    assign("nHHSampled", nHHSampled, globalenv())
    HHIs = clusterDatLong[, {
      .("HHIs"=sort(sample(HHI, get("nHHSampled", globalenv()), FALSE)))
    }, by=eaIs]
    # cumulativeHouseholds = c(0, cumsum(eaDat[["nHH"]]))
    # HHIs[["HHIs"]] = HHIs[["HHIs"]] + cumulativeHouseholds[HHIs[["eaIs"]] - 1]
    HHIs[["HHIs"]]
  }
  HHIs = sapply(1:nsim, sampleHouseholds)
  
  # compute the probability of each enumeration area and household being in our sample
  probsEA = rep(0, nrow(eaDat))
  for(i in 1:nrow(thisclustpc)) {
    countyName = counties[i]
    countyI = eaDat[,as.character(admin1) == countyName]
    
    # number of rural and urban clusters to sample
    numUrban = thisclustpc$clustUrb[thisclustpc$County == countyName]
    numRural = thisclustpc$clustRur[thisclustpc$County == countyName]
    
    # urban samples
    thisCountyI = countyI & urban
    if(SRS)
      probs = thisCountyI
    else
      probs = thisCountyI * eaDat[["nHH"]]
    probs = probs * (1/sum(probs))
    zs = probs[thisCountyI] * numUrban # the inclusion probabilities
    if(any(zs >= 1-1e-6)) {
      warning("inclusion probabilities greater than one, modifying to be .999")
      zs[zs >= 1-1e-6] = .999
    }
    probsEA[thisCountyI] = zs
    
    if(SRS) {
      if(any(abs(unique(zs) - thisclustpc$clustUrb[i]/sum(eaDat$admin1 == countyName & eaDat$urban == TRUE)) > .0001)) {
        print("woah there")
      }
    }
    
    # rural samples (if the county has any)
    if((countyName != "Mombasa") && (countyName != "Nairobi")) {
      thisCountyI = countyI & !urban
      # probs = thisCountyI * eaDat$popOrig
      if(SRS)
        probs = thisCountyI
      else
        probs = thisCountyI * eaDat[["nHH"]]
      probs = probs * (1/sum(probs))
      zs = probs[thisCountyI] * numRural # the inclusion probabilities
      if(any(zs >= 1-1e-6)) {
        warning("inclusion probabilities greater than one, modifying to be .999")
        zs[zs >= 1-1e-6] = .999
      }
      probsEA[thisCountyI] = zs
      
      if(SRS) {
        if(any(abs(unique(zs) - thisclustpc$clustRur[i]/sum(eaDat$admin1 == countyName & eaDat$urban == FALSE)) > .0001)) {
          print("woah there")
        }
      }
    }
  }
  probsHH = probsEA * (25 / eaDat[["nHH"]])
  
  # calculate the sampling weights. Each child has probability of being sampled same as its household, 
  # calculate the child's weight then multiply by number in HH or EA
  sampleWeightsLong = eaDatLong[["numChildren"]] / probsHH[eaDatLong[["eaIs"]]]
  # sampleWeights = eaDat[["numChildren"]] / probsHH
  # sampleWeightsEA = sampleWeights
  # sampleWeights = matrix(sampleWeights[eaIs], ncol=nsim)
  getSampleWeightsCluster = function(i) {
    sampleWeightsClusterLong = sampleWeightsLong[HHIs[,i]]
    thisEAIs = eaDatLong[HHIs[,i], eaIs]
    aggregate(sampleWeightsClusterLong, by=list(thisEAIs), FUN="sum")$x
  }
  sampleWeights = sapply(1:nsim, getSampleWeightsCluster)
  
  list(eaIs=eaIs, HHIs=HHIs, sampleWeights=sampleWeights)
}

expit = function(x) { exp(x)/(1+exp(x)) }

# function for simulating data given enumeration areas and their info.  Simulate using 
# LatticeKrig model on logit scale
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
simDat = function(eaDat, clustDat=NULL, urbanProps=NULL, counties=as.character(clustpc$County), 
                  nsim=1, margVar=4, nu=1.5, NC=5, effRange=300, nLayer=3, beta0=0, gamma=-1, 
                  tausq=1, normalize=TRUE, nBuffer=5, womenFrac=1/3, 
                  numClusters=423, seed=NULL, fixedWomenPerClust=TRUE, HHoldVar=0) {
  if(!is.null(seed))
    set.seed(seed)
  
  if(nsim != 1)
    stop("multiple simulations not yet implemented")
  
  ### first generate Binomial probabilities from transformed logit scale GP
  # generate Lattice Krig simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  LKArgs = list(coords=eaCoords, nsim=nsim, NC=NC, margVar=margVar, effRange=effRange, nu=nu, 
                     nLayer=nLayer, normalize=normalize, nBuffer=nBuffer)
  simVals = do.call("LKSimulator2", LKArgs)
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=nsim)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  if(is.null(clustDat))
    clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  
  # collect indices corresponding to EA and dataframe
  clustI = clustDat$sampleI
  clustDat = clustDat$clustDat
  
  # if fewer than 25 women in the EA, sample exactly all of them (only matters if fixedWomenPerClust == FALSE)
  numWomenClust = apply(cbind(numWomen, 25), 1, min)
  numChildrenClust = numWomenClust
  
  # helper function: for any EA index, get the number of mortalities sampled
  getClustMort = function(ind) {
    # convert from binomial data to multiple Bernouli trials
    thisDied = c(rep(1, eaDied[ind]), rep(0, numChildren[ind] - eaDied[ind]))
    
    # sum of sampled mortalities
    sum(sample(thisDied, numChildrenClust[ind], FALSE))
  }
  
  # get number of mortalities within each cluster out of 25 sampled households in cluster
  clustDied = sapply(clustI, getClustMort)
  
  # return simulated data
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  clustDat$died = clustDied
  clustDat$numWomen = numWomenClust[clustI]
  clustDat$numChildren = numChildrenClust[clustI]
  clustDat$trueProbDeath = probs[clustI]
  clustDat$trueProbDeathNoNug = probsNoNug[clustI]
  
  list(eaDat=eaDat, clustDat=clustDat, eaI=clustI)
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# SPDE or LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
# useSPDE: simulate the datasets with SPDE model or LatticeKrig model
simDat2 = function(eaDat, clustDat=NULL, nsim=1, margVar=4, nu=1, NC=5, effRange=300, 
                   nLayer=3, beta0=0, gamma=-1, tausq=1, normalize=TRUE, nBuffer=5, 
                   womenFrac=1/3, urbanOverSample=1, numClusters=423, seed=NULL, 
                   HHoldVar=0, fixedWomenPerClust=TRUE, useSPDE=TRUE) {
  if(!is.null(seed))
    set.seed(seed)
  
  ### first generate Binomial probabilities from transformed logit scale GP
  # generate Lattice Krig simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  print("Simulating nationwide mortality rates and data")
  if(!useSPDE) {
  LKArgs = list(coords=eaCoords, nsim=1, NC=NC, margVar=margVar, effRange=effRange, nu=nu, 
                nLayer=nLayer, normalize=normalize, nBuffer=nBuffer)
  simVals = do.call("LKSimulator2", LKArgs)
  }
  else {
    if(nu != 1)
      stop("SPDE model only supports nu=1")
    
    SPDEArgs = list(coords=eaCoords, nsim=1, margVar=margVar, effRange=effRange)
    simVals = do.call("simSPDE", SPDEArgs)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
  }
  
  # return simulated data
  print("finishing up...")
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIs(eaDat, clustDat$eaIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
simDatLK = function(eaDat, clustDat=NULL, nsim=1, margVar=4, nu=1.5, NC=5, effRange=300, 
                      nLayer=3, beta0=0, gamma=-1, tausq=1, normalize=TRUE, nBuffer=5, 
                      womenFrac=1/3, urbanOverSample=1, numClusters=423, seed=NULL, fullEADat=NULL, 
                      HHoldVar=0, fixedWomenPerClust=TRUE) {
  if(!is.null(seed))
    set.seed(seed)
  
  ### first generate Binomial probabilities from transformed logit scale GP
  # generate Lattice Krig simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  print("Simulating nationwide mortality rates and data")
  LKArgs = list(coords=eaCoords, nsim=1, NC=NC, margVar=margVar, effRange=effRange, nu=nu, 
                nLayer=nLayer, normalize=normalize, nBuffer=nBuffer)
  simVals = do.call("LKSimulator2", LKArgs)
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
  }
  
  # return simulated data
  print("finishing up...")
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIs(eaDat, clustDat$eaIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}

# function for simulating data given enumeration areas and their info.  Simulate using 
# SPDE or LatticeKrig model on logit scale.  Note that this function has some notable differences with simDat:
# 1. allows for multiple simulations, although each of them share the same simulate EA data
# 2. uses simClusters3 instead of simClusters2, which enables oversampling in urban or rural areas
# empiricalDistributions: the set of empirical distributions for households per cluster, 
#                         mothers per household, and children per mother stratified by urban
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDatShort: enumeration area data set in "short form" were each row is a cluster. 
#             "long form" is where each row is a household.
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the spatial process, excluding household end cluster effects. 
#          If 0, no spatial component is included
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
simDatEmpirical = function(empiricalDistributions, eaDat, clustDat=NULL, nsim=1, 
                           margVar=4, effRange=150, beta0=0, gamma=-1, tausq=1, 
                           urbanOverSamplefrac=0, seed=NULL, 
                           fullEADat=NULL, HHoldVar=0, nHHSampled=25) {
  if(!is.null(seed))
    set.seed(seed)
  
  ### simulate binomial data for each enumeration area
  # first generate the number of households
  numHouseholds = rep(0, nrow(eaDat))
  numHouseholds[eaDat$urban] = recdf(sum(eaDat$urban), empiricalDistributions$householdsUrban)
  numHouseholds[!eaDat$urban] = recdf(sum(!eaDat$urban), empiricalDistributions$householdsRural)
  eaDat$nHH = numHouseholds
  
  # now expand the eaDat table to be in the long format, were each row is a house
  rowsLong = rep(1:nrow(eaDat), numHouseholds)
  eaDatLong = eaDat[rowsLong, ]
  eaDatLong$eaIs = rowsLong
  eaUrbanLong = eaDatLong$urban
  
  # generate how many women and childen are in each cluster
  numWomenLong = rep(0, length(eaUrbanLong))
  numWomenLong[eaUrbanLong] = recdf(sum(eaUrbanLong), empiricalDistributions$mothersUrban)
  numWomenLong[!eaUrbanLong] = recdf(sum(!eaUrbanLong), empiricalDistributions$mothersRural)
  eaDatLong$numWomen = numWomenLong
  
  # now expand the eaDat table to be in the extra long format, were each row is a mother
  rowsExtraLong = rep(1:length(eaUrbanLong), numWomenLong)
  eaUrbanExtraLong = eaUrbanLong[rowsExtraLong]
  
  # generate the number of children per mother
  numChildrenExtraLong = rep(0, length(eaUrbanExtraLong))
  numChildrenExtraLong[eaUrbanExtraLong] = recdf(sum(eaUrbanExtraLong), empiricalDistributions$childrenUrban)
  numChildrenExtraLong[!eaUrbanExtraLong] = recdf(sum(!eaUrbanExtraLong), empiricalDistributions$childrenRural)
  
  # aggregate up from mother to household level
  numChildrenLong = rep(0, length(eaUrbanLong))
  numChildrenLong[numWomenLong!=0] = aggregate(numChildrenExtraLong, list(rowsExtraLong), sum)$x
  eaDatLong$numChildren = numChildrenLong
  
  ### generate Binomial probabilities from transformed logit scale GP
  # generate Lattice Krig simulations
  eaCoords = cbind(eaDat$east, eaDat$north)
  print("Simulating nationwide mortality rates and data")
  if(margVar != 0) {
    SPDEArgs = list(coords=eaCoords, nsim=1, margVar=margVar, effRange=effRange)
    simVals = do.call("simSPDE", SPDEArgs)
  } else {
    simVals = matrix(rep(0, nrow(eaCoords)), ncol=1)
  }
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  simValsNugLong = simValsNug[rowsLong]
  
  if(HHoldVar != 0) {
    simValsFinalLong = simValsNugLong + matrix(rnorm(length(simValsNugLong), sd=sqrt(HHoldVar)), ncol=1)
    probsFinalLong = expit(simValsFinalLong)
    eaDiedLong = matrix(rbinom(nrow(eaDatLong)*nsim, eaDatLong$numChildren, probsFinalLong), ncol=nsim)
  }
  else {
    probsFinalLong = expit(simValsNugLong)
    eaDiedLong = matrix(rbinom(nrow(eaDatLong)*nsim, eaDatLong$numChildren, rep(probsFinalLong, nsim)), ncol=nsim)
  }
  # transform back to original scale for the cluster level probabilities
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  probsNoNugLong = probsNoNug[rowsLong]
  probsNugLong = probs[rowsLong]
  
  eaDatLong$died = eaDiedLong
  eaDatLong$probsHH = probsFinalLong
  eaDatLong$trueProbDeathNoNug = probsNoNugLong
  eaDatLong$trueProbDeath = probsNugLong # technically, probsHH is the true probability at the household level
  
  # aggregate up from household to cluster level
  numChildren = aggregate(numChildrenLong, list(rowsLong), sum)$x
  numWomen = aggregate(numWomenLong, list(rowsLong), sum)$x
  eaDied = aggregate(c(eaDiedLong), list(rowsLong), sum)$x
  
  eaDat = as.data.frame(eaDat)
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numChildren
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    # clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
    clustDat = simClustersEmpirical(eaDat, eaDatLong, nsim, NULL, urbanOverSamplefrac, nHHSampled)
  }
  
  # return simulated data
  print("finishing up...")
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIsLong(eaDat, clustDat$eaIs, eaDatLong, clustDat$HHIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}

##### generate predictions for UM5R in counties (regions?), EAs
# datForDirect: 
fitMercerMod = function(datForDirect) {
  ### the following code that is commented out will theoretically run the modified 
  ### version of SUMMER in order to get the direct estimates
  # # Use SUMMER package.  See vignette("summer-vignette", "SUMMER")
  # require("SUMMER")
  # 
  # # We need: 
  # # clustid
  # # id
  # # region
  # # time (a factor)
  # # age (this will be 0 for each observations)
  # # weights
  # # strata (region * urban)
  # 
  # # NOTE: only use ageMonth as timeVar because it's the same for each row
  # datIn = datForDirect
  # datIn$ageMonth="0"
  # datIn$strata = datIn$regionRural
  # datIn$died = datIn$died.x
  # dat = myCountrySummary(births=datIn, years = "0", 
  #                        regionVar="admin1", timeVar="ageMonth", weightsVar="weight")
  
  
  # instead of using SUMMER, use Andrea's precomputed direct estimates and Mercer et al. code:
  load("resultsDirect.RData")
  
}

# Same as simDatLK, but uses SPDE model to simulate data set instead of LatticeKrig.
# urbanOverSample: ratio of probabilities passed to the sample function between urban and rural sampling probabilities
#                  (see simClusters3 for more details)
# eaDat: enumeration area data set
# clustDat: cluster samples data set
# urbanProps: see simClusters2
# counties: a vector of county names giving the order with which to simulate the data
# nsim: number of simulations
# margVar: marginal variance of the LatticeKrig process, excluding household end cluster effects
# nu: matern smoothness perimeter
# NC: number of latticed basis element for coarsest data grid along longest data dimension
# effRange: effective range of the latticeKrig process
# nLayer: number of layers in the LatticeKrig process
# beta0: intercept of logit model for mortality rate
# gamma: effect of urban on logit scale for logit model for mortality rate
# tausq: cluster effect variance for logit model of mortality rate
# normalize: whether or knot to normalize the LatticeKrig process
# nBuffer: number of buffer basis elements 4 LatticeKrig process
# womenFrac: proportion of total population in each cluster that is a woman
# numClusters: number of clusters in the faux cluster data set
# seed: random number generator seed
# fixedWomenPerClust: set the same number of women per cluster
simDatSPDE = function(eaDat, clustDat=NULL, nsim=1, margVar=1, effRange=300, 
                      beta0=0, gamma=-1, tausq=1, 
                      womenFrac=1/3, urbanOverSample=1, numClusters=423, seed=NULL, fullEADat=NULL, 
                      HHoldVar=0, fixedWomenPerClust=TRUE, plotProbs=FALSE, 
                      savePlots=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  
  ### first generate Binomial probabilities from transformed logit scale GP
  ## generate SPDE simulations
  # generate mesh grid
  eaCoords = cbind(eaDat$east, eaDat$north)
  mesh = getSPDEMeshGrid(eaCoords, doPlot = FALSE)
  
  print("Simulating nationwide mortality rates and data")
  SPDEArgs = list(coords=eaCoords, nsim=1, margVar=margVar, effRange=effRange)
  simVals = do.call("simSPDE", SPDEArgs)
  
  # add in intercept
  simVals = simVals + beta0
  
  # add in urban effect
  simVals = sweep(simVals, 1, gamma*eaDat$urban, "+")
  
  # add in nugget/cluster effect
  simValsNug = simVals + matrix(rnorm(length(simVals), sd=sqrt(tausq)), ncol=1)
  
  # transform back to original scale
  probs = expit(simValsNug)
  probsNoNug = expit(simVals)
  
  # plot first simulation if requested by user on both logit and linear scales
  if(plotProbs) {
    if(savePlots)
      pdf("figures/spdeSimTest.pdf", width=8, height=5)
    
    par(mfrow=c(1, 2))
    quilt.plot(eaCoords, simValsNug[, 1], main ="Simulated Logit Mortality Rates", 
               xlab="Easting", ylab="Northing")
    plotMapDat(project=TRUE)
    
    quilt.plot(eaCoords, probs[, 1], main ="Simulated Mortality Rates", 
               xlab="Easting", ylab="Northing")
    plotMapDat(project=TRUE)
    
    if(savePlots)
      dev.off()
  }
  
  ### simulate binomial data for each enumeration area
  # generate how many women and childen are in each cluster (assume ~1/3 of population by default, and 
  # 1 child per mother)
  if(! fixedWomenPerClust)
    numWomen = round(womenFrac * eaDat$pop)
  else
    numWomen = rep(25, nrow(eaDat))
  numChildren = numWomen
  
  # simulate mortalities
  if(HHoldVar != 0)
    eaDied = matrix(rLogisticNormBin(nrow(eaDat)*nsim, numChildren, rep(simValsNug, nsim), HHoldVar), ncol=nsim)
  else
    eaDied = matrix(rbinom(nrow(eaDat)*nsim, numChildren, rep(probs, nsim)), ncol=nsim)
  
  ### sample clusters and households within EAs
  # first generate clusters
  # if(is.null(clustDat))
  #   clustDat = simClusters2(eaDat, numClusters, urbanProps, counties, seed=NULL)
  if(is.null(clustDat)) {
    print("simulation cluster locations:")
    clustDat = simClusters3(eaDat, numClusters, urbanOverSample, nsim)
  }
  
  # return simulated data
  print("finishing up...")
  eaDat$died = eaDied
  eaDat$numWomen = numWomen
  eaDat$numChildren = numWomen
  eaDat$trueProbDeath = probs
  eaDat$trueProbDeathNoNug = probsNoNug
  
  # return cluster data in Andrea's format:
  clustList = genAndreaFormatFromEAIs(eaDat, clustDat$eaIs, clustDat$sampleWeights)
  
  list(eaDat=eaDat, clustDat=clustList)
}