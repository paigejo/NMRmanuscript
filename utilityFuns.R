# Kenya utility scripts

# for plotting administration data
# project: if FALSE, plot with lon/lat coordinates.  Otherwise, plot with projected coords 
#          using projKenya function.  This can be used when plotting the projected `east' 
#          and `north' variables in kenyaEAs for instance.
# ...: arguments to polygon function
plotMapDat = function(mapDat = adm0, plotVar=NULL, varCounties=sort(as.character(unique(mort$admin1))), zlim=NULL, project=FALSE, cols=tim.colors(), 
                      legend.mar=7, new=FALSE, plotArgs=NULL, main=NULL, xlim=NULL, xlab=NULL, scaleFun = function(x) {x}, scaleFunInverse = function(x) {x}, 
                      ylim=NULL, ylab=NULL, n.ticks=5, min.n=5, ticks=NULL, tickLabels=NULL, asp=1, legend.width=1.2, ...) {
  # do setup for ploting data by county if necessary
  if(!is.null(plotVar)) {
    if(is.null(zlim)) {
      zlim = range(plotVar)
    }
    
    # get region names from map data
    regionNames = mapDat@data$NAME_1
    
    # make sure county names are consistent for mapDat == adm1
    regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
    regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
  }
  
  # generate new plot if necessary
  if(new) {
    # set graphical parameters so the legend won't overlap with plot
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
    
    if(is.null(plotArgs)) {
      if(project) {
        if(is.null(xlab))
          xlab = "East (km)"
        if(is.null(xlim))
          xlim = eastLim
        if(is.null(ylab))
          ylab = "North (km)"
        if(is.null(ylim))
          ylim = northLim
      }
      else {
        if(is.null(xlab))
          xlab = "Longitude"
        if(is.null(xlim))
          xlim = kenyaLonRange
        if(is.null(ylab))
          ylab = "Latitude"
        if(is.null(ylim))
          ylim = kenyaLatRange
      }
      if(is.null(main))
        main = ""
      plotArgs = list(main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, asp=asp)
    }
    # par( oma=c( 0,0,0,6)) # leave room for the legend
    do.call("plot", c(list(1, 2, type="n"), plotArgs))
  }
  
  # add polygons to plot
  polys = mapDat@polygons
  plotCounty = function(i) {
    countyPolys = polys[[i]]@Polygons
    
    if(is.null(plotVar)) {
      if(!project)
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords), list(...)))})
      else
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(projKenya(countyPolys[[x]]@coords)), list(...)))})
    }
    else {
      # get index of plotVar corresponding to this county
      thisI = which(varCounties == regionNames[i])
      
      # get color to plot
      vals = c(zlim, scaleFun(plotVar[thisI]))
      vals = vals-vals[1]
      vals = vals/(vals[2] - vals[1])
      col = cols[round(vals[3]*(length(cols)-1))+1]
      
      if(!project)
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(countyPolys[[x]]@coords, col=col), list(...)))})
      else
        sapply(1:length(countyPolys), function(x) {do.call("polygon", c(list(projKenya(countyPolys[[x]]@coords), col=col), list(...)))})
    }
    
  }
  sapply(1:length(polys), plotCounty)
  
  if(!is.null(plotVar)) {
    # add legend
    # par( oma=c(0,0,0,2))
    if(is.null(ticks))
      ticks = scaleFun(pretty(scaleFunInverse(zlim), n=n.ticks, min.n=min.n))
    else
      ticks = scaleFun(ticks)
    if(is.null(tickLabels))
      tickLabels = scaleFunInverse(ticks)
    # par( oma=c( 0,0,0,3))
    image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
               col=cols, add = TRUE, axis.args=list(at=ticks, labels=tickLabels), 
               legend.mar=legend.mar, legend.width=legend.width)
    
    # image.plot(zlim=zlim, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
    #            col=cols, add = TRUE)
  }
  invisible(NULL)
}

# for computing what administrative regions the given points are in
# project: project to longitude/latitude coordinates
getRegion = function(points, mapDat = adm1, project=FALSE) {
  # project points to lon/lat coordinate system if user specifies
  if(project)
    points = projKenya(points, inverse=TRUE)
  
  regionNames = mapDat@data$NAME_1
  
  # make sure county names are consistent for mapDat == adm1
  regionNames[regionNames == "Elgeyo-Marakwet"] = "Elgeyo Marakwet"
  regionNames[regionNames == "Trans Nzoia"] = "Trans-Nzoia"
  
  # get region map polygons and set helper function for testing if points are in the regions
  polys = mapDat@polygons
  inRegion = function(i) {
    countyPolys = polys[[i]]@Polygons
    inside = sapply(1:length(countyPolys), function(x) {in.poly(points, countyPolys[[x]]@coords, inflation=0)})
    insideAny = apply(inside, 1, any)
    return(insideAny*i)
  }
  out = sapply(1:length(polys), inRegion)
  multipleRegs = apply(out, 1, function(vals) {sum(vals != 0) > 1})
  regionID = apply(out, 1, function(vals) {match(1, vals != 0)})
  regionNameVec = regionNames[regionID]
  list(regionID=regionID, regionNames=regionNameVec, multipleRegs=multipleRegs)
}

# convert county to region
countyToRegion = function(countyNames) {
  regionIs = match(as.character(countyNames), as.character(ctp[,1]))
  as.character(ctp[regionIs,2])
}

# do stratified sampling
# stratumVec: vector of stratum values, and each row being a sample
# sizes: vector with same length as unique(stratumVec) indicating the number of indices to sample
# probs: within each unique stratum level sample with probability proportional to probs (e.g. population density)
sampleByStratum = function(stratumVec, sizes, probs=1/length(stratumVec)) {
  uniqueStratumVec = unique(stratumVec)
  probMat = sapply(1:length(sizes), function(x) {probs*(stratumVec == uniqueStratumVec[x])/sum(probs*(stratumVec == uniqueStratumVec[x]))})
  unlist(sapply(1:length(sizes), function(x) {sample(1:length(stratumVec), size=sizes[x], prob=probMat[,x])}))
}

# do multilayer stratified sampling
# strataMat: each column is a stratum and each row is a sample
# sizes: vector of sizes of length equal to the number of unique 
#        combinations of values of levels in strataMat.  Number of 
#        elements to sample from that combination of strata levels
# probs: sample proportionally to probs within each strata level 
#        combination
sampleByStrata = function(strataMat, sizes, probs=1/nrow(strataMat)) {
  if(ncol(strataMat) == 1) 
    stop("only one stratum was input.  Use sampleByStratum in this case.")
  
  strataVec = strataMat[,1]
  for(i in 2:ncol(strataMat)) {
    strataVec = interaction(strataVec, strataMat[,i], drop=TRUE)
  }
  
  samples = sampleByStratum(strataVec, sizes, probs)
  
  list(strataVec=strataVec, samples=samples)
}

# generate a grid with a fixed latitude and longitude resolution.
# either set degrees per cell (res), number of grid cells on longest axis (nc), 
# or number of grid cells on both axes (nx and ny).
getKenyaGrid = function(res=.25, nc=NULL, nx=NULL, ny=NULL) {
  kenyaLonRange = c(33.5, 42)
  kenyaLatRange = c(-5,5.5)
  lonLen = kenyaLonRange[2] - kenyaLonRange[1]
  latLen = kenyaLatRange[2] - kenyaLatRange[1]
  
  # set individual axis grids
  if(!is.null(nc))
    res = latLen/nc
  if(is.null(nx) && is.null(ny)) {
    # we must go off of res
    lons = seq(kenyaLonRange[1], kenyaLonLength[2], by=res)
    lats = seq(kenyaLatRange[1], kenyaLatLength[2], by=res)
  }
  else {
    lons = seq(kenyaLonRange[1], kenyaLonLength[2], l=nx)
    lats = seq(kenyaLatRange[1], kenyaLatLength[2], l=ny)
  }
  
  # get full grid
  fullGrid = make.surface.grid(list(lon=lons, lat=lats))
  
  # return subset of points inside Kenya polygon
  polys = adm0@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  inKenya = in.poly(fullGrid, kenyaPoly)
  fullGrid[inKenya,]
}

##### put Kenya population density on a grid of the chosen resolution
# generate just the population density surface
makeKenyaPop = function(kmRes=5) {
  # load population density data
  require(raster)
  
  # pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
  load("Kenya2014Pop/pop.RData")
  
  # get a rectangular grid
  eastGrid = seq(eastLim[1], eastLim[2], by=kmRes)
  northGrid = seq(northLim[1], northLim[2], by=kmRes)
  utmGrid = make.surface.grid(list(east=eastGrid, north=northGrid))
  
  # project coordinates into lat/lon
  lonLatGrid = projKenya(utmGrid, inverse=TRUE)
  
  # subset grid so it's in Kenya
  polys = adm0@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  inKenya = in.poly(lonLatGrid, kenyaPoly)
  utmGrid = utmGrid[inKenya,]
  lonLatGrid = lonLatGrid[inKenya,]
  
  # get population density at those coordinates
  interpPopVals = extract(pop, SpatialPoints(lonLatGrid),method="bilinear")
  
  # compute counties associated with locations
  counties = getRegion(lonLatGrid, adm1)$regionNames
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], popOrig=interpPopVals, admin1=counties))
  
  newPop$east = utmGrid[,1]
  newPop$north = utmGrid[,2]
  
  newPop
}

# generate the population density surface along with urbanicity estimates
makeInterpPopGrid = function(kmRes=5, adjustPopSurface=FALSE, targetPop=c("children", "women")) {
  # load population density data
  require(raster)
  
  # pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
  load("Kenya2014Pop/pop.RData")
  
  # get a rectangular grid
  eastGrid = seq(eastLim[1], eastLim[2], by=kmRes)
  northGrid = seq(northLim[1], northLim[2], by=kmRes)
  utmGrid = make.surface.grid(list(east=eastGrid, north=northGrid))
  
  # project coordinates into lat/lon
  lonLatGrid = projKenya(utmGrid, inverse=TRUE)
  
  # subset grid so it's in Kenya
  polys = adm0@polygons
  kenyaPoly = polys[[1]]@Polygons[[77]]@coords
  inKenya = in.poly(lonLatGrid, kenyaPoly)
  utmGrid = utmGrid[inKenya,]
  lonLatGrid = lonLatGrid[inKenya,]
  
  # get population density at those coordinates
  interpPopVals = extract(pop, SpatialPoints(lonLatGrid),method="bilinear")
  
  # compute counties associated with locations
  counties = getRegion(lonLatGrid, adm1)$regionNames
  
  # determine which points are urban
  newPop = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2], popOrig=interpPopVals, admin1=counties))
  threshes = setThresholds2()
  popThreshes = sapply(1:nrow(newPop), function(i) {threshes$threshes[threshes$counties == newPop$admin1[i]]})
  urban = newPop$popOrig > popThreshes
  newPop$urban = urban
  
  newPop$east = utmGrid[,1]
  newPop$north = utmGrid[,2]
  
  # if necessary, adjust the population surface so that it better represents the the child population density 
  # rather than the total population density
  if(adjustPopSurface) {
    targetPop = match.arg(targetPop)
    
    # sort easpc by county name alphabetically
    counties=sort(unique(poppc$County))
    sortI = sort(easpc$County, index.return=TRUE)$ix
    temp = easpc[sortI,]
    
    # calculate the number of children per stratum using true total eas and empirical children per ea from census data
    load("empiricalDistributions.RData")
    if(targetPop == "children") {
      targetPopPerStratumUrban = temp$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
        ecdfExpectation(empiricalDistributions$childrenUrban)
      targetPopPerStratumRural = temp$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
        ecdfExpectation(empiricalDistributions$childrenRural)
    }
    else {
      targetPopPerStratumUrban = temp$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$womenUrban)
      targetPopPerStratumRural = temp$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$womenRural)
    }
    
    # generate 2 47 x nPixels matrices for urban and rural strata integrating pixels with respect to population density to get county estimates
    getCountyStratumIntegrationMatrix = function(getUrban=TRUE) {
      counties = as.character(counties)
      
      mat = t(sapply(counties, function(countyName) {newPop$admin1 == countyName}))
      mat = sweep(mat, 2, newPop$popOrig, "*")
      sweep(mat, 2, newPop$urban == getUrban, "*")
    }
    urbanIntegrationMat = getCountyStratumIntegrationMatrix()
    ruralIntegrationMat = getCountyStratumIntegrationMatrix(FALSE)
    
    # calculate number of people per stratum by integrating the population density surface
    urbanPopulations = rowSums(urbanIntegrationMat)
    ruralPopulations = rowSums(ruralIntegrationMat)
    
    # adjust each row of the integration matrices to get the correct expected number of children per stratum
    urbanIntegrationMat = sweep(urbanIntegrationMat, 1, targetPopPerStratumUrban / urbanPopulations, "*")
    ruralIntegrationMat = sweep(ruralIntegrationMat, 1, targetPopPerStratumRural / ruralPopulations, "*")
    ruralIntegrationMat[ruralPopulations == 0,] = 0
    
    # the column sums of the matrices give the correct modified population densities
    newPop$popOrig = colSums(urbanIntegrationMat) + colSums(ruralIntegrationMat)
  }
  
  newPop
}

testPopSurfaceAdjustment = function() {
  adjustedGrid = makeInterpPopGrid(adjustPopSurface = TRUE)
  
  # sort easpc by county name alphabetically
  counties=sort(unique(poppc$County))
  sortI = sort(easpc$County, index.return=TRUE)$ix
  temp = easpc[sortI,]
  
  # calculate the number of children per stratum using true total eas and empirical children per ea from census data
  load("empiricalDistributions.RData")
  householdsPerStratumUrban = temp$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban)
  householdsPerStratumRural = temp$EARur * ecdfExpectation(empiricalDistributions$householdsRural)
  childrenPerStratumUrban = householdsPerStratumUrban * ecdfExpectation(empiricalDistributions$mothersUrban) * 
    ecdfExpectation(empiricalDistributions$childrenUrban)
  childrenPerStratumRural = householdsPerStratumRural * ecdfExpectation(empiricalDistributions$mothersRural) * 
    ecdfExpectation(empiricalDistributions$childrenRural)
  
  # generate 2 47 x nPixels matrices for urban and rural strata integrating pixels with respect to population density to get county estimates
  getCountyStratumIntegrationMatrix = function(getUrban=TRUE) {
    counties = as.character(counties)
    
    mat = t(sapply(counties, function(countyName) {adjustedGrid$admin1 == countyName}))
    mat = sweep(mat, 2, adjustedGrid$popOrig, "*")
    sweep(mat, 2, adjustedGrid$urban == getUrban, "*")
  }
  urbanIntegrationMat = getCountyStratumIntegrationMatrix()
  ruralIntegrationMat = getCountyStratumIntegrationMatrix(FALSE)
  
  # calculate number of people per stratum by integrating the population density surface
  urbanPopulations = rowSums(urbanIntegrationMat)
  ruralPopulations = rowSums(ruralIntegrationMat)
  
  # load in the true enumeration areas
  margVar=.15^2
  gamma=-1
  beta0=-1.75
  tausq=0.1^2
  load(paste0("simDataMultiBeta", round(beta0, 4), "margVar", round(margVar, 4), "tausq", 
              round(tausq, 4), "gamma", round(gamma, 4), "HHoldVar0urbanOverSamplefrac0.RData"))
  
  # calculate the true number of households and children per stratum
  nHHUrban = aggregate(overSampDat$eaDat$nHH[overSampDat$eaDat$urban], by=list(overSampDat$eaDat$admin1[overSampDat$eaDat$urban]), FUN=sum)
  nHHRural = aggregate(overSampDat$eaDat$nHH[!overSampDat$eaDat$urban], by=list(overSampDat$eaDat$admin1[!overSampDat$eaDat$urban]), FUN=sum)
  nChildrenUrban = aggregate(overSampDat$eaDat$numChildren[overSampDat$eaDat$urban], by=list(overSampDat$eaDat$admin1[overSampDat$eaDat$urban]), FUN=sum)
  nChildrenRural = aggregate(overSampDat$eaDat$numChildren[!overSampDat$eaDat$urban], by=list(overSampDat$eaDat$admin1[!overSampDat$eaDat$urban]), FUN=sum)
  
  print(cbind(Expected=urbanPopulations, observed=nChildrenUrban$x, eaUrb=temp$EAUrb))
  
  genNumHouseholds = function(nSamplesPerEA, nEAs, urban=TRUE) {
    helperFun = function() {
      if(urban) {
        nHHs = sapply(nEAs, function(n) {sum(recdf(n, distribution=empiricalDistributions$householdsUrban))})
      }
      else {
        nHHs = sapply(nEAs, function(n) {sum(recdf(n, distribution=empiricalDistributions$householdsRural))})
      }
    }
    
    replicate(nSamplesPerEA, {helperFun()})
  }
  
  genNumChildren = function(nHouseholdsMat, urban=TRUE) {
    if(urban) {
      nMothersMat = matrix(sapply(nHouseholdsMat, function(n) {sum(recdf(n, distribution=empiricalDistributions$mothersUrban))}), ncol=ncol(nHouseholdsMat))
      matrix(sapply(nMothersMat, function(n) {sum(recdf(n, distribution=empiricalDistributions$childrenUrban))}), ncol=ncol(nHouseholdsMat))
    }
    else {
      nMothersMat = matrix(sapply(nHouseholdsMat, function(n) {sum(recdf(n, distribution=empiricalDistributions$mothersRural))}), ncol=ncol(nHouseholdsMat))
      matrix(sapply(nMothersMat, function(n) {sum(recdf(n, distribution=empiricalDistributions$childrenRural))}), ncol=ncol(nHouseholdsMat))
    }
  }
  
  ## calculate the number of households and the number of children per stratum
  # in urban areas
  testHouseholdsUrban = genNumHouseholds(2, temp$EAUrb)
  testChildrenUrban = genNumChildren(testHouseholdsUrban)
  
  # in rural areas
  testHouseholdsRural = genNumHouseholds(2, temp$EARur, FALSE)
  testHouseholdsRural = testHouseholdsRural[apply(testHouseholdsRural, 1, function(x) {any(x != 0)}),]
  testChildrenRural = genNumChildren(testHouseholdsRural, FALSE)
  
  pctEAUrban = temp$EAUrb/temp$EATotal
  
  print(cbind(datasetChildren=nChildrenUrban, simulatedChildren=testChildrenUrban, datasetHH=nHHUrban[,2], simulatedHH=testHouseholdsUrban, pctEAUrb=pctEAUrban, pctPopUrb=poppc$pctUrb[sortI]))
  print(cbind(datasetChildren=nChildrenRural, simulatedChildren=testChildrenRural, datasetHH=nHHRural[,2], simulatedHH=testHouseholdsRural, pctEARur=(1 - pctEAUrban)[pctEAUrban != 1], pctPopRur=(100-poppc$pctUrb[sortI])[pctEAUrban != 1]))
}

# takes the poppc table, containing the proportion of population that is urban and rural in each stratum, and
# adjusts it to be representative of the children in urban and rural areas per stratum based on census data
adjustPopulationPerCountyTable = function(dataType=c("children", "women")) {
  dataType = match.arg(dataType)
  
  # calculate the number of childrenor women per stratum using true total eas and empirical children per ea from census data
  load("empiricalDistributions.RData")
  if(dataType == "children") {
    targetPopPerStratumUrban = easpc$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
      ecdfExpectation(empiricalDistributions$childrenUrban)
    targetPopPerStratumRural = easpc$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
      ecdfExpectation(empiricalDistributions$childrenRural)
  }
  else {
    targetPopPerStratumUrban = easpc$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$womenUrban)
    targetPopPerStratumRural = easpc$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$womenRural)
  }
  
  
  # adjust poppc table to be representative of the number of children per stratum
  newPopTable = poppc
  targetPopPerCounty = targetPopPerStratumUrban + targetPopPerStratumRural
  newPopTable$popUrb = targetPopPerStratumUrban
  newPopTable$popRur = targetPopPerStratumRural
  newPopTable$popTotal = targetPopPerCounty
  newPopTable$pctUrb = newPopTable$popUrb / targetPopPerCounty  * 100
  newPopTable$pctTotal = newPopTable$popTotal/sum(newPopTable$popTotal) * 100
  
  # return results
  newPopTable
}

##### some random convenient function for breaking down number of urban counties for simulated 
##### and real datasets
numPerCounty = function(datSet, byUrban=TRUE, counties=clustpc$County) {
  if(!byUrban)
    results = aggregate(datSet$urban, list(datSet$admin1), length)
  else
    results = aggregate(datSet$urban, list(datSet$admin1), sum)
  sortI = match(counties, results$Group.1)
  out = matrix(results$x[sortI], ncol=1)
  rownames(out) = results$Group.1[sortI]
  out
}

propUrbPerCounty = function(datSet, counties=clustpc$County) {
  numPerCounty(datSet)/numPerCounty(datSet, FALSE, counties)
}

##### project from lat/lon to UTM northing/easting in kilometers.  Use epsg=21097
# either pass lon/east and lat/north, or a matrix with 2 columns: first being lon/east, second being lat/north
# inverse: if FALSE, projects from lon/lat to easting/northing.  Else from easting/northing to lon/lat
projKenya = function(lon, lat=NULL, inverse=FALSE) {
  if(is.null(lat)) {
    lat = lon[,2]
    lon = lon[,1]
  }
  
  if(!inverse) {
    # from lon/lat coords to easting/northing
    lonLatCoords = SpatialPoints(cbind(lon, lat), proj4string=CRS("+proj=longlat"))
    coordsUTM = spTransform(lonLatCoords, CRS("+init=epsg:21097 +units=km"))
    out = attr(coordsUTM, "coords")
  }
  else {
    # from easting/northing coords to lon/lat
    east = lon
    north = lat
    coordsUTM = SpatialPoints(cbind(east, north), proj4string=CRS("+init=epsg:21097 +units=km"))
    lonLatCoords = spTransform(coordsUTM, CRS("+proj=longlat"))
    out = attr(lonLatCoords, "coords")
  }
  
  out
}

##### convert Binomial data to Bernoulli
# inputs:
# binDat: a binomial data frame
# outputs:
# each row of binDat is expanded into multiple Bernoulli trials.  Additional 
# variables are added to the data
binToBern = function(binDat, nVar="numChildren", yVar="died", weightVar="samplingWeight") {
  tryCatch({weights = binDat[,weightVar]})
  # dataForDirect.Rdata
}

##### simulate from a logistic-normal-binomial distribution
rLogisticNormBin = function(nsim, n=1, logitProb, logitProbVar=0) {
  # first make the parameter vectors the correct length (they should each be nsim long)
  makeRightLength = function(vals) {
    currentLength = length(vals)
    numReps = nsim/currentLength
    
    # make sure the user didn't mess up
    if(round(numReps) != numReps)
      stop("nsim not multiple of all lengths of argument vectors")
    else if(numReps < 1)
      stop("nsim has longer length than other arguments")
    
    rep(vals, numReps)
  }
  n = makeRightLength(n)
  logitProb = makeRightLength(logitProb)
  logitProbVar = makeRightLength(logitProbVar)
  
  # simulate matrix of all bernoulli values in each binomial random variable
  N = sum(n)
  bernVar = rep(logitProbVar, n)
  bernProb = expit(rep(logitProb, n) + rnorm(N, 0, sqrt(bernVar)))
  bernVals = rbinom(sum(n), rep(1, N), bernProb)
  
  # function for simulating individual logistic normal binomial random variable from bernoulli values
  binIDs = as.factor(rep(1:length(n), n))
  # aggregate(bernVals, list(binIDs=binIDs), sum)$x # apparently the tapply function is wayyyyy faster than aggregate...
  tapply(bernVals, binIDs, sum)
}

# parallel matrix multiply
# leftMat %*% rightMat
# cl: the cluster, as initialized via:
# cores=detectCores()
# cl <- makeCluster(cores[1]-1) #not to overload your computer
# registerDoParallel(cl)
parMatMult = function(leftMat, rightMat, cl) {
  if (ncol(leftMat) != nrow(rightMat)) stop("Matrices do not conforme")
  idx   <- splitIndices(nrow(leftMat), length(cl))
  leftMatlist <- lapply(idx, function(ii) leftMat[ii,,drop=FALSE])
  ## ans   <- clusterApply(cl, leftMatlist, function(aa, rightMat) aa %*% rightMat, rightMat)
  ## Same as above, but faster:
  ans   <- clusterApply(cl, leftMatlist, get("%*%"), rightMat)
  do.call(rbind, ans)
}

# draw random numbers from an ecdf object
recdf = function(n, distribution) {
  probs = runif(n)
  quantile(distribution, probs, type=1) # type=1 signifies inverse ecdf
}

# get expectation of an ecdf object
ecdfExpectation = function(distribution) {
  distributionKnots = knots(distribution)
  distributionKnots = c(distributionKnots[1] - 1, distributionKnots)
  probs = distribution(distributionKnots[2:length(distributionKnots)]) - distribution(distributionKnots[1:(length(distributionKnots) - 1)])
  sum(distributionKnots[2:length(distributionKnots)] * probs)
}

# this function generates the sampling 
# probabilities of EAs within this strata
# i: the strata index, 1-47, corresponding to the row of easpc
# Nsamples: the total number of EAs drawn among strata with this urban/rural value
# m: the minimum number of draws in each strata
# urban: whether or not we are sampling from urban or rural strata
# log: return log probability or just probability
getStratumSamplingProb = function(i, nHHi, nHHTotal, Nsamples=617, m=3, urban=TRUE, dolog=FALSE) {
  # the number of households in this stratum (the commented out line should be the same as the one below)
  # Ni = ifelse(urban, easpc[i, "HHUrb"], easpc[i, "HHRur"])
  Ni = sum(nHHi)
  
  # the total number of households in Kenya that are either urban or rural
  Ntotal = ifelse(urban, sum(easpc$HHUrb), sum(easpc$HHRur))
  
  # the total number of strata with the given urban/rural categorization
  Nstrata = ifelse(urban, 47, 45)
  
  # the probability that you pick a particular ea in the first m draws from this stratum
  # pickEAinFirstM = dhyper(1, 1, Ni - 1, m)
  pickEAinFirstM = m * nHHi / Ni
  
  # the probability that you pick a particular ea in k additional draws from the stratum, 
  # given it wasn't drawn in the first m draws
  pickEAInK = function(k) {
    # successStates = Ni - m
    # totalStates = Ntotal - m * Nstrata
    # failureStates = totalStates - successStates
    totalDrawsLeft = Nsamples - m * Nstrata
    # probDrawStrataKTimes = dhyper(k, successStates, failureStates, totalDrawsLeft, log=TRUE)
    probDrawStrataKTimes = dbinom(k, totalDrawsLeft, Ni / Ntotal, log=TRUE)
    
    # countySuccessStates = 1
    # countyTotalStates = Ni - m
    # countyFailureStates = countyTotalStates - countySuccessStates
    # probDrawEaOutOfK = dhyper(1, countySuccessStates, countyFailureStates, k, log=TRUE)
    probDrawEaOutOfK = log(k * nHHi / Ni)
    
    # multiply the probabilities by summing on the log scale
    exp(probDrawStrataKTimes + probDrawEaOutOfK)
  }
  
  # the probability of picking this ea after the first m are drawn from the stratum 
  # is the probability it is drawn after the first m given that it was not drawn in 
  # the first m times the probability that it was not drawn in the first m
  maxDraws = min(Ni - m, Nsamples - m*Nstrata)
  pickEAafterFirstM = sum(sapply(1:maxDraws, pickEAInK)) * (1 - pickEAinFirstM)
  
  # probabilities of mutually exclusive events are summed
  pickEA = pickEAinFirstM + pickEAafterFirstM
  
  ifelse(dolog, log(pickEA), pickEA)
}

# this function generates the sampling probabilities of all households
# i: the strata index, 1-47, corresponding to the row of easpc
# Nsamples: the total number of EAs drawn among strata with this urban/rural value
# m: the minimum number of draws in each strata
# urban: whether or not we are sampling from urban or rural strata
# log: return log probability or just probability
getStratumSamplingProb2 = function(eaDat, eaDatLong, nSamplesUrban=617, nSamplesTotal=1612, m=3) {
  # the number of households in each stratum
  Nis = eaDat[, .(nHH=sum(nHH)), by=.(admin1, urban)]
  Ntotal = sum(Nis[["nHH"]])
  
  # the total number of households in Kenya that are either urban or rural
  Ntotal = ifelse(urban, sum(easpc$HHUrb), sum(easpc$HHRur))
  
  # the total number of strata with the given urban/rural categorization
  Nstrata = ifelse(urban, 47, 45)
  
  # the probability that you pick a particular ea in the first m draws from this stratum
  # pickEAinFirstM = dhyper(1, 1, Ni - 1, m)
  pickEAinFirstM = m * nHHi / Ni
  
  # the probability that you pick a particular ea in k additional draws from the stratum, 
  # given it wasn't drawn in the first m draws
  pickEAInK = function(k) {
    # successStates = Ni - m
    # totalStates = Ntotal - m * Nstrata
    # failureStates = totalStates - successStates
    totalDrawsLeft = Nsamples - m * Nstrata
    # probDrawStrataKTimes = dhyper(k, successStates, failureStates, totalDrawsLeft, log=TRUE)
    probDrawStrataKTimes = dbinom(k, totalDrawsLeft, Ni / Ntotal, log=TRUE)
    
    # countySuccessStates = 1
    # countyTotalStates = Ni - m
    # countyFailureStates = countyTotalStates - countySuccessStates
    # probDrawEaOutOfK = dhyper(1, countySuccessStates, countyFailureStates, k, log=TRUE)
    probDrawEaOutOfK = log(k * nHHi / Ni)
    
    # multiply the probabilities by summing on the log scale
    exp(probDrawStrataKTimes + probDrawEaOutOfK)
  }
  
  # the probability of picking this ea after the first m are drawn from the stratum 
  # is the probability it is drawn after the first m given that it was not drawn in 
  # the first m times the probability that it was not drawn in the first m
  maxDraws = min(Ni - m, Nsamples - m*Nstrata)
  pickEAafterFirstM = sum(sapply(1:maxDraws, pickEAInK)) * (1 - pickEAinFirstM)
  
  # probabilities of mutually exclusive events are summed
  pickEA = pickEAinFirstM + pickEAafterFirstM
  
  ifelse(dolog, log(pickEA), pickEA)
}

# get the probability of drawing a given household out of an enumeration area
getHHprob = function(nHH, nDraws=25, log=FALSE) {
  # dhyper(1, 1, nHH - 1, nDraws, log=log)
  nDraws/nHH
}

# match first occurrence of x in fun(values) using binary search assuming monotonicity of values
binarySearchMatch = function(x, values, fun=function(v, i) {v[i]}) {
  
  leftI = 1
  rightI = length(values)
  while(leftI < rightI) {
    matchI = floor((leftI + rightI) / 2)
    if(x <= fun(values, matchI))
      rightI = matchI
    else
      leftI = matchI + 1
  }
  
  if(x <= fun(values, leftI))
    leftI
  else
    NA
}

matchMultiple = function(vals, x) {
  x.map <- split(1:length(x), match(x, vals))
  # write a wrapper function that does a look-up on the unique list. and then returns all matches using the map.
  matchMake <- function(x) { x.map[[ fmatch(x, vals) ]] }
  sapply(vals, function(val) {matchMake(val)})
}

matchMultiple2 = function(vals, x) {
  dt_sol <- data.table(b, row = 1:length(b))[data.table(a), .(a, idx = row), on = .(b = a)][, .(idxs = list(idx)), by = a];
  dt_sol2 <- setNames(dt_sol$idxs, dt_sol$a)
  dt_sol2
}

# get expected children per cluster assuming independence of children per mother and mothers 
# per household conditional on urbanicity
getExpectedChildrenPerEa = function(distributionFile="empiricalDistributions.RData") {
  # load empirical distributions
  load(distributionFile)
  
  childrenUrban = ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
    ecdfExpectation(empiricalDistributions$childrenUrban) # 42.06381
  childrenRural = ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
    ecdfExpectation(empiricalDistributions$childrenRural) # 61.15194
  list(childrenPerClusterUrban=childrenUrban, childrenPerClusterRural=childrenRural)
}

# based on fields::discretize.image. Generates grid breaks from irregularly spaced data for producing an image
get2dCuts = function (x, m = 64, n = 64, grid = NULL, expand = c(1 + 1e-08, 1 + 1e-08), 
                    boundary.grid = FALSE, na.rm = TRUE) {
  if (length(expand) == 1) {
    expand <- rep(expand, 2)
  }
  if (is.null(grid)) {
    xr <- range(x[, 1], na.rm = na.rm)
    deltemp <- (xr[2] - xr[1]) * (expand[1] - 1) * 0.5
    gridX <- seq(xr[1] - deltemp, xr[2] + deltemp, , m)
    yr <- range(x[, 2], na.rm = na.rm)
    deltemp <- (yr[2] - yr[1]) * (expand[2] - 1) * 0.5
    gridY <- seq(yr[1] - deltemp, yr[2] + deltemp, , n)
    grid <- list(x = gridX, y = gridY)
    if (boundary.grid) {
      grid$x <- fields.convert.grid(grid$x)
      grid$y <- fields.convert.grid(grid$y)
    }
  }
  if (!boundary.grid) {
    xcut <- fields.convert.grid(grid$x)
    ycut <- fields.convert.grid(grid$y)
  }
  else {
    xcut <- grid$x
    ycut <- grid$y
  }
  list(x=grid$x, y=grid$y, xbreaks=xcut, ybreaks=ycut)
}

# faster version of fields::as.image. Converts irregularly spaced data into an image. 
# This implementation uses data.table
# Z: a vector of values FUN will be operated on
# x: a two column matrix of, for example, spatial coordinates. 2D bins will be constructed based on this, 
#    and FUN will be operated on the Z values within each of these bins
# weights: the weight of each observation, Z. This is currently not influencing FUN(Z) in any way, 
#          but the weights are summed within each bin, sort the result is the number of observations in each bin
my.as.image = function (Z, ind = NULL, grid = NULL, x = NULL, weights = rep(1, length(Z)), 
                        na.rm = FALSE, nx = 64, ny = 64, boundary.grid = FALSE, 
                        nrow = NULL, ncol = NULL, FUN = NULL, getUniqueN=getUniqueN) {
  Z <- c(Z)
  if (!is.null(ind)) {
    x <- ind
  }
  if (!is.null(nrow) & !is.null(ncol)) {
    nx <- nrow
    ny <- ncol
  }
  if (any(is.na(weights)) | any(is.na(c(x)))) {
    stop("missing values in weights or x")
  }
  # construct the grid of bins in which to evaluate FUN
  gridInfo <- get2dCuts(x, m = nx, n = ny, grid = grid, boundary.grid = boundary.grid)
  xbreaks = gridInfo$xbreaks
  ybreaks = gridInfo$ybreaks
  xcenters = gridInfo$x
  ycenters = gridInfo$y
  
  # group data by the breaks
  # xBinI = cut(x[,1], xbreaks, labels=FALSE)
  xBinI = findInterval(x[,1], xbreaks, left.open = FALSE, rightmost.closed = TRUE)
  # yBinI = cut(x[,2], ybreaks, labels=FALSE)
  yBinI = findInterval(x[,2], ybreaks, left.open = FALSE, rightmost.closed = TRUE)
  xvalues = xcenters[xBinI]
  yvalues = ycenters[yBinI]
  
  # construct the data table
  thisTable = data.table(x=xvalues, y=yvalues, z=Z, xBinI=xBinI, yBinI=yBinI, weights=weights)
  if(!getUniqueN)
    out = thisTable[,.(aggZ=FUN(z), weightSum=sum(weights)), by=.(xBinI, yBinI)]
  else
    out = thisTable[,.(aggZ=uniqueN(z), weightSum=sum(weights)), by=.(xBinI, yBinI)]
  zMat = matrix(nrow=length(xcenters), ncol=length(ycenters))
  zMat[cbind(out$xBinI, out$yBinI)] = out$aggZ
  weightMat = matrix(nrow=length(xcenters), ncol=length(ycenters))
  weightMat[cbind(out$xBinI, out$yBinI)] = out$weightSum
  # return results
  list(x = xcenters, y = ycenters, z = zMat, weights = weightMat)
}

my.quilt.plot = function (x, y, z, nx = 64, ny = 64, grid = NULL, add.legend = TRUE, 
                          add = FALSE, nlevel = 64, col = tim.colors(nlevel), nrow = NULL, 
                          ncol = NULL, FUN = NULL, plot = TRUE, na.rm = FALSE, getUniqueN=FALSE, ...) 
{
  if(is.null(FUN) && !getUniqueN)
    FUN = mean
  if (!is.null(nrow) | !is.null(nrow)) {
    nx <- nrow
    ny <- ncol
  }
  x <- as.matrix(x)
  if (ncol(x) == 2) {
    z <- y
  }
  if (ncol(x) == 1) {
    x <- cbind(x, y)
  }
  if (ncol(x) == 3) {
    z <- x[, 3]
    x <- x[, 1:2]
  }
  out.p <- my.as.image(z, x = x, nx = nx, ny = ny, grid = grid, 
                    FUN = FUN, na.rm = na.rm, getUniqueN=getUniqueN)
  if (plot) {
    if (add.legend) {
      image.plot(out.p, nlevel = nlevel, col = col, add = add, 
                 ...)
    }
    else {
      image(out.p, col = col, add = add, ...)
    }
  }
  invisible(out.p)
}

##### MVN log likelihood for zero mean multivariate normal
## inputs:
# dat: observations.  If multiple realizations, dat is a matrix with each column being a 
#      different set of observations for one realization
# SimgaU: upper triangular Cholesky decomp of covariance matrix
logLikGP = function(dat, SigmaU) {
  if(!is.matrix(dat))
    dat = matrix(dat, nrow=length(dat))
  n = nrow(dat)
  
  #calcualte GP log likelihood (log likelihood of S0j)
  # (-1/2) t(y) %*% Sigma0^-1 y - (1/2) log det Sigma0 - (n/2) log(2pi)
  # define z = U^-1 %*% y, and x = L^-1 %*% U^-1 %*% y
  # then:
  z = backsolve(SigmaU, dat, transpose=TRUE)
  x = backsolve(SigmaU, z) #make sure x is a column
  # get full loglik.  First sum() only used for multiple realizations
  log.likGP = sum(-(1/2) * colSums(dat * x) - sum(log(diag(SigmaU))) - (n/2)*log(2*pi))
  
  return(log.likGP)
}

# round a set of numbers to a certain number of digits passed the first significant digit 
# out of all the numbers
roundToFirstSigFigs = function(numbers, digits=3) {
  # count the number of zeros after the decimal point
  count0 <- function(x, tol = .Machine$double.eps ^ 0.5) { 
    x <- abs(x)
    y <- -log10(x - floor(x))
    floor(y) - (y %% 1 < tol)
  }
  
  # count the number of zeros pass the decimal point for each number, determine the number 
  # of significant figures to round to for each number
  zeros = count0(numbers)
  minNumber = min(zeros)
  extraZeros = zeros - minNumber
  finalDigits = sapply(extraZeros, function(x) {max(c(0, digits-x))})
  signif(numbers, finalDigits)
}

# same as the format function operated on tables, except in forces that some rows of numbers 
# are shown in scientific format if any of those numbers have zeros to the right of the decimal 
# point than the input maxZeros
formatSomeScientific = function(numberTable, maxZeros=4, digits=3, ...) {
  # count the number of zeros after the decimal point
  count0 <- function(x, tol = .Machine$double.eps ^ 0.5) { 
    x <- abs(x)
    y <- -log10(x - floor(x))
    floor(y) - (y %% 1 < tol)
  }
  
  # count the number of zeros pass the decimal point for each number, determine the number 
  # of significant figures to round to for each number
  zeros = count0(numberTable)
  rowZeros = apply(zeros, 1, max)
  doScientific = rowZeros > maxZeros
  temp = numberTable
  numberTable[!doScientific,] = format(temp[!doScientific,], scientific=FALSE, digits=digits, ...)
  numberTable[doScientific,] = format(temp[doScientific,], scientific=TRUE, digits=digits+1, ...)
  numberTable
}

getRadius = function(areaLevel=c("Region", "County")) {
  require(geosphere)
  require(shapefiles)
  areaLevel = match.arg(areaLevel)
  
  # load shape files
  require(maptools)
  if(areaLevel == "Region"){
    thisMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  } else if(areaLevel == "County"){
    out = load("../U5MR/adminMapData.RData")
    thisMap = adm1
  } else {
    stop(paste0("Unrecognized area level: ", areaLevel))
  }
  
  # thisMap@polygons[[1]]@Polygons[[1]]@coords
  # thisMap@polygons[[1]]@Polygons[[1]]@area
  getOneRadius = function(poly) {
    areas = sapply(poly@Polygons, function(x) {x@area})
    correctPolyI = which.max(areas)
    poly = poly@Polygons[[correctPolyI]]
    # thisCentroid = centroid(poly@coords)
    # allPoints = rbind(thisCentroid, 
    #                   poly@coords)
    allPoints = poly@coords
    allProjected = projKenya(allPoints[,1], allPoints[,2], inverse=FALSE)
    
    # downsample spatial polygons as necessary
    if(nrow(allProjected) >= 6000) {
      # simplify the polygon
      newProjected = allProjected
      tolerance = .05
      while(nrow(newProjected) >= 6000) {
        newProjected = do.call("cbind", dp(list(x=allProjected[,1], y=allProjected[,2]), tolerance=tolerance))
        tolerance = tolerance * 2
      }
      allProjected = newProjected
    }
    max(rdist(allProjected)) / 2
  }
  # calculate radii in km
  radii = sapply(thisMap@polygons, getOneRadius)
  
  # sort results by area name
  if(areaLevel == "Region") {
    areaNames = as.character(thisMap@data$name)
  } else if(areaLevel == "County") {
    areaNames = as.character(thisMap@data$NAME_1)
  }
  sortI = sort(areaNames, index.return=TRUE)$ix
  radii = radii[sortI]
  names(radii) = areaNames[sortI]
  radii
}

getMeanRadius = function(areaLevel=c("Region", "County")) {
  require(geosphere)
  require(shapefiles)
  areaLevel = match.arg(areaLevel)
  
  # load shape files
  require(maptools)
  if(areaLevel == "Region"){
    thisMap = readShapePoly("../U5MR/mapData/kenya_region_shapefile/kenya_region_shapefile.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
  } else if(areaLevel == "County"){
    out = load("../U5MR/adminMapData.RData")
    thisMap = adm1
  } else {
    stop(paste0("Unrecognized area level: ", areaLevel))
  }
  
  # thisMap@polygons[[1]]@Polygons[[1]]@coords
  # thisMap@polygons[[1]]@Polygons[[1]]@area
  getOneRadius = function(poly) {
    areas = sapply(poly@Polygons, function(x) {x@area})
    correctPolyI = which.max(areas)
    poly = poly@Polygons[[correctPolyI]]
    thisCentroid = centroid(poly@coords)
    allPoints = rbind(thisCentroid,
                      poly@coords)
    # allPoints = poly@coords
    allProjected = projKenya(allPoints[,1], allPoints[,2], inverse=FALSE)
    
    mean(rdist.vec(matrix(allProjected[1,], nrow=nrow(allProjected)-1, ncol=2, byrow=TRUE), allProjected[-1,]))
  }
  # calculate radii in km
  radii = sapply(thisMap@polygons, getOneRadius)
  
  # sort results by area name
  if(areaLevel == "Region") {
    areaNames = as.character(thisMap@data$name)
  } else if(areaLevel == "County") {
    areaNames = as.character(thisMap@data$NAME_1)
  }
  sortI = sort(areaNames, index.return=TRUE)$ix
  radii = radii[sortI]
  names(radii) = areaNames[sortI]
  radii
}