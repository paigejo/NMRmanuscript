# functions for matching clusters to EA or pixel locations
matchPoints = function(clustPoints, otherPoints, doRound=FALSE) {
  if(doRound) {
    seqX = sort(unique(otherPoints[,1]))
    seqY = sort(unique(otherPoints[,2]))
    xs = roundToNearest(clustPoints[,1], seqX)
    ys = roundToNearest(clustPoints[,2], seqY)
    clustPoints = cbind(xs, ys)
  }
  
  indices = findIndex(clustPoints, otherPoints)
  if(any(is.na(indices))) {
    # just take nearest one for these.  Maybe be due to rounding, irregular map border
    naClustPoints = matrix(clustPoints[is.na(indices),], ncol=2)
    distMat = rdist(naClustPoints, otherPoints)
    naIndices = apply(distMat, 1, which.min)
    indices[is.na(indices)] = naIndices
  }
  
  return(as.numeric(indices))
}

# for each element of vals, rounds it to the nearest value in seqVals, which 
# are assumed to be equidistant
roundToNearest = function(vals, seqVals, returnIndices=FALSE) {
  # get basic info about the grid
  width = seqVals[2] - seqVals[1]
  minVal = seqVals[1]
  maxVal = seqVals[length(seqVals)]
  
  # threshold values
  if(any(vals < maxVal)) {
    warning("some coordinates larger than maximum of grid")
    vals[vals > maxVal] = maxVal
  }
  if(any(vals < minVal)) {
    warning("some coordinates smaller than minimum of grid")
    vals[vals < minVal] = minVal
  }
  
  # round the values to the grid
  indices = round((vals - minVal)/width) + 1
  result = seqVals[indices]
  
  if(returnIndices)
    return(list(values=result, indices=indices))
  else
    return(result)
}

# fast function for getting indices of coords1 in coords2
findIndex = function(coords1, coords2) {
  match(as.list(data.frame(t(coords1))), as.list(data.frame(t(coords2))))
}

# add EAIndex to saved datasets
addEAPixelIndex = function(clustToEAI = FALSE) {
  # load dataset, set fixed coordinates of points to match
  out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  eaDat = SRSDat$eaDat
  eaCoords = cbind(eaDat$east, eaDat$north)
  pixelCoords = cbind(popGrid$east, popGrid$north)
  clustDatSRS = SRSDat$clustDat
  clustDatOverSamp = overSampDat$clustDat
  
  if(clustToEAI) {
    # now match the cluster coordinates with EA coordinates and pixel coordinates
    for(i in 1:100) {
      print(i)
      srsCoords = cbind(clustDatSRS[[i]]$east, clustDatSRS[[i]]$north)
      overSampCoords = cbind(clustDatOverSamp[[i]]$east, clustDatOverSamp[[i]]$north)
      clustDatSRS[[i]]$eaI = matchPoints(srsCoords, eaCoords)
      clustDatSRS[[i]]$pixelI = matchPoints(srsCoords, pixelCoords, doRound=TRUE)
      clustDatOverSamp[[i]]$eaI = matchPoints(overSampCoords, eaCoords)
      clustDatOverSamp[[i]]$pixelI = matchPoints(overSampCoords, pixelCoords, doRound=TRUE)
    }
  }
  
  # match ea coordinates with pixel grid:
  print("matching EA coordinates with pixel grid")
  srsI = matchPoints(cbind(SRSDat$eaDat$east, SRSDat$eaDat$north), pixelCoords, doRound=TRUE)
  SRSDat$eaDat$pixelI = srsI
  overI = matchPoints(cbind(overSampDat$eaDat$east, overSampDat$eaDat$north), pixelCoords, doRound=TRUE)
  overSampDat$eaDat$pixelI = overI
  
  # save results
  SRSDat$clustDat = clustDatSRS
  overSampDat$clustDat = clustDatOverSamp
  print("Done with cluster case")
  save(SRSDat, overSampDat, file="simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  
  # OR: cluster variance = 0 (no cluster effect)
  out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
  eaDat = SRSDat$eaDat
  eaCoords = cbind(eaDat$east, eaDat$north)
  pixelCoords = cbind(popGrid$east, popGrid$north)
  clustDatSRS = SRSDat$clustDat
  clustDatOverSamp = overSampDat$clustDat
  
  if(clustToEAI) {
    # now match the cluster coordinates with EA coordinates and pixel coordinates
    for(i in 1:100) {
      print(i)
      srsCoords = cbind(clustDatSRS[[i]]$east, clustDatSRS[[i]]$north)
      overSampCoords = cbind(clustDatOverSamp[[i]]$east, clustDatOverSamp[[i]]$north)
      clustDatSRS[[i]]$eaI = matchPoints(srsCoords, eaCoords)
      clustDatSRS[[i]]$pixelI = matchPoints(srsCoords, pixelCoords, doRound=TRUE)
      clustDatOverSamp[[i]]$eaI = matchPoints(overSampCoords, eaCoords)
      clustDatOverSamp[[i]]$pixelI = matchPoints(overSampCoords, pixelCoords, doRound=TRUE)
    }
  }
  
  srsI = matchPoints(cbind(SRSDat$eaDat$east, SRSDat$eaDat$north), pixelCoords, doRound=TRUE)
  SRSDat$eaDat$pixelI = srsI
  overI = matchPoints(cbind(overSampDat$eaDat$east, overSampDat$eaDat$north), pixelCoords, doRound=TRUE)
  overSampDat$eaDat$pixelI = overI
  
  # save results
  SRSDat$clustDat = clustDatSRS
  overSampDat$clustDat = clustDatOverSamp
  print("Done with no cluster case")
  save(SRSDat, overSampDat, file="simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
}

# for testing to make sure indices are correct:
# for(i in 1:100) {
#   # print(length(overSampDat$clustDat[[i]]$pixelI))
#   print(sum(overSampDat$eaDat[overSampDat$clustDat[[i]]$eaI,1:15] != overSampDat$clustDat[[i]][,1:15]))
#   print(sum(SRSDat$eaDat[SRSDat$clustDat[[i]]$eaI,1:15] != SRSDat$clustDat[[i]][,1:15]))
# }

matchEAToPixel = function(eaDat, maxRows = 100, kmres=5) {
  if(kmres == 5)
    load("popGrid.RData")
  else
    popGrid = makeInterpPopGrid(kmres)
  
  eaCoords = cbind(eaDat$east, eaDat$north)
  popCoords = cbind(popGrid$east, popGrid$north)
  
  matchToClosest = function(i) {
    startI = maxRows * (i - 1) + 1
    endI = min(startI + maxRows - 1, nrow(eaCoords))
    print(paste0("Matching up to row ", endI, "/", nrow(eaCoords)))
    
    eaCols = eaCoords[startI:endI,]
    distMat = rdist(popCoords, eaCols)
    pixelIs = apply(distMat, 2, which.min)
    
    pixelIs
  }
  
  c(unlist(sapply(1:ceiling(nrow(eaCoords) / maxRows), matchToClosest)))
}

addPixelIToEaDat = function(kmres=5) {
  # load the ea and cluster data
  out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  eaDat = SRSDat$eaDat
  clustDatSRS = SRSDat$clustDat
  clustDatOverSamp = overSampDat$clustDat
  
  # generate pixel indices
  pixelI = matchEAToPixel(eaDat, kmres=kmres)
  SRSDat$eaDat$pixelI = pixelI
  overSampDat$eaDat$pixelI = pixelI
  
  # add pixel indices to the cluster data
  for(i in 1:100) {
    SRSDat$clustDat[[i]]$pixelI = SRSDat$eaDat[SRSDat$clustDat[[i]]$eaIs,]$pixelI
    overSampDat$clustDat[[i]]$pixelI = overSampDat$eaDat[overSampDat$clustDat[[i]]$eaIs,]$pixelI
  }
  
  # save results
  print("Done with SRS case")
  save(SRSDat, overSampDat, file="simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  
  # OR: cluster variance = 0 (no cluster effect)
  # load the ea and cluster data
  out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
  eaDat = SRSDat$eaDat
  clustDatSRS = SRSDat$clustDat
  clustDatOverSamp = overSampDat$clustDat
  
  # generate pixel indices
  pixelI = matchEAToPixel(eaDat)
  SRSDat$eaDat$pixelI = pixelI
  overSampDat$eaDat$pixelI = pixelI
  
  # add pixel indices to the cluster data
  for(i in 1:100) {
    SRSDat$clustDat[[i]]$pixelI = SRSDat$eaDat[SRSDat$clustDat[[i]]$eaIs,]$pixelI
    overSampDat$clustDat[[i]]$pixelI = overSampDat$eaDat[overSampDat$clustDat[[i]]$eaIs,]$pixelI
  }
  
  # save results
  print("Done with SRS case")
  save(SRSDat, overSampDat, file="simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
}

# same as addPixelIToEaDat2, but returns modified clustDat instead of overwriting files
addPixelIToEaDat2 = function(clustDat, kmres=5) {
  # load the ea and cluster data
  out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
  eaDat = clustDat$eaDat
  clustDat = clustDat$clustDat
  
  # generate pixel indices
  pixelI = matchEAToPixel(eaDat, kmres=kmres)
  eaDat$pixelI = pixelI
  
  # add pixel indices to the cluster data
  for(i in 1:length(clustDat)) {
    clustDat[[i]]$pixelI = eaDat[clustDat[[i]]$eaIs,]$pixelI
  }
  
  # return results
  list(eaDat=eaDat, clustDat=clustDat)
}

