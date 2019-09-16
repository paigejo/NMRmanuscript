# Code for INLA SPDE model:

# Make the mesh grid for INLA's SPDE spatial model.  Accounts for spherical coordinate system.  
# Note that max.n and globe don't seem to affect the grid, and by default 3779 vertices are 
# generated.
# doPlot: If doPlot is true, plots the generated mesh
# offset: argument to inla.mesh.create
# doGlobe: if TRUE, makes grid on the sphere, assumes coordinates are lat/lon.  
#          Otherwise, assumes coordinates are esting/northing.  If FALSE, a good default 
#          for locs is locs=projKenya(ed$lon, ed$lat)
# cutoff: minimum leg size of trianglation
getSPDEMeshGrid = function(locs=NULL, max.n=100, doPlot=TRUE, 
                           doGlobe=FALSE, offset=-.18, cutoff=ifelse(doGlobe, .15, 15)) {
  # base mesh off of actual DHS dataset unless user specifies otherwise
  if(is.null(locs)) {
    locs = cbind(ed$lon, ed$lat)
    if(!doGlobe) {
      # locs = projKenya(locs)
      locs = cbind(ed$east, ed$north)
    }
  }
  
  if(doGlobe) {
    # generate the mesh on the globe
    mesh = inla.mesh.create(locs, globe=10, refine=list(max.n=max.n), extend=list(offset=offset), 
                            cutoff=cutoff)
  }
  else {
    # generate mesh on R2
    mesh = inla.mesh.create(loc=locs, refine=list(max.n=max.n), extend=list(offset=offset), 
                            cutoff=cutoff)
  }
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
    plotMapDat(project=!doGlobe, border="blue")
  }
  
  mesh
}

getSPDEMeshGrid2 = function(locs=cbind(ed$east, ed$north), n=1500, max.n=2000, doPlot=TRUE, max.edge=c(25, 250), 
                           offset=-.08, cutoff=15) {
  
  
  # generate mesh on R2
  mesh = inla.mesh.2d(loc=locs, n=n, max.n=max.n, offset=offset, cutoff=cutoff, max.edge=max.edge)
  
  # plot the mesh if user wants
  if(doPlot) {
    plot(mesh)
    plotMapDat(project=TRUE, border="blue")
  }
  
  mesh
}

# generate default priors for SPDE model
# from Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
# sigma0: field standard deviation
# NOTE: by default, this constructs a spde prior with unit median marginal variance 
#       and median effective range equal to a fifth of the spatial range 
# or use inla.spde2.pcmatern (possibly allow (1/4,4) variance here rather than (1/2,2))
getSPDEPrior = function(mesh, sigma0=1, strictPrior=FALSE) {
  size <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) # 1277.237
  # range0 <- size/5
  # kappa0 <- sqrt(8)/range0
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # spde <- inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),
  #                           B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0),
  #                           theta.prior.prec = c(0.1, 1))
  
  range0 <- size/5
  if(!strictPrior)
    spde = inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(1, 0.01))
  else
    spde = inla.spde2.pcmatern(mesh, prior.range=c(range0, 0.5), prior.sigma = c(.15, 0.01))
  spde
}

# fit the SPDE model given some data and some prediction locations
# Inputs:
# obsCoords: coordinates of the observations
# obsNs: the binomial n of the observations
# obsCounts: the response, education counts
# obsUrban: whether the observations were urban or knot
# predCoords: the spatial coordinates of the predictions: first the EAs, then the pixels
# predNs: the binomial n of the predictions
# predUrban: whether the prediction locations are urban or knot
# prior: the prior for the spde. Defaults to getSPDEPrior with the spatial mesh
# mesh: the spatial mesh. Defaults to something dat should be reasonable
# int.strategy: the introduction strategy argument for inla
# strategy: the strategy argument for inla
# genCountyLevel: whether or knot to generate county level predictions
# popGrid: population densities over the prediction grid for integration
# nPostSamples: number of samples from the joint posterior to take
# kmRes: the kilometers resolution of the pixels
# counties: the county names, determining the order of the county predictions
# clusterEffect: whether or not to include an iid gaussian error term at the EA level
# verbose: verbose argument to inla
# genRegionLevel: whether or not to generate regional level predictions. Not currently supported
# regions: the region names, determining the order of the region predictions
# keepPixelPreds: whether or not to return the pixel level predictions
# genEALevel: whether or not to generate ea level predictions
# eaIndices: indices giving which of the prediction coordinates in predCoords are EAs. NOTE: 
#            these must be the first indices of the predictions
# urbanEffect: if in urban fixed effect is included
# link: the link argument for inla. 1 means logit scale, 2 means probability scale (?)
# predictionType: either median or mean predictions are returned
# eaPixelIndices: the pixel indices each ea belongs to. This is of length equal to the number of EAs
# link: either one or two. If one, uses the identity link, otherwise uses
# exactAggregation: aggregate using the known ea locations rather than integrating over
#                   the aggregation area with respect to population density
# genCountLevel: whether or not to generate predictions at the count level versus logistic. DEPRECATED
# nSamplePixel: fewer samples are required for good approximation of the posterior at the pixel level 
#               than county, so we only take this many of the posterior samples for pixel level estimation
fitSPDEModel = function(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, predNs = rep(1, nrow(predCoords)), 
                        predUrban, prior=NULL, mesh=NULL, int.strategy="eb", strategy="gaussian", 
                        genCountyLevel=FALSE, popGrid=NULL, nPostSamples=100, kmRes=5, counties=sort(unique(eaDat$admin1)), 
                        verbose=TRUE, genRegionLevel=FALSE, regions=sort(unique(eaDat$region)), 
                        keepPixelPreds=FALSE, genEALevel=FALSE, eaIndices=1:nrow(kenyaEAs), 
                        urbanEffect=TRUE, link=1, predictionType=c("median", "mean"), eaDat=NULL, 
                        exactAggregation=FALSE, nSamplePixel=10, truthByPixel=NULL, 
                        truthByCounty=NULL, truthByRegion=NULL, truthByEa=NULL, clusterEffect=FALSE) {
  
  # match the prediction type
  predictionType = match.arg(predictionType)
  
  # for enumeration area level predictions, set enumeration area end pixel level prediction indices
  eaMat = NULL
  if(genEALevel) {
    pixelIndices = (1:nrow(predCoords))[-eaIndices]
  }
  else {
    pixelIndices = 1:nrow(predCoords)
  }
  
  # make default mesh
  if(is.null(mesh)) {
    mesh = getSPDEMeshGrid(doPlot = FALSE)
  }
  
  # make default prior
  if(is.null(prior)) {
    prior = getSPDEPrior(mesh)
  }
  
  # make covariate matrices (intercept plus urban/rural depending on whether urban effect included)
  if(urbanEffect) {
    X = cbind(1, obsUrban)
    XPred = cbind(1, predUrban)
  }
  else {
    X = matrix(1, ncol=1, nrow=nrow(obsCoords))
    XPred = matrix(1, ncol=1, nrow=nrow(predCoords))
  }
  
  
  # construct A matrix for observations
  m = nrow(obsCoords)
  AEst = inla.spde.make.A(mesh, loc = obsCoords)
  
  # construct A matrix for predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  
  # make inla stack
  ys = obsCounts
  n = ncol(AEst) # number of basis elements
  nObs = length(ys) # number of observations
  nPreds = nrow(predCoords)
  latticeInds = 1:n
  
  # do not include cluster effect at the pixel level
  stack.pred = inla.stack(A =list(APred[,pixelIndices], 1),
                          effects =list(field=pixelIndices, X=XPred[,pixelIndices]),
                          data =list(y=NA, Ntrials = predNs[pixelIndices], link=1),
                          tag ="pred",
                          remove.unused=FALSE)
  
  if(clusterEffect) {
    # clustIObs = (n+1):(n+nObs)
    # clustIPreds = (n+1):(n+nPreds)
    # clustIPreds = (n+nObs+1):(n+nObs+nPreds)
    clustIObs = 1:nObs
  }
  
  if(genEALevel) {
    
    if(clusterEffect) {
      
      # only include clusters in the enumeration area predictions
      clustIPreds = 1:nrow(eaDat)
      stack.predEA = inla.stack(A =list(APred[,eaIndices], 1, 1),
                                effects =list(field=eaIndices, clust=clustIPreds, X=XPred[eaIndices,]),
                                data =list(y=NA, Ntrials = predNs[eaIndices], link=1),
                                tag ="pred",
                                remove.unused=FALSE)
    } else {
      stack.predEA = inla.stack(A =list(APred[,eaIndices], 1),
                                effects =list(field=eaIndices, X=XPred[eaIndices,]),
                                data =list(y=NA, Ntrials = predNs[eaIndices], link=1),
                                tag ="pred",
                                remove.unused=FALSE)
    }
  }
  
  # construct the observation stack 
  if(clusterEffect) {
    stack.est = inla.stack(A =list(AEst, 1, 1),
                           effects =list(field=latticeInds, clust=clustIObs, X=X),
                           data =list(y=ys, Ntrials = obsNs, link=1),
                           tag ="est",
                           remove.unused=FALSE)
    
  } else {
    stack.est = inla.stack(A =list(AEst, 1), 
                           effects =list(field=latticeInds, X=X), 
                           data =list(y=ys, Ntrials = obsNs, link=link), 
                           tag ="est", 
                           remove.unused=FALSE)
  }
  
  # # whether or not cluster effect is included in the model, it should not be included in the pixel level predictions 
  # stack.pred = inla.stack(A =list(APred, 1),
  #                         effects =list(field=latticeInds, X=XPred),
  #                         data =list(y=NA, Ntrials = predNs, link=link),
  #                         tag ="pred",
  #                         remove.unused=FALSE)
  
  # construct the full prediction stack if necessary
  if(genEALevel) {
    stack.pred = inla.stack(stack.pred, stack.predEA)
  }
  
  # make mesh index
  mesh.index <- inla.spde.make.index(name = "field", n.spde = prior$n.spde)
  
  # fit model
  control.inla = list(strategy=strategy, int.strategy=int.strategy) 
  allNs = c(obsNs, predNs)
  stack.full = inla.stack(stack.est, stack.pred, remove.unused=FALSE)
  stackDat = inla.stack.data(stack.full, spde=prior)
  if(clusterEffect) {
    mod = inla(y ~ - 1 + X + f(field, model=prior) + 
                 f(clust, model="iid", 
                   hyper = list(prec = list(prior = "pc.prec", param = c(3,0.01)))),
               data = stackDat, Ntrials=stackDat$Ntrials,
               control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link),
               family="binomial", verbose=TRUE, control.inla=control.inla,
               control.compute=list(config=TRUE))
  } else {
    mod = inla(y ~ - 1 + X + f(field, model=prior), 
               data = stackDat, Ntrials=stackDat$Ntrials, 
               control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link), 
               family="binomial", verbose=verbose, control.inla=control.inla, 
               control.compute=list(config=TRUE))
  }
  
  # get predictive surface, SD, and data
  n = nrow(obsCoords)
  # obsInds = 1:n
  obsInds = inla.stack.index(stack.full, "est")$data
  predInds = inla.stack.index(stack.full, "pred")$data
  index = inla.stack.index(stack.full, "pred")$data
  if(predictionType == "mean")
    linpreds = mod[["summary.fitted.values"]]$mean
  else
    linpreds = mod[["summary.fitted.values"]]$`0.5quant`
  linpred.sd = mod[["summary.fitted.values"]]$sd
  
  # predInds = (n+1):(n+nrow(predCoords))
  preds = linpreds
  predSDs = linpred.sd
  
  countyPredMat = NULL
  if(genCountyLevel) {
    ## generate county level predictions
    
    # if not supplied, get grid of population densities for pop-weighted integration
    if(is.null(popGrid))
      popGrid = makeInterpPopGrid(kmRes=kmRes)
    
    # make sure prediction coordinates correspond to population density grid coordinates
    if(length(pixelIndices) != nrow(popGrid)) {
      stop("for county level predictions, prediction points must be on same grid as popGrid")
    }
    
    # generate samples from posterior
    postSamples = inla.posterior.sample(nPostSamples, mod)
    latentMat = sapply(postSamples, function(x) x$latent)
    
    # get logit education rates at the prediction locations
    index.pred <- inla.stack.index(stack.full,tag="pred")$data
    predMat <- latentMat[index.pred,]
    
    eaMarginals = NULL
    pixelMarginals = NULL
    if(genEALevel) {
      # make sure to separate enumeration area and pixel level predictions
      eaMat = expit(predMat[eaIndices, ])
      predMat = predMat[pixelIndices, ]
      
      # get the marginals of the enumeration areas
      eaMarginals = mod$marginals.linear.predictor[eaIndices]
    }
    if(keepPixelPreds)
      pixelMarginals = mod$marginals.linear.predictor[pixelIndices]
    
    # the below code is commented since we removed the cluster effect from the model. The commented code is not yet tested
    # Also, if the cluster effect is included, it should only be included 4 the enumeration area predictions, 
    # not the areal level predictions
    # # remove cluster effect for integration
    # nPreds = nrow(predCoords)
    # clusterIndices = inla.stack.index(stack.full,tag="pred")$effects[(1 + mesh$n):(mesh$n + nPreds)]
    # clusterEffects = latentMat[clusterIndices, ]
    # noClustMat = predMat - clusterEffects
    
    # county level predictions
    if(!exactAggregation) {
      if(genCountLevel)
        stop("not using exact aggregation is not supported if genCountLevel is set to TRUE")
      
      # for each county, integrate predictions with respect to population density
      predCounties = popGrid$admin1
      pops = popGrid$popOrig
      
      integratePredsByCounty = function(countyName) {
        # subset data by county of interest
        countyI = predCounties == countyName
        countyPredMat = expit(predMat[countyI,])
        countyPops = pops[countyI]
        
        # compute average integral of predictions
        countyPops = countyPops/sum(countyPops)
        integrals = t(countyPredMat) %*% countyPops
        integrals
      }
      countyPredMat <- t(sapply(counties, integratePredsByCounty))
    } else {
      
      # for each county, integrate predictions over EAs
      if(is.null(eaDat))
        stop("eaDat is null, but a non-null eaDat must be provided for account level predictions for EAs.")
      predCounties = eaDat$admin1
      
      integratePredsByCounty = function(countyName) {
        # subset data by county of interest
        countyI = predCounties == countyName
        countyProbMat = matrix(eaMat[eaIndices[countyI],], nrow=sum(countyI))
        
        # combine results by EA
        # numClusters = sum(countyI)
        numChildren = eaDat$numChildren[countyI]
        distribution = dSumBinomRandom(0:sum(numChildren), numChildren, countyProbMat)
        distribution
      }
      
      # in this case, countyPredMat is a data frame with predictions and confidence intervals, with rejection probabilities
      distributions <- lapply(counties, integratePredsByCounty)
      countyPredMatExact <- matrix(sapply(distributions, function(masses) {sum(masses * (0:(length(masses) - 1)))}), ncol=1)
      countyVarMat <- matrix(sapply(distributions, function(masses) {sum(masses * (0:(length(masses) - 1))^2) - sum(masses * (0:(length(masses) - 1)))^2}), ncol=1)
      # intervals = matrix(as.numeric(sapply(distributions, generateBinomialInterval)), nrow=4)
      intervals = matrix(as.numeric(sapply(distributions, getHDI)), nrow=4)
      
      # calculate crps
      calcCrpsByI = function(i) {
        crpsCounts(truthByCounty$truth[i] * (length(distributions[[i]]) - 1), distributions[[i]])
      }
      countyCrps = sapply(1:nrow(truthByCounty), calcCrpsByI)
      
      countyPredMatExact <- data.frame(preds=countyPredMatExact, vars=countyVarMat, 
                                  lower=intervals[1,], upper=intervals[2,], 
                                  leftRejectProb=intervals[3,], 
                                  rightRejectProb=intervals[4,], 
                                  crps=countyCrps)
    }
  }
  
  regionPredMat = NULL
  if(genRegionLevel) {
    ## generate region level predictions
    
    # if not supplied, get grid of population densities for pop-weighted integration
    if(is.null(popGrid))
      popGrid = makeInterpPopGrid(kmRes=kmRes)
    
    # make sure prediction coordinates correspond to population density grid coordinates
    if(length(pixelIndices) != nrow(popGrid)) {
      stop("for county level predictions, prediction points must be on same grid as popGrid")
    }
    
    # generate samples from posterior
    if(is.null(predMat)) {
      postSamples = inla.posterior.sample(nPostSamples, mod)
      latentMat = sapply(postSamples, function(x) x$latent)
      
      # get logit education rates at the prediction locations
      index.pred <- inla.stack.index(stack.full,tag="pred")$data
      predMat <- latentMat[index.pred,]
      
      if(genEALevel) {
        # make sure to separate enumeration area and pixel level predictions
        eaMat = expit(predMat[eaIndices, ])
        predMat = predMat[pixelIndices, ]
      }
    }
    
    # the below code is commented since we removed the cluster effect from the model. The commented code is not yet tested
    # # remove cluster effect for integration
    # nPreds = nrow(predCoords)
    # clusterIndices = inla.stack.index(stack.full,tag="pred")$effects[(1 + mesh$n):(mesh$n + nPreds)]
    # clusterEffects = latentMat[clusterIndices, ]
    # noClustMat = predMat - clusterEffects
    
    # for each region, integrate predictions with respect to population density
    if(!exactAggregation) {
      if(genCountLevel)
        stop("not using exact aggregation is not supported if genCountLevel is set to TRUE")
      
      predCounties = popGrid$admin1
      predRegions = countyToRegion(predCounties)
      pops = popGrid$popOrig
      integratePredsByRegion = function(regionName) {
        # subset data by region of interest
        regionI = predRegions == regionName
        regionPredMat = expit(predMat[regionI,])
        regionPops = pops[regionI]
        
        # compute average integral of predictions
        regionPops = regionPops/sum(regionPops)
        integrals = t(regionPredMat) %*% regionPops
      }
      regionPredMat <- t(sapply(regions, integratePredsByRegion))
    } else {
      # for each region, integrate predictions over EAs
      predCounties = eaDat$admin1
      predRegions = countyToRegion(predCounties)
      
      integratePredsByRegion = function(regionName) {
        # subset data by region of interest
        regionI = predRegions == regionName
        regionProbMat = matrix(eaMat[regionI,], nrow=sum(regionI))
        
        # combine results by EA
        # numClusters = sum(regionI)
        numChildren = eaDat$numChildren[regionI]
        distribution = dSumBinomRandom(0:sum(numChildren), numChildren, regionProbMat)
        distribution
      }
      
      # in this case, countyPredMat is a data frame with predictions and confidence intervals, with rejection probabilities
      distributions <- lapply(regions, integratePredsByRegion)
      regionPredMat <- matrix(sapply(distributions, function(masses) {sum(masses * (0:(length(masses) - 1)))}), ncol=1)
      regionVarMat <- matrix(sapply(distributions, function(masses) {sum(masses * (0:(length(masses) - 1))^2) - sum(masses * (0:(length(masses) - 1)))^2}), ncol=1)
      # intervals = matrix(as.numeric(sapply(distributions, generateBinomialInterval)), nrow=4)
      intervals = matrix(as.numeric(sapply(distributions, getHDI)), nrow=4)
      
      # calculate crps
      calcCrpsByI = function(i) {
        crpsCounts(truthByRegion$truth[i] * (length(distributions[[i]]) - 1), distributions[[i]])
      }
      regionCrps = sapply(1:nrow(truthByRegion), calcCrpsByI)
      
      regionPredMat <- data.frame(preds=regionPredMat, vars=regionVarMat, 
                                  lower=intervals[1,], upper=intervals[2,], 
                                  leftRejectProb=intervals[3,], 
                                  rightRejectProb=intervals[4,], 
                                  crps=regionCrps)
    }
  }
  
  # generate pixel level predictions if necessary
  pixelPreds = NULL
  if(keepPixelPreds) {
    if(!exactAggregation)
      pixelPreds = expit(predMat)
    else {
      # for each region, integrate predictions over EAs
      eaToPixel = eaDat$pixelI
      pixelsWithData = sort(unique(eaToPixel))
      thisEaMat = eaMat[, 1:nSamplePixel]
      
      integratePredsByRegion = function(pixelIndex) {
        # subset data by region of interest
        pixelI = eaToPixel == pixelIndex
        
        # if no data in the pixel, return NA
        if(sum(pixelI) == 0) {
          return(NA)
        } else if(sum(pixelI) == 1) {
          # otherwise, if the pixel only has one ea, use the marginal of the pixel
          # get the marginal "binomial" densities at each location
          
          # this function evaluates the "binomial" probability mass for a fixed n, k, and a marginal
          n = eaDat$numChildren[pixelI]
          binomProb = function(k) {
            inla.emarginal(function(logitP) {dbinom(k, n, expit(logitP))}, pixelMarginals[[pixelIndex]])
          }
          
          ## make highest density coverage interval on count scale
          # generate the "binomial" pmfs for this marginal
          probs = sapply(0:n, binomProb)
          probs
        } else {
          # if we have more than one EAs in the pixel, approximate the convolution of binomials with 
          # the pearson distribution
          
          # get matrix of EA simulated joint distribution education probabilities for this pixel. 
          # Only use a small amount of posterior samples for efficient computation
          pixelProbMat = matrix(thisEaMat[pixelI,], nrow=sum(pixelI))
          
          # combine results by pixel
          # numClusters = sum(pixelI)
          numChildren = eaDat$numChildren[pixelI]
          n = sum(numChildren)
          distribution = dSumBinomRandom(0:n, numChildren, pixelProbMat)
          distribution
        }
      }
      
      # in this case, pixelPredMat is a data frame with predictions and confidence intervals, with rejection probabilities
      distributions <- lapply(1:nrow(predMat), integratePredsByRegion)
      pixelPredMat <- matrix(sapply(distributions, function(masses) {sum(masses * (0:(length(masses) - 1)))}), ncol=1)
      pixelVarMat <- matrix(sapply(distributions, function(masses) {sum(masses * (0:(length(masses) - 1))^2) - sum(masses * (0:(length(masses) - 1)))^2}), ncol=1)
      # intervals = matrix(as.numeric(sapply(distributions, generateBinomialInterval)), nrow=4)
      intervals = matrix(as.numeric(sapply(distributions, getHDI)), nrow=4)
      
      # calculate crps
      calcCrpsByI = function(i) {
        crpsCounts(truthByPixel$truth[i] * (length(distributions[[i]]) - 1), distributions[[i]])
      }
      pixelCrps = sapply(1:nrow(truthByPixel), calcCrpsByI)
      
      pixelPredMat <- data.frame(preds=pixelPredMat[pixelsWithData], vars=pixelVarMat[pixelsWithData], 
                                 lower=intervals[1,pixelsWithData], upper=intervals[2,pixelsWithData], 
                                 leftRejectProb=intervals[3,pixelsWithData], 
                                 rightRejectProb=intervals[4,pixelsWithData], 
                                 crps=pixelCrps)
    }
  }
  
  list(mod=mod, preds=preds, SDs=predSDs, obsInds=obsInds, predInds=index, mesh=mesh, 
       prior=prior, stack=stack.full, countyPredMat=countyPredMat, regionPredMat=regionPredMat, 
       pixelPredMat=pixelPredMat, pixelMarginals=pixelMarginals, eaPredMat=eaMat, eaMarginals=eaMarginals)
}

# fit the SPDE model given some data and some prediction locations
# Inputs:
# obsCoords: coordinates of the observations
# obsNs: the binomial n of the observations
# obsCounts: the response, education counts
# obsUrban: whether the observations were urban or knot
# predCoords: the spatial coordinates of the predictions: first the EAs, then the pixels
# predNs: the binomial n of the predictions
# predUrban: whether the prediction locations are urban or not
# clusterIndices: the row indices in eaDat corresponding to the observations
# prior: the prior for the spde. Defaults to getSPDEPrior with the spatial mesh
# mesh: the spatial mesh. Defaults to something dat should be reasonable
# int.strategy: the introduction strategy argument for inla
# strategy: the strategy argument for inla
# genCountyLevel: whether or knot to generate county level predictions
# popGrid: population densities over the prediction grid for integration
# nPostSamples: number of samples from the joint posterior to take
# kmRes: the kilometers resolution of the pixels
# counties: the county names, determining the order of the county predictions
# clusterEffect: whether or not to include an iid gaussian error term at the EA level
# verbose: verbose argument to inla
# genRegionLevel: whether or not to generate regional level predictions. Not currently supported
# regions: the region names, determining the order of the region predictions
# keepPixelPreds: whether or not to return the pixel level predictions
# genEALevel: whether or not to generate ea level predictions
# eaIndices: indices giving which of the prediction coordinates in predCoords are EAs. NOTE: 
#            these must be the first indices of the predictions
# urbanEffect: if in urban fixed effect is included
# link: the link argument for inla. 1 means logit scale, 2 means probability scale (?)
# predictionType: either median or mean predictions are returned
# eaPixelIndices: the pixel indices each ea belongs to. This is of length equal to the number of EAs
# link: either one or two. If one, uses the identity link, otherwise uses
# exactAggregation: aggregate using the known ea locations rather than integrating over
#                   the aggregation area with respect to population density
# genCountLevel: whether or not to generate predictions at the count level versus logistic. DEPRECATED
# nSamplePixel: fewer samples are required for good approximation of the posterior at the pixel level 
#               than county, so we only take this many of the posterior samples for pixel level estimation
fitSPDEModel2 = function(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, 
                         predNs = rep(1, nrow(predCoords)), predUrban, clusterIndices, prior=NULL, 
                         mesh=NULL, int.strategy="eb", strategy="gaussian", 
                         genCountyLevel=FALSE, popGrid=NULL, nPostSamples=100, kmRes=5, 
                         counties=sort(unique(eaDat$admin1)), verbose=TRUE, genRegionLevel=FALSE, 
                         regions=sort(unique(eaDat$region)), keepPixelPreds=FALSE, genEALevel=FALSE, 
                         eaIndices=1:nrow(kenyaEAs), urbanEffect=TRUE, link=1, 
                         predictionType=c("median", "mean"), eaDat=NULL, nSamplePixel=10, 
                         truthByPixel=NULL, truthByCounty=NULL, truthByRegion=NULL, 
                         truthByEa=NULL, clusterEffect=FALSE, significance=.8, adjustPopSurface=FALSE) {
  
  # match the prediction type
  predictionType = match.arg(predictionType)
  
  # for enumeration area level predictions, set enumeration area end pixel level prediction indices
  eaMat = NULL
  if(genEALevel) {
    pixelIndices = (1:nrow(predCoords))[-eaIndices]
  }
  else {
    pixelIndices = 1:nrow(predCoords)
  }
  
  # make default mesh
  if(is.null(mesh)) {
    mesh = getSPDEMeshGrid(doPlot = FALSE)
  }
  
  # make default prior
  if(is.null(prior)) {
    prior = getSPDEPrior(mesh)
  }
  
  # make covariate matrices (intercept plus urban/rural depending on whether urban effect included)
  if(urbanEffect) {
    X = cbind(1, obsUrban)
    XPred = cbind(1, predUrban)
  }
  else {
    X = matrix(1, ncol=1, nrow=nrow(obsCoords))
    XPred = matrix(1, ncol=1, nrow=nrow(predCoords))
  }
  
  
  # construct A matrix for observations
  m = nrow(obsCoords)
  AEst = inla.spde.make.A(mesh, loc = obsCoords)
  
  # construct A matrix for predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  
  # make inla stack
  ys = obsCounts
  n = ncol(AEst) # number of basis elements
  nObs = length(ys) # number of observations
  nPreds = nrow(predCoords)
  latticeInds = 1:n
  
  # do not include cluster effect at the pixel level
  stack.pred = inla.stack(A =list(APred[pixelIndices,], 1),
                          effects =list(field=latticeInds, X=XPred[pixelIndices,]),
                          data =list(y=NA, Ntrials = predNs[pixelIndices], link=1),
                          tag ="pred",
                          remove.unused=FALSE)
  
  if(clusterEffect) {
    # clustIObs = (n+1):(n+nObs)
    # clustIPreds = (n+1):(n+nPreds)
    # clustIPreds = (n+nObs+1):(n+nObs+nPreds)
    clustIObs = clusterIndices
  }
  
  if(genEALevel) {
    
    if(clusterEffect) {
      
      # only include clusters in the enumeration area predictions
      clustIPreds = 1:nrow(eaDat)
      stack.predEA = inla.stack(A =list(APred[eaIndices,], 1, 1),
                                effects =list(field=latticeInds, clust=clustIPreds, X=XPred[eaIndices,]),
                                data =list(y=NA, Ntrials = predNs[eaIndices], link=1),
                                tag ="pred",
                                remove.unused=FALSE)
    } else {
      stack.predEA = inla.stack(A =list(APred[eaIndices,], 1),
                                effects =list(field=latticeInds, X=XPred[eaIndices,]),
                                data =list(y=NA, Ntrials = predNs[eaIndices], link=1),
                                tag ="pred",
                                remove.unused=FALSE)
    }
  }
  
  # construct the observation stack 
  if(clusterEffect) {
    stack.est = inla.stack(A =list(AEst, 1, 1),
                           effects =list(field=latticeInds, clust=clustIObs, X=X),
                           data =list(y=ys, Ntrials = obsNs, link=1),
                           tag ="est",
                           remove.unused=FALSE)
    
  } else {
    stack.est = inla.stack(A =list(AEst, 1), 
                           effects =list(field=latticeInds, X=X), 
                           data =list(y=ys, Ntrials = obsNs, link=link), 
                           tag ="est", 
                           remove.unused=FALSE)
  }
  
  # # whether or not cluster effect is included in the model, it should not be included in the pixel level predictions 
  # stack.pred = inla.stack(A =list(APred, 1),
  #                         effects =list(field=latticeInds, X=XPred),
  #                         data =list(y=NA, Ntrials = predNs, link=link),
  #                         tag ="pred",
  #                         remove.unused=FALSE)
  
  # construct the full prediction stack if necessary
  if(genEALevel) {
    stack.pred = inla.stack(stack.predEA, stack.pred)
  }
  
  # make mesh index
  mesh.index <- inla.spde.make.index(name = "field", n.spde = prior$n.spde)
  
  # fit model
  control.inla = list(strategy=strategy, int.strategy=int.strategy) 
  allNs = c(obsNs, predNs)
  stack.full = inla.stack(stack.est, stack.pred, remove.unused=FALSE)
  stackDat = inla.stack.data(stack.full, spde=prior)
  allQuantiles = c((1-significance) / 2, 1 - (1-significance) / 2, seq(0, 1 - 1 / nSamplePixel, by = 1/nSamplePixel) + 1 / (2 * nSamplePixel))
  if(clusterEffect) {
    mod = inla(y ~ - 1 + X + f(field, model=prior) + 
                 f(clust, model="iid", 
                   hyper = list(prec = list(prior = "pc.prec", param = c(3,0.01)))),
               data = stackDat, Ntrials=stackDat$Ntrials,
               control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles),
               family="binomial", verbose=TRUE, control.inla=control.inla,
               control.compute=list(config=TRUE))
  } else {
    mod = inla(y ~ - 1 + X + f(field, model=prior), 
               data = stackDat, Ntrials=stackDat$Ntrials, 
               control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
               family="binomial", verbose=verbose, control.inla=control.inla, 
               control.compute=list(config=TRUE))
  }
  
  # get predictive surface, SD, and data
  n = nrow(obsCoords)
  # obsInds = 1:n
  obsInds = inla.stack.index(stack.full, "est")$data
  predInds = inla.stack.index(stack.full, "pred")$data
  index = inla.stack.index(stack.full, "pred")$data
  if(predictionType == "mean")
    linpreds = mod[["summary.fitted.values"]]$mean
  else
    linpreds = mod[["summary.fitted.values"]]$`0.5quant`
  linpred.sd = mod[["summary.fitted.values"]]$sd
  
  # predInds = (n+1):(n+nrow(predCoords))
  preds = linpreds
  predSDs = linpred.sd
  
  # if not supplied, get grid of population densities for pop-weighted integration
  if(is.null(popGrid))
    popGrid = makeInterpPopGrid(kmRes=kmRes)
  
  # make sure prediction coordinates correspond to population density grid coordinates
  if(length(pixelIndices) != nrow(popGrid)) {
    stop("for county level predictions, prediction points must be on same grid as popGrid")
  }
  
  # generate samples from posterior
  postSamples = inla.posterior.sample(nPostSamples, mod)
  latentMat = sapply(postSamples, function(x) x$latent)
  
  # get logit education rates at the prediction locations
  index.pred <- inla.stack.index(stack.full,tag="pred")$data
  predMat <- latentMat[index.pred,]
  
  eaMarginals = NULL
  pixelMarginals = NULL
  if(genEALevel) {
    # make sure to separate enumeration area and pixel level predictions
    eaMat = expit(predMat[eaIndices, ])
    predMat = predMat[pixelIndices, ]
    
    # get the marginals of the enumeration areas
    eaMarginals = mod$marginals.linear.predictor[index.pred[eaIndices]]
  }
  if(keepPixelPreds)
    pixelMarginals = mod$marginals.linear.predictor[index.pred[pixelIndices]]
  
  # The below code shouldn't be necessary since the matrix should already be produced:
  # # generate samples from posterior to use for prediction and aggregation
  # if(is.null(predMat)) {
  #   postSamples = inla.posterior.sample(nPostSamples, mod)
  #   latentMat = sapply(postSamples, function(x) x$latent)
  #   
  #   # get logit education rates at the prediction locations
  #   index.pred <- inla.stack.index(stack.full,tag="pred")$data
  #   predMat <- latentMat[index.pred,]
  #   
  #   if(genEALevel) {
  #     # make sure to separate enumeration area and pixel level predictions
  #     eaMat = expit(predMat[eaIndices, ])
  #     predMat = predMat[pixelIndices, ]
  #   }
  # }
  
  pops = popGrid$popOrig
  
  # regionPixelI: the pixel indices (from 1 to the number of pixels) over which to aggregate predictions
  inexactIntegration = function(regionPixelI, nSample=nPostSamples) {
    ## aggregate predictions over the region by integrating with respect to population density
    
    # subset pixel predictions by region of interest
    regionPredMat = expit(predMat[regionPixelI,1:nSample])
    regionPops = pops[regionPixelI]
    
    # compute integrals of predictions with respect to population density
    regionPops = regionPops/sum(regionPops)
    integrals = t(regionPredMat) %*% regionPops
    
    # the result is a column, where each element is a aggregated probability from a different simulation
    integrals
  }
  
  # regionEAI: the indices of enumeration areas (from 1 to ~96k) within the region to aggregate over
  # withBinom: whether or not to include binomial variation
  exactIntegration = function(regionEAI, withBinom, nSample=nPostSamples) {
    ## aggregate predictions over the region by integrating with respect to number of children per EA
    
    # subset data by county of interest
    regionProbMat = matrix(eaMat[regionEAI,], nrow=sum(regionEAI))
    
    # combine results by EA
    numChildren = eaDat$numChildren[regionEAI]
    if(!withBinom) {
      # same as the inexact aggregation, but integrate with respect to number of children per EA
      numChildren = numChildren/sum(numChildren)
      integrals = t(regionProbMat[,1:nSample]) %*% numChildren
      return(integrals)
    } else {
      # Use Pearson approximation to account for number of children per EA while accounting for binomial variation. 
      # Return the resulting discrete probability mass function
      distribution = dSumBinomRandom(0:sum(numChildren), numChildren, regionProbMat[,1:nSample])
      return(distribution)
    }
  }
  
  countyPreds = NULL
  if(genCountyLevel) {
    ### county level predictions
    print("Aggregating over counties")
    
    predCountiesPixel = popGrid$admin1
    predCountiesEa = eaDat$admin1
    integratePredsByCounty = function(countyName, exact, withBinom=NULL) {
      if(exact) {
        countyI = predCountiesEa == countyName
        exactIntegration(countyI, withBinom)
      }
      else {
        countyI = predCountiesPixel == countyName
        inexactIntegration(countyI)
      }
    }
    
    # integrate over the pixels to get aggregated predictions
    countyPredMatInexact <- t(sapply(counties, integratePredsByCounty, exact=FALSE))
    cat(".")
    
    # integrate over the EAs to get aggregated predictions
    if(genEALevel) {
      # the following are respectively a matrix of probability draws from the county 
      # distributions and a list of county probability distributions
      if(is.null(eaDat))
        stop("eaDat is null, but a non-null eaDat must be provided for account level predictions for EAs.")
      countyPredMatExact = t(sapply(counties, integratePredsByCounty, exact=TRUE, withBinom=FALSE))
      cat(".")
      countyPredMatExactB <- lapply(counties, integratePredsByCounty, exact=TRUE, withBinom=TRUE)
      cat(".")
    } else {
      countyPredMatExact = NULL
      countyPredMatExactB = NULL
    }
    
    countyPreds <- list(countyPredMatInexact=countyPredMatInexact, countyPredMatExact=countyPredMatExact, 
                        countyPredMatExactB=countyPredMatExactB)
  }
  
  regionPreds = NULL
  if(genRegionLevel) {
    print("Aggregating over regions")
    
    ## generate region level predictions
    predCountiesPixel = popGrid$admin1
    predCountiesEa = eaDat$admin1
    predRegionsPixel = countyToRegion(predCountiesPixel)
    predRegionsEa = countyToRegion(predCountiesEa)
    integratePredsByRegion = function(regionName, exact, withBinom=NULL) {
      if(exact) {
        regionI = predRegionsEa == regionName
        exactIntegration(regionI, withBinom)
      }
      else {
        regionI = predRegionsPixel == regionName
        inexactIntegration(regionI)
      }
    }
    
    # integrate over the pixels to get aggregated predictions
    regionPredMatInexact <- t(sapply(regions, integratePredsByRegion, exact=FALSE))
    cat(".")
    
    # integrate over the EAs to get aggregated predictions
    if(genEALevel) {
      # the following are respectively a matrix of probability draws from the county 
      # distributions and a list of county probability distributions
      if(is.null(eaDat))
        stop("eaDat is null, but a non-null eaDat must be provided for account level predictions for EAs.")
      regionPredMatExact = t(sapply(regions, integratePredsByRegion, exact=TRUE, withBinom=FALSE))
      cat(".")
      regionPredMatExactB <- lapply(regions, integratePredsByRegion, exact=TRUE, withBinom=TRUE)
      cat(".")
    } else {
      regionPredMatExact = NULL
      regionPredMatExactB = NULL
    }
    
    regionPreds <- list(regionPredMatInexact=regionPredMatInexact, regionPredMatExact=regionPredMatExact, 
                        regionPredMatExactB=regionPredMatExactB)
  }
  
  # generate pixel level predictions if necessary
  eaToPixel = eaDat$pixelI
  # pixelsWithData = sort(unique(eaToPixel))
  pixelsWithData = unique(eaToPixel) # don't sort so that these indices matchup with the indices in truth
  pixelPredIndices = index.pred[pixelIndices[pixelsWithData]]
  pixelPreds = NULL
  if(keepPixelPreds) {
    print("Aggregating over pixels")
    
    ##### redefine the integration functions, since care must be taken for pixel level predictions 
    ##### due to the fineness of the spatial resolution
    # pixelI: the pixel indices (from 1 to the number of pixels) over which to aggregate predictions
    inexactIntegrationPixel = function(pixelI) {
      ## aggregate predictions over the pixel by taking the quantiles of the marginal for the pixel
      
      thisMarginal = pixelMarginals[[pixelI]]
      quantiles = seq(0, 1 - 1 / nSamplePixel, by = 1/nSamplePixel) + 1 / (2 * nSamplePixel)
      inla.qmarginal(quantiles, thisMarginal)
    }
    
    # pixelI: the pixel indices (from 1 to the number of pixels) over which to aggregate predictions
    # pixelEAI: logical vector of enumeration areas (from 1 to ~96k) within the region to aggregate over
    # withBinom: whether or not to include binomial variation
    exactIntegrationPixel = function(pixelI, pixelEAI, withBinom) {
      ## aggregate predictions over the region by integrating with respect to number of children per EA
      nEAs = sum(pixelEAI)
      numChildren = eaDat$numChildren[pixelEAI]
      
      # subset data by county of interest
      regionProbMat = matrix(eaMat[pixelEAI,1:nSamplePixel], nrow=length(numChildren))
      
      # combine results by EA
      if(!withBinom) {
        if(nEAs == 1) {
          # just use the quantiles of the EA marginal
          expit(matrix(unlist(mod$summary.linear.predictor[pixelIndices[pixelI], 5:(4 + nSamplePixel)]), nrow=1))
        } else {
          # take the weighted average of joint probability draws of EAs in this pixel
          # apply(thisEaMat[pixelEAI], 2, weighted.mean, w=numChildren)
          colSums(regionProbMat * (numChildren / sum(numChildren)))
        }
      } else {
        # Use Pearson approximation to account for number of children per EA while accounting for binomial variation. 
        # Return the resulting discrete probability mass function. Aside from that, separate into the
        # cases from without binomial variation: 
        
        if(nEAs == 1) {
          # just use the quantiles of the EA marginal
          if(length(numChildren) > 1)
            stop("Too many children")
          thisProbMat = expit(matrix(unlist(mod$summary.linear.predictor[eaIndices[pixelEAI], 5:(4 + nSamplePixel)]), nrow=1))
          dSumBinomRandom(0:sum(numChildren), numChildren, thisProbMat)
        } else {
          # use the joint probability draws of EAs in this pixel
          dSumBinomRandom(0:sum(numChildren), numChildren, regionProbMat)
        }
      }
    }
    
    ## generate pixel level predictions
    # for each region, integrate predictions over EAs
    pixelToEa = match(pixelsWithData, eaToPixel)
    predPixelPixel = 1:nrow(popGrid)
    predPixelEa = pixelToEa
    integratePredsByPixel = function(pixelID, withBinom=NULL) {
      # if(exact) {
      pixelI = pixelID
      pixelEaI = eaToPixel == pixelID
      exactIntegrationPixel(pixelI, pixelEaI, withBinom)
      # }
      # else {
      #   pixelI = predRegionsPixel == pixelID
      #   inexactIntegrationPixel(pixelI)
      # }
    }
    
    # use the precomputed logit scale mean, variance, CIs, and equally spaced pixel marginal 
    # quantiles for integration
    pixelPredIndices = index.pred[pixelIndices[pixelsWithData]]
    logitLower = mod$summary.linear.predictor[pixelPredIndices, 3]
    logitUpper = mod$summary.linear.predictor[pixelPredIndices, 4]
    logitEst = mod$summary.linear.predictor[pixelPredIndices, 1]
    logitVar = mod$summary.linear.predictor[pixelPredIndices, 2]
    pixelPredMatInexact = expit(as.matrix(unlist(mod$summary.linear.predictor[pixelPredIndices, 5:(4 + nSamplePixel)])))
    est = rowMeans(pixelPredMatInexact)
    cat(".")
    
    # integrate over the EAs to get aggregated predictions
    if(genEALevel) {
      # the following are respectively a matrix of probability draws from the county 
      # distributions and a list of county probability distributions
      if(is.null(eaDat))
        stop("eaDat is null, but a non-null eaDat must be provided for account level predictions for EAs.")
      pixelPredMatExact = t(sapply(pixelsWithData, integratePredsByPixel, withBinom=FALSE))
      cat(".")
      pixelPredMatExactB <- lapply(pixelsWithData, integratePredsByPixel, withBinom=TRUE)
      cat(".")
    } else {
      pixelPredMatExact = NULL
      pixelPredMatExactB = NULL
    }
    
    pixelPreds <- list(pixelPredMatInexact=pixelPredMatInexact, pixelPredMatExact=pixelPredMatExact, 
                       pixelPredMatExactB=pixelPredMatExactB, pixelsWithData=pixelsWithData, 
                       logitLower=logitLower, logitUpper=logitUpper, logitEst=logitEst, logitVar=logitVar, 
                       est=est)
  }
  
  if(genEALevel) {
    # use the precomputed logit scale mean, variance, CIs, and equally spaced pixel marginal 
    # quantiles for integration
    eaPredIndices = index.pred[eaIndices]
    logitLower = mod$summary.linear.predictor[eaPredIndices, 3]
    logitUpper = mod$summary.linear.predictor[eaPredIndices, 4]
    logitEst = mod$summary.linear.predictor[eaPredIndices, 1]
    logitVar = mod$summary.linear.predictor[eaPredIndices, 2]
    eaPredMat = expit(unlist(as.matrix(mod$summary.linear.predictor[eaPredIndices, 5:(4 + nSamplePixel)])))
    est = rowMeans(eaPredMat)
    
    eaPreds = list(eaPredMat=eaPredMat, logitLower=logitLower, logitUpper=logitUpper, 
                   logitEst=logitEst, logitVar=logitVar, est=est)
  } else {
    eaPreds = NULL
  }
  
  list(mod=mod, preds=preds, SDs=predSDs, obsInds=obsInds, predInds=index, pixelInds=pixelIndices, 
       eaInds=eaIndices, mesh=mesh, prior=prior, stack=stack.full, 
       countyPreds=countyPreds, regionPreds=regionPreds, pixelPreds=pixelPreds, eaPreds=eaPreds)
}

# fit the SPDE model given some data and some prediction locations
# Inputs:
# obsCoords: coordinates of the observations
# obsNs: the binomial n of the observations
# obsCounts: the response, secondary education completion counts
# obsUrban: whether the observations were urban or knot
# predCoords: the spatial coordinates of the predictions: first the EAs, then the pixels
# predNs: the binomial n of the predictions
# predUrban: whether the prediction locations are urban or not
# clusterIndices: the row indices in eaDat corresponding to the observations
# prior: the prior for the spde. Defaults to getSPDEPrior with the spatial mesh
# mesh: the spatial mesh. Defaults to something dat should be reasonable
# int.strategy: the introduction strategy argument for inla
# strategy: the strategy argument for inla
# genCountyLevel: whether or knot to generate county level predictions
# popGrid: population densities over the prediction grid for integration
# nPostSamples: number of samples from the joint posterior to take
# kmRes: the kilometers resolution of the pixels
# counties: the county names, determining the order of the county predictions
# clusterEffect: whether or not to include an iid gaussian error term at the EA level
# verbose: verbose argument to inla
# genRegionLevel: whether or not to generate regional level predictions. Not currently supported
# regions: the region names, determining the order of the region predictions
# keepPixelPreds: whether or not to return the pixel level predictions
# genEALevel: whether or not to generate ea level predictions
# eaIndices: indices giving which of the prediction coordinates in predCoords are EAs. NOTE: 
#            these must be the first indices of the predictions
# urbanEffect: if in urban fixed effect is included
# link: the link argument for inla. 1 means logit scale, 2 means probability scale (?)
# predictionType: either median or mean predictions are returned
# eaPixelIndices: the pixel indices each ea belongs to. This is of length equal to the number of EAs
# link: either one or two. If one, uses the identity link, otherwise uses
# exactAggregation: aggregate using the known ea locations rather than integrating over
#                   the aggregation area with respect to population density
# genCountLevel: whether or not to generate predictions at the count level versus logistic. DEPRECATED
# nSamplePixel: fewer samples are required for good approximation of the posterior at the pixel level 
#               than county, so we only take this many of the posterior samples for pixel level estimation
# allPixels: if TRUE, predict at all pixels. If FALSE, predict only at pixels with data
# newMesh: this should definitely be set to TRUE, or the mesh will probably not be fine enough
# doValidation: whether or not to perform validation
# previousResult: the output of the inla call of a previous call to this function. Used for initialization
fitSPDEModel3 = function(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, 
                         predNs = rep(1, nrow(predCoords)), predUrban, clusterIndices, prior=NULL, 
                         mesh=NULL, int.strategy="grid", strategy="laplace", 
                         genCountyLevel=FALSE, popGrid=NULL, nPostSamples=100, kmRes=5, 
                         counties=sort(unique(eaDat$admin1)), verbose=TRUE, genRegionLevel=FALSE, 
                         regions=sort(unique(eaDat$region)), keepPixelPreds=FALSE, genEALevel=FALSE, 
                         eaIndices=1:nrow(kenyaEAs), urbanEffect=TRUE, link=1, 
                         predictionType=c("mean", "median"), eaDat=NULL, nSamplePixel=10, 
                         clusterEffect=FALSE, significance=.8, 
                         onlyInexact=FALSE, allPixels=FALSE, newMesh=TRUE, doValidation=FALSE, 
                         previousResult=NULL, predCountyI=NULL, continuousOnly=FALSE, strictPrior=FALSE, 
                         integrateOutCluster=TRUE, returnUnintegratedResults=TRUE, adjustPopSurface=TRUE, 
                         seed=NULL, popGridAdjusted=NULL, targetPop=c("children", "women"), 
                         strictSpatialPrior=FALSE) {
  if(!is.null(seed))
    set.seed(seed)
  
  if(!is.null(predCountyI) && !onlyInexact)
    stop("If generating predictions for a fixed county (i.e. predCountyI is not NULL) then onlyInexact currently must be set to TRUE")
  
  # match the prediction type
  predictionType = match.arg(predictionType)
  
  # for enumeration area level predictions, set enumeration area end pixel level prediction indices
  eaMat = NULL
  if(genEALevel) {
    pixelIndices = (1:nrow(predCoords))[-eaIndices]
  }
  else {
    pixelIndices = 1:nrow(predCoords)
  }
  
  # make default mesh
  if(is.null(mesh)) {
    if(newMesh)
      mesh = getSPDEMeshGrid2(doPlot = FALSE)
    else
      mesh = getSPDEMeshGrid(doPlot = FALSE)
  }
  
  # make default prior
  if(is.null(prior)) {
    prior = getSPDEPrior(mesh, strictPrior=strictPrior||strictSpatialPrior)
  }
  
  # make covariate matrices (intercept plus urban/rural depending on whether urban effect included)
  if(urbanEffect) {
    X = cbind(1, obsUrban)
    XPred = cbind(1, predUrban)
  }
  else {
    X = matrix(1, ncol=1, nrow=nrow(obsCoords))
    XPred = matrix(1, ncol=1, nrow=nrow(predCoords))
  }
  
  
  # construct A matrix for observations
  m = nrow(obsCoords)
  AEst = inla.spde.make.A(mesh, loc = obsCoords)
  
  # construct A matrix for predictions
  APred = inla.spde.make.A(mesh, loc = predCoords)
  
  # make inla stack
  ys = obsCounts
  n = ncol(AEst) # number of basis elements
  nObs = length(ys) # number of observations
  nPreds = nrow(predCoords)
  latticeInds = 1:n
  
  if(clusterEffect) {
    # clustIObs = (n+1):(n+nObs)
    # clustIPreds = (n+1):(n+nPreds)
    # clustIPreds = (n+nObs+1):(n+nObs+nPreds)
    clustIObs = clusterIndices
  }
    
  # construct the observation stack 
  if(clusterEffect) {
    stack.est = inla.stack(A =list(AEst, 1, 1),
                           effects =list(field=latticeInds, clust=clustIObs, X=X),
                           data =list(y=ys, Ntrials = obsNs, link=1),
                           tag ="est",
                           remove.unused=FALSE)
    
  } else {
    stack.est = inla.stack(A =list(AEst, 1), 
                           effects =list(field=latticeInds, X=X), 
                           data =list(y=ys, Ntrials = obsNs, link=link), 
                           tag ="est", 
                           remove.unused=FALSE)
  }
  
  # make mesh index
  mesh.index <- inla.spde.make.index(name = "field", n.spde = prior$n.spde)
  
  # fit model
  if(doValidation) {
    # as recommended by http://www.r-inla.org/faq#TOC-How-can-I-compute-cross-validation-or-predictive-measures-of-fit-
    if(strategy != "laplace" || int.strategy != "grid")
      warning("Since doValidation has been set to TRUE, changing integration strategy to 'grid', and strategy to 'laplace'")
    control.inla = list(strategy="laplace", int.strategy="grid", diff.logdens=4, npoints=21) 
  }
  else {
    control.inla = list(strategy=strategy, int.strategy=int.strategy)
  }
  modeControl = inla.set.control.mode.default()
  if(!is.null(previousResult)) {
    # initialize the fitting process based on a previous optimum
    # modeControl$result = previousResult
    modeControl$theta = previousResult$mode$theta
    modeControl$x = previousResult$mode$x
    modeControl$restart = TRUE
  }
  
  allNs = obsNs
  stack.full = stack.est
  stackDat = inla.stack.data(stack.full, spde=prior)
  allQuantiles = c(0.5, (1-significance) / 2, 1 - (1-significance) / 2)
  if(!strictPrior)
    clusterList = list(param=c(1, 0.01), prior="pc.prec")
  else
    clusterList = list(param=c(.15, 0.01), prior="pc.prec")
  if(clusterEffect) {
    mod = inla(y ~ - 1 + X + f(field, model=prior) + 
                 f(clust, model="iid", hyper = list(prec = clusterList)),
               data = stackDat, Ntrials=stackDat$Ntrials,
               control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles),
               family="binomial", verbose=verbose, control.inla=control.inla,
               control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
               control.mode=modeControl)
  } else {
    mod = inla(y ~ - 1 + X + f(field, model=prior), 
               data = stackDat, Ntrials=stackDat$Ntrials, 
               control.predictor=list(A=inla.stack.A(stack.full), compute=TRUE, link=stackDat$link, quantiles=allQuantiles), 
               family="binomial", verbose=verbose, control.inla=control.inla, 
               control.compute=list(config=TRUE, cpo=doValidation, dic=doValidation, waic=doValidation), 
               control.mode=modeControl)
  }
  print("Finish fitting model")
  
  # get predictive surface, SD, and data
  n = nrow(obsCoords)
  # obsInds = 1:n
  obsInds = inla.stack.index(stack.full, "est")$data
  predInds = inla.stack.index(stack.full, "pred")$data
  index = inla.stack.index(stack.full, "pred")$data
  if(predictionType == "mean")
    linpreds = mod[["summary.fitted.values"]]$mean
  else
    linpreds = mod[["summary.fitted.values"]]$`0.5quant`
  linpred.sd = mod[["summary.fitted.values"]]$sd
  
  # predInds = (n+1):(n+nrow(predCoords))
  preds = linpreds
  predSDs = linpred.sd
  
  # if not supplied, get grid of population densities for pop-weighted integration
  if(kmRes == 5) {
    if(is.null(popGridAdjusted) && adjustPopSurface) {
      if(targetPop == "children") {
        load("popGridAdjusted.RData")
        popGridAdjusted = popGrid
      }
      else {
        load("popGridAdjustedWomen.RData")
        popGridAdjusted = popGrid
      }
    }
    if(is.null(popGrid)) {
      load("popGrid.RData")
    }
  }
  else {
    popGrid = makeInterpPopGrid(kmRes, adjustPopSurface, targetPop)
  }
  
  # make sure prediction coordinates correspond to population density grid coordinates
  if(length(pixelIndices) != nrow(popGrid)) {
    stop("for county level predictions, prediction points must be on same grid as popGrid")
  }
  
  # generate samples from posterior
  ("Sampling from the posterior...")
  postSamples = inla.posterior.sample(nPostSamples, mod)
  latentMat = sapply(postSamples, function(x) {x$latent})
  if(clusterEffect)
    clusterVars = sapply(postSamples, function(x) {1 / x$hyperpar[3]})
  latentVarNames = rownames(postSamples[[1]]$latent)
  fieldIndices = which(grepl("field", latentVarNames))
  fixedIndices = which(grepl("X", latentVarNames))
  if(clusterEffect)
    clustIndices = grepl("clust", latentVarNames)
  
  # generate logit predictions (first without cluster effect then add the cluster effect in)
  paste("Computing prediction surface and EA predictions for integration...")
  predMatNoClust = XPred  %*% latentMat[fixedIndices,] + APred %*% latentMat[fieldIndices,]
  predMat <- predMatNoClust
  
  if(clusterEffect) {
    # add in estimated cluster effects at the sampled enumeration areas
    predMat[clusterIndices,] = predMat[clusterIndices,] + latentMat[clustIndices,]
    
    # addend cluster effects we have no information about at the on sampled enumeration areas
    predMat[eaIndices[-clusterIndices],] = predMat[eaIndices[-clusterIndices],] + 
      rnorm(length(eaIndices[-clusterIndices]) * nPostSamples, sd = rep(sqrt(clusterVars), each=length(eaIndices[-clusterIndices])))
  }
  
  if(genEALevel) {
    # make sure to separate enumeration area and pixel level predictions
    eaMat = expit(predMat[eaIndices, ])
    predMat = predMat[pixelIndices, ]
  }
  
  unintegratedMat = predMat
  if(integrateOutCluster && clusterEffect) {
    # In this case, we shift each pixel level prediction by the the amount that will remove the 
    # bias induced by not accounting for the cluster effect (we integrate out the cluster effect 
    # for each simulated pixel value)
    predMat = matrix(logit(logitNormMean(cbind(c(as.matrix(predMat)), rep(sqrt(clusterVars), each=nrow(predMat))))), nrow=nrow(predMat))
  }
  
  pops = popGrid$popOrig
  
  # regionPixelI: the pixel indices (from 1 to the number of pixels) over which to aggregate predictions
  # integrateByPixel: either integrate over the pixels or enumeration areas
  noBinomialIntegration = function(integrationMat, integrateByPixel=TRUE, nSample=nPostSamples, 
                                   pixelMatrix=predMat, eaMatrix=eaMat) {
    ## aggregate predictions over the region by integrating with respect to population density
    
    # either integrate predictions of the pixels or enumeration areas
    if(integrateByPixel)
      regionPredMat = expit(pixelMatrix[,1:nSample])
    else
      regionPredMat = eaMatrix[,1:nSample]
    
    # compute integrals of predictions with respect to population density or the number of children
    integrationMat %*% regionPredMat
  }
  
  # regionEAI: the indices of enumeration areas (from 1 to ~96k) within the region to aggregate over
  binomialIntegration = function(regionEAI, nSample=nPostSamples) {
    ## aggregate predictions over the region by integrating with respect to number of children per EA
    
    # subset data by county of interest
    regionProbMat = matrix(eaMat[regionEAI,1:nSample], nrow=sum(regionEAI))
    
    # combine results by EA
    numChildren = eaDat$numChildren[regionEAI]
    
    # Use Pearson approximation to account for number of children per EA while accounting for binomial variation. 
    # Return the resulting discrete probability mass function
    distribution = dSumBinomRandom2(0:sum(numChildren), numChildren, regionProbMat)
    distribution
  }
  
  countyPreds = NULL
  if(genCountyLevel) {
    ### county level predictions
    print("Aggregating over counties")
    
    predCountiesPixel = popGrid$admin1
    
    getCountyIntegrationMatrix = function(integrateByPixel=TRUE, adjustDensities=FALSE) {
      counties = as.character(counties)
      
      if(integrateByPixel) {
        mat = t(sapply(counties, function(countyName) {popGrid$admin1 == countyName}))
        if(!adjustDensities)
          mat = sweep(mat, 2, popGrid$popOrig, "*")
        else
          mat = sweep(mat, 2, popGridAdjusted$popOrig, "*")
      }
      else {
        mat = t(sapply(counties, function(countyName) {eaDat$admin1 == countyName}))
        mat = sweep(mat, 2, eaDat$numChildren, "*")
      }
      sweep(mat, 1, 1/rowSums(mat), "*")
    }
    
    # integrate over the pixels to get aggregated predictions. We consider three cases from most to least complex: 
    # integrate out cluster and adjust population density surface to account for different number of children in urban/rural areas
    # adjust population density surface to account for different number of children in urban/rural areas (Unintegrated)
    # do no integration or adjustments (Unadjusted)
    countyPredMatInexact <- noBinomialIntegration(getCountyIntegrationMatrix(TRUE, adjustPopSurface), TRUE, pixelMatrix=predMat)
    if(integrateOutCluster && returnUnintegratedResults) {
      countyPredMatInexactUnintegrated <- noBinomialIntegration(getCountyIntegrationMatrix(TRUE, adjustPopSurface), TRUE, pixelMatrix=unintegratedMat)
      
      if(adjustPopSurface)
        countyPredMatInexactUnadjusted <- noBinomialIntegration(getCountyIntegrationMatrix(TRUE, FALSE), TRUE, pixelMatrix=unintegratedMat)
      else
        countyPredMatInexactUnadjusted = NULL
    }
    else {
      countyPredMatInexactUnintegrated = NULL
    }
    cat(".")
    
    # integrate over the EAs to get aggregated predictions
    if(genEALevel && !onlyInexact && !continuousOnly) {
      predCountiesEa = eaDat$admin1
      
      countyBinomialIntegration = function(countyName) {
        countyI = predCountiesEa == countyName
        binomialIntegration(countyI)
      }
      
      # the following are respectively a matrix of probability draws from the county 
      # distributions and a list of county probability distributions
      if(is.null(eaDat))
        stop("eaDat is null, but a non-null eaDat must be provided for account level predictions for EAs.")
      countyPredMatExact = noBinomialIntegration(getCountyIntegrationMatrix(FALSE), FALSE)
      cat(".")
      countyPredMatExactB <- lapply(counties, countyBinomialIntegration)
      cat(".")
    } else {
      countyPredMatExact = NULL
      countyPredMatExactB = NULL
    }
    
    countyPreds <- list(countyPredMatInexact=countyPredMatInexact, countyPredMatExact=countyPredMatExact, 
                        countyPredMatExactB=countyPredMatExactB, countyPredMatInexactUnintegrated=countyPredMatInexactUnintegrated, 
                        countyPredMatInexactUnadjusted=countyPredMatInexactUnadjusted)
  }
  
  regionPreds = NULL
  if(genRegionLevel) {
    print("Aggregating over regions")
    
    ## generate region level predictions
    predCountiesPixel = popGrid$admin1
    predRegionsPixel = countyToRegion(predCountiesPixel)
    
    
    getRegionIntegrationMatrix = function(integrateByPixel=TRUE) {
      
      if(integrateByPixel) {
        mat = t(sapply(regions, function(regionName) {predRegionsPixel == regionName}))
        mat = sweep(mat, 2, popGrid$popOrig, "*")
      }
      else {
        mat = t(sapply(regions, function(regionName) {predRegionsEa == regionName}))
        mat = sweep(mat, 2, eaDat$numChildren, "*")
      }
      sweep(mat, 1, 1/rowSums(mat), "*")
    }
    
    regionBinomialIntegration = function(regionName) {
      regionI = predRegionsEa == regionName
      binomialIntegration(regionI)
    }
    
    # integrate over the pixels to get aggregated predictions
    regionPredMatInexact <- noBinomialIntegration(getRegionIntegrationMatrix(TRUE), TRUE)
    cat(".")
    
    # integrate over the EAs to get aggregated predictions
    if(genEALevel && !onlyInexact && !continuousOnly) {
      
      predCountiesEa = eaDat$admin1
      predRegionsEa = countyToRegion(predCountiesEa)
      
      # the following are respectively a matrix of probability draws from the county 
      # distributions and a list of county probability distributions
      if(is.null(eaDat))
        stop("eaDat is null, but a non-null eaDat must be provided for account level predictions for EAs.")
      regionPredMatExact = noBinomialIntegration(getRegionIntegrationMatrix(FALSE), FALSE)
      cat(".")
      regionPredMatExactB <- lapply(regions, regionBinomialIntegration)
      cat(".")
    } else {
      regionPredMatExact = NULL
      regionPredMatExactB = NULL
    }
    
    regionPreds <- list(regionPredMatInexact=regionPredMatInexact, regionPredMatExact=regionPredMatExact, 
                        regionPredMatExactB=regionPredMatExactB)
  }
  
  # generate pixel level predictions if necessary
  # pixelsWithData = sort(unique(eaToPixel))
  if(!onlyInexact || continuousOnly) {
    eaToPixel = eaDat$pixelI
    pixelsWithData = unique(eaToPixel) # don't sort so that these indices matchup with the indices in truth
  }
  else
    pixelsWithData = NULL
  pixelPreds = NULL
  if(keepPixelPreds) {
    print("Aggregating over pixels")
    
    # this integration function will use a sparse matrix due to memory constraints
    getPixelIntegrationMatrix = function(integrateByPixel=TRUE) {
      # row is pixel, in the same order as pixelsWithData
      # column is either pixel, indexed in the same way as popGrid, or enumeration area, 
      # indexed in the same way as eaMat/eaDat
      
      if(integrateByPixel) {
        stop("Do not call getPixelIntegrationMatrix with integrateByPixel set to TRUE. 
             Instead take the indices of predMat directly.")
      }
      else {
        allcols = matchMultiple(pixelsWithData, eaToPixel)
        allrows = lapply(1:length(allcols), function(i) {rep(i, l=length(allcols[[i]]))})
        vals = lapply(allcols, function(eaIs) {eaDat$numChildren[eaIs] / sum(eaDat$numChildren[eaIs])})
        sparseMatrix(unlist(allrows), unlist(allcols), x=unlist(vals))
      }
    }
    
    pixelBinomialIntegration = function(thisPixelToEA) {
      logicalIndices = rep(FALSE, nrow(eaDat))
      logicalIndices[thisPixelToEA] = TRUE
      binomialIntegration(logicalIndices, nSamplePixel)
    }
    
    # integrate over the pixels to get aggregated predictions
    if(!allPixels && !onlyInexact)
      pixelPredMatInexact <- expit(predMat[pixelsWithData,1:nSamplePixel])
    else
      pixelPredMatInexact <- expit(predMat[,1:nSamplePixel])
    cat(".")
    
    # integrate over the EAs to get aggregated predictions
    if(genEALevel && !onlyInexact && !continuousOnly) {
      # the following are respectively a matrix of probability draws from the county 
      # distributions and a list of county probability distributions
      if(is.null(eaDat))
        stop("eaDat is null, but a non-null eaDat must be provided for account level predictions for EAs.")
      pixelPredMatExact = noBinomialIntegration(getPixelIntegrationMatrix(FALSE), FALSE, nSamplePixel)
      cat(".")
      pixelToEA = matchMultiple(pixelsWithData, eaToPixel)
      pixelPredMatExactB <- lapply(pixelToEA, pixelBinomialIntegration)
      cat(".")
    } else {
      pixelPredMatExact = NULL
      pixelPredMatExactB = NULL
    }
    
    pixelPreds <- list(pixelPredMatInexact=pixelPredMatInexact, pixelPredMatExact=pixelPredMatExact, 
                       pixelPredMatExactB=pixelPredMatExactB, pixelsWithData=pixelsWithData)
  }
  
  if(genEALevel) {
    eaPreds = list(eaPredMat=eaMat)
  } else {
    eaPreds = NULL
  }
  
  list(mod=mod, preds=preds, SDs=predSDs, obsInds=obsInds, predInds=index, pixelInds=pixelIndices, 
       eaInds=eaIndices, mesh=mesh, prior=prior, stack=stack.full, 
       countyPreds=countyPreds, regionPreds=regionPreds, pixelPreds=pixelPreds, eaPreds=eaPreds)
}

testFitSPDEModel = function(predCoords=NULL, nPredPts=NULL, predUrban=NULL, seed=1234, numClusters=1000) {
  # generate simulated dataset
  sigmasq = .5^2
  effRange = 100
  beta0 = -1.5
  gamma = -.75
  tausq = .25^2
  out = simDat(kenyaEAs, beta0=beta0, margVar=sigmasq, tausq=tausq, gamma=gamma, effRange=effRange, 
               seed=seed, NC=40, nLayer=2, numClusters=numClusters)
  eaDat = out$eaDat
  clustDat = out$clustDat
  
  # load in simulated dataset
  # load("simDat.RData") # "eaDat"    "clustDat" "eaID"
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = rep(25, nrow(obsCoords))
  obsCounts = clustDat$died
  obsUrban = clustDat$urban
  
  # set default prediction locations and urbanicity
  if(is.null(predCoords)) {
    popGrid = makeInterpPopGrid()
    predCoords = cbind(popGrid$east, popGrid$north)
    nPredPts = nrow(predCoords)
    predUrban = popGrid$urban
  }
  
  sampleIs = sample(1:nrow(predCoords), nPredPts, FALSE)
  predCoords = predCoords[sampleIs,]
  predNs = rep(25, nPredPts)
  predUrban = predUrban[sampleIs]
  
  # fit model
  out = fitSPDEModel(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, predNs = rep(25, nrow(predCoords)), 
                     predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=100)
  mod = out$mod
  preds = out$preds
  sds = out$SDs
  obsInds = out$obsInds
  predInds = out$predInds
  mesh = out$mesh
  prior = out$prior
  stack = out$stack
  countyPredMat = out$countyPredMat
  
  # plot things
  # true datasets for all enerumeration areas and the cluster sample
  pdf("figures/spdeDataset.pdf", width=8, height=8)
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
  
  # predictions
  pdf("figures/spdeTestPredictions.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  zlim = range(c(0, preds[obsInds]))
  quilt.plot(obsCoords, obsCounts/obsNs, main="Observed Mortality", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, preds[obsInds], main="Posterior mean at observation locations", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, preds[predInds], main="Posterior Mean", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, sds[predInds], main="Prediction SD", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim)
  plotMapDat(project=TRUE)
  dev.off()
  
  # get difference between urban and rural education
  urbanDiff = mean(preds[obsInds][!obsUrban]) - mean(preds[obsInds][obsUrban])
  print(paste0("mean urban prediction difference: ", urbanDiff))
  
  # plot SPDE parameter posteriors
  result <- inla.spde2.result(mod, "field", prior)
  pdf("figures/spdeTestPosterior1.pdf", width=9, height=5)
  par(mfrow=c(1,2))
  plot(result[["marginals.range.nominal"]][[1]], type = "l",
       main = "Nominal range, posterior density")
  abline(v=effRange, col="green")
  
  plot(result$marginals.variance.nominal[[1]], type = "l",
       main = "Variance, posterior density")
  abline(v=sigmasq, col="green")
  dev.off()
  
  pdf("figures/spdeTestPosterior2.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  plot(mod$marginals.fixed[[1]], type="l", main="Intercept posterior (logit scale)")
  abline(v=beta0, col="green")
  plot(mod$marginals.fixed[[2]], type="l", main="Urban effect posterior (logit scale)")
  abline(v=gamma, col="green")
  # plot(mod$, type="l", main="Urban effect posterior (logit scale)")
  # abline(v=tausq, col="green")
  clustMarg = inla.tmarginal(function(x) {1/exp(x)}, mod$internal.marginals.hyperpar$`Log precision for clust`)
  plot(clustMarg, type="l", main="Cluster variance posterior (logit scale)")
  abline(v=tausq, col="green")
  dev.off()
  
  countyMeans = rowMeans(countyPredMat)
  pdf("figures/spdeTestCountyPreds.pdf", width=5, height=5)
  plotMapDat(plotVar=countyMeans, adm1, project=TRUE, new=TRUE, main="Predicted secondary education completion rate")
  dev.off()
  
  # print summary
  print(summary(mod))
  
  out
}

# same as testFitSPDEModel, accept generates data using spde simulation instead of
# LatticeKrig simulation
testFitSPDEModel2 = function(predCoords=NULL, nPredPts=NULL, predUrban=NULL, seed=1234, 
                             numClusters=1000, popGrid=NULL) {
  # generate simulated dataset
  sigmasq = .5^2
  effRange = 100
  beta0 = -1.5
  gamma = -.75
  tausq = .25^2
  out = simDat(kenyaEAs, beta0=beta0, margVar=sigmasq, tausq=tausq, gamma=gamma, effRange=effRange, 
               seed=seed, NC=40, nLayer=2, numClusters=numClusters)
  eaDat = out$eaDat
  clustDat = out$clustDat
  
  # load in simulated dataset
  # load("simDat.RData") # "eaDat"    "clustDat" "eaID"
  obsCoords = cbind(clustDat$east, clustDat$north)
  obsNs = rep(25, nrow(obsCoords))
  obsCounts = clustDat$died
  obsUrban = clustDat$urban
  
  # set default prediction locations and urbanicity
  if(is.null(predCoords)) {
    if(is.null(popGrid)) {
      popGrid = makeInterpPopGrid()
    }
    predCoords = cbind(popGrid$east, popGrid$north)
    nPredPts = nrow(predCoords)
    predUrban = popGrid$urban
  }
  
  sampleIs = sample(1:nrow(predCoords), nPredPts, FALSE)
  predCoords = predCoords[sampleIs,]
  predNs = rep(25, nPredPts)
  predUrban = predUrban[sampleIs]
  
  # fit model
  out = fitSPDEModel(obsCoords, obsNs=rep(25, nrow(obsCoords)), obsCounts, obsUrban, predCoords, predNs = rep(25, nrow(predCoords)), 
                     predUrban, genCountyLevel=TRUE, popGrid=popGrid, nPostSamples=100)
  mod = out$mod
  preds = out$preds
  sds = out$SDs
  obsInds = out$obsInds
  predInds = out$predInds
  mesh = out$mesh
  prior = out$prior
  stack = out$stack
  countyPredMat = out$countyPredMat
  
  # plot things
  # true datasets for all enerumeration areas and the cluster sample
  pdf("figures/spdeDataset.pdf", width=8, height=8)
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
  
  # predictions
  pdf("figures/spdeTestPredictions.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  zlim = range(c(0, preds[obsInds]))
  quilt.plot(obsCoords, obsCounts/obsNs, main="Observed Mortality", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(obsCoords, preds[obsInds], main="Posterior mean at observation locations", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, preds[predInds], main="Posterior Mean", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim, zlim=zlim)
  plotMapDat(project=TRUE)
  quilt.plot(predCoords, sds[predInds], main="Prediction SD", 
             xlab="Easting", ylab="Northing", xlim=eastLim, ylim=northLim)
  plotMapDat(project=TRUE)
  dev.off()
  
  # get difference between urban and rural education
  urbanDiff = mean(preds[obsInds][!obsUrban]) - mean(preds[obsInds][obsUrban])
  print(paste0("mean urban prediction difference: ", urbanDiff))
  
  # plot SPDE parameter posteriors
  result <- inla.spde2.result(mod, "field", prior)
  pdf("figures/spdeTestPosterior1.pdf", width=9, height=5)
  par(mfrow=c(1,2))
  plot(result[["marginals.range.nominal"]][[1]], type = "l",
       main = "Nominal range, posterior density")
  abline(v=effRange, col="green")
  
  plot(result$marginals.variance.nominal[[1]], type = "l",
       main = "Variance, posterior density")
  abline(v=sigmasq, col="green")
  dev.off()
  
  pdf("figures/spdeTestPosterior2.pdf", width=8, height=8)
  par(mfrow=c(2,2))
  plot(mod$marginals.fixed[[1]], type="l", main="Intercept posterior (logit scale)")
  abline(v=beta0, col="green")
  plot(mod$marginals.fixed[[2]], type="l", main="Urban effect posterior (logit scale)")
  abline(v=gamma, col="green")
  # plot(mod$, type="l", main="Urban effect posterior (logit scale)")
  # abline(v=tausq, col="green")
  clustMarg = inla.tmarginal(function(x) {1/exp(x)}, mod$internal.marginals.hyperpar$`Log precision for clust`)
  plot(clustMarg, type="l", main="Cluster variance posterior (logit scale)")
  abline(v=tausq, col="green")
  dev.off()
  
  countyMeans = rowMeans(countyPredMat)
  pdf("figures/spdeTestCountyPreds.pdf", width=5, height=5)
  plotMapDat(plotVar=countyMeans, adm1, project=TRUE, new=TRUE, main="Predicted secondary education completion rate")
  dev.off()
  
  # print summary
  print(summary(mod))
  
  out
}

simSPDE = function(coords, nsim=1, mesh=NULL, effRange=(max(coords[,1])-min(coords[,1]))/3, margVar=1) {
  # generate mesh grid if necessary
  if(is.null(mesh)) {
    mesh = getSPDEMeshGrid(coords, doPlot = FALSE)
  }
  
  # calculate SPDE model parameters based on Lindgren Rue (2015) "Bayesian Spatial Modelling with R-INLA"
  meshSize <- min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2]))))
  # it is easier to use theta and set sigma0 to 1 then to set sigma0 and the effective range directly
  # kappa0 <- sqrt(8)/effRange * meshSize # since nu = 1
  # kappa0 <- sqrt(8)/effRange # since nu = 1
  # kappa0 = sqrt(8) / 5
  # logKappa = log(kappa0)
  sigma0 = 1
  # tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
  # logTau = log(tau0)
  
  # from page 5 of the paper listed above:
  logKappa = 0.5 * log(8)
  logTau = 0.5 * (lgamma(1) - (lgamma(2) + log(4*pi))) - logKappa
  theta = c(log(sqrt(margVar)), log(effRange))
  spde <- inla.spde2.matern(mesh, B.tau = cbind(logTau, -1, +1),
                            B.kappa = cbind(logKappa, 0, -1), theta.prior.mean = theta,
                            theta.prior.prec = c(0.1, 1))
  
  # generate A and Q precision matrix
  Q = inla.spde2.precision(spde, theta = theta)
  A = inla.spde.make.A(mesh, coords)
  
  # generate simulations
  simField = inla.qsample(nsim, Q)
  simDat = as.matrix(A %*% simField)
  
  simDat
}





