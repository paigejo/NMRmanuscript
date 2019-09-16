# Kenya U5MR

setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")


# read in the STATA data:
# DHS recode manual:
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf

library(haven)
library(fields)

##### first plot our actual dataset

pdf(file="figures/mortDat.pdf", width=5, height=5)
par(mfrow=c(1,1))
quilt.plot(mort$lon, mort$lat, mort$y/mort$n, main=TeX("Survey Empirical U5MR (2003-2007)"), 
           xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, ylim=kenyaLatRange, nx=120, ny=120)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

pdf(file="figures/EAsUrban.pdf", width=5, height=5)
urban = kenyaEAs$urban
plot(kenyaEAs$lon[!urban], kenyaEAs$lat[!urban], pch=".", col="green", main=TeX("Urban vs. rural enumeration areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(kenyaEAs$lon[urban], kenyaEAs$lat[urban], pch=".", col="blue")
world(add=TRUE)
plotMapDat(adm1)
dev.off()

# lines represent births
data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))

length(unique(data$v001)) # 1593 unqiue clusters
urbRur = by(as.numeric(data$v025), factor(data$v001), mean)
sum(urbRur==1)
# [1] 617
sum(urbRur==2)
# [1] 976

# load Kenya data
out = load("kenyaData.RData")
rm(y)
attach(mort)
out = load("gpsDat.RData")
# out

library(glmm)
mod = glmm(y ~ region + urban, list(~region, ~clusterID), c("region", "clustID"), family=binomial.glmm)

# test urbanicity Data set
library(raster)
library(rgdal)
urban = raster("GHS_SMOD_POP2015_GLOBE_R2016A_54009_1k_v1_0/GHS_SMOD_POP2015_GLOBE_R2016A_54009_1k_v1_0.tif")

# first project to be on same scale as pop dataset
urban = projectRaster(urban, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", method="ngb")
lonRange=c(-180, 180)
latRange=c(-90,90)
kenyaLonRange = c(33.5, 42)
kenyaLatRange = c(-5,5.5)
kenyaLonLength = kenyaLonRange[2] - kenyaLonRange[1]
kenyaLatLength = kenyaLatRange[2] - kenyaLatRange[1]
kenyaExtent = extent(c(xmin=kenyaLonRange[1], xmax=kenyaLonRange[2],
                       ymin=kenyaLatRange[1], ymax=kenyaLatRange[2]))
numEAs = 96251
#  Number of rows and columns at the original resolution
# >  range(rI)
# [1]  805 1008
# >  range(cI)
# [1] 4082 4560
totalRows = 15236 # latitude
totalCols = 35497 # longitude
origLonRes = 360/totalCols # ~98.60/degree
origLatRes = 180/totalRows # ~84.64/degree
resPerDegLon = 1/origLonRes
resPerDegLat = 1/origLatRes
extentCols = round(kenyaLonLength*resPerDegLon)
extentRows = round(kenyaLatLength*resPerDegLat)
increaseFac = 2
lonsInterp = seq(kenyaLonRange[1], kenyaLonRange[2], l=round(extentCols*increaseFac))
latsInterp = seq(kenyaLatRange[1], kenyaLatRange[2], l=round(extentRows*increaseFac))
locsInterp = make.surface.grid(list(x=lonsInterp, y=latsInterp))
# test = extract(urban, kenyaExtent)
kenyaUrbanVals = extract(urban, SpatialPoints(locsInterp),method="simple") # all of these methods have
# kenyaUrbanVals2 = extract(urban, data.frame(lon=locsInterp[,1], lat=locsInterp[,2]), method="simple")
lonRes = lonsInterp[2] - lonsInterp[1]
latRes = latsInterp[2] - latsInterp[1]

sum(!is.na(kenyaUrbanVals))
sum(kenyaUrbanVals != round(kenyaUrbanVals)) # should be 0
kenyaUrban = data.frame(list(lon=locsInterp[,1], lat=locsInterp[,2], urban=kenyaUrbanVals))

# read in administrative map areas
# https://stackoverflow.com/questions/17723822/administrative-regions-map-of-a-country-with-ggmap-and-ggplot2
library(ggplot2)
library(rgdal)
library(sp)
# out = load("mapData/KEN_adm1.rds")
# pakistan.adm2.spdf <- get("gadm")
adm1 = readRDS("mapData/KEN_adm1.rds")
plot(adm1)
adm0 = readRDS("mapData/KEN_adm0.rds")
plot(adm0)
names(adm0)

# subset population density to be within Kenya
polys = adm0@polygons
kenyaPoly = polys[[1]]@Polygons[[77]]@coords
plot(kenyaPoly, type="l")
inKenya = in.poly(cbind(kenyaUrban$lon, kenyaUrban$lat), kenyaPoly)
kenyaUrban = kenyaUrban[inKenya,]
dim(kenyaUrban) # make sure we have much more than 96,000 points
# [1] 396516     3 for increaseFac=1
# [1] 1588080       3 for increaseFac=2

# renormalize Kenya population
kenyaUrban$urbanOrig = kenyaUrban$urban

save(kenyaUrban, file="kenyaUrban.RData")

varRange = range(kenyaUrban$urban)
cols = tim.colors(3)
plotVar = kenyaUrban$urban
varRange=range(plotVar)
#############################################
kenyaEAs = simEAs(kenyaPop)
png("figures/EAsAndUrbanicity.png", width=1000, height=600)
par(mfrow=c(1,2))
# set.panel(1,2)
par( oma=c( 0,0,0,5))
plot(kenyaEAs$lon, kenyaEAs$lat, pch=".", col="blue", main=TeX("Enumeration Areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
quilt.plot(kenyaUrban[,1:2], kenyaUrban$urban, 
           nx=400, ny=400, main=TeX("Kenya urbanicity"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
dev.off()

###################################
numClusters = 1593
clustersPerCounty = aggregate(rep(1, nrow(mort)), list(mort$admin1), sum)
clusters = kenyaEAs[sample(1:numEAs, numClusters, replace= FALSE),]
clusters2 = kenyaEAs[sampleByStratum(kenyaEAs$admin1, clustersPerCounty$x),]
png("figures/EAsPopClusters.png", width=1000, height=1000)
par(mfrow=c(2,2))
# set.panel(1,2) I
par( oma=c( 0,0,0,5))
plot(clusters2$lon, clusters2$lat, typ="n", main=TeX("Simulated Clusters (by County)"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
plotMapDat(adm1, lwd=.5)
points(clusters2$lon, clusters2$lat, pch=".", col="blue")
plot(mort$lon, mort$lat, type="n", main=TeX("True Clusters"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
plotMapDat(adm1, lwd=.5)
points(mort$lon, mort$lat, pch=".", col="blue")
plot(clusters$lon, clusters$lat, type="n", main=TeX("Simulated Clusters"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
plotMapDat(adm1, lwd=.5)
points(clusters$lon, clusters$lat, pch=".", col="blue")
plot(kenyaPop[,1:2], type="n", main=TeX("Kenya Population Density (people/mi$^2$)"), ylim=kenyaLatRange, 
     xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
quilt.plot(kenyaPop[,1:2], kenyaPop$pop/cellArea, FUN = function(x) {x = x[x != 0]; log10(mean(x))}, 
           nx=400, ny=400, add.legend=FALSE, add=TRUE)
plotMapDat(adm1, lwd=.5)
world(add=TRUE)
par( oma=c(0,0,0,2))
image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
           col=cols, add = TRUE, axis.args=list(at=log10(ticks), labels=ticks))
dev.off()


### determine the threshold for urban/rural
## first get population density at cluster locations
# clusterPop = extract(pop, SpatialPoints(cbind(mort$lon, mort$lat)),method="bilinear")
clusterUrban = extract(urban, SpatialPoints(cbind(mort$lon, mort$lat)),method="simple")
range(clusterUrban[ mort$urban])
range(clusterUrban[ !mort$urban])
pdf(file="figures/urbanicityByUrbanHist.pdf", width=8, height=5)
par(mfrow=c(1,2))
breaks = seq(-.5, 3.5, by=1)
hist(clusterUrban[ mort$urban], main="Urban Classifications", breaks=breaks, freq=F, xlab="Urbanicity")
hist(clusterUrban[ !mort$urban], main="Rural Classifications", breaks=breaks, freq=F, xlab="Urbanicity")
dev.off()

freqUrban = table(clusterUrban[ mort$urban])
freqRural = table(clusterUrban[ !mort$urban])
total = freqUrban + freqRural
probUrban = freqUrban/total
pdf("figures/empiricalProbUrbanVsUrbanicity.pdf", width=5, height=5)
barplot(probUrban, names.arg=c("0", "1", "2", "3"), main="Empirical probability of urban", 
        xlab="Urbanicity", ylim=c(0,1))
abline(h=.5, lty=2, col="blue")
dev.off()
probUrban
# 0         1         2         3 
# 0.2376238 0.2605042 0.4300000 0.9850746 

png("figures/UrbanicityVsUrban.png", width=1000, height=600)
par(mfrow=c(1,2))
par( oma=c( 0,0,0,5))
plot( mort$lon[ mort$urban], mort$lat[ mort$urban], pch=".",col="blue", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("DHS Cluster Urban Classification"))
points( mort$lon[ !mort$urban], mort$lat[ !mort$urban], pch=".",col="green")
world(add=TRUE)
quilt.plot(kenyaUrban[,1:2], kenyaUrban$urban, 
           nx=400, ny=400, add.legend=FALSE, main=TeX("Kenya Urbanicity"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
dev.off()
# par( oma=c(0,0,0,2))
# image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
#            col=cols, add = TRUE, axis.args=list(at=log10(ticks), labels=ticks))

par(mfrow=c(1,3))
urbanicity = factor(clusterUrban)
modUrban = mort$urban
modReg = mort$region
modAd = mort$admin1
mod = glm(modUrban ~ clusterUrban, data=mort, family=binomial("probit"))
mod = glm(modUrban ~ clusterUrban + modReg, data=mort, family=binomial("probit"))
mod = glm(modUrban ~ clusterUrban + modAd, data=mort, family=binomial("probit"))
summary( mod)
urbanLevels = c(0, 1, 2, 3)
probs = predict(mod, list(clusterUrban=urbanLevels), type="response")
plot(log10(clusterPop), urban, pch="+")
lines(clustLog10s, probs, col="blue", lwd=2)
thetas = coef(mod)
abline(v=-thetas[1]/thetas[2], col="green") # threshold: logit,probit,cauchit: 4.260,4.269,4.250
sum((modUrban - fitted(mod))^2) # cauchit works the best
sum((modUrban - round(fitted(mod)))^2) # 402 for all of them, or 372 including region, 279 including county


###### plot thresholded population and see if it matches up reasonably with the data
thresh = 18583.71
threshes = c(8000, 10000, 12500, 15000, 17500, 18583)
for(i in 1:length(threshes)) {
  thresh = threshes[i]
  png(paste0("figures/PopThresh", thresh, ".png"), width=1000, height=600)
  par(mfrow=c(1,2))
  par( oma=c( 0,0,0,5))
  plot( mort$lon[ mort$urban], mort$lat[ mort$urban], pch=".",col="blue", ylim=kenyaLatRange, 
        xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("DHS Cluster Urban Classification"))
  points( mort$lon[ !mort$urban], mort$lat[ !mort$urban], pch=".",col="green")
  world(add=TRUE)
  quilt.plot(kenyaPop[,1:2], kenyaPop$popOrig > thresh, FUN=max, col=c("green", "blue"),
             nx=400, ny=400, main=TeX(paste0("Urban Population Threshold (", thresh, ")")), ylim=kenyaLatRange, 
             xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
  world(add=TRUE)
  dev.off()
}

##### test sampling clusters by region and urban/rural
par(mfrow=c(1,1))
plot(easpc[,2]/easpc[,4], clustpc[,2]/clustpc[,4], xlim=c(0,1), ylim=c(0,1), col="blue", 
     xlab="Percent urban clusters in county", ylab="Percent urban EAs in county", 
     main="Percent urban EAs vs clusters by county")
abline(0,1)

# threshes = setThresholds2()
# eaThreshes = sapply(1:nrow(kenyaEAs), function(i) {threshes$threshes[threshes$counties == kenyaEAs$admin1[i]]})
# urban = kenyaEAs$popOrig > eaThreshes
pdf(file="figures/EAsClustsUrban.pdf", width=8, height=5)
par(mfrow=c(1,2))
urban = kenyaEAs$urban
plot(kenyaEAs$lon[!urban], kenyaEAs$lat[!urban], pch=".", col="green", main=TeX("Enumeration Areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(kenyaEAs$lon[urban], kenyaEAs$lat[urban], pch=".", col="blue")
world(add=TRUE)
plotMapDat(adm1)

kenyaClusts = simClusters2(kenyaEAs)$clustDat
urban = kenyaClusts$urban
plot(kenyaClusts$lon[!urban], kenyaClusts$lat[!urban], pch=".", col="green", main=TeX("Enumeration Areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(kenyaClusts$lon[urban], kenyaClusts$lat[urban], pch=".", col="blue")
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### test new thresholds based on population:
threshes = setThresholds2()
clusterPop = extract(pop, SpatialPoints(cbind(mort$lon, mort$lat)),method="bilinear")
urban = setUrbanByThreshes(clusterPop, mort$admin1, threshes)
pdf("figures/urbanicityThreshTest.pdf", width=8, height=5)
par(mfrow=c(1,2))
plot( mort$lon[urban], mort$lat[urban], pch=".",col="blue", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("Our Cluster Urban Classification"))
points( mort$lon[ !urban], mort$lat[ !urban], pch=".",col="green")
world(add=TRUE)
plotMapDat(adm1)

plot( mort$lon[mort$urban], mort$lat[mort$urban], pch=".",col="blue", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("DHS Cluster Urban Classification"))
points( mort$lon[ !mort$urban], mort$lat[ !mort$urban], pch=".",col="green")
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### test simulating data (on logit scale not including nugget should vary between -2.5 and -2)
set.seed(580252)
beta0 = -2
margVar = .15^2
# tausq = .1^2
tausq = .2^2
gamma = -.5
HHoldVar = 0
# HHoldVar = .3^2
HHoldVar = .4^2
urbanOverSample = 2
effRange = 150
out = simDat2(kenyaEAs, clustDat=NULL, nsim=1, urbanOverSample=2,
              beta0=beta0, margVar=margVar, tausq=tausq, gamma=gamma, HHoldVar=HHoldVar, 
              effRange=effRange)
eaDat = out$eaDat
clustDat = out$clustDat

# compute empirical mortality rates
eaPHat = eaDat$died/eaDat$numChildren
clustPHat = clustDat$died/clustDat$numChildren
allDat = c(eaPHat, clustPHat, eaDat$trueProbDeath, eaDat$trueProbDeathNoNug)
allDat = allDat[is.finite(allDat)]
colRange = quantile(allDat, probs=c(0.025, 0.975))

# plot simulated enumeration areas, clusters, and urbanicity
par(mfrow=c(1,2))
plot( eaDat$lon[ !eaDat$urban], eaDat$lat[ !eaDat$urban], pch=".",col="green", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("EAs and urban classification"))
points( eaDat$lon[eaDat$urban], eaDat$lat[eaDat$urban], pch=".",col="blue")
world(add=TRUE)
plotMapDat(adm1)

plot( clustDat$lon[ !clustDat$urban], clustDat$lat[ !clustDat$urban], pch=".",col="green", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude", main=TeX("Clusters and urban classification"))
points( clustDat$lon[clustDat$urban], clustDat$lat[clustDat$urban], pch=".",col="blue")
world(add=TRUE)
plotMapDat(adm1)

# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb2.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

# true probabilities without nugget/cluster effect
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaDat$trueProbDeathNoNug, zlim=colRange, nx=200, ny=200, FUN = min, 
           main="Median true probabilities", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# true probabilities with nugget/cluster effect
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaDat$trueProbDeath, zlim=colRange, nx=200, ny=200, FUN = min, 
           main="True mortality rate (cluster/nugget RE)", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=50, ny=50, 
           main="EA empirical mortality rate", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# empirical probabilities in clusters
quilt.plot(cbind(clustDat$lon, clustDat$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster empirical mortality rate", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### how much do we oversample each urban clusters in each county
pdf(file="figures/urbanOversampling.pdf", width=5, height=5)
plot(1,1, main="County urban oversampling", col="blue", type="n", 
     xlab="Proportion or urban EAs", ylab="Proportion of urban clusters", xlim=c(0,1), ylim=c(0,1))
abline(0,1)
points(propUrbPerCounty(eaDat), propUrbPerCounty(clustDat), col="blue")
dev.off()

##### test new cluster sampling scheme (simClusters3) and simDat2
overSampDat = simDat2(kenyaEAs, clustDat=NULL, nsim=4, urbanOverSample=2, 
                      beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5)
eaDat = overSampDat$eaDat
clustDat1 = overSampDat$clustDat[[1]]
clustDat2 = overSampDat$clustDat[[2]]
clustDat3 = overSampDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3OverUrban.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=50, ny=50, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### do the same except without oversampling
tmpClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=4)
clustList = genAndreaFormatFromEAIs(overSampDat$eaDat, tmpClustDat$eaIs, tmpClustDat$sampleWeights)
SRSDat = list(eaDat=overSampDat$eaDat, clustDat=clustList)
eaDat = SRSDat$eaDat
clustDat1 = SRSDat$clustDat[[1]]
clustDat2 = SRSDat$clustDat[[2]]
clustDat3 = SRSDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3SRS.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=50, ny=50, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=50, ny=50, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()


##### high-res versions:

eaDat = overSampDat$eaDat
clustDat1 = overSampDat$clustDat[[1]]
clustDat2 = overSampDat$clustDat[[2]]
clustDat3 = overSampDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3OverUrbanHR.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=200, ny=200, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### do the same except without oversampling
eaDat = SRSDat$eaDat
clustDat1 = SRSDat$clustDat[[1]]
clustDat2 = SRSDat$clustDat[[2]]
clustDat3 = SRSDat$clustDat[[3]]

# plot 4 simulations
# plot true and empirical probability of death in EAs and clusters along with urbanicity
pdf(file="figures/simDatEmpiricalProb3SRSHR.pdf", width=10, height=10)
par(mfrow=c(2,2))
par( oma=c( 0,0,0,3))

## empirical probabilities in enumeration areas
quilt.plot(cbind(eaDat$lon, eaDat$lat), eaPHat, zlim=colRange, nx=200, ny=200, 
           main="EA empirical mortality rate", xlab="", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

## empirical probabilities in clusters
# sim 1
quilt.plot(cbind(clustDat1$lon, clustDat1$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 1", xlab="", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 2
quilt.plot(cbind(clustDat2$lon, clustDat2$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 2", xlab="Longitude", ylab="Latitude", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)

# sim 3
quilt.plot(cbind(clustDat3$lon, clustDat3$lat), clustPHat, zlim=colRange, nx=70, ny=70, 
           main="Cluster simulation 3", xlab="Longitude", ylab="", xlim=kenyaLonRange, 
           ylim=kenyaLatRange)
world(add=TRUE)
plotMapDat(adm1)
dev.off()

##### show urban oversampling versus SRS variation
overSampDat = simDat2(kenyaEAs, clustDat=NULL, nsim=4, urbanOverSample=2, 
                      beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5)
eaDat = overSampDat$eaDat
clustDat = overSampDat$clustDat[[1]]
clustDat2 = overSampDat$clustDat[[2]]
clustDat3 = overSampDat$clustDat[[3]]

tmpClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=4)
clustList = genAndreaFormatFromEAIs(overSampDat$eaDat, tmpClustDat$eaIs, tmpClustDat$sampleWeights)
SRSDat = list(eaDat=overSampDat$eaDat, clustDat=clustList)

pdf(file="figures/clusterSimsAll.pdf", width=8, height=12)
par(mfrow=c(3,2))
clustDat = SRSDat$clustDat[[1]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("SRS urban vs. rural 1"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = overSampDat$clustDat[[1]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("Stratified urban vs. rural 1"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = SRSDat$clustDat[[2]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("SRS urban vs. rural 2"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = overSampDat$clustDat[[2]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("Stratified urban vs. rural 2"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = SRSDat$clustDat[[3]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("SRS urban vs. rural 3"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

clustDat = overSampDat$clustDat[[3]]
urban = clustDat$urban
plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX("Stratified urban vs. rural 3"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
world(add=TRUE)
plotMapDat(adm1)

dev.off()

overSampDat = simDat2(kenyaEAs, clustDat=NULL, nsim=9, urbanOverSample=2, 
                      beta0=-2, margVar=.15^2, tausq=.1^2, gamma=-.5)
eaDat = overSampDat$eaDat
pdf(file="figures/clusterSimsUrbanOver.pdf", width=12, height=12)
par(mfrow=c(3,3))

for(i in 1:9) {
  clustDat = overSampDat$clustDat[[i]]
  urban = clustDat$urban
  plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX(paste0("Stratified urban vs. rural ", i)), xlim=kenyaLonRange, 
       ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
  points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
  world(add=TRUE)
  plotMapDat(adm1)
}

dev.off()

tmpClustDat = simClusters3(kenyaEAs, urbanOverSample=1, nsim=9)
clustList = genAndreaFormatFromEAIs(overSampDat$eaDat, tmpClustDat$eaIs, tmpClustDat$sampleWeights)
SRSDat = list(eaDat=overSampDat$eaDat, clustDat=clustList)
pdf(file="figures/clusterSimSRS.pdf", width=12, height=12)
par(mfrow=c(3,3))

for(i in 1:9) {
  clustDat = SRSDat$clustDat[[i]]
  urban = clustDat$urban
  plot(clustDat$lon[!urban], clustDat$lat[!urban], pch=19, col="green", main=TeX(paste0("SRS urban vs. rural ", i)), xlim=kenyaLonRange, 
       ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude", cex=.3)
  points(clustDat$lon[urban], clustDat$lat[urban], pch=19, col="blue", cex=.3)
  world(add=TRUE)
  plotMapDat(adm1)
}
dev.off()

# examine census and surveys to try and find empirical distribution of households per cluster
library(haven)
library(fields)
library(zoo)
library(latex2exp)
library(maptools)
setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")


data <- data.frame(read_dta("popSurvey2009/Housing_2009KPHC_10PCT_STATA.dta"))

names(data)
test = data.frame(RecreatedEANO=data$RecreatedEANO, HHNO=data$HHNO)
test = test[!is.na(data$HHNO), ]
head(test)
getNumHHolds = function(x) {length(unique(x))}
nHHolds = aggregate(data$HHNO, data.frame(list(x=data$RecreatedEANO)), getNumHHolds)
head(nHHolds) # these seem way too large
hist(nHHolds[,2])
mean(nHHolds[,2] >= 150) # 0.152. Too large

getNumHHolds = function(x) {max(x)}
nHHolds = aggregate(test$HHNO, data.frame(list(x=test$RecreatedEANO)), getNumHHolds)
head(nHHolds)
hist(nHHolds[,2]) # this is wayyy too large
mean(nHHolds[,2] >= 150) # 0.459. Way too large

data <- data.frame(read_dta("popSurvey2009/Population_2009KPHC_10PCT_STATA.dta"))
childrenBorn = data[!is.na(data$P24), c("P24", "P25")]
totalChildrenBorn = rowSums(childrenBorn)
hist(totalChildrenBorn, breaks=seq(-.5, max(totalChildrenBorn) + .5, by=1), col="skyblue", 
     freq=FALSE, main="Histogram of Total Children Born", xlab="Total Children Born")

pdf("figures/totalChildrenBorn.pdf", width=5, height=8)
par(mfrow=c(2, 1))
childrenBorn = data[(!is.na(data$P24)) & (data$EATYPE == 2), c("P24", "P25")]
totalChildrenBorn = rowSums(childrenBorn)
hist(totalChildrenBorn, breaks=seq(-.5, max(totalChildrenBorn) + .5, by=1), col="skyblue", 
     freq=FALSE, main="Histogram of Total Children Born (Urban)", xlab="Total Children Born", 
     ylim=c(0, .4))

childrenBorn = data[(!is.na(data$P24)) & (data$EATYPE == 1), c("P24", "P25")]
totalChildrenBorn = rowSums(childrenBorn)
hist(totalChildrenBorn, breaks=seq(-.5, max(totalChildrenBorn) + .5, by=1), col="skyblue", 
     freq=FALSE, main="Histogram of Total Children Born (Rural)", xlab="Total Children Born", 
     ylim=c(0, .4))
dev.off()

library(foreign)
owners = read.spss("housingSurvey2012/Owners_Household data.sav", to.data.frame=TRUE)
renters = read.spss("housingSurvey2012/Renters_Household data.sav", to.data.frame=TRUE)
ownersi = read.spss("housingSurvey2012/Owners_Individual data.sav", to.data.frame=TRUE)
rentersi = read.spss("housingSurvey2012/Renters_Individual data.sav", to.data.frame=TRUE)
names(owners)
length(unique(owners$CLUSTER))
length(unique(renters$CLUSTER))
length(unique(ownersi$CLUSTER))
length(unique(rentersi$CLUSTER))

length(unique(owners$CLUSTER)) + length(unique(renters$CLUSTER))
test = c(owners$CLUSTER, renters$CLUSTER, ownersi$CLUSTER, rentersi$CLUSTER)

length(unique(test)) # 1263 (due to three counties not being included yet? Or nonresponse?)
5360 / 4 # 1340

test = data.frame(CLUSTER=c(owners$CLUSTER, renters$CLUSTER), 
                  HOUSEHOLD_NUMBER=c(owners$HOUSEHOLD_NUMBER, renters$HOUSEHOLD_NUMBER))
naHouseholds = is.na(test$HOUSEHOLD_NUMBER)
test = test[!naHouseholds, ]
getNumHHolds = function(x) {length(unique(x))}
nHHolds = aggregate(test$HOUSEHOLD_NUMBER, data.frame(list(x=test$CLUSTER)), getNumHHolds)
head(nHHolds)
hist(nHHolds[,2])
sum(nHHolds[,2]) # 18423 (goal was to complete 19,140 according to page 12)
15 *nrow(nHHolds) # 18945 (this seems to indicate there were some clusters missing)
(19140 - 15 *nrow(nHHolds)) / 15 # 13 clusters missing apparently

getNumHHolds = function(x) {max(x)}
nHHolds = aggregate(test$HOUSEHOLD_NUMBER, data.frame(list(x=test$CLUSTER)), getNumHHolds)
head(nHHolds)
hist(nHHolds[,2], breaks=30, col="skyblue", main="2012/13 Kenya Household Survey\nMax Household ID Per Cluster", 
     xlab="Max Household ID")
# the following proportions should be 0... (could this be due to people moving in and out of EAs 
# and the extra houses being renumbered?)
mean(nHHolds[,2] >= 150)
mean(nHHolds[,2] <= 49)
nHHolds[nHHolds[,2] >= 150,]

# now to find empirical distribution of mothers per household
data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))

getNumMothers = function(x) {max(x)}
nHHolds = aggregate(data$v002, data.frame(list(v001=data$v001)), getNumHHolds)
head(nHHolds)
hist(nHHolds[,2], xlab="Max Household ID", main="2014 Kenya DHS\nMax Household ID Per Cluster", 
     breaks=100, col="skyblue")

##### Plot calculated empirical distributions stratified by urban/rural
load("empiricalDistributions.RData")
pdf("figures/empiricalHouseholdDistributionsECDF.pdf", width=5, height=8)
par(mfrow=c(3, 1))
xRange = range(c(knots(empiricalDistributions$households), knots(empiricalDistributions$householdsUrban), 
                 knots(empiricalDistributions$householdsRural)))
plot(empiricalDistributions$households, main="Households Per Cluster", xlab="Households Per Cluster", 
     ylab="Empirical CDF", xlim=xRange)
plot(empiricalDistributions$householdsUrban, main="Households Per Cluster (Urban)", xlab="Households Per Cluster", 
     ylab="Empirical CDF", xlim=xRange)
plot(empiricalDistributions$householdsRural, main="Households Per Cluster (Rural)", xlab="Households Per Cluster", 
     ylab="Empirical CDF", xlim=xRange)
dev.off()

pdf("figures/empiricalMothersDistributionsECDF.pdf", width=5, height=8)
par(mfrow=c(3, 1))
xRange = range(c(knots(empiricalDistributions$mothers), knots(empiricalDistributions$mothersUrban), 
                 knots(empiricalDistributions$mothersRural)))
plot(empiricalDistributions$mothers, main="Mothers Per Household", xlab="Mothers Per Household", 
     ylab="Empirical CDF", xlim=xRange)
plot(empiricalDistributions$mothersUrban, main="Mothers Per Household (Urban)", xlab="Mothers Per Household", 
     ylab="Empirical CDF", xlim=xRange)
plot(empiricalDistributions$mothersRural, main="Mothers Per Household (Rural)", xlab="Mothers Per Household", 
     ylab="Empirical CDF", xlim=xRange)
dev.off()

pdf("figures/empiricalChildrenDistributionsECDF.pdf", width=5, height=8)
par(mfrow=c(3, 1))
xRange = range(c(knots(empiricalDistributions$children), knots(empiricalDistributions$childrenUrban), 
                 knots(empiricalDistributions$childrenRural)))
plot(empiricalDistributions$children, main="Children Per Mother", xlab="Children Per Mother", 
     ylab="Empirical CDF", xlim=xRange)
plot(empiricalDistributions$childrenUrban, main="Children Per Mother (Urban)", xlab="Children Per Mother", 
     ylab="Empirical CDF", xlim=xRange)
plot(empiricalDistributions$childrenRural, main="Children Per Mother (Rural)", xlab="Children Per Mother", 
     ylab="Empirical CDF", xlim=xRange)
dev.off()

##### now plot the empirical pdfs
# household per cluster
by = 5
householdKnots = knots(empiricalDistributions$households)
householdKnotsUrban = knots(empiricalDistributions$householdsUrban)
householdKnotsRural = knots(empiricalDistributions$householdsRural)
xRange = range(c(knots(empiricalDistributions$households), knots(empiricalDistributions$householdsUrban), 
                 knots(empiricalDistributions$householdsRural)))
breaks = seq(-1, ceiling(xRange[2]/by)*by, by=by)
householdVals = empiricalDistributions$households(breaks[2:length(breaks)]) - 
  empiricalDistributions$households(breaks[1:(length(breaks) - 1)])
householdValsUrban = empiricalDistributions$householdsUrban(breaks[2:length(breaks)]) - 
  empiricalDistributions$householdsUrban(breaks[1:(length(breaks) - 1)])
householdValsRural = empiricalDistributions$householdsRural(breaks[2:length(breaks)]) - 
  empiricalDistributions$householdsRural(breaks[1:(length(breaks) - 1)])

xRange[2] = xRange[2] + .5
yRange = c(0, max(c(householdVals, householdValsUrban, householdValsRural))) / by
pdf("figures/empiricalHouseholdDistributionsPMF.pdf", width=5, height=8)
nBuffer = match(1, breaks[2:length(breaks)] >= xRange[1]) - 1
nRest = length(householdVals) - nBuffer
cVec = c(rep(rgb(1,1,1,0), nBuffer), rep("black", nRest))
par(mfrow=c(3, 1), xpd=FALSE)
barplot(householdVals / by, main="Households Per Cluster", xlab="Households Per Cluster", 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue", space=0, 
        axes=FALSE, border=cVec)
axis(1,at=seq(50, 150, by=50),labels=seq(50, 150, by=50))
axis(2,at=c(0, .005, .01, .015),labels=c("", "0.005", "", "0.015"))
barplot(householdValsUrban / by, main="Households Per Cluster (Urban)", xlab="Households Per Cluster", 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue", space=0, 
        axes=FALSE, border=cVec)
axis(1,at=seq(50, 150, by=50),labels=seq(50, 150, by=50))
axis(2,at=c(0, .005, .01, .015),labels=c("", "0.005", "", "0.015"))
barplot(householdValsRural / by, main="Households Per Cluster (Rural)", xlab="Households Per Cluster", 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue", space=0, 
        axes=FALSE, border=cVec)
axis(1,at=seq(50, 150, by=50),labels=seq(50, 150, by=50))
axis(2,at=c(0, .005, .01, .015),labels=c("", "0.005", "", "0.015"))
dev.off()

# mothers per household
by = 1
motherKnots = knots(empiricalDistributions$mothers)
motherKnotsUrban = knots(empiricalDistributions$mothersUrban)
motherKnotsRural = knots(empiricalDistributions$mothersRural)
xRange = range(c(knots(empiricalDistributions$mothers), knots(empiricalDistributions$mothersUrban), 
                 knots(empiricalDistributions$mothersRural)))
breaks = seq(ceiling(xRange[1]/by)*by - by, ceiling(xRange[2]/by)*by, by=by)
motherVals = empiricalDistributions$mothers(breaks[2:length(breaks)]) - 
  empiricalDistributions$mothers(breaks[1:(length(breaks) - 1)])
motherValsUrban = empiricalDistributions$mothersUrban(breaks[2:length(breaks)]) - 
  empiricalDistributions$mothersUrban(breaks[1:(length(breaks) - 1)])
motherValsRural = empiricalDistributions$mothersRural(breaks[2:length(breaks)]) - 
  empiricalDistributions$mothersRural(breaks[1:(length(breaks) - 1)])
xRange[2] = xRange[2] + .5
yRange = c(0, max(c(motherVals, motherValsUrban, motherValsRural))) / by
pdf("figures/empiricalMothersDistributionsPMF.pdf", width=5, height=8)
par(mfrow=c(3, 1))
barplot(motherVals / by, main="Mothers Per Household", xlab="Mothers Per Household", space=0, 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue")
axis(1,at=0:max(xRange)+.5,labels=0:max(xRange))
barplot(motherValsUrban / by, main="Mothers Per Household (Urban)", xlab="Mothers Per Household", space=0, 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue")
axis(1,at=0:max(xRange)+.5,labels=0:max(xRange))
barplot(motherValsRural / by, main="Mothers Per Household (Rural)", xlab="Mothers Per Household", space=0, 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue")
axis(1,at=0:max(xRange)+.5,labels=0:max(xRange))
dev.off()

# children per mother
by = 1
childrenKnots = knots(empiricalDistributions$children)
childrenKnotsUrban = knots(empiricalDistributions$childrenUrban)
childrenKnotsRural = knots(empiricalDistributions$childrenRural)
xRange = range(c(knots(empiricalDistributions$children), knots(empiricalDistributions$childrenUrban), 
                 knots(empiricalDistributions$childrenRural)))
breaks = seq(ceiling(xRange[1]/by)*by - by, ceiling(xRange[2]/by)*by, by=by)
childrenVals = empiricalDistributions$children(breaks[2:length(breaks)]) - 
  empiricalDistributions$children(breaks[1:(length(breaks) - 1)])
childrenValsUrban = empiricalDistributions$childrenUrban(breaks[2:length(breaks)]) - 
  empiricalDistributions$childrenUrban(breaks[1:(length(breaks) - 1)])
childrenValsRural = empiricalDistributions$childrenRural(breaks[2:length(breaks)]) - 
  empiricalDistributions$childrenRural(breaks[1:(length(breaks) - 1)])
xRange[1] = 0
xRange[2] = xRange[2] + .5
yRange = c(0, max(c(childrenVals, childrenValsUrban, childrenValsRural))) / by
pdf("figures/empiricalChildrenDistributionsPMF.pdf", width=5, height=8)
par(mfrow=c(3, 1))
barplot(childrenVals / by, main="Children Per Mother", xlab="Children Per Mother", space=0, 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue")
axis(1,at=0:(max(xRange)-1)+.5,labels=1:max(xRange))
barplot(childrenValsUrban / by, main="Children Per Mother (Urban)", xlab="Children Per Mother", space=0, 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue")
axis(1,at=0:(max(xRange)-1)+.5,labels=1:max(xRange))
barplot(childrenValsRural / by, main="Children Per Mother (Rural)", xlab="Children Per Mother", space=0, 
        ylab="Binned Empirical PMF", xlim=xRange, ylim=yRange, width=by, col="skyblue")
axis(1,at=0:(max(xRange)-1)+.5,labels=1:max(xRange))
dev.off()

# expected households per cluster
ecdfExpectation(empiricalDistributions$householdsUrban) # 92.81475
ecdfExpectation(empiricalDistributions$householdsRural) # 87.7717

# expected children per cluster (assuming independence of children per mother and mothers per household conditional on urbanicity)
ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
  ecdfExpectation(empiricalDistributions$childrenUrban) # 42.06381
ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
  ecdfExpectation(empiricalDistributions$childrenRural) # 61.15194

# expected children per household (assuming independence of children per mother and mothers per household conditional on urbanicity)
ecdfExpectation(empiricalDistributions$mothersUrban) * ecdfExpectation(empiricalDistributions$childrenUrban) # 0.4532017 urban
ecdfExpectation(empiricalDistributions$mothersRural) * ecdfExpectation(empiricalDistributions$childrenRural) # 0.6967159 rural

# Integrate population density within each county and stratum
popGrid = makeInterpPopGrid(kmRes=5)
counties=sort(unique(poppc$County))
getCountyStratumIntegrationMatrix = function(getUrban=TRUE) {
  counties = as.character(counties)
  
  mat = t(sapply(counties, function(countyName) {popGrid$admin1 == countyName}))
  mat = sweep(mat, 2, popGrid$popOrig, "*")
  sweep(mat, 2, popGrid$urban == getUrban, "*")
}
urbanPopulations = rowSums(getCountyStratumIntegrationMatrix())
ruralPopulations = rowSums(getCountyStratumIntegrationMatrix(FALSE))
popTable = cbind(urban=urbanPopulations, rural=ruralPopulations, pctUrban=urbanPopulations / (urbanPopulations + ruralPopulations))

# compare that to the average proportion of children that are in urban and rural strata per county
load("~/git/U5MR/empiricalDistributions.RData")
sortI = sort(easpc$County, index.return=TRUE)$ix
temp = easpc[sortI,]
childrenPerStratumUrban = temp$EAUrb * ecdfExpectation(empiricalDistributions$householdsUrban) * ecdfExpectation(empiricalDistributions$mothersUrban) * 
  ecdfExpectation(empiricalDistributions$childrenUrban)
childrenPerStratumRural = temp$EARur * ecdfExpectation(empiricalDistributions$householdsRural) * ecdfExpectation(empiricalDistributions$mothersRural) * 
  ecdfExpectation(empiricalDistributions$childrenRural)
childTable = cbind(urban=childrenPerStratumUrban, rural=childrenPerStratumRural, pctUrban=childrenPerStratumUrban / (childrenPerStratumUrban + childrenPerStratumRural))

compareTable = cbind(pctUrbanPop=popTable[,3], pctUrbanChild=childTable[,3])
format(compareTable, digits=1)
colMeans(compareTable)
# pctUrbanPop pctUrbanChild 
# 0.2536728     0.2308785 

popEstimateUnintegrated = 0.2536728 * mean(expit(-1.777 -1.003)) + (1 - 0.2536728) * mean(expit(-1.777 + rnorm(100000, sd=.103)))
childEstimate = 0.2308785 * mean(expit(-1.777 -1.003+ rnorm(100000, sd=.103))) + (1 - 0.2308785) * mean(expit(-1.777+ rnorm(100000, sd=.103)))
childEstimate - popEstimateUnintegrated
# [1] 0.002076907

# test sampling weights:
getWeightedMean = function(clusterDat) {
  weights = clusterDat$samplingWeight
  weights = weights / sum(weights)
  props = clusterDat$died / clusterDat$numChildren
  sum(props * weights)
}
getWeightedSum = function(clusterDat, N=3906862) {
  weights = clusterDat$samplingWeight
  weights = weights / sum(weights)
  died = clusterDat$died
  sum(died * weights) * (N / sum(clusterDat$numChildren))
}
# if(tausq == .1^2) {
#   # out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
#   if(test) {
#     out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
#   } else {
    out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOverSamplefrac0.25.RData")
#   }
# }
# else {
#   if(tausq != 0)
#     stop("tausq can only be equal to .1^2 or 0")
#   
#   # out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOver2.RData")
#   if(test) {
#     out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0.25Test.RData")
#   } else {
    out = load("simDataMultiBeta-1.75margVar0.0225tausq0gamma-1HHoldVar0urbanOverSamplefrac0.25.RData")
#   }
# }
mean(sapply(SRSDat$clustDat, getWeightedMean))
mean(SRSDat$eaDat$died / SRSDat$eaDat$numChildren)
sum(SRSDat$eaDat$died) / sum(SRSDat$eaDat$numChildren)

mean(sapply(SRSDat$clustDat, function(x) {sum(x$samplingWeight)}))

sum(SRSDat$eaDat$numChildren)
mean(sapply(SRSDat$clustDat, getWeightedSum))
sum(SRSDat$eaDat$died)

## are urban areas overrepresented by the DHS?
# see if urban clusters are overrepresented/oversampled (yes they are when averaging over county, but total urban versus total rural is good)
mean((clustpc$clustUrb / clustpc$clustRur)  / (easpc$EAUrb / easpc$EARur), na.rm=TRUE)
mean((clustpc$clustUrb / clustpc$clustRur)  - (easpc$EAUrb / easpc$EARur), na.rm=TRUE)
sum(clustpc$clustUrb) / sum(clustpc$clustTotal)
sum(easpc$EAUrb) / sum(easpc$EATotal)
cbind(clustpc$clustTotal, (clustpc$clustUrb / clustpc$clustRur)  / (easpc$EAUrb / easpc$EARur))
plot(clustpc$clustTotal, (clustpc$clustUrb / clustpc$clustRur)  / (easpc$EAUrb / easpc$EARur)) # as the total number of clusters increases, urban areas are over5sampled less

# see if urban children are overrepresented/oversampled (yes they are when averaging over county, but total urban versus total rural is good)
mean(((clustpc$clustUrb * 0.4532017) / (clustpc$clustRur * 0.6967159))  / ((easpc$EAUrb * 0.4532017) / (easpc$EARur * 0.6967159)), na.rm=TRUE)
mean(((clustpc$clustUrb * 0.4532017) / (clustpc$clustRur * 0.6967159)) - ((easpc$EAUrb * 0.4532017) / (easpc$EARur * 0.6967159)), na.rm=TRUE)
sum(clustpc$clustUrb * 0.4532017) / sum(clustpc$clustRur * 0.6967159)
sum(easpc$EAUrb * 0.4532017) / sum(easpc$EARur * 0.6967159)

# expected children per household (assuming independence of children per mother and mothers per household conditional on urbanicity)
ecdfExpectation(empiricalDistributions$mothersUrban) * ecdfExpectation(empiricalDistributions$childrenUrban) # 0.4532017 urban
ecdfExpectation(empiricalDistributions$mothersRural) * ecdfExpectation(empiricalDistributions$childrenRural) # 0.6967159 rural


##### test how well the pearson distribution works for approximating sums of binomials
# in the first case, assume the probabilities are normally distributed on the logit scale
# also, assume exchangeable correlation structure
logitMuFirst = -2.6
logitMuSecond = -1.7
logitVar = .1^2

# calculate some "scoring rules" (mean, variance, crps)
# to calculate crps, take a fake draw or calculate at the mean
allN = c(1, 10, 1000)
allNSims = c(1, 10, 100, 1000) # 1 signals to take mean
scoringRulesPearson = matrix(nrow = length(allN) * length(allNSims) * length(allNSims), ncol = 11)
scoringRulesBinom = matrix(nrow = length(allN) * length(allNSims) * length(allNSims), ncol = 11)
thisRow = 1
for(i in 1:length(allN)) {
  nFirst = allN[i]
  nSecond = allN[i]
  
  # generate logit mean vector
  logitMu = c(rep(logitMuFirst, nFirst), rep(logitMuSecond, nSecond))
  
  # generate covariance matrix and Cholesky decomposition for simulations
  rho = .75 # rho should be high for pixel, and low for county
  Sigma = matrix(rho, nrow = nFirst + nSecond, ncol = nFirst + nSecond)
  diag(Sigma) = 1
  Sigma = logitVar * Sigma
  L = t(chol(Sigma))
  
  for(j in 1:length(allNSims)) {
    nSimProbs = allNSims[j]
    
    # simulate values of ps
    if(nSimProbs != 1) {
      logitProbs = L %*% matrix(rnorm((nFirst + nSecond) * nSimProbs), nrow = nFirst + nSecond)
      logitProbs = sweep(logitProbs, 1, logitMu, "+")
      probMat = expit(logitProbs)
    }
    else {
      probMat = matrix(expit(logitMu), ncol=1)
    }
    
    # use a pearson distribution approximation
    xs = 0:round(25 * (nFirst + nSecond))
    pearsonTime = system.time(pearsonApproximation <- dSumBinomRandom(xs, 25, probMat))[3]
    
    # generate false data based on the mode of the pearson approximation
    index = which.max(pearsonApproximation)
    falseData = xs[index]
    
    # calculate the "scoring rules" for the pearson approximation
    pearsonMean = sum(xs * pearsonApproximation)
    pearsonVar = sum(xs^2 * pearsonApproximation) - pearsonMean^2
    pearsonCrps = crpsCounts(falseData, pearsonApproximation)
    pearsonTime = pearsonTime + system.time(pearsonCi <- generateBinomialInterval(pearsonApproximation))[3]
    
    for(k in 1:length(allNSims)) {
      
      # generate simulation from a binomial for approximation
      nSimBinom = allNSims[k]
      simTime = system.time(simApproximation <- dSumBinomRandomSim(xs, 25, probMat, nSim=nSimBinom))[3]
      
      # calculate the "scoring rules" for the binomial approximation
      simMean = sum(xs * simApproximation)
      simVar = sum(xs^2 * simApproximation) - simMean^2
      simCrps = crpsCounts(falseData, simApproximation)
      simTime = simTime + system.time(simCi <- generateBinomialInterval(simApproximation))[3]
      
      # update scoring rules tables
      scoringRulesPearson[thisRow, ] = c(nFirst * 2, nSimProbs, nSimBinom, 
                                         pearsonMean, pearsonVar, pearsonCrps, pearsonTime, 
                                         pearsonCi)
      scoringRulesBinom[thisRow, ] = c(nFirst * 2, nSimProbs, nSimBinom, 
                                       simMean, simVar, simCrps, simTime, 
                                       simCi)
      
      # plot the distributions
      maxMass = max(c(pearsonApproximation, simApproximation))
      
      pdf(paste0("figures/pearsonTest_", nFirst, "_", nSecond, "_nSimProb", nSimProbs, "_nSimBinom", nSimBinom, ".pdf"), 
          width=5, height=5)
      plot(xs / (25 * (nFirst + nSecond)), pearsonApproximation, pch=19, ylim=c(0, maxMass), main="Pearson Versus Binomial Simulation With 80% CIs", 
           ylab="Probability Mass", xlab="Proportion Mortality", cex=.1, xlim=c(0, 1/4))
      points(xs / (25 * (nFirst + nSecond)), simApproximation, pch=19, cex=.1, col="blue")
      legend("topright", c("Pearson", "Simulation"), pch=19, col=c("black", "blue"))
      abline(v=pearsonCi[1] / (25 * (nFirst + nSecond)), lty=2)
      abline(v=pearsonCi[2] / (25 * (nFirst + nSecond)), lty=2)
      abline(v=simCi[1] / (25 * (nFirst + nSecond)), lty=2, col="blue")
      abline(v=simCi[2] / (25 * (nFirst + nSecond)), lty=2, col="blue")
      dev.off()
      
      thisRow = thisRow + 1
    }
  }
}

colnames(scoringRulesPearson) = c("totClusters", "nSimProb", "nSimBinom", 
                                  "mean", "var", "crps", "time", "lower80", "upper80", 
                                  "leftReject80", "rightReject80")
colnames(scoringRulesBinom) = c("totClusters", "nSimProb", "nSimBinom", 
                                "mean", "var", "crps", "time", "lower80", "upper80", 
                                "leftReject80", "rightReject80")
scoringRulesPearson = data.frame(scoringRulesPearson)
scoringRulesBinom = data.frame(scoringRulesBinom)
scoringRulesAll = cbind(scoringRulesPearson[, 1:4], scoringRulesBinom$mean, scoringRulesPearson$var, 
                        scoringRulesBinom$var, scoringRulesPearson$crps, scoringRulesBinom$crps, 
                        scoringRulesPearson$time, scoringRulesBinom$time, 
                        paste0("(", round(scoringRulesPearson$lower80, 4), ", ", round(scoringRulesPearson$upper80, 4), ")"), 
                        paste0("(", round(scoringRulesBinom$lower80, 4), ", ", round(scoringRulesBinom$upper80, 4), ")"), 
                        paste0("(", round(scoringRulesPearson$leftReject80, 4), ", ", round(scoringRulesPearson$rightReject80, 4), ")"), 
                        paste0("(", round(scoringRulesBinom$leftReject80, 4), ", ", round(scoringRulesBinom$rightReject80, 4), ")"))
names(scoringRulesAll) = c(names(scoringRulesPearson)[1:3], "meanPearson", "meanBinom", 
                           "varPearson", "varBinom", "crpsPearson", "crpsBinom", 
                           "timePearson", "timeBinom", 
                           "CI80Pearson", "CI80Binom", "reject80Pearson", "reject80Binom")
cbind(round(scoringRulesAll[,1:11], 4), scoringRulesAll[,12:15])
save(scoringRulesAll, "pearsonTests.RData")

##### test parallelization
n = 5000
N = 25
truth = rep(12, n) / N
preds = expit(rnorm(n))
system.time(out <- crpsBinomial(truth, preds, N))[3]
system.time(out <- crpsBinomial(truth, preds, N, parClust = NULL))[3]

n = 5000
m = 100
ns = rep(25, n)
pMat = matrix(.5, nrow=n, ncol=m)
k = 0:sum(ns)
system.time(out <- dSumBinomRandom(k, ns, pMat))[3]
system.time(out <- dSumBinomRandom(k, ns, pMat, parClust=NULL))[3]

##### test intra-cluster correlation
logitMu = seq(-3, -1, length=30)
logitVarEps = seq(.1^2, 2, length=30)

# calculate some "scoring rules" (mean, variance, crps)
# to calculate crps, take a fake draw
nSim = 100
n = 25
tausq = .1^2
correlationMat = matrix(nrow=length(logitMu), ncol=length(logitVarEps))
thisRow = 1
for(i in 1:length(logitMu)) {
  thisLogitMu = logitMu[i]
  
  for(j in 1:length(logitVarEps)) {
    thisLogitVarEps = logitVarEps[j]
    
    correlationMat[i, j] = logitNormCor(muSigmaMat = matrix(c(thisLogitMu, sqrt(thisLogitVarEps)), ncol=2))
    
    if(is.na(correlationMat[i, j])) {
      debug(logitNormCor)
      logitNormCor(muSigmaMat = matrix(c(thisLogitMu, sqrt(thisLogitVarEps)), ncol=2))
    }
  }
}

image.plot(logitMu, logitVarEps, correlationMat, xlab="Logit Mean", ylab="Nugget Variance", 
           main="Within Cluster Correlation")
abline(v=-1.75, col="green", lwd=2)
abline(v=-2.75, col="blue", lwd=2)

out = load("simDataMultiBeta-1.75margVar0.0225tausq0.01gamma-1HHoldVar0urbanOver2.RData")
pdf("~/Google Drive/UW/Wakefield/WakefieldShared/UM5R/figures/mortalityCountHist.pdf", width=5, height=15)
par(mfrow=c(3,1))
hist(mort$y, freq=FALSE, main="Mortality Count Histogram (Empirical)", xlab="Mortality Counts", breaks=seq(0, 11) - .5)
hist(SRSDat$clustDat[[1]]$died, freq=FALSE, main="Mortality Count Histogram (SRS)", xlab="Mortality Counts", breaks=seq(0, 11)-.5)
hist(overSampDat$clustDat[[1]]$died, freq=FALSE, main="Mortality Count Histogram (Urban Over Sampled)", xlab="Mortality Counts", breaks=seq(0, 11) - .5)
dev.off()