# Kenya U5MR
# setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")
# library(foreign)

# read in the STATA data:
# DHS recode manual:
# https://dhsprogram.com/pubs/pdf/DHSG4/Recode6_DHS_22March2013_DHSG4.pdf
# ?:
# https://dhsprogram.com/pubs/pdf/SAR8/SAR8.pdf
# final report:
# https://dhsprogram.com/pubs/pdf/FR308/FR308.pdf
library(haven)
library(fields)
library(zoo)
library(latex2exp)
library(maptools)
library(data.table)

wd = getwd()
setwd("~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/")

# lines represent births
data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))

# extract the columns of interest
# b1 - month of birth of child,
# b2 - year of birth of child, 
# b4 - gender of child, # TODO: important?
# b5 - child alive at time of interview, 
# b7 - age of child at death in months completed, 
# v024 - region of residence
# v025 - type of place of residence urban rural (32.3% of pop. is urban in 2009 and 2014, see FR p.30/2)
# v001 - cluster number (in theory there are 1,612 total.  In actuality there are 1,593)
# v002 - household number
# v005 - sample weight out of 1000000. normalized weights sum to number of households
# v006 - month of interview
# V007 - year of interview
# V136 - total number of household members
# V138 - number of women in the household aged 15 to 49
# V137 - number of children in the household aged 5 or under
subdata <- data.frame(data[,c('b2', 'b1', 'b5', 'b7', 'v024', 'v025', 'v001', 'v002', 'v005', 'v006', 'v007')])

# extract births in the range 2005 to 2010 (most recent five years)
lowYear <- 2005
highYear <- 2009
subdata <- subdata[(subdata[,'b2'] >= lowYear & subdata[,'b2'] <= highYear),]

# only consider children that died within their first year
died = !is.na(subdata[,'b7'])
totalChildren = nrow(subdata)
subdata = subdata[died,]
subdata <- subdata[subdata[,'b7'] < 12,]
totalFirstYearDied = nrow(subdata)

averageFirstYearMortality = totalFirstYearDied / totalChildren
print(paste0("Average first year mortality: ", averageFirstYearMortality))

# do the same thing for first month
subdata <- data.frame(data[,c('b2', 'b1', 'b5', 'b7', 'v024', 'v025', 'v001', 'v002', 'v005', 'v006', 'v007')])

# extract births in the range 2005 to 2010 (most recent five years)
lowYear <- 2010
highYear <- 2014
subdata <- subdata[(subdata[,'b2'] >= lowYear & subdata[,'b2'] <= highYear),]

# only consider children that died within their first month
died = !is.na(subdata[,'b7'])
totalChildren = nrow(subdata)
subdata = subdata[died,]
subdata <- subdata[subdata[,'b7'] < 1,]
totalFirstMonthDied = nrow(subdata)

averageFirstMonthMortality = totalFirstMonthDied / totalChildren # 0.03909
print(paste0("Average first month mortality: ", averageFirstMonthMortality)) # 0.02159
print(paste0("Average first month mortality rate: ", averageFirstMonthMortality * 12)) # 0.259

# subset by urban/rural
subdata <- data.frame(data[,c('b2', 'b1', 'b5', 'b7', 'v024', 'v025', 'v001', 'v002', 'v005', 'v006', 'v007')])

# extract births in the range 2005 to 2010 (most recent five years)
lowYear <- 2005
highYear <- 2009
subdata <- subdata[(subdata[,'b2'] >= lowYear & subdata[,'b2'] <= highYear),]

# separate urban and rural births
urban = subdata$v025 == 1
urbanData = subdata[urban,]
ruralData = subdata[!urban,]

totalUrbanChildren = nrow(urbanData)
totalRuralChildren = nrow(ruralData)
urbanDied = urbanData[!is.na(urbanData$b7),]
ruralDied = ruralData[!is.na(ruralData$b7),]
firstYearUrban = nrow(urbanDied[urbanDied$b7 < 12,])
firstYearRural = nrow(ruralDied[ruralDied$b7 < 12,])
print(paste0("Average first year mortality urban: ", firstYearUrban / totalUrbanChildren)) # 0.03819
print(paste0("Average first year mortality rural: ", firstYearRural / totalRuralChildren)) # 0.03949

# add a column for the stratification variable as an interaction between
# the urban/rural indicator 'v025' (1: urban, 2:rural) and the region indicator 'v024'
subdata$regionUral <- with(subdata, interaction(v024, v025), drop=TRUE)

# add a column for the unique households with interaction between
# the household indicator 'v002' and the cluster indicator 'v001'
subdata$hhold <- with(subdata, interaction(v001, v002), drop=TRUE)

## find for each cluster the regionUral indicator
clStrat = subdata[,c("v001", "regionUral", "v024", "v025")]
clStrat = clStrat[!duplicated(clStrat), ]
colnames(clStrat) = c("clusterid", "regionRural", "region", "urban-rural")

# get age of child  at the end of 2003 to 2007 period in months
subdata$cage =  12*(2008 - subdata[,"b2"]) + (1 - subdata[,"b1"])

# if child didn't die during time period, pretend they didn't die
ageofdeath = subdata$b7
ageofdeath[is.na(ageofdeath)] = -1
cage = subdata$cage
lived = is.na(subdata$b7) | (ageofdeath > subdata$cage)

# calculate number of trials in binomial by child, age group, and household
getAgeGroup = function(age) {
  if(age == 0)
    return(1)
  else if(age < 12)
    return(2)
  else if(age < 24)
    return(3)
  else if(age < 36)
    return(4)
  else if(age < 48)
    return(5)
  else if(age < 60)
    return(6)
  else
    return(-1)
}
getAgeGroups = function(ages) {
  sapply(ages, getAgeGroup)
}
getTrialsPerAgeGroup = function(tabRow) {
  thisCage = tabRow$cage
  c(1,
    (thisCage >= 12)*11 + (thisCage < 12)*thisCage*(thisCage > 0), 
    (thisCage >= 24)*11 + (thisCage < 24)*(thisCage-11)*(thisCage > 11), 
    (thisCage >= 36)*11 + (thisCage < 36)*(thisCage-23)*(thisCage > 23), 
    (thisCage >= 48)*11 + (thisCage < 48)*(thisCage-35)*(thisCage > 35), 
    (thisCage >= 60)*11 + (thisCage < 60)*(thisCage-47)*(thisCage > 47))
}
trials = (!lived) * (ageofdeath + 1) + lived * cage # Number of trials for each row of subdata in binomial
subdata$trial
# get the number of birth by cluster
n <- table(subdata[,'v001'])
clusterid <- dimnames(n)[[1]]
n.data = data.frame(n=as.vector(n), clusterid=clusterid)

# remove births which are alive or where there is no information on age at death
noinfo = is.na(subdata[,'b7'])
cat("We have information that", sum(!noinfo)/nrow(subdata)*100, "% of the children died")
subdata <- subdata[!noinfo,]

# remove kids that died at an age older than 5 years
older5 <- subdata[,'b7'] > 60
cat(sum(older5)/nrow(subdata)*100, "% of the kids who died, died with an age older than 5 years")
subdata <- subdata[!older5,]

# get number of death below 5 by cluster/village
y <- table(subdata[,'v001'])
clusterid <- dimnames(y)[[1]]

y.data = data.frame(y=as.vector(y), clusterid=clusterid)

# generate one dataset
my.data <- merge(y.data, n.data, by='clusterid', all=TRUE, sort=TRUE)
my.data <- merge(my.data, clStrat, by='clusterid', all=TRUE, sort=TRUE)

# Set NAs in y to zero
my.data$y[is.na(my.data$y)] = 0

# Read geographical information
library(rgdal)
spObj = readOGR(dsn = "Kenya2014gps/", layer = "KEGE71FL")

# Extract (lon, lat) coordinates of all clusters
geoObj = data.frame(cId = spObj$DHSCLUST, lon = spObj$LONGNUM, lat = spObj$LATNUM)

# Extract coordinates of clusters with data
idx = match(my.data$clusterid, geoObj$cId)
my.data$lon = geoObj$lon[idx]
my.data$lat = geoObj$lat[idx]

# Missing geographical information is assigned value (0,0)
# Remove these
missIdx = which(my.data$lon == 0)
my.data = my.data[-missIdx,]

mort = my.data
newNames = names(mort)
newNames[1] = "clusterID"
newNames[6] = "urban"
mort[,6] = mort[,6] == 1
names(mort)=newNames
# save(file = "kenyaData.RData", mort) # don't save until we've added admin1 data from gps dataset





# (kends00ag.asc, admToCounty.txt)
# Read population density map
# popMap = as.matrix(read.table(file = "KenyaPopDens/kends00ag.asc", skip = 6, na.strings = "-9999"))
# lonMap = matrix(rep(33 + (0:239)*0.0416667, each = 288), ncol = 240)
# latMap = matrix(rep(5.9583 - (0:287)*0.0416667, 240), ncol = 240)
# image.plot(lonMap, latMap, popMap, asp = 1)
# pop.data = list(dens = popMap, lonG = lonMap, latG = latMap)
# save(file = "Data/kenyaPop.RData", pop.data)

# Administrative region to county
# admToCounty = read.table('Data/admToCounty.txt', header = TRUE)
# save(file = "Data/kenyaRegionCounty.RData", admToCounty)

# v001 is supposedly cluster number, while v002 is supposedly houshold/dwelling number
# (there can be multiple dwellings per household).  Household numbers restart for each
# cluster.
length(unique(subdata[,'v001']))
# 1593
length(unique(subdata[,'v002']))
# 171

# why are there fewer births after 2008? (for the original subdata, not the modified one)
hist(subdata[,'b2'], xlab="Birth year", freq=F, main="Histogram of birth year", 
     breaks=seq(1975.5, 2014.5, by=1))
abline(v=2008, col="blue")

# for each unique cluster, calculate the number of unique households/dwellings in it, and 
# make a histogram.  Why isn't it always exactly 25? Should this be used in the sample 
# weights instead of 1/25?
test = data.frame(data[,c('b2', 'b5', 'b7', 'v024', 'v025', 'v001', 'v002')])
getNumHHolds = function(v002s) {length(unique(v002s))}
nHHolds = aggregate(test$v002, data.frame(list(v001=test$v001)), getNumHHolds)
head(nHHolds)
hist(nHHolds[,2], xlab="Households sampled per cluster", freq=F, main="Households sampled per cluster")

getNumHHolds = function(v002s) {max(v002s)}
nHHolds = aggregate(test$v002, data.frame(list(v001=test$v001)), getNumHHolds)
head(nHHolds)
hist(nHHolds[,2], xlab="Approximate households", freq=F, main="Approximate households per cluster", 
     breaks=100, col="skyblue")


getNumHHolds = function(v002s) {max(v002s)}
nHHolds = aggregate(data$v002, data.frame(list(v001=data$v001)), getNumHHolds)
head(nHHolds)
hist(nHHolds[,2], xlab="Max Household ID", main="2014 Kenya DHS\nMax Household ID Per Cluster", 
     breaks=100, col="skyblue")

# plot data spatially
quilt.plot(my.data$lon, my.data$lat, my.data$y/my.data$n, xlim=c(33.75, 42.1), 
           ylim=c(-4.9, 5.5), main="Kenya 2003-2008 U5MR")
world(add=TRUE)




###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
# now read gps data
library(SUMMER)
library(foreign)
# gpsDat = read_shape("Kenya2014gps/KEGE71FL.shp")
# out = readShapePoly("Kenya2014gps/KEGE71FL.shp", delete_null_obj=TRUE, force_ring=TRUE, repair=TRUE)
# out = readShapePoly("allDat/KEGE71FL.shp", delete_null_obj=TRUE, force_ring=TRUE)
gpsDat = readShapePoints("Kenya2014gps/KEGE71FL.shp")
# tmp = readShapeLines("Kenya2014gps/KEGE71FL.shp")
plot(gpsDat)
world(add=TRUE)
names(attributes(gpsDat))
coords = attr(gpsDat, "coords")
plot(coords)
world(add=TRUE)
test = coords[,1] < 20 #  these are the observations  whose source is missing.  remove these
sum(test)
gpsDat[test,]

#  remove bad observations
gpsDat=gpsDat[!test,]
plot(gpsDat)
world(add=TRUE)
dim(gpsDat)
names(gpsDat)=c("ID","countryID", "year", "clustID", "countryIDFIPS", "countryIDAdminID", 
                  "admin1FIPS","admin1IDSALB", 
                  "admin1SALB", "admin1ID", "admin1", "regionID", "region", 
                  "source", "urban", "lat", "lon", "altGPS", "altRadar", "coordRef")
gpsDat$altGPS[gpsDat$altGPS == 9999] = NA
gpsDat$altGPS[gpsDat$altRadar == 9999] = NA
gpsDat$urban = gpsDat$urban == "U"
gpsDat$countryIDFIPS[gpsDat$countryIDFIPS == "NULL"] = NA
gpsDat$admin1FIPS[gpsDat$admin1FIPS == "NULL"] = NA
gpsDat$admin1IDSALB[gpsDat$admin1IDSALB == "NULL"] = NA
gpsDat$admin1SALB[gpsDat$admin1SALB == "NULL"] = NA
gpsDat$countryIDAdminID[gpsDat$countryIDAdminID == "NULL"] = NA
head(gpsDat)

save(gpsDat, file="gpsDat.RData")

# get region and admin data from gps data, add to clusters in mort dataset
gpsI = match(data.frame(rbind(mort$lon, mort$lat)), data.frame(rbind(gpsDat$lon, gpsDat$lat)))
mortAdmin1 = gpsDat$admin1[gpsI]
mort$admin1 = mortAdmin1
mortReg = gpsDat$region[gpsI]
mort$region = mortReg
save(mort, file="kenyaData.RData")

#  load in the world population density data
# library(tiff)
library(raster)
# pop = readTIFF("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", convert= TRUE)
pop = raster("Kenya2014Pop/worldpop_total_1y_2014_00_00.tif", values= TRUE)
plot(pop)
names(attributes(pop))
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
totalRows = 4320 # latitude
totalCols = 8640 # longitude
origLonRes = 0.04166667 # 24/degree
origLatRes = 0.04166667
resPerDeg = 24
extentCols = round(kenyaLonLength*resPerDeg)
extentRows = round(kenyaLatLength*resPerDeg)
increaseFac = 1
lonsInterp = seq(kenyaLonRange[1], kenyaLonRange[2], l=round(extentCols*increaseFac))
latsInterp = seq(kenyaLatRange[1], kenyaLatRange[2], l=round(extentRows*increaseFac))
locsInterp = make.surface.grid(list(x=lonsInterp, y=latsInterp))
# test = extract(pop, kenyaExtent)
kenyaPopVals = extract(pop, SpatialPoints(locsInterp),method="bilinear") # this may take quite a while
lonRes = lonsInterp[2] - lonsInterp[1]
latRes = latsInterp[2] - latsInterp[1]

sum(!is.na(kenyaPop))
kenyaPop = data.frame(list(lon=locsInterp[,1], lat=locsInterp[,2], pop=kenyaPopVals))

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
inKenya = in.poly(cbind(kenyaPop$lon, kenyaPop$lat), kenyaPoly) # takes a very long time
kenyaPop = kenyaPop[inKenya,]
dim(kenyaPop) # make sure we have enough points
# [1] >98628     3 for increaseFac=1.9
# [1] 985983      3 for increaseFac=6 (we need a sqrt(10) x increase to get enough points in Nairobi)
# [1] 10966143        3 for increaseFac=20

# renormalize Kenya population
totalKenyaPop = 43.0 * 10^6 # in 2014 from DHS documents
decreaseFac = totalKenyaPop/sum(kenyaPop$pop) # instead of rescaling population using resolution, rescale population exactly for same total
kenyaPop$popOrig = kenyaPop$pop
kenyaPop$pop = kenyaPop$popOrig*decreaseFac
kenyaArea = 224081 # in miles^2
cellArea = kenyaArea/nrow(kenyaPop)
kenyaPop$popDens = kenyaPop$pop/cellArea # density in people/mi^2

# add county and region data
countyDat = getRegion(cbind(kenyaPop$lon, kenyaPop$lat), adm1) # takes a few minutes
# sum(is.na(countyDat$regionNames)) # test to make sure every grid cell has exactly 1 associated county
kenyaPop$admin1 = countyDat$regionNames
# kenyaPop$admin1MultAdmins = countyDat$multipleRegs # this is 0 when increaseFac=1
regions = countyToRegion(kenyaPop$admin1)
kenyaPop$region = regions

save(kenyaPop, file="kenyaPop.RData")

varRange = range(kenyaPop$pop[kenyaPop$pop!= 0])
cols = tim.colors()
# ticks = axisTicks(varRange, log=TRUE)
ticks=c(10^seq(-4, 4, by=2))
plotVar = log10(kenyaPop$pop)
varRange=range(plotVar[kenyaPop$pop != 0])

kenyaEAs = simEAs(kenyaPop)

png("figures/EAsAndPop.png", width=1000, height=600)
par(mfrow=c(1,2))
# set.panel(1,2)
par( oma=c( 0,0,0,5))
plot(kenyaEAs$lon, kenyaEAs$lat, pch=".", col="blue", main=TeX("Enumeration Areas"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
quilt.plot(kenyaPop[,1:2], kenyaPop$pop/cellArea, FUN = function(x) {x = x[x != 0]; log10(mean(x))}, 
           nx=400, ny=400, add.legend=FALSE, main=TeX("Kenya Population Density (people/mi$^2$)"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
par( oma=c(0,0,0,2))
image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
           col=cols, add = TRUE, axis.args=list(at=log10(ticks), labels=ticks))
dev.off()

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
clusterPop = extract(pop, SpatialPoints(cbind(mort$lon, mort$lat)),method="bilinear")
range(clusterPop[ mort$urban])
range(clusterPop[ !mort$urban])
par(mfrow=c(1,2))
hist(clusterPop[ mort$urban], main="Urban population", breaks=30)
hist(clusterPop[ !mort$urban], main="Rural population", breaks=30)

par(mfrow=c(1,2))
par( oma=c( 0,0,0,5))
plot( mort$lon[ mort$urban], mort$lat[ mort$urban], pch=".",col="blue", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
points( mort$lon[ !mort$urban], mort$lat[ !mort$urban], pch=".",col="green")
world(add=TRUE)
quilt.plot(kenyaPop[,1:2], kenyaPop$pop/cellArea, FUN = function(x) {x = x[x != 0]; log10(mean(x))}, 
           nx=400, ny=400, add.legend=FALSE, main=TeX("Kenya Population Density (people/mi$^2$)"), ylim=kenyaLatRange, 
           xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
par( oma=c(0,0,0,2))
image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE,
           col=cols, add = TRUE, axis.args=list(at=log10(ticks), labels=ticks))

par(mfrow=c(1,3))
clusterLog10Pop = log10(clusterPop)
doPred = (clusterPop != 0) & (clusterPop < 40000) & (clusterPop > 250)
modUrban = mort$urban[doPred]
modPop = log10(clusterPop[doPred])
modReg = mort$region[doPred]
modAd = mort$admin1[doPred]
mod = glm(modUrban ~ modPop, data=mort, family=binomial("probit"))
# summary( mod)
clustPopLog10s = seq(1, 6, l=500)
clustPops = 10^clustPopLog10s
probs = predict(mod, list(modPop=clustPopLog10s), type="response")
plot(clusterPop, mort$urban, pch="+", log="x", xlab="Population", ylab="Probability Urban", 
     main="Probability Urban")
lines(clustPops, probs, col="blue", lwd=2)
thetas = coef(mod)
popThresh = 10^(-thetas[1]/thetas[2])  # ~ 18583.71 
abline(v=popThresh, col="green") # threshold: logit,probit,cauchit: 4.260,4.269,4.250
sum((modUrban - fitted(mod))^2) # cauchit works the best
sum((modUrban - round(fitted(mod)))^2) # 402 for all of them, or 372 including region, 279 including county

misclassified = (modUrban - round(fitted(mod))) != 0
par(mfrow=c(1,1))
plot( mort$lon[misclassified], mort$lat[misclassified], type="n", ylim=kenyaLatRange, 
      xlim=kenyaLonRange, xlab="Longitude", ylab="Latitude")
world(add=TRUE)
plotMapDat(adm1)
points( mort$lon[misclassified], mort$lat[misclassified], pch=".",col="red")
points( mort$lon[!misclassified], mort$lat[!misclassified], pch=".",col="green")

# numClusters = table(mort$admin1)
meanPop = aggregate(clusterPop, list(mort$admin1), mean)
numClusters = aggregate(rep(1, nrow(mort)), list(mort$admin1), sum)
par(mfrow=c(1,1))
pdf(file="figures/ClusterByCountyPop.pdf", width=5, height=5)
plot(meanPop$x, numClusters$x, log="x", col="blue", xlab="Mean Cluster Population", ylab="Number of Clusters", 
     main="Sampled Clusters per County")
dev.off()

# do the same but for all EAs
meanPopEA = aggregate(kenyaEAs$pop, list(kenyaEAs$admin1), mean)
numClustersEA = aggregate(rep(1, nrow(kenyaEAs)), list(kenyaEAs$admin1), sum)
par(mfrow=c(1,1))
pdf(file="figures/ClusterByCountyPop.pdf", width=5, height=5)
plot(meanPopEA$x, numClustersEA$x, log="x", col="blue", xlab="Mean Cluster Population", ylab="Number of Clusters", 
     main="Sampled Clusters per County")
dev.off()

# check predictions of the county specific threshold model
clusterLog10Pop = log10(clusterPop)
doPred = (clusterPop != 0) & (clusterPop < 40000) & (clusterPop > 250)
modUrban = mort$urban[doPred]
modPop = log10(clusterPop[doPred])
modReg = mort$region[doPred]
modAd = mort$admin1[doPred]
mod = glm(modUrban ~ modPop + modAd, data=mort, family=binomial("probit"))
preds = mort$urban
preds[doPred] = round(fitted(mod))
pdf(file="figures/clusterFittedUrban.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(mort$lon[preds == 1], mort$lat[preds == 1], pch=".", col="blue", main=TeX("Predicted urbanicity"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(mort$lon[preds != 1], mort$lat[preds != 1], pch=".", col="green")
world(add=TRUE)
plotMapDat(adm1)

plot(mort$lon[mort$urban], mort$lat[mort$urban], pch=".", col="blue", main=TeX("True urbanicity"), xlim=kenyaLonRange, 
     ylim=kenyaLatRange, xlab="Longitude", ylab="Latitude")
points(mort$lon[!mort$urban], mort$lat[!mort$urban], pch=".", col="green")
world(add=TRUE)
plotMapDat(adm1)
dev.off()

# check to see if any county has no urban predictions
aggregate(preds, by=list(mort$admin1), sum) # this is very bad:
# 1          Baringo 10
# 2            Bomet  0
# 3          Bungoma  4
# 4            Busia  5
# 5  Elgeyo Marakwet  2
# 6             Embu 11
# 7          Garissa 11
# 8         Homa Bay  6
# 9           Isiolo 11
# 10         Kajiado 19
# 11        Kakamega  4
# 12         Kericho 15
# 13          Kiambu 23
# 14          Kilifi 12
# 15       Kirinyaga  5
# 16           Kisii  5
# 17          Kisumu 14
# 18           Kitui  9
# 19           Kwale 10
# 20        Laikipia 10
# 21            Lamu 14
# 22        Machakos 20
# 23         Makueni  4
# 24         Mandera 10
# 25        Marsabit 10
# 26            Meru  6
# 27          Migori 12
# 28         Mombasa 36
# 29        Murang'a  1
# 30         Nairobi 56
# 31          Nakuru 20
# 32           Nandi  3
# 33           Narok  5
# 34         Nyamira  0
# 35       Nyandarua  7
# 36           Nyeri 11
# 37         Samburu  9
# 38           Siaya  2
# 39    Taita Taveta 12
# 40      Tana River 10
# 41   Tharaka-Nithi 11
# 42     Trans-Nzoia  8
# 43         Turkana  9
# 44     Uasin Gishu 16
# 45          Vihiga  4
# 46           Wajir 11
# 47      West Pokot  8

##### set number of EAs per county/region
easpc = read.csv("EAdata/tableEA.csv")
easpc$County = as.character(easpc$County)
for(i in 2:ncol(easpc)) {
  easpc[,i] = as.numeric(gsub(",","",as.character(easpc[,i])))
}

# save the dataset
save(easpc, file="easpc.RData")

# compute region totals in each column of easpc
regs = countyToRegion(easpc$County)
getRegTot = function(colDat) {
  aggregate(colDat, list(regs), sum)$x
}
easpr = data.frame(c(list(sort(unique(regs))), list(apply(easpc[,2:ncol(easpc)], 2, getRegTot))))
names(easpr)[1] = "Region"

save(easpr, file="easpr.RData")

##### do the same except with clusters
clustpc = read.csv("EAdata/tableClust.csv")
clustpc$County = as.character(clustpc$County)
clustpc = clustpc[1:47,] # delete the last row, which is blank (not sure why an extra row was included...)
for(i in 2:ncol(clustpc)) {
  clustpc[,i] = as.numeric(gsub(",","",as.character(clustpc[,i])))
}

# save the dataset
save(clustpc, file="clustpc.RData")

# compute region totals in each column of clustpc
regs = countyToRegion(clustpc$County)
getRegTot = function(colDat) {
  aggregate(colDat, list(regs), sum)$x
}
clustpr = data.frame(c(list(sort(unique(regs))), list(apply(clustpc[,2:ncol(clustpc)], 2, getRegTot))))
names(clustpr)[1] = "Region"

save(clustpr, file="clustpr.RData")

##### do the same except for population per county by urban/rural (in 2009.  Inflate to have correct total population for 2014)
poppc = read.csv("EAdata/tablePop.csv")
poppc$County = as.character(poppc$County)
poppc = poppc[1:47,] # delete the last row, which is blank (not sure why an extra row was included...)
for(i in 2:4) {
  poppc[,i] = as.numeric(gsub(",","",as.character(poppc[,i])))
}
increaseFac = totalKenyaPop/sum(poppc[,4])
poppc[,2:4] = round(increaseFac * poppc[,2:4])

# save the dataset
save(poppc, file="poppc.RData")

# get regional population data
regs = countyToRegion(poppc$County)
getRegTot = function(colDat) {
  aggregate(colDat, list(regs), sum)$x
}
poppr = data.frame(c(list(sort(unique(regs))), list(apply(poppc[,2:4], 2, getRegTot))))
names(poppr)[1] = "Region"
poppr$pctTotal = poppr$popTotal/totalKenyaPop
poppr$pctUrb = poppr$popUrb/poppr$popTotal

save(poppr, file="poppr.RData")

## compute the neighbourhood structure as needed by INLA
spP <- SpatialPolygons(adm1@polygons, proj4string = adm1@proj4string)
neighb <- poly2nb(spP)

## save the neighbourhood graph
nb2INLA(file = 'Kenyaadm1.graph', neighb)

setwd(wd)
