# load packages, data sets, and scripts

library(zoo)
library(latex2exp)
# library(haven)
library(fields)
library(ggplot2)
# library(rgdal)
library(sp)
library(data.table)
# library(raster)
library(sampling)
library(fastmatch)

setwd("~/git/U5MR/")

source('~/git/U5MR/utilityFuns.R')
source('~/git/U5MR/matchPoints.R')
source('~/git/U5MR/simStudy.R')
source('~/git/U5MR/getMercer.R')
source('~/git/U5MR/designBased.R')
source('~/git/U5MR/compareModels.R')
source('~/git/U5MR/validation.R')
source('~/git/U5MR/scores.R')
source('~/git/U5MR/spdeMod.R')
source('~/git/U5MR/spdeResults.R')
source("~/git/U5MR/neonatalSimStudyWeighted.R")
source('~/git/LK-INLA/LKinla_rgeneric.R')
source('~/git/LK-INLA/LKinla.R')

inf = sessionInfo()
if(inf$platform != "x86_64-apple-darwin15.6.0 (64-bit)") {
  INLA:::inla.dynload.workaround()
  # avoid setting too many threads and thereby using too much memory
  inla.setOption(num.threads=1)
  options(error=traceback)
} else {
  options(error=recover)
}

setwd("~/git/U5MR/")

# load Kenya data
out = load("kenyaData.RData")
out = load("kenyaDataEd.RData")
# out = load("kenyaPop.RData")
# tmp = projKenya(kenyaPop$lon, kenyaPop$lat)
# kenyaPop$east = tmp[,1]
# kenyaPop$north = tmp[,2]
# kenyaPop=makeKenyaPop(kmres=5)
# save(kenyaPop, file="kenyaPopProj.RData")
load("kenyaPopProj.RData")

# map data
# adm1 = readRDS("mapData/KEN_adm1.rds")
# adm0 = readRDS("mapData/KEN_adm0.rds")
# save(adm1, adm0, file="adminMapData.RData")
load("adminMapData.RData")

# county to region mapping
ctp = read.csv("mapData/kenya-prov-county-map.csv")

# number of EAs and clusters per county/region
numEAs = 96251
load("easpc.RData")
load("easpr.RData")
load("clustpc.RData")
load("clustpr.RData")
load("poppc.RData")
load("poppr.RData")

# set some final parameters to hold constant throughout analysis
kenyaLonRange = c(33.5, 42)
kenyaLatRange = c(-5, 5.5)
kenyaLonLength = kenyaLonRange[2] - kenyaLonRange[1]
kenyaLatLength = kenyaLatRange[2] - kenyaLatRange[1]
totalKenyaPop = 43*10^6 # from DHS 2014 survey final report page 2
totalRows = 4320 # latitude
totalCols = 8640 # longitude
increaseFac=1
resPerDeg = 24
extentCols = round(kenyaLonLength*resPerDeg)
extentRows = round(kenyaLatLength*resPerDeg)
lonsInterp = seq(kenyaLonRange[1], kenyaLonRange[2], l=round(extentCols*increaseFac))
latsInterp = seq(kenyaLatRange[1], kenyaLatRange[2], l=round(extentRows*increaseFac))
lonRes = lonsInterp[2] - lonsInterp[1]
latRes = latsInterp[2] - latsInterp[1]

# get limits of easting/northing for plotting
# tmp = projKenya(kenyaLonRange, kenyaLatRange)
# eastLim = tmp[,1]
# northLim = tmp[,2]
# save(eastLim, northLim, file="lims.RData")
load("lims.RData")

# generate 5km population density grid over Kenya
# popGrid = makeInterpPopGrid(kmRes=5)
# save(popGrid, file="popGrid.RData")
# popGrid = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE)
# save(popGrid, file="popGridAdjusted.RData")
# popGrid = makeInterpPopGrid(kmRes=5, adjustPopSurface=TRUE, "women")
# save(popGrid, file="popGridAdjustedWomen.RData")

# set enumeration areas
# kenyaEAs = simEAs2(popGrid, numEAs, totalKenyaPop, fixNumUrbanAtTruth=TRUE)
# save(kenyaEAs, file="kenyaEAs.RData")
load("kenyaEAs.RData")

# load 5km population density grid over Kenya
# load("popGrid.RData")
load("popGridAdjusted.RData")

# project mort dataset lon/lat coords to easting/westing in km
# tmp = projKenya(mort$lon, mort$lat)
# mort$east = tmp[,1]
# mort$north = tmp[,2]
# save(mort, file="mortProj.RData")
# load("mortProj.RData")

# Bernoulli datasets
out = load("data4direct.RData")

# empirical distributions
# empiricalDistributions = getSurveyEmpiricalDistributions2(maxAge=0)
# empiricalDistributions = c(empiricalDistributions, getSurveyEmpiricalDistributionsWomen())
# save(empiricalDistributions, file="empiricalDistributions.RData")

# parallelization
if(!exists("doParallel") || (exists("doParallel") && doParallel == FALSE)) {
  assign("doParallel", FALSE, envir=.GlobalEnv)
  assign("cores", NULL, envir=.GlobalEnv)
  assign("cl", NULL, envir=.GlobalEnv)
}