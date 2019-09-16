# an update of the readDat.R script with cleaner code for educational attainment


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

wd1 = getwd()
wd2 = "~/Google Drive/UW/Wakefield/WakefieldShared/U5MR/"
setwd(wd2)

# lines represent births
# data <- data.frame(read_dta("Kenya2014BirthRecode/KEBR70FL.DTA"))
data <- data.frame(read_dta("Kenya2014IndividualRecode/KEIR71FL.DTA"))

## variables from the birth recode
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

## variables from this (individual/women) recode
# v024 - region of residence
# v025 - type of place of residence urban rural (32.3% of pop. is urban in 2009 and 2014, see FR p.30/2)
# v001 - cluster number (in theory there are 1,612 total.  In actuality there are 1,593)
# v002 - household number
# v005 - sample weight out of 1000000. normalized weights sum to number of households
# v006 - month of interview
# V007 - year of interview
# v012 - current age in completed years
# v106 - highest education level attended in order of: no educ, primary, secondary, higher
# v107 - highest year comleted (years completed at level given in v106)
# v149 - educational achievement: none (0), incomplete primary (1), completed primary (2), incomplete secondary (3), 
#        complete secondary (4), higher education (5)
subdata <- data.frame(data[,c('v024', 'v025', 'v001', 'v002', 'v005', 'v012', 'v106', 'v107', 'v149')])

# extract women with the given age range
lowAge = 20
highAge = 29
subdata <- subdata[(subdata[,'v012'] >= lowAge & subdata[,'v012'] <= highAge),]

# add a column for the stratification variable as an interaction between
# the urban/rural indicator 'v025' (1: urban, 2:rural) and the region indicator 'v024'
subdata$regionUral <- with(subdata, interaction(v024, v025), drop=TRUE)

# add a column for the unique households with interaction between
# the household indicator 'v002' and the cluster indicator 'v001'
subdata$hhold <- with(subdata, interaction(v001, v002), drop=TRUE)

# find for each cluster the regionUral indicator
clStrat = subdata[,c("v001", "regionUral", "v024", "v025", "v005")]
clStrat = clStrat[!duplicated(clStrat), ]
colnames(clStrat) = c("clusterID", "regionRural", "region", "urban", "samplingWeight")
clStrat$urban =  clStrat$urban == 1

# determine whether each woman completed their secondary education
completed = subdata$v149 >= 4
completed[is.na(completed)] = TRUE

# get the number of women by cluster
n <- table(subdata[,'v001'])
clusterid <- dimnames(n)[[1]]
n.data = data.frame(clusterID=clusterid, n=as.vector(n))

# get the number of women who completed their secondary education by cluster
# y <- table(subdata[,'v001'])
y = aggregate(completed, list(subdata[,'v001']), sum)
n.data$y = y$x

# add in strata
ed <- merge(n.data, clStrat, by='clusterID', all=TRUE, sort=TRUE)

# Read geographical information
library(rgdal)
spObj = readOGR(dsn = "Kenya2014gps/", layer = "KEGE71FL")

# Extract (lon, lat) coordinates of all clusters
geoObj = data.frame(cId = spObj$DHSCLUST, lon = spObj$LONGNUM, lat = spObj$LATNUM)

# Extract coordinates of clusters with data
idx = match(ed$clusterID, geoObj$cId)
ed$lon = geoObj$lon[idx]
ed$lat = geoObj$lat[idx]

# Missing geographical information is assigned value (0,0)
# Remove these
missIdx = which(ed$lon == 0)
ed = ed[-missIdx,]

library(SUMMER)
library(foreign)

gpsDat = readShapePoints("Kenya2014gps/KEGE71FL.shp")
coords = attr(gpsDat, "coords")
plot(coords)
world(add=TRUE)
test = coords[,1] < 20 #  these are the observations  whose source is missing.  remove these
sum(test)

#  remove observations with unknown locations and set unknown data to NA
gpsDat=gpsDat[!test,]
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

setwd(wd1)
save(gpsDat, file="gpsDatEd.RData")
setwd(wd2)

# get region and admin data from gps data, add to clusters in ed dataset
gpsI = match(data.frame(rbind(ed$lon, ed$lat)), data.frame(rbind(gpsDat$lon, gpsDat$lat)))
ed$admin1 = gpsDat$admin1[gpsI]
ed$region = gpsDat$region[gpsI]

# get easting and northing using projection
tmp = projKenya(ed$lon, ed$lat)
ed$east = tmp[,1]
ed$north = tmp[,2]
setwd(wd1)
save(ed, file="kenyaDataEd.RData")
setwd(wd2)

aggregate(ed$y/ed$n, by=list(ed$urban), mean)


setwd(wd1)
