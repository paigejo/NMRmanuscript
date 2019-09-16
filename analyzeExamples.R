source("setup.R")
source("neonatalSimStudyWeighted.R")
source('generateExampleResults.R')

# generate example results for the education attainment dataset
generateExampleResults(ed, resultNameRoot="Ed")

# generate example results for the neonatal mortality dataset
generateExampleResults(mort, resultNameRoot="Mort")

