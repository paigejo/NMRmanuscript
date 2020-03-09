source("setup.R")
source("getDirect.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("directCommandArgsNew.RData")
argList = directCommandArgs[[index]]
do.call("getDirectNaive", argList)