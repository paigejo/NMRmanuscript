source("setup.R")
source("getDirect.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("directCommandArgs.RData")
argList = directCommandArgs[[index]]
do.call("getDirectNaive", argList)