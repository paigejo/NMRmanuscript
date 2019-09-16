source("setup.R")
source("compareModels.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("compareModelCommandArgs.RData")
argList = compareModelCommandArgs[[index]]
do.call("runCompareModels2", argList)