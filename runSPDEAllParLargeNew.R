source("setupParallelLarge.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("spdeCommandArgsNew.RData")
argList = spdeCommandArgs[[index]]
do.call("resultsSPDE", argList)
stopCluster(cl)