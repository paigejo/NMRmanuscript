source("setup.R")
source("getMercer.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("mercerCommandArgs.RData")
argList = mercerCommandArgs[[index]]
do.call("getMercer", argList)