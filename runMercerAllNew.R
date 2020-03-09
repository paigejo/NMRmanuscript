source("setup.R")
source("getMercer.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("mercerCommandArgsNew.RData")
argList = mercerCommandArgs[[index]]
do.call("getMercer", argList)