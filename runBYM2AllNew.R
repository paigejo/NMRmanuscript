source("setup.R")
source("designBased.R")
index = as.numeric(commandArgs(trailingOnly = TRUE))
load("bym2CommandArgsNew.RData")
argList = bym2CommandArgs[[index]]
do.call("runBYM2", argList)