source("setup.R")
source("validation.R")
startFrom = as.numeric(commandArgs(trailingOnly = TRUE))
do.call("validateExample", list(startFrom=startFrom, dat=mort))