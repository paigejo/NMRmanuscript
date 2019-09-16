# sets up global variables and libraries to run code in parallel
assign("doParallel", TRUE, envir=.GlobalEnv)

library(parallel)
library(doParallel)
library(foreach)
setwd("~/git/NMRmanuscript//")
source("setup.R")

# assign("cl", detectCores() - 1, envir=.GlobalEnv)
assign("cl", makeCluster(3, outfile="log.txt"), envir=.GlobalEnv)
clusterEvalQ(cl, {setwd("~/git/NMRmanuscript/"); source("setup.R")})
registerDoParallel(cl)
# options(error=recover)
# stopCluster(cl)