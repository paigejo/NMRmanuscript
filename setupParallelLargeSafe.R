# sets up global variables and libraries to run code in parallel
assign("doParallel", TRUE, envir=.GlobalEnv)

library(parallel)
library(doParallel)
library(foreach)
setwd("~/git/U5MR/")
source("setup.R")

# assign("cl", detectCores() - 1, envir=.GlobalEnv)
assign("cl", makeCluster(5, outfile="log.txt"), envir=.GlobalEnv)
clusterEvalQ(cl, {setwd("~/git/U5MR/"); source("setup.R")})
registerDoParallel(cl)
options(error=recover)
# stopCluster(cl)