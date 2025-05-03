# Set working directory to 'vim-mc-rf/simulation':

setwd("~/vim-mc-rf/simulation")



# Number of parallel cores to use:

# NOTE: You can choose a smaller or larger number here, depending on the resources available on your system.

ncores_used <- 100


# Load necessary libraries

library(parallel)
library(doParallel)



# Make table of settings (each row in the table contains the information for one iteration):

n <- c(100, 500, 1000, 2000)
K <- c(4, 6, 10)
itind <- 1:1000

scenariogrid <- expand.grid(itind=itind, 
                            n=n, K=K, stringsAsFactors = FALSE)
scenariogrid <- scenariogrid[,ncol(scenariogrid):1, drop=FALSE]


# Add a specific seed to each row, such that the individual iterations are reproducible:
set.seed(1234)
seeds <- sample(1000:10000000, size=length(n)*length(K)*length(itind))

scenariogrid$seed <- rep(seeds, each=1)

nrow(scenariogrid)



# Randomly permute the rows of "scenariogrid" so that the computational burden
# is ensured to be distributed evenly across the parallel nodes:

set.seed(1234)
reorderind <- sample(1:nrow(scenariogrid))
scenariogrid <- scenariogrid[reorderind,,drop=FALSE]
rownames(scenariogrid) <- NULL

scenariogrid$settingid <- 1:nrow(scenariogrid)
rownames(scenariogrid) <- 1:nrow(scenariogrid)



# Save scenariogrid, needed in evaluation of the results:

save(scenariogrid, file="./intermediate_results/scenariogrid_simulation.Rda")




# Source the functions that are used in performing the calculations 
# on the cluster:

source("./simulation_functions.R")




# Start the cluster:

cl <- makeCluster(ncores_used, type = "PSOCK")



# Register the parallel backend:

registerDoParallel(cl)

mywd <- getwd()



# Export the objects in the workspace to the
# parallel jobs:

clusterExport(cl, varlist = ls())




# Perform the calculations:

Results <- parLapply(cl, 1:nrow(scenariogrid), function(z) {
  try({
    evaluatesetting(z, mywd=mywd)
  })
})



# Save the results:

save(Results, file="./intermediate_results/results_simulation.Rda")



# Stop the cluster:

stopCluster(cl)
