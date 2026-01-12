####################################################################################

# NOTE: Before the code can be executed, the R working directory *MUST* 
# be set to the directory 'vim-mc-rf/simulation' (that is
# the directory in which this R script is contained):

# Remove '#' from the line below and replace 'here/is/my/path/' by the path
# to the directory 'vim-mc-rf/simulation':

# setwd("here/is/my/path/")

####################################################################################


# Simulate datasets and obtain runtime estimates of the different VIMs for these:
#################################################################################

source("./simulation_functions.R")

n <- 500
Ks <- c(4, 6, 10)

seed <- 1234

perm_times_sim <- gini_times_sim <- muobj_times_sim <- 0

for (i in seq(along=Ks)) {
  
  K <- Ks[i]
  
  # Simulate dataset:
  
  set.seed(seed)
  dataset <- simdata(n=n, K=K)
  
  
  # Apply the methods:
  
  set.seed(seed)
  perm_times_sim[i] <- system.time({
    perm <- ranger::ranger(dependent.variable.name = "cl", data=dataset, importance="permutation", 
                           num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
                           min.node.size=5)$variable.importance
  })["elapsed"]
  
  set.seed(seed)
  gini_times_sim[i] <- system.time({
    gini_corr <- ranger::ranger(dependent.variable.name = "cl", data=dataset, importance="impurity_corrected",
                                num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                                probability=TRUE, min.node.size=5)$variable.importance
  })["elapsed"]
  
  set.seed(seed)
  muobj_times_sim[i] <- system.time({
    muobj <- diversityForest::multifor(dependent.variable.name = "cl", data=dataset, importance="both", 
                                       npervar = 5000, num.trees=5, replace = FALSE, sample.fraction = 0.7,
                                       probability=TRUE, min.node.size=5)
  })["elapsed"]
  
}



# Save the obtained runtimes:

save(perm_times_sim, gini_times_sim, muobj_times_sim, file="./intermediate_results/runtime.Rda")
