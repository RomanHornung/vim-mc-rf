# This function performs one repetition of the simulation. It simulates
# one dataset and estimates the four different variable importance metrics
# for this dataset.
#
# It takes the whole number 'iter', which corresponds to the iter-th line 
# of 'scenariogrid', which contains the necessary information
# on the iter-th repetition.

# Input parameters:

# iter  - iter-th line of 'scenariogrid', which corresponds to the iter-th line 
#         of 'scenariogrid', which contains the necessary information
#         on the iter-th repetition.
# mywd  - working directory. This needs to be the path to 
#         the directory "vim-mc-rf/simulation".

# Output:

# A list of length four, where each element contains the variable importance values
# for a specific variable importance measure

evaluatesetting <- function(iter, mywd) {
  
  setwd(mywd)
  
  # Obtain information for the iter-th setting:
  
  n <- scenariogrid$n[iter]
  K <- scenariogrid$K[iter]
  it <- scenariogrid$itind[iter]
  settingid <- scenariogrid$settingid[iter]
  seed <- scenariogrid$seed[iter]
  
  # Simulate dataset:
  
  set.seed(seed)
  dataset <- simdata(n=n, K=K)
  
  
  # Apply the methods:
  
  set.seed(seed)
  perm <- ranger::ranger(dependent.variable.name = "cl", data=dataset, importance="permutation", 
                         num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
                         min.node.size=5, num.threads=1)$variable.importance
  
  set.seed(seed)
  gini_corr <- ranger::ranger(dependent.variable.name = "cl", data=dataset, importance="impurity_corrected",
                              num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                              probability=TRUE, min.node.size=5, num.threads=1)$variable.importance
  
  set.seed(seed)
  muobj <- diversityForest::multifor(dependent.variable.name = "cl", data=dataset, importance="both", 
                                      npervar = 5000, num.trees=5, replace = FALSE, sample.fraction = 0.7,
                                      probability=TRUE, min.node.size=5, num.threads=1)
  class_foc <- muobj$class_foc_vim
  discr <- muobj$discr_vim
  rm(muobj); gc()
  
  
  # Save the results in a list:
  
  res <- list(perm=perm,
              gini_corr=gini_corr,
              class_foc=class_foc,
              discr=discr)


  # Save results:

  # save(res, file=paste0("./intermediate_results/res_simulation_study_", iter, ".Rda"))

  
  # Return results:
  
  return(res)
  
}







# Function for simulating a dataset:

# Input parameters:

# n - sample size
# K - number of outcome classes, one of: 4, 6, 10

# Output:

# Simulated dataset

simdata <- function(n=500, K=6) {
  
  if (K==4)
    data <- simdataK4(n=n)
  
  if (K==6)
    data <- simdataK6(n=n)
  
  if (K==10)
    data <- simdataK10(n=n)
  
  return(data)
  
}






# Function for simulating a dataset with 4 outcome classes:

# Input parameters:

# n - sample size

# Output:

# Simulated dataset

simdataK4 <- function(n=500) {
  
  K <- 4
  
  # Construct the outcome:
  
  ncls <- rep(floor(n/K), K)
  sum_ncls <- sum(ncls)
  if (n > sum_ncls)
    ncls[1:(n - sum_ncls)] <- ncls[1:(n - sum_ncls)] + 1
  
  cl <- factor(rep(1:K, times=sample(ncls)))
  
  
  # Construct the means of the variables in the
  # different outcome classes:
  
  means_bin <- c(rep(0, K/2), rep(1.5, K/2))
  
  tempval <- 1
  means_1_cl <- c(rep(0, 3), tempval)
  
  means_2_cl <- c(rep(0, 2), tempval, 2*tempval)
  
  tempval <- 0.75
  means_3_cl <- c(0, tempval, 2*tempval, 3*tempval)
  
  
  
  # Simulate the variable values:
  
  X_noise <- matrix(nrow=n, ncol=50, data=rnorm(n*50))
  
  nblock <- 3
  
  X_bin <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_bin[cl==i,j] <- rnorm(sum(cl==i), mean=means_bin[i], sd=1)
    }
  }
  
  
  X_1_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_1_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_1_cl[i], sd=1)
    }
  }
  
  
  X_2_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_2_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_2_cl[i], sd=1)
    }
  }
  
  
  X_3_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_3_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_3_cl[i], sd=1)
    }
  }
  
  
  data <- data.frame(cbind(X_noise, X_bin, X_1_cl, X_2_cl, X_3_cl))
  data$cl <- cl
  
  return(data)
  
  
}







# Function for simulating a dataset with 6 outcome classes:

# Input parameters:

# n - sample size

# Output:

# Simulated dataset

simdataK6 <- function(n=500) {
  
  K <- 6
  
  ncls <- rep(floor(n/K), K)
  sum_ncls <- sum(ncls)
  if (n > sum_ncls)
    ncls[1:(n - sum_ncls)] <- ncls[1:(n - sum_ncls)] + 1
  
  cl <- factor(rep(1:K, times=sample(ncls))) 
  
  
  means_bin <- c(rep(0, K/2), rep(1.5, K/2))
  
  tempval <- 1
  means_interm <- c(rep(0, 2), rep(tempval, 2), rep(2*tempval, 2))
  
  tempval <- 1
  means_1_cl <- c(rep(0, 5), tempval)
  
  means_2_cl <- c(rep(0, 4), tempval, 2*tempval)
  
  tempval <- 0.75
  means_3_cl <- c(rep(0, 3), tempval, 2*tempval, 3*tempval)
  
  
  X_noise <- matrix(nrow=n, ncol=50, data=rnorm(n*50))
  
  nblock <- 3
  
  
  X_bin <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_bin[cl==i,j] <- rnorm(sum(cl==i), mean=means_bin[i], sd=1)
    }
  }
  
  
  X_interm <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_interm[cl==i,j] <- rnorm(sum(cl==i), mean=means_interm[i], sd=1)
    }
  }
  
  
  X_1_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_1_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_1_cl[i], sd=1)
    }
  }
  
  
  X_2_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_2_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_2_cl[i], sd=1)
    }
  }
  
  
  X_3_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_3_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_3_cl[i], sd=1)
    }
  }
  
  
  data <- data.frame(cbind(X_noise, X_bin, X_interm, X_1_cl, X_2_cl, X_3_cl))
  data$cl <- cl
  
  return(data)
  
}








# Function for simulating a dataset with 10 outcome classes:

# Input parameters:

# n - sample size

# Output:

# Simulated dataset

simdataK10 <- function(n=500) {
  
  K <- 10
  
  ncls <- rep(floor(n/K), K)
  sum_ncls <- sum(ncls)
  if (n > sum_ncls)
    ncls[1:(n - sum_ncls)] <- ncls[1:(n - sum_ncls)] + 1
  
  cl <- factor(rep(1:K, times=sample(ncls))) 
  
  
  means_bin <- c(rep(0, K/2), rep(1.5, K/2))
  
  tempval <- 1
  means_interm <- c(rep(0, 4), rep(tempval, 3), rep(2*tempval, 3))
  
  
  tempval <- 1
  means_1_cl <- c(rep(0, 9), tempval)
  
  means_2_cl <- c(rep(0, 8), tempval, 2*tempval)
  
  
  tempval <- 0.75
  means_3_cl <- c(rep(0, 6), rep(tempval, 2), 2*tempval, 3*tempval)
  
  means_4_cl <- c(rep(0, 4), rep(tempval, 2), rep(2*tempval, 2), 3*tempval, 4*tempval)
  
  
  X_noise <- matrix(nrow=n, ncol=50, data=rnorm(n*50))
  
  nblock <- 3
  
  
  
  
  X_bin <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_bin[cl==i,j] <- rnorm(sum(cl==i), mean=means_bin[i], sd=1)
    }
  }
  
  
  X_interm <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_interm[cl==i,j] <- rnorm(sum(cl==i), mean=means_interm[i], sd=1)
    }
  }
  
  
  X_1_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_1_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_1_cl[i], sd=1)
    }
  }
  
  
  X_2_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_2_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_2_cl[i], sd=1)
    }
  }
  
  
  X_3_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_3_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_3_cl[i], sd=1)
    }
  }
  
  
  X_4_cl <- matrix(nrow=n, ncol=nblock)
  
  for(j in 1:nblock) {
    for(i in 1:K) {
      X_4_cl[cl==i,j] <- rnorm(sum(cl==i), mean=means_4_cl[i], sd=1)
    }
  }
  
  data <- data.frame(cbind(X_noise, X_bin, X_interm, X_1_cl, X_2_cl, X_3_cl, X_4_cl))
  data$cl <- cl
  
  return(data)
  
}
