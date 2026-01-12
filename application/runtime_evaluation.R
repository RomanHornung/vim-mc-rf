####################################################################################

# NOTE: Before the code can be executed, the R working directory *MUST* 
# be set to the directory 'vim-mc-rf/application' (that is
# the directory in which this R script is contained):

# Remove '#' from the line below and replace 'here/is/my/path/' by the path
# to the directory 'vim-mc-rf/application':

# setwd("here/is/my/path/")

####################################################################################


# Vectors that will contain the runtimes:

perm_times_appl <- gini_times_appl <- muobj_times_appl <- 0



# hars dataset:
###############

# Compute the VIMs:

library("diversityForest")
data(hars)

library("ranger")

seed <- 1234


# Compute the VIMs:

set.seed(seed)
perm_times_appl[1] <- system.time({
  perm <- ranger(dependent.variable.name = "Activity", data=hars, importance="permutation", 
                 num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
                 min.node.size=5)$variable.importance
})["elapsed"]


set.seed(seed)
gini_times_appl[1] <- system.time({
  gini_corr <- ranger(dependent.variable.name = "Activity", data=hars, importance="impurity_corrected",
                      num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                      probability=TRUE, min.node.size=5)$variable.importance
})["elapsed"]

set.seed(seed)
muobj_times_appl[1] <- system.time({
  muobj <- multifor(dependent.variable.name = "Activity", data=hars, importance="both", 
                    npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, min.node.size=5)
})["elapsed"]





#  ctg dataset:
###############

# Load the dataset:

data(ctg)


# Compute the VIMs:

set.seed(seed)
perm_times_appl[2] <- system.time({
perm <- ranger(dependent.variable.name = "CLASS", data=ctg, importance="permutation", 
               num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
               min.node.size=5)$variable.importance
})["elapsed"]


set.seed(seed)
gini_times_appl[2] <- system.time({
gini_corr <- ranger(dependent.variable.name = "CLASS", data=ctg, importance="impurity_corrected",
                    num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, min.node.size=5)$variable.importance
})["elapsed"]


set.seed(seed)
muobj_times_appl[2] <- system.time({
muobj <- multifor(dependent.variable.name = "CLASS", data=ctg, importance="both", 
                  npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                  probability=TRUE, min.node.size=5)
})["elapsed"]





# gas-drift dataset:
######################

# Load dataset:

load("gas-drift.Rda")


# Compute the VIMs:

set.seed(seed)
perm_times_appl[3] <- system.time({
perm <- ranger(dependent.variable.name = "Class", data=gasdrift, importance="permutation",
               num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE,
               min.node.size=5)$variable.importance
})["elapsed"]


set.seed(seed)
gini_times_appl[3] <- system.time({
gini_corr <- ranger(dependent.variable.name = "Class", data=gasdrift, importance="impurity_corrected",
                    num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, min.node.size=5)$variable.importance
})["elapsed"]


set.seed(seed)
muobj_times_appl[3] <- system.time({
muobj <- multifor(dependent.variable.name = "Class", data=gasdrift, importance="both",
                  npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                  probability=TRUE, min.node.size=5)
})["elapsed"]

class_foc_gas <- muobj$class_foc_vim
discr_gas <- muobj$discr_vim




# Reduce the maximum depth to 10 to save computation time:

set.seed(seed)
muobj_times_gas_maxdepth10 <- system.time({
  muobj <- multifor(dependent.variable.name = "Class", data=gasdrift, importance="both",
                    npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, max.depth=10)
})["elapsed"]

class_foc_gas_maxdepth10 <- muobj$class_foc_vim
discr_gas_maxdepth10 <- muobj$discr_vim




# Reduce the maximum depth further to 7 to save computation time:

set.seed(seed)
muobj_times_gas_maxdepth7 <- system.time({
  muobj <- multifor(dependent.variable.name = "Class", data=gasdrift, importance="both",
                    npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, max.depth=7)
})["elapsed"]

class_foc_gas_maxdepth7 <- muobj$class_foc_vim
discr_gas_maxdepth7 <- muobj$discr_vim



# Reduce the maximum depth further to 5 to save computation time:

set.seed(seed)
muobj_times_gas_maxdepth5 <- system.time({
  muobj <- multifor(dependent.variable.name = "Class", data=gasdrift, importance="both",
                    npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, max.depth=5)
})["elapsed"]

class_foc_gas_maxdepth5 <- muobj$class_foc_vim
discr_gas_maxdepth5 <- muobj$discr_vim




# Save the obtained runtimes:

save(perm_times_appl, gini_times_appl, muobj_times_appl,
     muobj_times_gas_maxdepth10, muobj_times_gas_maxdepth7,
     muobj_times_gas_maxdepth5,
     class_foc_gas, discr_gas, 
     class_foc_gas_maxdepth10, discr_gas_maxdepth10,
     class_foc_gas_maxdepth7, discr_gas_maxdepth7,
     class_foc_gas_maxdepth5, discr_gas_maxdepth5,
     file="./runtime.Rda")

load("./runtime.Rda")



# Load the runtimes for the simulated datasets:

load("../simulation/intermediate_results/runtime.Rda")


# Function for formatting the runtimes for the LaTeX table:

formatit <- function(x) paste(format(round(x/60, 2), nsmall=2), collapse = " & ")




# Table S11: Runtime analysis (in minutes).
###########################################

# Output file name
outfile <- "../tables/TabS11.tex"

# Table header
table_lines <- c(
  "\\begin{table}[ht]",
  "\\scriptsize",
  "\\centering",
  "\\caption{Runtime analysis (in minutes). Perm: permutation importance; Gini\\_corr: corrected Gini importance; Class-foc\\_Discr\\_def: class-focused and discriminatory VIMs with default hyperparameters; Class-foc\\_Discr\\_md10, Class-foc\\_Discr\\_md7, and Class-foc\\_Discr\\_md5: class-focused and discriminatory VIMs computed using a maximum tree depth of 10, 7, and 5, respectively}",
  "\\label{tab:runtime_analysis}",
  "\\begin{tabular}{p{1.2cm} p{1.3cm} p{1.3cm} p{1.9cm} p{1.9cm} p{1.9cm} p{1.9cm}}",
  "\\hline",
  " & Perm & Gini\\_corr & Class-foc\\_Discr\\_def & Class-foc\\_Discr\\_md10 & Class-foc\\_Discr\\_md7 & Class-foc\\_Discr\\_md5 \\\\",
  "\\hline"
)

# Table body
table_body <- c(
  paste0("Sim.\\ data $C=4$ & ",
         formatit(c(perm_times_sim[1], gini_times_sim[1], muobj_times_sim[1])),
         " & -- & -- & -- \\\\[8pt]"),
  
  paste0("Sim.\\ data $C=6$ & ",
         formatit(c(perm_times_sim[2], gini_times_sim[2], muobj_times_sim[2])),
         " & -- & -- & -- \\\\[8pt]"),
  
  paste0("Sim.\\ data $C=10$ & ",
         formatit(c(perm_times_sim[3], gini_times_sim[3], muobj_times_sim[3])),
         " & -- & -- & -- \\\\[8pt]"),
  
  paste0("hars & ",
         formatit(c(perm_times_appl[1], gini_times_appl[1], muobj_times_appl[1])),
         " & -- & -- & -- \\\\[8pt]"),
  
  paste0("ctg & ",
         formatit(c(perm_times_appl[2], gini_times_appl[2], muobj_times_appl[2])),
         " & -- & -- & -- \\\\[8pt]"),
  
  paste0("gas-drift & ",
         formatit(c(perm_times_appl[3], gini_times_appl[3], muobj_times_appl[3])),
         " & ",
         formatit(c(muobj_times_gas_maxdepth10,
                    muobj_times_gas_maxdepth7,
                    muobj_times_gas_maxdepth5)),
         " \\\\[8pt]")
)

# Table footer
table_footer <- c(
  "\\end{tabular}",
  "\\end{table}"
)

# Combine all parts:
latex_table <- c(table_lines, table_body, table_footer)


# Table S11:

writeLines(latex_table, outfile)





# Correlation between multi-class/discriminatory VIM values when using unrestricted and
# restricted tree depth:

# Function used to format the estimated calculation coefficients:
formatit2 <- function(x) format(round(x, 4), nsmall=4)

formatit2(cor(class_foc_gas, class_foc_gas_maxdepth10))
formatit2(cor(class_foc_gas, class_foc_gas_maxdepth7))
formatit2(cor(class_foc_gas, class_foc_gas_maxdepth5))

formatit2(cor(discr_gas, discr_gas_maxdepth10))
formatit2(cor(discr_gas, discr_gas_maxdepth7))
formatit2(cor(discr_gas, discr_gas_maxdepth5))
