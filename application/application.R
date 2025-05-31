####################################################################################

# NOTE: Before the code can be executed, the R working directory *MUST* 
# be set to the directory 'vim-mc-rf/application' (that is
# the directory in which this R script is contained):

# Remove '#' from the line below and replace 'here/is/my/path/' by the path
# to the directory 'vim-mc-rf/application':

# setwd("here/is/my/path/")

####################################################################################


# Analysis of the hars dataset:
###############################
###############################


# Load data:

library("diversityForest")
data(hars)
levels(hars$Activity) <- c("Laying", "Sitting", "Standing", "Walking", "Walking downstairs", "Walking upstairs")



# Compute the VIMs:

seed <- 1234

library("ranger")

set.seed(seed)
perm <- ranger(dependent.variable.name = "Activity", data=hars, importance="permutation", 
               num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
               min.node.size=5)$variable.importance

set.seed(seed)
gini_corr <- ranger(dependent.variable.name = "Activity", data=hars, importance="impurity_corrected",
                    num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, min.node.size=5)$variable.importance

set.seed(seed)
muobj <- multifor(dependent.variable.name = "Activity", data=hars, importance="both", 
                  npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                  probability=TRUE, min.node.size=5)
class_foc <- muobj$class_foc_vim
discr <- muobj$discr_vim
rm(muobj); gc()





# Figure 6: VIM values for all covariates in the hars dataset.
###############################################################

library("ggplot2")

# Assuming perm, gini_corr, discr, and class_foc are numeric vectors
df_perm <- data.frame(x = seq_along(perm), y = perm)
df_gini_corr <- data.frame(x = seq_along(gini_corr), y = gini_corr)
df_discr <- data.frame(x = seq_along(discr), y = discr)
df_class_foc <- data.frame(x = seq_along(class_foc), y = class_foc)

titles <- c("Permutation VIM", "Corrected Gini importance", "Discriminatory VIM", "Class-focused VIM")

# Assuming perm, gini_corr, discr, and class_foc are numeric vectors
df_perm <- data.frame(x = seq_along(perm), y = perm, method = titles[1])
df_gini_corr <- data.frame(x = seq_along(gini_corr), y = gini_corr, method = titles[2])
df_discr <- data.frame(x = seq_along(discr), y = discr, method = titles[3])
df_class_foc <- data.frame(x = seq_along(class_foc), y = class_foc, method = titles[4])

# Combine all into a single dataframe
df_combined <- rbind(df_perm, df_gini_corr, df_discr, df_class_foc)

df_combined$method <- factor(df_combined$method, levels=titles)

# Facet plot
p <- ggplot(df_combined, aes(x = x, y = y)) +
  geom_point(size=1) +
  facet_wrap(~method, scales = "free_y") +  # Separate y-scales per method
  xlab("Variable index") +
  ylab("VIM values") +
  theme_bw() + theme(axis.text.x = element_text(color="black"))
p

# Figure 6:

ggsave("../figures/Fig6.eps", width=10*0.8, height=7*0.8)




# --> There is obviously a block of covariates starting with the 29th covariates,
# the VIMs of which have a very similar pattern, except for the class-focused VIM.
# They all have similar values above zero (about 0.01 in the case of the permuation VIM).

# --> Investigate these covariates further.


# Identify the atypical covariates:

plot(perm)
abline(v=28.5)
upper_limit <- 0.0115
lower_limit <- 0.008
abline(h=upper_limit)
abline(h=lower_limit)


# Select the atypical covariates and calculate the correlations between them:

X_atypical_block <- hars[,names(perm)[29:length(perm)][perm[29:length(perm)] > lower_limit & perm[29:length(perm)] < upper_limit]]
cor_atypical_block <- cor(X_atypical_block)
cors_atypical_block <- cor_atypical_block[upper.tri(cor_atypical_block)]


# For comparison, calculate the correlations between the first 28 covariates,
# which do not show such an atypical pattern in the VIMs:

X_regular_block <- hars[,names(perm)[1:28]]
cor_regular_block <- cor(X_regular_block)
cors_regular_block <- cor_regular_block[upper.tri(cor_regular_block)]


# Compare the correlations:

boxplot(cors_atypical_block, cors_regular_block)

mean(cors_atypical_block)
length(cors_atypical_block)

quantile(cors_atypical_block)

# --> Obviously the correlations between the covariates in the atypical
# block are strongly correlated, which is not at all the case for the 
# covariates in the regular block.



# Investigate the influence of the covariates in the atypical block on the outcome_

# plotMcl(data = hars, yvarname = "Activity", varnames = names(X_atypical_block))

# --> All of these covariates strongly separate the two groups of classes 
# laying, sitting, standing and walking, walking downstairs, walking upstairs.





# Load slightly modified plot functions, which were modified to
# produce figures which can also be printed as black and white figures:

source("plot_funs.R")






# Figure 7: Activity-specific distributions of the five covariates ranked highest by the
# class-focused VIM in the hars dataset.
#########################################################################################

tempnames <- names(class_foc)[order(class_foc, decreasing=TRUE)][1:5]

yvarname <- "Activity"
ps <- list()

i <- 1
plotobj <- myplotVar(hars[,tempnames[i]], y=hars[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.73, 0.63), legend.title = element_text(size = 16), legend.text = element_text(size = 12), axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

library("patchwork")
library("purrr")
library("gridExtra")
library("grid")

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- myplotVar(hars[,tempnames[i]], y=hars[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}


combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1) #+ plot_layout(guides = "collect")


# Display the combined plot
combined_plot

# Figure 7

ggsave("../figures/Fig7.eps", plot = combined_plot, width = 10*1.3, height = 14*1.3)






# Figure 8: Activity-specific distributions of the five covariates ranked highest by the
# permutation VIM in the hars dataset.
#########################################################################################

tempnames <- names(perm)[order(perm, decreasing=TRUE)][1:5]

yvarname <- "Activity"
ps <- list()

i <- 1
plotobj <- myplotVar(hars[,tempnames[i]], y=hars[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.73, 0.63), legend.title = element_text(size = 16), legend.text = element_text(size = 12), axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- myplotVar(hars[,tempnames[i]], y=hars[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}



combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1) #+ plot_layout(guides = "collect")


# Display the combined plot
combined_plot

# Figure 8

ggsave("../figures/Fig8.eps", plot = combined_plot, width = 10*1.3, height = 14*1.3)













# Analysis of the ctg dataset:
###############################
###############################

# Load the dataset:

data(ctg)




# Compute the VIMs:

seed <- 1234

set.seed(seed)
perm <- ranger(dependent.variable.name = "CLASS", data=ctg, importance="permutation", 
               num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE, 
               min.node.size=5)$variable.importance

set.seed(seed)
gini_corr <- ranger(dependent.variable.name = "CLASS", data=ctg, importance="impurity_corrected",
                    num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, min.node.size=5)$variable.importance

set.seed(seed)
muobj <- multifor(dependent.variable.name = "CLASS", data=ctg, importance="both", 
                  npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                  probability=TRUE, min.node.size=5)
class_foc <- muobj$class_foc_vim
discr <- muobj$discr_vim
rm(muobj); gc()





# Fig. S5: VIM values for all covariates in the ctg dataset.
#############################################################

# Assuming perm, gini_corr, discr, and class_foc are numeric vectors
df_perm <- data.frame(x = seq_along(perm), y = perm)
df_gini_corr <- data.frame(x = seq_along(gini_corr), y = gini_corr)
df_discr <- data.frame(x = seq_along(discr), y = discr)
df_class_foc <- data.frame(x = seq_along(class_foc), y = class_foc)

titles <- c("Permutation VIM", "Corrected Gini importance", "Discriminatory VIM", "Class-focused VIM")

# Assuming perm, gini_corr, discr, and class_foc are numeric vectors
df_perm <- data.frame(x = seq_along(perm), y = perm, method = titles[1])
df_gini_corr <- data.frame(x = seq_along(gini_corr), y = gini_corr, method = titles[2])
df_discr <- data.frame(x = seq_along(discr), y = discr, method = titles[3])
df_class_foc <- data.frame(x = seq_along(class_foc), y = class_foc, method = titles[4])

# Combine all into a single dataframe
df_combined <- rbind(df_perm, df_gini_corr, df_discr, df_class_foc)

df_combined$method <- factor(df_combined$method, levels=titles)

# Facet plot
p <- ggplot(df_combined, aes(x = x, y = y)) +
  geom_point(size=1) +
  facet_wrap(~method, scales = "free_y") +  # Separate y-scales per method
  xlab("Variable index") +
  ylab("VIM values") +
  theme_bw() + theme(axis.text.x = element_text(color="black"))
p

# Figure S5:

ggsave("../figures/FigS5.pdf", width=10*0.8, height=7*0.8)








# Table S9: Frequencies of heart rate patterns.
###############################################

ctg_2 <- ctg
ctg_2

names(ctg_2)[names(ctg_2)=="CLASS"] <- "Heart rate pattern"


# Create the table
tab <- table(ctg_2$`Heart rate pattern`)

# Extract names and values
patterns <- names(tab)
counts <- as.integer(tab)

# Collapse into LaTeX lines
latex_line1 <- paste(patterns, collapse = " & ")
latex_line2 <- paste(counts, collapse = " & ")

# Add LaTeX line endings
latex_table <- paste0(latex_line1, " \\\\\n\\hline\n", latex_line2, " \\\\")

# Optional: add surrounding LaTeX table environment
latex_full <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Frequencies of heart rate patterns.}\n",
  "\\label{tab:heart_rate_pattern}\n",
  "\\begin{tabular}{", paste(rep("c", length(tab)), collapse = ""), "}\n",
  "\\hline\n",
  latex_table, "\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\end{table}"
)

# Table S9:

# Save to .tex file
writeLines(latex_full, "../tables/TabS9.tex")







# Fig. S6: Heart rate pattern-specific distributions of the five covariates ranked highest by the
# class-focused VIM in the ctg dataset.
##################################################################################################

tempnames <- names(class_foc)[order(class_foc, decreasing=TRUE)][1:5]

yvarname <- "Heart rate pattern"
ps <- list()

i <- 1
plotobj <- plotVar(ctg_2[,tempnames[i]], y=ctg_2[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.73, 0.58), legend.title = element_text(size = 16), legend.text = element_text(size = 12), axis.title = element_text(size=16), axis.text = element_text(size=12), legend.key.height = unit(0.3, "lines")) + ylab("(Scaled) density") 
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- plotVar(ctg_2[,tempnames[i]], y=ctg_2[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}



combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1) #+ plot_layout(guides = "collect")


# Display the combined plot
combined_plot

# Figure S6:

ggsave("../figures/FigS6.pdf", plot = combined_plot, width = 10*1.3, height = 14*1.3)







# Fig. S7: Heart rate pattern-specific distributions of the five covariates ranked highest by the
# permutation VIM in the ctg dataset.
##################################################################################################

tempnames <- names(perm)[order(perm, decreasing=TRUE)][1:5]

yvarname <- "Heart rate pattern"
ps <- list()

i <- 1
plotobj <- plotVar(ctg_2[,tempnames[i]], y=ctg_2[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.73, 0.58), legend.title = element_text(size = 16), legend.text = element_text(size = 12), axis.title = element_text(size=16), axis.text = element_text(size=12), legend.key.height = unit(0.3, "lines")) + ylab("(Scaled) density")
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- plotVar(ctg_2[,tempnames[i]], y=ctg_2[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density") 
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}



combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1) #+ plot_layout(guides = "collect")


# Display the combined plot
combined_plot

# Figure S7:

ggsave("../figures/FigS7.pdf", plot = combined_plot, width = 10*1.3, height = 14*1.3)














# Analysis of the gas-drift dataset:
####################################
####################################


# Load dataset and pre-process it:

# library("OpenML")

# datacompl <- getOMLDataSet(data.id = 1476)

# gasdrift <- datacompl$data

# save(gasdrift, file="gas-drift.Rda")

load("gas-drift.Rda")

names(gasdrift)[names(gasdrift)=="Class"] <- "Gas"

gasdrift$Gas <- as.character(gasdrift$Gas)

gasdrift$Gas[gasdrift$Gas=="1"] <- "Ethanol"
gasdrift$Gas[gasdrift$Gas=="2"] <- "Ethylene"
gasdrift$Gas[gasdrift$Gas=="3"] <- "Ammonia"
gasdrift$Gas[gasdrift$Gas=="4"] <- "Acetaldehyde"
gasdrift$Gas[gasdrift$Gas=="5"] <- "Acetone"
gasdrift$Gas[gasdrift$Gas=="6"] <- "Toluene"

gasdrift$Gas <- factor(gasdrift$Gas, levels=c("Ethanol", "Ethylene", "Ammonia", "Acetaldehyde", "Acetone", "Toluene"))




# Table S10: Frequencies of gases.
###################################

# Create the table
tab <- table(gasdrift$Gas)

# Extract names and values
gases <- names(tab)
counts <- as.integer(tab)

# Collapse into LaTeX lines
latex_line1 <- paste(gases, collapse = " & ")
latex_line2 <- paste(counts, collapse = " & ")

# Add LaTeX line endings
latex_table <- paste0(latex_line1, " \\\\\n\\hline\n", latex_line2, " \\\\")

# Optional: add surrounding LaTeX table environment
latex_full <- paste0(
  "\\begin{table}[ht]\n",
  "\\centering\n",
  "\\caption{Frequencies of gases.}\n",
  "\\label{tab:gases}\n",
  "\\begin{tabular}{", paste(rep("c", length(tab)), collapse = ""), "}\n",
  "\\hline\n",
  latex_table, "\n",
  "\\hline\n",
  "\\end{tabular}\n",
  "\\end{table}"
)

# Table S10:

# Save to .tex file
writeLines(latex_full, "../tables/TabS10.tex")






# Compute the VIMs:

seed <- 1234

set.seed(seed)
perm <- ranger(dependent.variable.name = "Gas", data=gasdrift, importance="permutation",
               num.trees=5000, replace = FALSE, sample.fraction = 0.7, probability=TRUE,
               min.node.size=5)$variable.importance

set.seed(seed)
gini_corr <- ranger(dependent.variable.name = "Gas", data=gasdrift, importance="impurity_corrected",
                    num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                    probability=TRUE, min.node.size=5)$variable.importance

set.seed(seed)

# The following can take a while...
muobj <- multifor(dependent.variable.name = "Gas", data=gasdrift, importance="both",
                  npervar = 5, num.trees=5000, replace = FALSE, sample.fraction = 0.7,
                  probability=TRUE, min.node.size=5)
class_foc <- muobj$class_foc_vim
discr <- muobj$discr_vim
rm(muobj); gc()





# Fig. S8: VIM values for all covariates in the gas-drift dataset.
###################################################################


# Assuming perm, gini_corr, discr, and class_foc are numeric vectors
df_perm <- data.frame(x = seq_along(perm), y = perm)
df_gini_corr <- data.frame(x = seq_along(gini_corr), y = gini_corr)
df_discr <- data.frame(x = seq_along(discr), y = discr)
df_class_foc <- data.frame(x = seq_along(class_foc), y = class_foc)

titles <- c("Permutation VIM", "Corrected Gini importance", "Discriminatory VIM", "Class-focused VIM")

# Assuming perm, gini_corr, discr, and class_foc are numeric vectors
df_perm <- data.frame(x = seq_along(perm), y = perm, method = titles[1])
df_gini_corr <- data.frame(x = seq_along(gini_corr), y = gini_corr, method = titles[2])
df_discr <- data.frame(x = seq_along(discr), y = discr, method = titles[3])
df_class_foc <- data.frame(x = seq_along(class_foc), y = class_foc, method = titles[4])

# Combine all into a single dataframe
df_combined <- rbind(df_perm, df_gini_corr, df_discr, df_class_foc)

df_combined$method <- factor(df_combined$method, levels=titles)

# Facet plot
p <- ggplot(df_combined, aes(x = x, y = y)) +
  geom_point(size=1) +
  facet_wrap(~method, scales = "free_y") +  # Separate y-scales per method
  xlab("Variable index") +
  ylab("VIM values") +
  theme_bw() + theme(axis.text.x = element_text(color="black"))
p

# Figure S8:

ggsave("../figures/FigS8.pdf", width=10*0.8, height=7*0.8)






# Fig. S9: Gas-specific distributions of the five covariates ranked highest by the class-focused VIM
# in the gas-drift dataset.
#####################################################################################################


tempnames <- names(class_foc)[order(class_foc, decreasing=TRUE)][1:5]

yvarname <- "Gas"
ps <- list()

i <- 1
plotobj <- plotVar(gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.73, 0.63), legend.title = element_text(size = 16), legend.text = element_text(size = 12), axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density") 
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- plotVar(gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density") 
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}



combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1) #+ plot_layout(guides = "collect")


# Display the combined plot
combined_plot

# Figure S9:

ggsave("../figures/FigS9.pdf", plot = combined_plot, width = 10*1.3, height = 14*1.3)






# Fig. S10: Gas-specific distributions of the five covariates ranked highest by the permutation VIM
# in the gas-drift dataset.
####################################################################################################


tempnames <- names(perm)[order(perm, decreasing=TRUE)][1:5]

yvarname <- "Gas"
ps <- list()

i <- 1
plotobj <- plotVar(gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.73, 0.63), legend.title = element_text(size = 16), legend.text = element_text(size = 12), axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density") 
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- plotVar(gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=tempnames[i], y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density") 
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}



combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1) #+ plot_layout(guides = "collect")


# Display the combined plot
combined_plot

# Figure S10:

ggsave("../figures/FigS10.pdf", plot = combined_plot, width = 10*1.3, height = 14*1.3)






# This function is a version of plotVar from the diversityForest package that allows to plot 
# the covariate distributions after logarithmization.

# Input parameters:

# x - Metric variable or ordered categorical variable that has at least as many unique values as y
# y - Factor variable with at least three categories.

# Output:

# A list returned invisibly containing:

# - Only the element dens_pl if plot_type = "density" in plotVar;

# - Only the element boxplot_pl if plot_type = "boxplot" in plotVar;

# - The elements dens_pl, boxplot_pl, and combined_pl if plot_type = "both" in plotVar.

# All returned plots are ggplot2 objects, with combined_pl being a patchwork object.

plotVartemp <- function(x, y, ...) {
  
  if (mean(x > 0) > 0.5) {
    
    cat(paste0("Proportion excluded: ", mean(mean(x <= 0))), "\n")
    
    bool_sel <- x > 0
    x <- x[bool_sel]
    y <- y[bool_sel]
    
    x <- log(x)
    
    return(plotVar(x=x, y=y, ...))
  }
  
  if (mean(x > 0) < 0.5) {
    
    cat(paste0("Proportion excluded: ", mean(mean(x >= 0))), "\n")
    
    bool_sel <- x < 0
    x <- x[bool_sel]
    y <- y[bool_sel]
    
    x <- -log(-x)
    
    return(plotVar(x=x, y=y, ...))
  }
  
}






# Fig. S11: Gas-specific distributions of the five logarithmized covariates ranked highest by the
# class-focused VIM in the gas-drift dataset.
#################################################################################################

tempnames <- names(class_foc)[order(class_foc, decreasing=TRUE)][1:5]

yvarname <- "Gas"
ps <- list()

i <- 1
plotobj <- plotVartemp(x=gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=ifelse(mean(gasdrift[,tempnames[i]] > 0) > 0.5, paste0("log(", tempnames[i], ")"), paste0("-log(-", tempnames[i], ")")), y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.13, 0.63), legend.title = element_text(size = 15), legend.text = element_text(size = 10), axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- plotVartemp(x=gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=ifelse(mean(gasdrift[,tempnames[i]] > 0) > 0.5, paste0("log(", tempnames[i], ")"), paste0("-log(-", tempnames[i], ")")), y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by class-focused VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}



combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1)


# Display the combined plot
combined_plot

# Figure S11:

ggsave("../figures/FigS11.pdf", plot = combined_plot, width = 10*1.3, height = 14*1.3)








# Fig. S12: Gas-specific distributions of the five logarithmized covariates ranked highest by the
# permutation VIM in the gas-drift dataset.
##################################################################################################

tempnames <- names(perm)[order(perm, decreasing=TRUE)][1:5]

yvarname <- "Gas"
ps <- list()

i <- 1
plotobj <- plotVartemp(gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=ifelse(mean(gasdrift[,tempnames[i]] > 0) > 0.5, paste0("log(", tempnames[i], ")"), paste0("-log(-", tempnames[i], ")")), y_label=yvarname, plotit=FALSE)
dens_pl <- plotobj$dens_pl + theme(legend.position = "inside", legend.position.inside = c(0.17, 0.63), legend.title = element_text(size = 16), legend.text = element_text(size = 12), axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
boxplot_pl

# Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):

ps <- list()

# Create a title with increased font size
title_grob <- textGrob(
  paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
  gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
)

# Arrange plots with a larger title
combined_pl <- grid.arrange(
  dens_pl, boxplot_pl, 
  top = title_grob,  # Use custom title
  nrow = 1
)

ps[[i]] <- combined_pl

for(i in 2:5) {
  
  plotobj <- plotVartemp(gasdrift[,tempnames[i]], y=gasdrift[,yvarname], x_label=ifelse(mean(gasdrift[,tempnames[i]] > 0) > 0.5, paste0("log(", tempnames[i], ")"), paste0("-log(-", tempnames[i], ")")), y_label=yvarname, plotit=FALSE)
  dens_pl <- plotobj$dens_pl + theme(legend.position = "none", axis.title = element_text(size=16), axis.text = element_text(size=12)) + ylab("(Scaled) density")
  boxplot_pl <- plotobj$boxplot_pl + theme(axis.text.x = element_text(color = "transparent"), axis.ticks.x = element_line(color = "transparent"), axis.title = element_text(size=16), axis.text = element_text(size=12))
  
  # Arrange them using grid.arrange()
  # Create a title with increased font size
  title_grob <- textGrob(
    paste0("Ranked #", i, " by permutation VIM:  ", tempnames[i]), 
    gp = gpar(fontsize = 18)  # Increase size & make it bold ### , fontface = "bold"
  )
  
  # Arrange plots with a larger title
  combined_pl <- grid.arrange(
    dens_pl, boxplot_pl, 
    top = title_grob,  # Use custom title
    nrow = 1
  )
  
  ps[[i]] <- combined_pl
}


combined_plot <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ps[[5]], ncol=1) #+ plot_layout(guides = "collect")


# Display the combined plot
combined_plot

# Figure S12:

ggsave("../figures/FigS12.pdf", plot = combined_plot, width = 10*1.3, height = 14*1.3)




# It seems like there are two clusters of informative covariates:
# The first one, in the region of the first covariates, seems to contain
# mostly class-relevant covariates.
# In contrast, the second one, in a region in the middle of the
# covariates, seems to contain mostly class-group-differentiating covariates.

# --> Investigate whether this is actually true.


# Determine the two regions:

par(mfrow=c(1,2))
plot(class_foc)
abline(v=c(8.5, 16.5))
abline(v=c(64.5, 80.5))
abline(v=c(99.5, 125.5))
plot(perm)
abline(v=c(8.5, 16.5))
abline(v=c(64.5, 80.5))
abline(v=c(99.5, 125.5))
par(mfrow=c(1,1))


# Plot the class-specific ditributions of the covariates in the first region:

class_foc2 <- class_foc[9:16]
tempnames <- names(class_foc2)[order(class_foc2, decreasing=TRUE)][1:5]
tempnames

for(i in seq(along=tempnames)) {
  Sys.sleep(1)
  plotVartemp(x=gasdrift[,tempnames[i]], y=gasdrift$Gas)
  Sys.sleep(1)
}

# --> Yes, these covariates seem to be rather class-relevant for class "acetone".



# There also a third cluster of covariates:

class_foc2 <- class_foc[100:125]
tempnames <- names(class_foc2)[order(class_foc2, decreasing=TRUE)][1:5]
tempnames

for(i in seq(along=tempnames)) {
  Sys.sleep(1)
  plotVartemp(x=gasdrift[,tempnames[i]], y=gasdrift$Gas)
  Sys.sleep(1)
}

# --> Yes, these covariates seem to be rather class-relevant for class "ethanol".






# Plot the class-specific ditributions of the covariates in the second region:

class_foc2 <- class_foc[65:80]
tempnames <- names(class_foc2)[order(class_foc2, decreasing=TRUE)][1:5]
tempnames


for(i in seq(along=tempnames)) {
  Sys.sleep(1)
  plotVartemp(x=gasdrift[,tempnames[i]], y=gasdrift$Gas)
  Sys.sleep(1)
}

# --> Yes, these feature seem to be rather class-group-differentiating, where they
# discriminate the classes ethylene and ammonia from the other classes.
