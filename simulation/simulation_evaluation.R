####################################################################################

# NOTE: Before the code can be executed, the R working directory *MUST* 
# be set to the directory 'vim-mc-rf/simulation' (that is
# the directory in which this R script is contained):

# Remove '#' from the line below and replace 'here/is/my/path/' by the path
# to the directory 'vim-mc-rf/simulation':

# setwd("here/is/my/path/")

####################################################################################

# Load and pre-process the results:
#####################################

load("./intermediate_results/scenariogrid_simulation.Rda")
load("./intermediate_results/results_simulation.Rda")


reorderind <- order(scenariogrid$n, scenariogrid$K, scenariogrid$itind)
scenariogrid <- scenariogrid[reorderind,]
Results <- Results[reorderind]

scenariogrid$seed <- scenariogrid$settingid <- NULL


results <- scenariogrid[rep(1:nrow(scenariogrid), each=length(Results[[1]])),]
results$method_all <- rep(names(Results[[1]]), times=nrow(scenariogrid))

Results <- unlist(Results, recursive = FALSE)

resultsall <- results[rep(1:nrow(results), times=sapply(Results, length)),]
resultsall$rank <- unlist(lapply(Results, function(x) rank(-x)))
resultsall$vim <- unlist(Results)
rownames(resultsall) <- 1:nrow(resultsall)



resultsall$method_all <- factor(resultsall$method_all, levels=c("perm", "gini_corr", "discr", "class_foc"))
levels(resultsall$method_all) <- c("Perm", "Gini_corr", "Discr", "Class-foc")





# Function used for saving a table produced in R as a LaTeX table:

# Input parameters:

# combined_stats - object containing the table (has to be of a very specific
#                  format, see below)
# filename       - name of the file to which the table should be saved, can include
#                  path to that file
# table_title    - title to include in the produced LaTeX table

# Output:

# NULL

prepare_and_save_latex <- function(combined_stats, filename, table_title = NULL) {
  combined_stats$n <- as.integer(combined_stats$n) # Convert n to integer to remove decimal points
  
  # Check if a table title is provided and set the caption and caption.placement accordingly
  if (!is.null(table_title)) {
    latex_table <- xtable(combined_stats, digits=0, auto=FALSE, include.rownames=FALSE, caption = table_title)
    caption_placement <- "top"
  } else {
    latex_table <- xtable(combined_stats, digits=0, auto=FALSE, include.rownames=FALSE)
    caption_placement <- NULL  # No title, so no caption placement needed
  }
  
  # Customize the LaTeX output to add horizontal lines after each sample size
  print(latex_table, 
        type = "latex", 
        file = filename,
        include.rownames=FALSE,  # Ensure row names are not included
        add.to.row = list(pos = list(which(diff(combined_stats$n) != 0)), 
                          command = "\\hline "),
        caption.placement = caption_placement) # Conditionally place the caption based on if a title was provided
}





# Tables S1 to S8: Mean AUC values with 95% confidence intervals
################################################################


library("dplyr")
library("xtable")
library("stringr")
library("forcats")


resultstemp <- resultsall[resultsall$K==4,]
resultstemp$K <- NULL


# Provided AUC calculation function
auroc <- function(score, bool) {
  n1 <- sum(!bool)
  n2 <- sum(bool)
  U  <- sum(rank(score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

calculate_mean_ci_l <- function(auc_values) {
  mean_auc <- mean(auc_values)
  sd_auc <- sd(auc_values)
  n <- length(auc_values)
  se_auc <- sd_auc / sqrt(n)
  error_margin <- qnorm(0.975) * se_auc
  lower_ci <- mean_auc - error_margin
  return(lower_ci)
}

calculate_mean_ci_u <- function(auc_values) {
  mean_auc <- mean(auc_values)
  sd_auc <- sd(auc_values)
  n <- length(auc_values)
  se_auc <- sd_auc / sqrt(n)
  error_margin <- qnorm(0.975) * se_auc
  upper_ci <- mean_auc + error_margin
  upper_ci
}


# Compute AUC for specified groups and aggregate results
resultsK4_vsnoise <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(1:50, 51:53)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_2 = auroc(vim[c(1:50, 54:56)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_3 = auroc(vim[ c(1:50, 57:59)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(1:50, 60:62)], c(rep(FALSE, 50), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    .groups = 'drop'
  ) 

widths <- c(resultsK4_vsnoise$ci_auc_1_u - resultsK4_vsnoise$ci_auc_1_l,
            resultsK4_vsnoise$ci_auc_2_u - resultsK4_vsnoise$ci_auc_2_l,
            resultsK4_vsnoise$ci_auc_3_u - resultsK4_vsnoise$ci_auc_3_l,
            resultsK4_vsnoise$ci_auc_4_u - resultsK4_vsnoise$ci_auc_4_l)

resultsK4_vsnoise <- resultsK4_vsnoise %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4) %>%
  rename_with(
    ~ c("two_gr", "cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_drop(method_all)) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc",
  ))

# Table S1:

prepare_and_save_latex(resultsK4_vsnoise, "../tables/TabS1.tex", table_title = "Influential vs noise variables, C = 4")








# Compute AUC for specified groups and aggregate results
resultsK4_vstwo_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 54:56)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    .groups = 'drop'
  )  

widths_2 <- c(resultsK4_vstwo_gr$ci_auc_1_u - resultsK4_vstwo_gr$ci_auc_1_l,
              resultsK4_vstwo_gr$ci_auc_2_u - resultsK4_vstwo_gr$ci_auc_2_l,
              resultsK4_vstwo_gr$ci_auc_3_u - resultsK4_vstwo_gr$ci_auc_3_l)

resultsK4_vstwo_gr <- resultsK4_vstwo_gr %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc"
  ))

# Table S4:

prepare_and_save_latex(resultsK4_vstwo_gr, "../tables/TabS4.tex", table_title = "cl\\_sp vs two\\_gr, C = 4")










resultstemp <- resultsall[resultsall$K==6,]
resultstemp$K <- NULL



# Compute AUC for specified groups and aggregate results
resultsK6_vsnoise <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(1:50, 51:53)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_2 = auroc(vim[c(1:50, 54:56)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_3 = auroc(vim[ c(1:50, 57:59)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(1:50, 60:62)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_5 = auroc(vim[c(1:50, 63:65)], c(rep(FALSE, 50), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    mean_auc_5 = mean(auc_5),
    ci_auc_5_l = calculate_mean_ci_l(auc_5),
    ci_auc_5_u = calculate_mean_ci_u(auc_5),
    .groups = 'drop'
  ) 


widths <- c(widths, resultsK6_vsnoise$ci_auc_1_u - resultsK6_vsnoise$ci_auc_1_l,
            resultsK6_vsnoise$ci_auc_2_u - resultsK6_vsnoise$ci_auc_2_l,
            resultsK6_vsnoise$ci_auc_3_u - resultsK6_vsnoise$ci_auc_3_l,
            resultsK6_vsnoise$ci_auc_4_u - resultsK6_vsnoise$ci_auc_4_l,
            resultsK6_vsnoise$ci_auc_5_u - resultsK6_vsnoise$ci_auc_5_l)

resultsK6_vsnoise <- resultsK6_vsnoise %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u),
    auc_5 = sprintf("%.2f [%.2f, %.2f]", mean_auc_5, ci_auc_5_l, ci_auc_5_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4, auc_5) %>%
  rename_with(
    ~ c("two_gr", "thr_gr", "cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_drop(method_all)) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc",
  ))

# Table S2:

prepare_and_save_latex(resultsK6_vsnoise, "../tables/TabS2.tex", table_title = "Influential vs noise variables, C = 6")





# Compute AUC for specified groups and aggregate results
resultsK6_vstwo_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    .groups = 'drop'
  )  

widths_2 <- c(widths_2, resultsK6_vstwo_gr$ci_auc_1_u - resultsK6_vstwo_gr$ci_auc_1_l,
              resultsK6_vstwo_gr$ci_auc_2_u - resultsK6_vstwo_gr$ci_auc_2_l,
              resultsK6_vstwo_gr$ci_auc_3_u - resultsK6_vstwo_gr$ci_auc_3_l)

resultsK6_vstwo_gr <- resultsK6_vstwo_gr %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc"
  ))

# Table S5:

prepare_and_save_latex(resultsK6_vstwo_gr, "../tables/TabS5.tex", table_title = "cl\\_sp vs two\\_gr, C = 6")




# Compute AUC for specified groups and aggregate results
resultsK6_vsthr_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    .groups = 'drop'
  )  

widths_2 <- c(widths_2, resultsK6_vsthr_gr$ci_auc_1_u - resultsK6_vsthr_gr$ci_auc_1_l,
              resultsK6_vsthr_gr$ci_auc_2_u - resultsK6_vsthr_gr$ci_auc_2_l,
              resultsK6_vsthr_gr$ci_auc_3_u - resultsK6_vsthr_gr$ci_auc_3_l)

resultsK6_vsthr_gr <- resultsK6_vsthr_gr %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc"
  ))

# Table S6:

prepare_and_save_latex(resultsK6_vsthr_gr, "../tables/TabS6.tex", table_title = "cl\\_sp vs thr\\_gr, C = 6")







resultstemp <- resultsall[resultsall$K==10,]
resultstemp$K <- NULL



# Compute AUC for specified groups and aggregate results
resultsK10_vsnoise <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(1:50, 51:53)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_2 = auroc(vim[c(1:50, 54:56)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_3 = auroc(vim[ c(1:50, 57:59)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(1:50, 60:62)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_5 = auroc(vim[c(1:50, 63:65)], c(rep(FALSE, 50), rep(TRUE, 3))),
    auc_6 = auroc(vim[c(1:50, 66:68)], c(rep(FALSE, 50), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    mean_auc_5 = mean(auc_5),
    ci_auc_5_l = calculate_mean_ci_l(auc_5),
    ci_auc_5_u = calculate_mean_ci_u(auc_5),
    mean_auc_6 = mean(auc_6),
    ci_auc_6_l = calculate_mean_ci_l(auc_6),
    ci_auc_6_u = calculate_mean_ci_u(auc_6),
    .groups = 'drop'
  ) %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u),
    auc_5 = sprintf("%.2f [%.2f, %.2f]", mean_auc_5, ci_auc_5_l, ci_auc_5_u),
    auc_6 = sprintf("%.2f [%.2f, %.2f]", mean_auc_6, ci_auc_6_l, ci_auc_6_u)
  )  


widths <- c(widths, resultsK10_vsnoise$ci_auc_1_u - resultsK10_vsnoise$ci_auc_1_l,
            resultsK10_vsnoise$ci_auc_2_u - resultsK10_vsnoise$ci_auc_2_l,
            resultsK10_vsnoise$ci_auc_3_u - resultsK10_vsnoise$ci_auc_3_l,
            resultsK10_vsnoise$ci_auc_4_u - resultsK10_vsnoise$ci_auc_4_l,
            resultsK10_vsnoise$ci_auc_5_u - resultsK10_vsnoise$ci_auc_5_l,
            resultsK10_vsnoise$ci_auc_6_u - resultsK10_vsnoise$ci_auc_6_l)

resultsK10_vsnoise <- resultsK10_vsnoise %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4, auc_5, auc_6) %>%
  rename_with(
    ~ c("two_gr", "thr_gr", "cl_rel_1", "cl_rel_2", "cl_rel_3", "cl_rel_4"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_drop(method_all)) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc",
  ))

# Table S3:

prepare_and_save_latex(resultsK10_vsnoise, "../tables/TabS3.tex", table_title =
                         "Influential vs noise variables, C = 10")







# Compute AUC for specified groups and aggregate results
resultsK10_vstwo_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(51:53, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    .groups = 'drop'
  )  

widths_2 <- c(widths_2, resultsK10_vstwo_gr$ci_auc_1_u - resultsK10_vstwo_gr$ci_auc_1_l,
              resultsK10_vstwo_gr$ci_auc_2_u - resultsK10_vstwo_gr$ci_auc_2_l,
              resultsK10_vstwo_gr$ci_auc_3_u - resultsK10_vstwo_gr$ci_auc_3_l,
              resultsK10_vstwo_gr$ci_auc_4_u - resultsK10_vstwo_gr$ci_auc_4_l)

resultsK10_vstwo_gr <- resultsK10_vstwo_gr %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3", "cl_rel_4"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc"
  ))
  
# Table S7:

prepare_and_save_latex(resultsK10_vstwo_gr, "../tables/TabS7.tex", table_title = "cl\\_sp vs two\\_gr, C = 10")





# Compute AUC for specified groups and aggregate results
resultsK10_vsthr_gr <- resultstemp %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(54:56, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    ci_auc_1_l = calculate_mean_ci_l(auc_1),
    ci_auc_1_u = calculate_mean_ci_u(auc_1),
    mean_auc_2 = mean(auc_2),
    ci_auc_2_l = calculate_mean_ci_l(auc_2),
    ci_auc_2_u = calculate_mean_ci_u(auc_2),
    mean_auc_3 = mean(auc_3),
    ci_auc_3_l = calculate_mean_ci_l(auc_3),
    ci_auc_3_u = calculate_mean_ci_u(auc_3),
    mean_auc_4 = mean(auc_4),
    ci_auc_4_l = calculate_mean_ci_l(auc_4),
    ci_auc_4_u = calculate_mean_ci_u(auc_4),
    .groups = 'drop'
  )  

widths_2 <- c(widths_2, resultsK10_vsthr_gr$ci_auc_1_u - resultsK10_vsthr_gr$ci_auc_1_l,
              resultsK10_vsthr_gr$ci_auc_2_u - resultsK10_vsthr_gr$ci_auc_2_l,
              resultsK10_vsthr_gr$ci_auc_3_u - resultsK10_vsthr_gr$ci_auc_3_l,
              resultsK10_vsthr_gr$ci_auc_4_u - resultsK10_vsthr_gr$ci_auc_4_l)

resultsK10_vsthr_gr <- resultsK10_vsthr_gr %>%
  mutate(
    auc_1 = sprintf("%.2f [%.2f, %.2f]", mean_auc_1, ci_auc_1_l, ci_auc_1_u),
    auc_2 = sprintf("%.2f [%.2f, %.2f]", mean_auc_2, ci_auc_2_l, ci_auc_2_u),
    auc_3 = sprintf("%.2f [%.2f, %.2f]", mean_auc_3, ci_auc_3_l, ci_auc_3_u),
    auc_4 = sprintf("%.2f [%.2f, %.2f]", mean_auc_4, ci_auc_4_l, ci_auc_4_u)
  ) %>%
  select(n, method_all, auc_1, auc_2, auc_3, auc_4) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3", "cl_rel_4"), 
    .cols = starts_with("auc_")
  ) %>%
  mutate(method_all = fct_recode(method_all,
                                 "Perm" = "Perm",
                                 "Gini_corr" = "Gini_corr",
                                 "Discr" = "Discr",
                                 "Class-foc" = "Class-foc"
  ))
  
# Table S8:

prepare_and_save_latex(resultsK10_vsthr_gr, "../tables/TabS8.tex", table_title = "cl\\_sp vs thr\\_gr, C = 10")





# Range of confidence interval widths for the comparison between the informative covariates
# versus the noise covariates:
range(widths)




# Range of confidence interval widths for the comparison between the class-related covariates
# versus the class-group-differentiating covariates:
range(widths_2)







resultsall$n <- factor(paste0("n = ", resultsall$n), levels=c("n = 100", "n = 500", "n = 1000", "n = 2000"))





# Figure 2: Mean AUC values per considered sample size and method for C = 4.
#############################################################################

library("tidyr")

resultstemp <- resultsall[resultsall$K==4,]
resultstemp$K <- NULL


# Compute AUC for specified groups and aggregate results
resultsK4_vstwo_gr <- resultstemp %>%
  filter(method_all != "Discr") %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 54:56)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK4_vstwo_gr %>%
  pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value")

results_long <- results_long %>%
  mutate(
    method_all = recode(method_all, class_foc = "Class-foc"),
    method_all = fct_relevel(method_all, "Class-foc", "Perm", "Gini_corr")
  )


library("ggplot2")

p <- ggplot(results_long, aes(x = cl_rel, y = value, group = method_all, color = method_all)) +
  geom_line() +
  geom_point(aes(shape = method_all), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=15),
        strip.text.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14)) + 
  labs(color="VIM type", shape="VIM type", y="AUC") +
  scale_shape_manual(values=c("Class-foc" = 1, "Perm" = 2, "Gini_corr" = 3)) +
  scale_x_discrete(labels = c("cl_rel_1" = expression(X[cl_rel_1]), "cl_rel_2" = 
                                expression(X[cl_rel_2]), "cl_rel_3" = expression(X[cl_rel_3]), 
                              "cl_rel_4" = expression(X[cl_rel_4]))) +
  theme(legend.position = "right") +
  facet_wrap(~ n) +  # Add facet_wrap for separate plots by "n"
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))

p

# Figure 2:

ggsave("../figures/Fig2.eps", width=10, height=6)








# Figure 3: Mean AUC values per considered sample size and method for C = 6.
#############################################################################

# Calculate for K=6 and save
resultstemp <- resultsall[resultsall$K==6,]
resultstemp$K <- NULL

# Compute AUC for specified groups and aggregate results
resultsK6_vstwo_gr <- resultstemp %>%
  filter(method_all != "Discr") %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK6_vstwo_gr %>%
  pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value")

results_long <- results_long %>%
  mutate(
    method_all = recode(method_all, class_foc = "Class-foc"),
    method_all = fct_relevel(method_all, "Class-foc", "Perm", "Gini_corr")
  )

p1 <- ggplot(results_long, aes(x = cl_rel, y = value, group = method_all, color = method_all)) +
  geom_line() +
  geom_point(aes(shape = method_all), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=17),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)
  ) + 
  ylab("AUC") +
  scale_shape_manual(values=c("Class-foc" = 1, "Perm" = 2, "Gini_corr" = 3)) +
  scale_x_discrete(labels = c("cl_rel_1" = expression(X[cl_rel_1]), "cl_rel_2" = expression(X[cl_rel_2]), 
                              "cl_rel_3" = expression(X[cl_rel_3]), "cl_rel_4" = expression(X[cl_rel_4]))) +
  theme(legend.position = "none") +
  ggtitle(expression("Comparison with" ~ X[two_gr])) + 
  facet_wrap(~ n)  # Add facet_wrap for separate plots by "n"



# Compute AUC for specified groups and aggregate results
resultsK6_vsthr_gr <- resultstemp %>%
  filter(method_all != "Discr") %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK6_vsthr_gr %>%
  pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value")

results_long <- results_long %>%
  mutate(
    method_all = recode(method_all, class_foc = "Class-foc"),
    method_all = fct_relevel(method_all, "Class-foc", "Perm", "Gini_corr")
  )

p2 <- ggplot(results_long, aes(x = cl_rel, y = value, group = method_all, color = method_all)) +
  geom_point(aes(shape = method_all), size=2.5) +
  geom_line() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=17),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16)
  ) + 
  labs(color="VIM type", shape="VIM type", y="AUC") +
  scale_shape_manual(values=c("Class-foc" = 1, "Perm" = 2, "Gini_corr" = 3)) +
  scale_x_discrete(labels = c("cl_rel_1" = expression(X[cl_rel_1]), "cl_rel_2" = 
                                expression(X[cl_rel_2]), "cl_rel_3" = expression(X[cl_rel_3]), 
                              "cl_rel_4" = expression(X[cl_rel_4]))) +
  theme(legend.position = "right") +
  ggtitle(expression("Comparison with" ~ X[thr_gr])) + 
  facet_wrap(~ n) + # Add facet_wrap for separate plots by "n"
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))


results_combined <- bind_rows(
  resultsK6_vstwo_gr %>%
    pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value"),
  resultsK6_vsthr_gr %>%
    pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value")
)

y_min <- min(results_combined$value, na.rm = TRUE)
y_max <- max(results_combined$value, na.rm = TRUE)

p1 <- p1 + coord_cartesian(ylim = c(y_min, y_max))
p2 <- p2 + coord_cartesian(ylim = c(y_min, y_max))


library("gridExtra")
library("grid") 

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, nullGrob(), p2, ncol = 3, widths = c(2, 0.1, 2.5)) # c(2, 0.1, 2.83))  # `nullGrob()` is an empty plot

# Figure 3:

ggsave("../figures/Fig3.eps", grid_plot, width=14, height=5.5)








# Figure 4: Mean AUC values per considered sample size and method for C = 10.
#############################################################################

resultstemp <- resultsall[resultsall$K==10,]
resultstemp$K <- NULL


# Compute AUC for specified groups and aggregate results
resultsK10_vstwo_gr <- resultstemp %>%
  filter(method_all != "Discr") %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(51:53, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(51:53, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(51:53, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(51:53, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    mean_auc_4 = mean(auc_4),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3, mean_auc_4) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3", "cl_rel_4"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK10_vstwo_gr %>%
  pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value")

results_long <- results_long %>%
  mutate(
    method_all = recode(method_all, class_foc = "Class-foc"),
    method_all = fct_relevel(method_all, "Class-foc", "Perm", "Gini_corr")
  )

p1 <- ggplot(results_long, aes(x = cl_rel, y = value, group = method_all, color = method_all)) +
  geom_line() +
  geom_point(aes(shape = method_all), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=16),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)
  ) +
  ylab("AUC") +
  scale_shape_manual(values=c("Class-foc" = 1, "Perm" = 2, "Gini_corr" = 3)) +
  scale_x_discrete(labels = c("cl_rel_1" = expression(X[cl_rel_1]), "cl_rel_2" = 
                                expression(X[cl_rel_2]), "cl_rel_3" = expression(X[cl_rel_3]), 
                              "cl_rel_4" = expression(X[cl_rel_4]))) +
  theme(legend.position = "none") +
  ggtitle(expression("Comparison with" ~ X[two_gr])) + 
  facet_wrap(~ n)  # Add facet_wrap for separate plots by "n"




# Compute AUC for specified groups and aggregate results
resultsK10_vsthr_gr <- resultstemp %>%
  filter(method_all != "Discr") %>%
  group_by(n, method_all, itind) %>%
  summarise(
    auc_1 = auroc(vim[c(54:56, 57:59)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_2 = auroc(vim[ c(54:56, 60:62)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_3 = auroc(vim[c(54:56, 63:65)], c(rep(FALSE, 3), rep(TRUE, 3))),
    auc_4 = auroc(vim[c(54:56, 66:68)], c(rep(FALSE, 3), rep(TRUE, 3))),
    .groups = 'drop'
  ) %>%
  group_by(n, method_all) %>%
  summarise(
    mean_auc_1 = mean(auc_1),
    mean_auc_2 = mean(auc_2),
    mean_auc_3 = mean(auc_3),
    mean_auc_4 = mean(auc_4),
    .groups = 'drop'
  ) %>%
  select(n, method_all, mean_auc_1, mean_auc_2, mean_auc_3, mean_auc_4) %>%
  rename_with(
    ~ c("cl_rel_1", "cl_rel_2", "cl_rel_3", "cl_rel_4"), 
    .cols = starts_with("mean_auc_")
  )

# Reshape data to long format
results_long <- resultsK10_vsthr_gr %>%
  pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value")

results_long <- results_long %>%
  mutate(
    method_all = recode(method_all, class_foc = "Class-foc"),
    method_all = fct_relevel(method_all, "Class-foc", "Perm", "Gini_corr")
  )

p2 <- ggplot(results_long, aes(x = cl_rel, y = value, group = method_all, color = method_all)) +
  geom_line() +
  geom_point(aes(shape = method_all), size=2.5) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(colour="black", size=16),
        strip.text.x=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=12),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16),
        plot.title = element_text(size = 16)
  ) +
  labs(color="VIM type", shape="VIM type", y="AUC") +
  scale_shape_manual(values=c("Class-foc" = 1, "Perm" = 2, "Gini_corr" = 3)) +
  scale_x_discrete(labels = c("cl_rel_1" = expression(X[cl_rel_1]), "cl_rel_2" = expression(X[cl_rel_2]), 
                              "cl_rel_3" = expression(X[cl_rel_3]), "cl_rel_4" = expression(X[cl_rel_4]))) +
  theme(legend.position = "right") +
  ggtitle(expression("Comparison with" ~ X[thr_gr])) + 
  facet_wrap(~ n) + # Add facet_wrap for separate plots by "n"
  guides(color = guide_legend(order = 1),
         shape = guide_legend(order = 1))


results_combined <- bind_rows(
  resultsK10_vstwo_gr %>%
    pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value"),
  resultsK10_vsthr_gr %>%
    pivot_longer(cols = starts_with("cl_rel_"), names_to = "cl_rel", values_to = "value")
)

y_min <- min(results_combined$value, na.rm = TRUE)
y_max <- max(results_combined$value, na.rm = TRUE)

p1 <- p1 + coord_cartesian(ylim = c(y_min, y_max))
p2 <- p2 + coord_cartesian(ylim = c(y_min, y_max))


# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, nullGrob(), p2, ncol = 3, widths = c(2, 0.1, 2.4)) # c(2, 0.1, 2.83))  # `nullGrob()` is an empty plot

# Figure 4:

ggsave("../figures/Fig4.eps", grid_plot, width=16.5, height=5.5)








# Figure 1: VIM values obtained for class-focused and discriminatory VIM and the permutation VIM (perm) 
# obtained for all simulated datasets with n = 500.
#################################################################################################


library("RColorBrewer")

# display.brewer.pal(n = 9, name = "YlGnBu")
colors <- brewer.pal(9, "YlGnBu")[c(2, 5)]

resi <- resultsall %>% filter(n=="n = 500", method_all %in% c("Discr", "Class-foc"))

resi <- resi %>% mutate(method = factor(method_all, levels = c("Discr", "Class-foc")))

resi <- resi %>%
  group_by(n, K, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p1 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(aes(fill = method), position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Class-focused and discriminatory VIM", fill="VIM type", x = "Covariates", y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("Discr", "Class-foc")) +
  theme(legend.position.inside = c(0.145, 0.94), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16),
        legend.background = element_rect(fill = NA, colour = NA),  # Transparent legend background and no border
        legend.key = element_rect(fill = NA, colour = NA),  # Transparent keys and no border
        legend.box.background = element_rect(fill = NA, colour = NA)) +
  guides(fill = guide_legend(position = "inside")) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)



resi <- resultsall %>% filter(n=="n = 500", method_all %in% c("Perm"))

resi <- resi %>%
  group_by(n, K, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), 
             levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p2 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Permutation VIM", x = "Covariates", y = "VIM values") +
  theme(legend.position = "none", axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, p2, ncol = 2)  # `nullGrob()` is an empty plot

# Figure 1:

ggsave("../figures/Fig1.eps", grid_plot, width = 12, height = 10)







# Figure S2: VIM values obtained for class-focused and discriminatory VIM and the permutation VIM (perm) 
# obtained for all simulated datasets with n = 100.
#################################################################################################

resi <- resultsall %>% filter(n=="n = 100", method_all %in% c("Discr", "Class-foc"))

resi <- resi %>% mutate(method = factor(method_all, levels = c("Discr", "Class-foc")))

resi <- resi %>%
  group_by(n, K, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p1 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(aes(fill = method), position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Class-focused and discriminatory VIM", fill="VIM type", x = "Covariates", y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("Discr", "Class-foc")) +
  theme(legend.position.inside = c(0.145, 0.94), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16),
        legend.background = element_rect(fill = NA, colour = NA),  # Transparent legend background and no border
        legend.key = element_rect(fill = NA, colour = NA),  # Transparent keys and no border
        legend.box.background = element_rect(fill = NA, colour = NA)) +
  guides(fill = guide_legend(position = "inside")) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)



resi <- resultsall %>% filter(n=="n = 500", method_all %in% c("Perm"))

resi <- resi %>%
  group_by(n, K, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), 
             levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p2 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Permutation VIM", x = "Covariates", y = "VIM values") +
  theme(legend.position = "none", axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, p2, ncol = 2)  # `nullGrob()` is an empty plot

# Figure S2:

ggsave("../figures/FigS2.pdf", grid_plot, width = 12, height = 10)








# Figure S3: VIM values obtained for class-focused and discriminatory VIM and the permutation VIM (perm) 
# obtained for all simulated datasets with n = 1000.
#################################################################################################

resi <- resultsall %>% filter(n=="n = 1000", method_all %in% c("Discr", "Class-foc"))

resi <- resi %>% mutate(method = factor(method_all, levels = c("Discr", "Class-foc")))

resi <- resi %>%
  group_by(n, K, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p1 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(aes(fill = method), position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Class-focused and discriminatory VIM", fill="VIM type", x = "Covariates", y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("Discr", "Class-foc")) +
  theme(legend.position.inside = c(0.145, 0.94), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16),
        legend.background = element_rect(fill = NA, colour = NA),  # Transparent legend background and no border
        legend.key = element_rect(fill = NA, colour = NA),  # Transparent keys and no border
        legend.box.background = element_rect(fill = NA, colour = NA)) +
  guides(fill = guide_legend(position = "inside")) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)



resi <- resultsall %>% filter(n=="n = 500", method_all %in% c("Perm"))

resi <- resi %>%
  group_by(n, K, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), 
             levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p2 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Permutation VIM", x = "Covariates", y = "VIM values") +
  theme(legend.position = "none", axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, p2, ncol = 2)  # `nullGrob()` is an empty plot

# Figure S3:

ggsave("../figures/FigS3.pdf", grid_plot, width = 12, height = 10)








# Figure S4: VIM values obtained for class-focused and discriminatory VIM and the permutation VIM (perm) 
# obtained for all simulated datasets with n = 2000.
#################################################################################################

resi <- resultsall %>% filter(n=="n = 2000", method_all %in% c("Discr", "Class-foc"))

resi <- resi %>% mutate(method = factor(method_all, levels = c("Discr", "Class-foc")))

resi <- resi %>%
  group_by(n, K, method_all, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p1 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(aes(fill = method), position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Class-focused and discriminatory VIM", fill="VIM type", x = "Covariates", y = "VIM values") +
  scale_fill_manual(values=colors, labels = c("Discr", "Class-foc")) +
  theme(legend.position.inside = c(0.145, 0.94), 
        legend.title = element_text(size=13),
        legend.text = element_text(size=12.5),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16),
        legend.background = element_rect(fill = NA, colour = NA),  # Transparent legend background and no border
        legend.key = element_rect(fill = NA, colour = NA),  # Transparent keys and no border
        legend.box.background = element_rect(fill = NA, colour = NA)) +
  guides(fill = guide_legend(position = "inside")) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)



resi <- resultsall %>% filter(n=="n = 500", method_all %in% c("Perm"))

resi <- resi %>%
  group_by(n, K, itind) %>%
  mutate(seq = row_number()) %>%
  ungroup()

resi <- resi %>% filter(!(seq %in% 1:45))

resi$K <- factor(paste0("C = ", resi$K), levels=c("C = 4", "C = 6", "C = 10"))

resi2 <- data.frame(
  K = factor(rep(c("C = 4", "C = 6", "C = 10"), times=c(5, 6, 7)), 
             levels=c("C = 4", "C = 6", "C = 10")),
  xpos = c(2.75, 7, 10, 13, 16, 2.75, 7, 10, 13, 16, 19, 2.75, 7, 10, 13, 16, 19, 22),
  labels = c("X[no]", "X[two_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]",
             "X[no]", "X[two_gr]", "X[thr_gr]", "X[cl_rel_1]", "X[cl_rel_2]", "X[cl_rel_3]", "X[cl_rel_4]")
)

# Create a data frame with x-intercept values specific for each K group
vline_data <- data.frame(
  K = factor(c(rep("C = 4", 4), rep("C = 6", 5), rep("C = 10", 6)),
             levels = c("C = 4", "C = 6", "C = 10")),
  xintercept = c(5.5, 8.5, 11.5, 14.5, 5.5, 8.5, 11.5, 14.5, 17.5, 5.5, 8.5, 11.5, 14.5, 17.5, 20.5)
)

# Calculate the minimum y-value for each K group and determine the y-position for labels
resi_label_pos <- resi %>%
  group_by(K) %>%
  summarize(min_y = min(vim), .groups = "drop") %>%
  mutate(label_y = min_y + (4 * min_y))

# Join this information back with resi2 to include label positions
resi2 <- resi2 %>%
  left_join(resi_label_pos, by = "K") %>%
  select(K, xpos, labels, label_y)

# Plot
p2 <- ggplot(resi, aes(x = factor(seq), y = vim)) +
  facet_wrap(~ K, ncol=1, scales="free") +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_vline(data = vline_data, aes(xintercept = xintercept)) +
  theme_bw() +
  labs(title = "Permutation VIM", x = "Covariates", y = "VIM values") +
  theme(legend.position = "none", axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.text.y=element_text(size=12),
        plot.title = element_text(size = 16)) +
  geom_text(data = resi2, aes(x = xpos, y = label_y, label = labels), parse = TRUE, size=5)

# Assuming p1 and p2 are your ggplot objects
# Create the plots with an additional empty column for spacing
grid_plot <- arrangeGrob(p1, p2, ncol = 2)  # `nullGrob()` is an empty plot

# Figure S4:

ggsave("../figures/FigS4.pdf", grid_plot, width = 12, height = 10)








# Figure S1: Class-specific distributions of the informative covariates.
######################################################################

color_palette <- brewer.pal(4, "Set1")

xseq <- seq(-4, 7, length.out = 200)

dtemp <- 0.015


plotdata <- data.frame(c=factor(rep(1:4, each=200)),
                       x=rep(xseq, 4),
                       y=c(dnorm(xseq, mean=0, sd=1), 
                           dnorm(xseq, mean=0, sd=1) + dtemp,
                           dnorm(xseq, mean=1.5, sd=1), 
                           dnorm(xseq, mean=1.5, sd=1) + dtemp))

p1 <- ggplot(data=plotdata, aes(x=x, y=y, color=c)) + theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[two_gr]))



plotdata <- data.frame(
  c=factor(rep(1:4, each=200)),
  x=rep(xseq, 4),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=1, sd=1))
)


p2 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_1]))





plotdata <- data.frame(
  c=factor(rep(1:4, each=200)),
  x=rep(xseq, 4),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=1, sd=1),
      dnorm(xseq, mean=2, sd=1))
)


p3 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_2]))





plotdata <- data.frame(
  c=factor(rep(1:4, each=200)),
  x=rep(xseq, 4),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=2.25, sd=1))
)


p4 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_3]))

# Arrange plots in a grid
C4 <- grid.arrange(p1, p2, p3, p4, ncol=3)

# Create a text grob for the title with increased font size
title_grob <- textGrob("C = 4", gp=gpar(fontsize=20))

# Use this grob as the top parameter in a new grid.arrange call
C4 <- grid.arrange(grobs = list(C4), top = title_grob)




color_palette <- brewer.pal(6, "Set1")

xseq <- seq(-4, 7, length.out = 200)


plotdata <- data.frame(c=factor(rep(1:6, each=200)),
                       x=rep(xseq, 6),
                       y=c(dnorm(xseq, mean=0, sd=1), 
                           dnorm(xseq, mean=0, sd=1) + dtemp,
                           dnorm(xseq, mean=0, sd=1) + 2*dtemp,
                           dnorm(xseq, mean=1.5, sd=1), 
                           dnorm(xseq, mean=1.5, sd=1) + dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 2*dtemp))

p1 <- ggplot(data=plotdata, aes(x=x, y=y, color=c)) + theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[two_gr]))



plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=1, sd=1), 
      dnorm(xseq, mean=1, sd=1) + dtemp,
      dnorm(xseq, mean=2, sd=1), 
      dnorm(xseq, mean=2, sd=1) + dtemp)
)


p2 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[thr_gr]))




plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=1, sd=1))
)


p3 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_1]))



plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=1, sd=1),
      dnorm(xseq, mean=2, sd=1))
)


p4 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_2]))





plotdata <- data.frame(
  c=factor(rep(1:6, each=200)),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=2.25, sd=1))
)


p5 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_3]))

# Arrange plots in a grid
C6 <- grid.arrange(p1, p2, p3, p4, p5, ncol=3)

# Create a text grob for the title with increased font size
title_grob <- textGrob("C = 6", gp=gpar(fontsize=20))

# Use this grob as the top parameter in a new grid.arrange call
C6 <- grid.arrange(grobs = list(C6), top = title_grob)




set1_original <- brewer.pal(8, "Set1")
color_palette <- colorRampPalette(set1_original)(10)

xseq <- seq(-4, 7, length.out = 200)


plotdata <- data.frame(c=factor(rep(1:10, each=200)),
                       x=rep(xseq, 10),
                       y=c(dnorm(xseq, mean=0, sd=1), 
                           dnorm(xseq, mean=0, sd=1) + dtemp,
                           dnorm(xseq, mean=0, sd=1) + 2*dtemp,
                           dnorm(xseq, mean=0, sd=1) + 3*dtemp,
                           dnorm(xseq, mean=0, sd=1) + 4*dtemp,
                           dnorm(xseq, mean=1.5, sd=1), 
                           dnorm(xseq, mean=1.5, sd=1) + dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 2*dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 3*dtemp,
                           dnorm(xseq, mean=1.5, sd=1) + 4*dtemp))

p1 <- ggplot(data=plotdata, aes(x=x, y=y, color=c)) + theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[two_gr]))



plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=1, sd=1), 
      dnorm(xseq, mean=1, sd=1) + dtemp,
      dnorm(xseq, mean=1, sd=1) + 2*dtemp,
      dnorm(xseq, mean=2, sd=1), 
      dnorm(xseq, mean=2, sd=1) + dtemp,
      dnorm(xseq, mean=2, sd=1) + 2*dtemp)
)


p2 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[thr_gr]))




plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=0, sd=1) + 5*dtemp,
      dnorm(xseq, mean=0, sd=1) + 6*dtemp,
      dnorm(xseq, mean=0, sd=1) + 7*dtemp,
      dnorm(xseq, mean=0, sd=1) + 8*dtemp,
      dnorm(xseq, mean=1, sd=1))
)


p3 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_1]))



plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=0, sd=1) + 5*dtemp,
      dnorm(xseq, mean=0, sd=1) + 6*dtemp,
      dnorm(xseq, mean=0, sd=1) + 7*dtemp,
      dnorm(xseq, mean=1, sd=1),
      dnorm(xseq, mean=2, sd=1))
)


p4 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_2]))




plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0, sd=1) + 4*dtemp,
      dnorm(xseq, mean=0, sd=1) + 5*dtemp,
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=0.75, sd=1) + dtemp,
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=2.25, sd=1))
)


p5 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_3]))





plotdata <- data.frame(
  c=factor(rep(1:10, each=200)),
  x=rep(xseq, 10),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=0.75, sd=1),
      dnorm(xseq, mean=0.75, sd=1) + dtemp,
      dnorm(xseq, mean=1.5, sd=1),
      dnorm(xseq, mean=1.5, sd=1) + dtemp,
      dnorm(xseq, mean=2.25, sd=1),
      dnorm(xseq, mean=3, sd=1))
)

p6 <- ggplot(plotdata, aes(x=x, y=y, color=c)) + 
  theme_bw() + 
  geom_line() +
  scale_color_manual(values=color_palette) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_blank(),
        legend.position = "none") +
  ylab("density") +
  ggtitle(expression(X[cl_rel_4]))

# Arrange plots in a grid
C10 <- grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)

# Create a text grob for the title with increased font size
title_grob <- textGrob("C = 10", gp=gpar(fontsize=20))

# Use this grob as the top parameter in a new grid.arrange call
C10 <- grid.arrange(grobs = list(C10), top = title_grob)




# Create a spacer plot
spacer <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank())


library("cowplot")

combined_plot <- plot_grid(
  C4,
  spacer,
  C6,
  spacer,
  C10,
  ncol = 1,
  rel_heights = c(1, 0.1, 1, 0.1, 1)
)


# Figure S1:

ggsave("../figures/FigS1.pdf", width=10*0.8, height=13*0.8)








# Figure 5: Class-specific distributions of the covariates X_cl_rel_2 and X_thr_gr for C = 6.
#############################################################################################

xseq <- seq(-4, 7, length.out = 200)

dtemp <- 0.015


plotdata1 <- data.frame(
  c=factor(rep(1:6, each=200)),
  cl_group=factor(rep(c("1, 2, 3, 4", "5", "6"), times=c(4, 1, 1)*200), levels=c("1, 2, 3, 4", "5", "6")),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=1, sd=1), 
      dnorm(xseq, mean=1, sd=1) + dtemp,
      dnorm(xseq, mean=2, sd=1), 
      dnorm(xseq, mean=2, sd=1) + dtemp),
  variable=1
)

plotdata2 <- data.frame(
  c=factor(rep(1:6, each=200)),
  cl_group=factor(rep(c("1, 2, 3, 4", "5", "6"), times=c(4, 1, 1)*200), levels=c("1, 2, 3, 4", "5", "6")),
  x=rep(xseq, 6),
  y=c(dnorm(xseq, mean=0, sd=1), 
      dnorm(xseq, mean=0, sd=1) + dtemp,
      dnorm(xseq, mean=0, sd=1) + 2*dtemp,
      dnorm(xseq, mean=0, sd=1) + 3*dtemp,
      dnorm(xseq, mean=1, sd=1),
      dnorm(xseq, mean=2, sd=1)),
  variable=2
)

plotdata <- rbind(plotdata1, plotdata2)
plotdata$variable <- factor(plotdata$variable, levels=c("2", "1"))


facet_labels <- c(
  `1` = "X[thr_gr]",
  `2` = "X[cl_rel_2]"
)

p <- ggplot(plotdata, aes(x = x, y = y, group = c, linetype = cl_group)) +
  theme_bw() +
  geom_line() +
  facet_wrap(~ variable, nrow = 1, labeller = labeller(variable = facet_labels, .default = label_parsed)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    strip.text.x=element_text(size=14),
    legend.title = element_text(size=14),
    legend.text = element_text(size=14),
    axis.title.y=element_text(size=14)
  ) +
  ylab("Density") +
  labs(linetype = "Class")
p

# Figure 5:

ggsave("../figures/Fig5.eps", width=10*0.8, height=4*0.8)
