#############################################################################################################
# This script contains modified versions of the visualization functions from the diversityForest R package.
# They were modified to produce figures which can also be printed as black and white figures.
############################################################################################################



# -------------------------------------------------------------------------------
#   This file is part of 'diversityForest'.
#
# 'diversityForest' is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# 'diversityForest' is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with 'diversityForest'. If not, see <http://www.gnu.org/licenses/>.
#
#  NOTE: 'diversityForest' is a fork of the popular R package 'ranger', written by Marvin N. Wright.
#  Most R and C++ code is identical with that of 'ranger'. The package 'diversityForest'
#  was written by taking the original 'ranger' code and making any
#  changes necessary to implement diversity forests.
#
# -------------------------------------------------------------------------------

##' Plot function for \code{multifor} objects that allows to obtain a first overview of the result of the
##' multi-class VIM analysis. This function visualises the distribution of the multi-class VIM values
##' together with that of the corresponding discriminatory VIM values and
##' the estimated dependency structures of the multi-class outcome on the variables 
##' with largest multi-class VIM values. These estimated dependency structures are visualised
##' using density plots and/or boxplots.
##' 
##' In the plot showing the distribution of the multi-class VIM values along with 
##' that of the discriminatory VIM values, the discriminatory VIM values are 
##' normalized to make them comparable to the multi-class VIM values. This is 
##' achieved by dividing the discriminatory VIM values by their mean and multiplying 
##' it by that of the multi-class VIM values. Although the discriminatory VIM 
##' values are computed for all variables, only those variables for which the 
##' multi-class VIM values were computed are included in this analysis (i.e., 
##' all variables that have at least as many unique values as there are classes
##' in the outcome variable).\cr
##' For details on the plots of the estimated dependency structures of the 
##' multi-class outcome on the variables, see \code{\link{plotMcl}}.
##' The latter function allows to visualise these estimated dependency structures
##' for arbitrary variables in the data.
##' 
##' @title Plot method for \code{multifor} objects
##' @param x Object of class \code{multifor}.
##' @param plot_type Plot type, one of the following: "both" (the default), "density", "boxplot".  If "density", \code{"density"} plots are produced, if "boxplot", \code{"boxplot"} plots are produced, and if "both", both \code{"density"} plots and \code{"boxplot"} plots are produced. See the 'Details' section of \code{\link{plotMcl}} for details.
##' @param num_best The number of variables with largest multi-class VIM values to plot. Default is 5.
##' @param ... Further arguments passed to or from other methods.
##' @return A ggplot2 plot.
##' @examples
##' \dontrun{
##' 
##' ## Load package:
##' 
##' library("diversityForest")
##' 
##' 
##' 
##' ## Set seed to make results reproducible:
##' 
##' set.seed(1234)
##' 
##' 
##' 
##' ## Construct multi forest and calculate multi-class and discriminatory VIM values:
##' 
##' data(hars)
##' model <- multifor(dependent.variable.name = "Activity", data = hars, 
##'                   num.trees = 100, probability=TRUE)
##' 
##' # NOTE: num.trees = 100 (in the above) would be likely too small for practical 
##' # purposes. This small number of trees was simply used to keep the
##' # runtime of the example short.
##' # The default number of trees is num.trees = 5000 for datasets with a maximum of
##' # 5000 observations and num.trees = 1000 for datasets larger than that.
##' 
##' 
##' 
##' ## By default the estimated class-specific distributions of the num_best=5
##' ## variables with the largest multi-class VIM values are plotted:
##' 
##' plot(model)
##' 
##' ## Consider only the 2 variables with the largest multi-class VIM values:
##' 
##' plot(model, num_best = 2)
##' 
##' ## Show only the density plots or only the boxplots:
##' 
##' plot(model, plot_type = "density", num_best = 2)
##' plot(model, plot_type = "boxplot", num_best = 2)
##' 
##' ## Show only the plot of the distributions of the multi-class and
##' ## discriminatory VIM values:
##' 
##' plot(model, num_best = 0)
##' 
##' }
##'
##' @author Roman Hornung
##' @references
##' \itemize{
##'   \item Hornung, R., Hapfelmeier, A. (2024). Multi forests: Variable importance for multi-class outcomes. arXiv:2409.08925, <\doi{10.48550/arXiv.2409.08925}>.
##'   \item Hornung, R. (2022). Diversity forests: Using split sampling to enable innovative complex split procedures in random forests. SN Computer Science 3(2):1, <\doi{10.1007/s42979-021-00920-1}>.
##'   }
##' @seealso \code{\link{plotMcl}}
##' @encoding UTF-8
##' @importFrom ggplot2 ggplot aes theme_bw geom_point scale_color_manual ylab theme element_blank
##' @importFrom rlang .data
##' @rdname plot.multifor
##' @export
plotmultifor <- function(x, plot_type=c("both", "density", "boxplot")[1], num_best=5, ...) {
  
  if (num_best < 0)
    stop("'num_best' must be an integer greater than or equal to zero.")
  
  # Extract the multi-class and discriminatory VIM values sort the values in
  # decreasing order according to the multi-clas VIM values:
  
  vim_multiclass <- x$var.imp.multiclass
  vim_discr <- x$var.imp.discr
  
  if (all(is.na(vim_multiclass)))
    stop("There are no (non-NA) multi-class VIM values.")
  
  vim_multiclass_noNA <- vim_multiclass[!is.na(vim_multiclass)]
  vim_discr_noNA <- vim_discr[!is.na(vim_multiclass)]
  
  reorderind <- order(vim_multiclass_noNA, decreasing=TRUE)
  vim_multiclass_order <- vim_multiclass_noNA[reorderind]
  vim_discr_order <- vim_discr_noNA[reorderind]
  
  
  # Rescale the discriminatory VIM values, so that they have the same
  # mean as the multi-class VIM values. This is done so that the two types of
  # VIM can be compared visually:
  
  vim_discr_order_resc <- vim_discr_order*mean(vim_multiclass_order)/mean(vim_discr_order)
  
  
  # Plot the multi-class and discriminatory VIM values:
  
  vim_multiclass_names <- names(vim_multiclass_order)
  
  datacov <- x$plotres$data[,vim_multiclass_names, drop=FALSE]
  y_outcome <- x$plotres$data[,x$plotres$yvarname]
  
  dataplot <- data.frame(x2=rep(1:length(vim_multiclass_order), 2), vim=c(vim_multiclass_order, vim_discr_order_resc),
                         type=factor(rep(c("multi-class", "discriminatory (normalized)"), each=length(vim_multiclass_order)), levels=c("multi-class", "discriminatory (normalized)")))
  
  p <- ggplot(data=dataplot, aes(x=.data$x2, y=.data$vim, color=.data$type)) + theme_bw() + geom_point() + 
    scale_color_manual(values=c("black", "grey"), name="VIM type") + ylab("VIM value") +
    theme(legend.position = c(0.95, 0.95),  # Coordinates for top-right corner
          legend.justification = c(1, 1),
          axis.title.x = element_blank())
  print(p)
  readline(prompt="Press [enter] for next plot.")
  
  if (num_best > 0) {
    
    if (num_best > ncol(datacov)) {
      warning(paste0("The value num_best=", num_best, " is larger than the number of covariates with multi-class VIM values. --> The value num_best was set to ", ncol(datacov), "."))
      num_best <- ncol(datacov)
    }
    
    
    # Make plots of the class-specific distributions of the values of the num_best
    # covariates with the largest multi-class VIM values:
    
    for(i in 1:num_best) {
      plot_title <- paste0(vim_multiclass_names[i], "  -  rank ", i, " according to the multi-class VIM")
      myplotVar(datacov[,i], y_outcome, x_label=vim_multiclass_names[i], y_label = x$plotres$yvarname, plot_title=plot_title, plot_type=plot_type, plotit=TRUE)
      if(i < num_best)
        readline(prompt="Press [enter] for next plot.")
    }
    
  }
  
}









# -------------------------------------------------------------------------------
#   This file is part of 'diversityForest'.
#
# 'diversityForest' is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# 'diversityForest' is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with 'diversityForest'. If not, see <http://www.gnu.org/licenses/>.
#
#  NOTE: 'diversityForest' is a fork of the popular R package 'ranger', written by Marvin N. Wright.
#  Most R and C++ code is identical with that of 'ranger'. The package 'diversityForest'
#  was written by taking the original 'ranger' code and making any
#  changes necessary to implement diversity forests.
#
# -------------------------------------------------------------------------------

##' This function allows to visualise the (estimated) distributions of one or several variables for each of the classes of the outcomes.
##' This allows to study how exactly variables of interest influence the outcome, which is crucial for interpretive purposes.
##' Two types of visualisations are available: density plots and boxplots. See the 'Details' section below for further explanation.
##'
##' For the \code{"density"} plots, kernel density estimates (obtained using the 
##' \code{density()} function from base R) of the within-class distributions are 
##' plotted in the same plot using different colors and, depending on the number 
##' of classes, different line types. To account for the different number of
##' observations per class, each density is multiplied by the proportion of 
##' observations from that class. The resulting scaled densities can be interpreted 
##' in terms of the local density of the observations from each class relative to 
##' those from the other classes. For example, if a scaled density has the largest 
##' value in a particular region, this can be interpreted as the respective class 
##' being the most frequent in that region. Another example: If the scaled density 
##' of class "A" is twice as large as the scaled density of class "B" in a particular 
##' region, this can be interpreted to mean that there are twice as many observations 
##' of class "A" as of class "B" in that region.
##' 
##' In the \code{"density"} plots, only classes represented by at least two 
##' observations are considered. If the number of classes is greater than 7, 
##' the different classes are distinguished using both colors and line styles. 
##' To indicate the absolute numbers of observations in the different regions, 
##' the locations of the observations from the different classes are visualized 
##' using a rug plot on the x-axis, using the same colors and line types as for 
##' the density plots. If the number of observations is greater than 1,000, a 
##' random subset of 1,000 observations is shown in the rug plot instead of all 
##' observations for visual clarity.
##' 
##' The \code{"boxplot"} plots show the (estimated) within-class distributions 
##' side by side using boxplots. All classes are considered, even those represented 
##' by only a single observation. For the \code{plot_type="both"} option, which 
##' displays both \code{"density"} and \code{"boxplot"} plots, the boxplots are 
##' displayed using the same colors and ( if applicable) line styles as the kernel 
##' density estimates, for clarity. Boxplots of classes for which no kernel density 
##' estimates were obtained (i.e., those of the classes represented by single 
##' observations) are shown in grey. 
##' 
##' Note that plots are only generated for those variables in \code{varnames} 
##' that have at least as many unique values as there are outcome classes. For 
##' categorical variables, the category labels are printed on the x- or y-axis 
##' of the \code{"density"} or \code{"boxplot"} plots, respectively. The rug plots 
##' of the \code{"density"} plots are produced only for numeric variables.
##' 
##' @title Plots of the (estimated) within-class distributions of variables
##' @param data Data frame containing the variables.
##' @param yvarname Name of outcome variable.
##' @param varnames Names of the variables for which plots should be created.
##' @param plot_type Plot type, one of the following: "both" (the default), "density", "boxplot".  If "density", \code{"density"} plot are produced, if "boxplot", \code{"boxplot"} plots are produced, and if "both", both \code{"density"} plots and \code{"boxplot"} plots are produced. See the 'Details' section below for details.
##' @param addtitles Set to \code{TRUE} (default) to add headings providing the names of the respective variables to the plots.
##' @param plotit This states whether the plots are actually plotted or merely returned as \code{ggplot} objects. Default is \code{TRUE}.
##' @return A list of ggplot2 plots returned invisibly.
##' @examples
##' \dontrun{
##'
##' ## Load package:
##' 
##' library("diversityForest")
##' 
##' 
##' 
##' ## Plot "density" and "boxplot" plots (default: plot_type = "both") for the 
##' ## first three variables in the "hars" dataset:
##' 
##' data(hars)
##' plotMcl(data = hars, yvarname = "Activity", varnames = c("tBodyAcc.mean...X", 
##'                                                          "tBodyAcc.mean...Y", 
##'                                                          "tBodyAcc.mean...Z"))
##' 
##' 
##' ## Plot only the "density" plots for these variables:
##' 
##' plotMcl(data = hars, yvarname = "Activity", 
##'         varnames = c("tBodyAcc.mean...X", "tBodyAcc.mean...Y", 
##'                      "tBodyAcc.mean...Z"), plot_type = "density")
##' 
##' ## Plot the "density" plots for these variables, but without titles of the
##' ## plots:
##' 
##' plotMcl(data = hars, yvarname = "Activity", varnames = 
##'           c("tBodyAcc.mean...X", "tBodyAcc.mean...Y", "tBodyAcc.mean...Z"), 
##'         plot_type = "density", addtitles = FALSE)
##' 
##' 
##' ## Make density plots for these variables, but only save them in a list "ps"
##' ## without plotting them ("plotit = FALSE"):
##' 
##' ps <- plotMcl(data = hars, yvarname = "Activity", varnames = 
##'                 c("tBodyAcc.mean...X", "tBodyAcc.mean...Y", 
##'                   "tBodyAcc.mean...Z"), plot_type = "density", 
##'               addtitles = FALSE, plotit = FALSE)
##' 
##' 
##' ## The plots can be manipulated later by using ggplot2 functionalities:
##' 
##' library("ggplot2")
##' p1 <- ps[[1]] + ggtitle("First variable in the dataset") + 
##'   labs(x="Variable values", y="my scaled density")
##' 
##' p2 <- ps[[3]] + ggtitle("Third variable in the dataset") + 
##'   labs(x="Variable values", y="my scaled density")
##' 
##' 
##' ## Combine both of the above plots:
##' 
##' library("ggpubr")
##' p <- ggarrange(p1, p2, ncol = 2)
##' p
##' 
##' ## # Save as PDF:
##' ## ggsave(file="mypathtofolder/FigureXY1.pdf", width=14, height=6)
##' 
##' }
##'
##' @author Roman Hornung
##' @references
##' \itemize{
##'   \item Hornung, R., Hapfelmeier, A. (2024). Multi forests: Variable importance for multi-class outcomes. arXiv:2409.08925, <\doi{10.48550/arXiv.2409.08925}>.
##'   \item Hornung, R. (2022). Diversity forests: Using split sampling to enable innovative complex split procedures in random forests. SN Computer Science 3(2):1, <\doi{10.1007/s42979-021-00920-1}>.
##'   }
##' @seealso \code{\link{plot.multifor}}, \code{\link{plotVar}}
##' @encoding UTF-8
##' @import stats 
##' @import utils
##' @export
myplotMcl <- function(data, yvarname, varnames, plot_type=c("both", "density", "boxplot")[1], addtitles = TRUE, plotit = TRUE) {
  
  if (!all(varnames %in% names(data)))
    stop("Not all entries of 'varnames' are found in 'data'.")
  
  datacov <- data[,varnames, drop=FALSE]
  y_outcome <- data[,yvarname]
  
  # Plots are created only for covariates that have at least as many unique values 
  # as there are outcome classes:
  suit_inds <- which(apply(datacov, 2, function(x) length(unique(x)) >= length(unique(y_outcome))))
  
  if (length(suit_inds)==0)
    stop("None of the variables in 'varnames' have at least as many unique values as there are outcome classes. --> Nothing to plot.")
  
  if (length(suit_inds) < ncol(datacov)) {
    varnames_notused <- paste0("\"", varnames[-suit_inds], "\"")
    warning(paste0("Not all variables in 'varnames' have at least as many unique values as there are outcome classes.\n--> For the following variables in 'varnames' no plots were generated:\n", paste(varnames_notused, collapse=", ")))
    datacov <- datacov[,suit_inds]
    varnames <- varnames[suit_inds]
  }
  
  
  # Create the plots and store them in a list:
  
  ps <- list()
  for(i in 1:ncol(datacov)) {
    plot_title <- ""
    if (addtitles)
      plot_title <- varnames[i]
    p <- myplotVar(datacov[,i], y_outcome, x_label=varnames[i], y_label=yvarname, plot_title=plot_title, plot_type=plot_type, plotit=plotit)
    if (plotit) {
      if(i < ncol(datacov))
        readline(prompt="Press [enter] for next plot.")
    }
    ps[[i]] <- p
  }
    
  names(ps) <- names(datacov)
  
  
  # Return the list of plots invisibly:
  
  invisible(ps)
  
}












# -------------------------------------------------------------------------------
#   This file is part of 'diversityForest'.
#
# 'diversityForest' is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# 'diversityForest' is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with 'diversityForest'. If not, see <http://www.gnu.org/licenses/>.
#
#  NOTE: 'diversityForest' is a fork of the popular R package 'ranger', written by Marvin N. Wright.
#  Most R and C++ code is identical with that of 'ranger'. The package 'diversityForest'
#  was written by taking the original 'ranger' code and making any
#  changes necessary to implement diversity forests.
#
# -------------------------------------------------------------------------------

##' This function allows to visualise the (estimated) distributions of a variable \code{x} for each of the categories of a categorical variable \code{y}.
##' This allows to study the dependency structure of \code{y} on \code{x}.
##' Two types of visualisations are available: density plots and boxplots.
##'
##' See the 'Details' section of \code{\link{plotMcl}}.
##' 
##' @title Plot of the (estimated) dependency structure of a variable \code{x} on a categorical variable \code{y}
##' @param x Metric variable or ordered categorical variable that has at least as many unique values as \code{y}
##' @param y Factor variable with at least three categories.
##' @param plot_type Plot type, one of the following: "both" (the default), "density", "boxplot".  If "density", a \code{"density"} plot is produced, if "boxplot", a \code{"boxplot"} is produced, and if "both", both a \code{"density"} plot and a \code{"boxplot"} are produced. See the 'Details' section of \code{\link{plotMcl}} for details.
##' @param x_label Optional. The label of the x-axis.
##' @param y_label Optional. The label (heading) of the legend that differentiates the categories of \code{y}.
##' @param plot_title Optional. The title of the plot.
##' @return A ggplot2 plot.
##' @examples
##' \dontrun{
##' 
##' ## Load package:
##' 
##' library("diversityForest")
##' 
##' 
##' 
##' ## Load the "ctg" dataset:
##' 
##' data(ctg)
##' 
##' 
##' ## Set seed to make results reproducible (this is necessary because
##' ## the rug plot produced by 'plotVar' does not show all observations, but
##' ## only a random subset of 1000 observations):
##' 
##' set.seed(1234)
##' 
##' 
##' ## Using a "density" plot and a "boxplot", visualise the (estimated) 
##' ## distributions of  the variable "Mean" for each of the categories of the 
##' # variable "Tendency":
##' 
##' plotVar(x = ctg$Mean, y = ctg$Tendency)
##' 
##' 
##' ## Re-create this plot with labels:
##' 
##' plotVar(x = ctg$Mean, y = ctg$Tendency, x_label = "Mean of the histogram ('Mean')",
##'         y_label = "Histogram tendency ('Tendency')", 
##'         plot_title = "Relationship between 'Mean' and 'Tendency'")
##' 
##' 
##' ## Re-create this plot, but only show the "density" plot:
##' 
##' plotVar(x = ctg$Mean, y = ctg$Tendency, plot_type = "density",
##'         x_label = "Mean of the histogram ('Mean')", 
##'         y_label = "Histogram tendency ('Tendency')", 
##'         plot_title = "Relationship between 'Mean' and 'Tendency'")
##' 
##' 
##' ## Use ggplot2 and RColorBrewer functionalities to change the line colors and
##' ## the labels of the categories of "Tendency":
##' 
##' library("ggplot2")
##' library("RColorBrewer")
##' p <- plotVar(x = ctg$Mean, y = ctg$Tendency, plot_type = "density",
##'              x_label = "Mean of the histogram ('Mean')", 
##'              y_label = "Histogram tendency ('Tendency')", 
##'              plot_title = "Relationship between 'Mean' and 'Tendency'") +
##'   scale_color_manual(values = brewer.pal(n = 3, name = "Set2"),
##'                      labels = c("left asymmetric", "symmetric", 
##'                                 "right asymmetric")) +
##'   scale_linetype_manual(values = rep(1, 3),
##'                         labels = c("left asymmetric", "symmetric", 
##'                                    "right asymmetric"))
##' 
##' p
##' 
##' ## # Save as PDF:
##' ## ggsave(file="mypathtofolder/FigureXY1.pdf", width=10, height=7)
##' 
##' }
##'
##' @author Roman Hornung
##' @references
##' \itemize{
##'   \item Hornung, R., Hapfelmeier, A. (2024). Multi forests: Variable importance for multi-class outcomes. arXiv:2409.08925, <\doi{10.48550/arXiv.2409.08925}>.
##'   \item Hornung, R. (2022). Diversity forests: Using split sampling to enable innovative complex split procedures in random forests. SN Computer Science 3(2):1, <\doi{10.1007/s42979-021-00920-1}>.
##'   }
##' @seealso \code{\link{plotMcl}}, \code{\link{plot.multifor}}
##' @encoding UTF-8
##' @importFrom ggplot2 ggplot aes geom_line geom_rug theme_bw theme labs ggtitle xlab ylab scale_x_continuous scale_y_continuous scale_color_manual scale_linetype_manual geom_boxplot element_text
##' @export
myplotVar <- function(x, y, plot_type=c("both", "density", "boxplot")[1], x_label="", y_label="", plot_title="", plotit=TRUE) {
  
  # If plot_type=="density", create a density plot:
  if (plot_type=="density") {
    dens_pl <- plotVarDensity(x=x, y=y, x_label=x_label, y_label=y_label, plot_title=plot_title)$p
	if (plotit) {
	  print(dens_pl)
	}
	res <- list(dens_pl=dens_pl)
	invisible(res)
  }
  
  # If plot_type=="boxplot", create a boxplot:
  if (plot_type=="boxplot") {
    boxplot_pl <- plotVarBoxplot(x=x, y=y, x_label=x_label, y_label=y_label, plot_title=plot_title)
	if (plotit) {
	  print(boxplot_pl)
	}
	res <- list(boxplot_pl=boxplot_pl)
	invisible(res)
  }
  
  # If plot_type=="both", create both a density plot a boxplot:
  if (plot_type=="both") {
    # Create the density plot:
    dens_res <- plotVarDensity(x=x, y=y, x_label=x_label, y_label=y_label, plot_title="")
	dens_pl <- dens_res$p
	boxplot_pl <- plotVarBoxplot(x=x, y=y, x_label=x_label, y_label=y_label, plot_title="", plotres=dens_res$plotres)
    # Add the boxplot using the same colors and line types as the density plot (through 'plotres=dens_res$plotres'):
    combined_pl <- patchwork::wrap_plots(dens_pl, boxplot_pl, ncol = 2)
    combined_pl <- combined_pl +
      patchwork::plot_annotation(
        title = plot_title,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5)
        )
      )
	if (plotit) {
	  print(combined_pl)
	}
	res <- list(combined_pl=combined_pl, dens_pl=dens_pl, boxplot_pl=boxplot_pl)
    invisible(res)
  }
  
}


plotVarDensity <- function(x, y, x_label="", y_label="", plot_title="") {

  classtab <- table(y)
  
  # The densities are plotted only for classes with at least two observations:
  levels_to_keep <- names(classtab[classtab >= 2])
  
  filterbool <- y %in% levels_to_keep
  
  x <- x[filterbool]
  y <- y[filterbool]
  
  if (length(unique(x)) < length(unique(y)))
    stop("The number of unique covariate values must be at least as large as the number of classes.")
  
  allclasses <- levels(y)[levels(y) %in% unique(y)]
  
  classtab <- classtab[classtab >= 2]
  classprob <- classtab/sum(classtab)
  
  # The maximum number of different colors used. If the number of classes is larger
  # than this, the different classes are differentiated visually using both
  # colors and line types:
  nmax <- min(c(length(allclasses), 7))
  
  colors <- scales::hue_pal()(nmax)
  
  #if (length(allclasses) == nmax) {
  #  colorsvec <- colors
  #  linetypesvec <- rep("solid", length=length(colorsvec))
  #} else {
  #  colorsvec <- rep(colors, length=length(allclasses))
  #  
  #  linetypesvec <- rep(c("solid", "longdash", "dotdash"), each=nmax)[1:length(colorsvec)]
  #  linetypesvec <- c(linetypesvec, rep("dotdash", times=length(colorsvec) - length(linetypesvec)))
  #}
  
  colorsvec <- colors
  
  linetypesvec <- c("solid", "longdash", "dotdash", "solid", "longdash", "dotdash")
  shapesvec <- c(16, 17, 18, 16, 17, 18)  # Different symbols for points (circle, triangle, diamond)
  
  
  # Create a density plot for a numeric covariate:
  
  if (inherits(x, "numeric")) {
    
    denstemps <- list()
    
    for(i in seq(along=allclasses)) {
      xtemp <- x[y==allclasses[i]]
      
      denstemp <- density(xtemp)
      denstemp <- data.frame(x=denstemp$x, y=denstemp$y)
      # The density values are scaled by the class sizes:
      denstemp$y <- denstemp$y*classprob[i]
      denstemps[[i]] <- denstemp
    }
    
    plotdata <- do.call("rbind", denstemps)
    plotdata$class <- factor(rep(allclasses, times=sapply(denstemps, nrow)), levels=allclasses)
    
	  npoints <- 4
	  gridpoints <- seq(0, 1, length=npoints+2)
	  gridpoints <- gridpoints[-c(1,length(gridpoints))]
	  
	  allclasses <- unique(plotdata$class)
	  
	  plotdatapoints <- plotdata[1, , drop=FALSE]
	  
	  for(i in seq(along=allclasses)) {
	  
	  plotdatatemp <- plotdata[plotdata$class==allclasses[i],]
	  plotdatapoints <- rbind(plotdatapoints, plotdatatemp[sapply(quantile(plotdatatemp$x, gridpoints), function(x) which.min(abs(plotdatatemp$x - x))),])
	  
	  }
	  plotdatapoints <- plotdatapoints[-1,]
	  
    pointdata <- data.frame(x=x, class=y)
    pointdata$class <- droplevels(pointdata$class)
    
    # If there are more than 1000 observations, the rug plot on the lower margin
    # only shows a random subset of 1000 observations:
    if (nrow(pointdata) > 1000) {
      pointdata <- pointdata[sample(1:nrow(pointdata), size=1000),]
    }
    
   # p <- ggplot(plotdata, aes(x=.data$x, color=.data$class, linetype=.data$class, shape = .data$class)) + theme_bw() + geom_line(aes(y=.data$y)) + 
#	  geom_point(data = plotdatapoints, aes(y=.data$y), size = 2) + 
      #scale_color_manual(values=colorsvec) + scale_linetype_manual(values = linetypesvec) +
	  #scale_shape_manual(values = shapesvec) + theme(legend.key.width = unit(2, "cm")) + 
      #ylab("(scaled) density") + geom_rug(data=pointdata, sides="b")
	      
	p <- ggplot(plotdata, aes(x = .data$x, color = .data$class, linetype = .data$class, shape = .data$class)) + 
  theme_bw() + 
  geom_line(aes(y = .data$y)) + 
  
  # Adjust size dynamically: shape 18 → size 3, others → size 2
  geom_point(data = plotdatapoints, aes(y = .data$y, size = ifelse(shapesvec[as.numeric(.data$class)] == 18, 3, 2)))+ 
  
  scale_color_manual(values = colorsvec) + 
  scale_linetype_manual(values = linetypesvec) + 
  scale_shape_manual(values = shapesvec) + 
  scale_size_identity() +  # Prevents ggplot from treating size as a factor
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  theme(legend.key.width = unit(2, "cm")) + 
  ylab("(scaled) density") + 
  geom_rug(data = pointdata, sides = "b") 
		  
		  
  }
  
  # Create a density plot for a factor covariate:
  
  if (inherits(x, "ordered") || inherits(x, "factor")) {
    
    if (inherits(x, "factor"))
      warning("The plot is likely not meaningful because the variable is an unordered factor..")
    
    x_levels <- levels(x)[levels(x) %in% unique(x)]
    
    # For plotting, the factor variable is transformed to a continuous variable:
    x <- as.numeric(x)
    
    denstemps <- list()
    
    for(i in seq(along=allclasses)) {
      xtemp <- x[y==allclasses[i]]
      
      denstemp <- density(xtemp)
      denstemp <- data.frame(x=denstemp$x, y=denstemp$y)
      denstemp$y <- denstemp$y*classprob[i]
      denstemps[[i]] <- denstemp
    }
    
    plotdata <- do.call("rbind", denstemps)
    plotdata$class <- factor(rep(allclasses, times=sapply(denstemps, nrow)), levels=allclasses)
    
    if (x_label=="")
      xlabadd <- theme(axis.title.x=element_blank())
    else
      xlabadd <- xlab(x_label)
    
    x_unique_sorted <- sort(unique(x))
    
    p <- ggplot(plotdata, aes(x=.data$x, y=.data$y, color=.data$class, linetype=.data$class)) + theme_bw() + geom_line() + 
      scale_color_manual(values=colorsvec) + scale_linetype_manual(values = linetypesvec) +
      # The labels of the categories of the covariate are added to the x-axis:
      scale_x_continuous(breaks=x_unique_sorted, labels=x_levels) +
      ylab("density") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
  }
  
  
  # Add labels to the plot if provided:
  
  if (x_label=="")
    p <- p + theme(axis.title.x=element_blank())
  else
    p <- p + xlab(x_label)
  
  if (y_label!="")
    p <- p + labs(colour=y_label, linetype=y_label, shape=y_label)
  
  if (plot_title!="")
    p <- p + ggtitle(plot_title)
	
  # The information on the colors and linetypes of the classes are returned too
  # because these are required by "plotVarBoxplot" in cases in which both the
  # densities and the boxplots are plotted:
  plotres <- list(allclasses=allclasses, colorsvec=colorsvec, linetypesvec=linetypesvec, shapesvec=shapesvec)
  reslist <- list(p=p, plotres=plotres)
  
  return(reslist)
  
}


plotVarBoxplot <- function(x, y, x_label="", y_label="", plot_title="", plotres=NULL) {
  
  # Create a boxplot for a numeric covariate:
  
  if (inherits(x, "numeric")) {
    
    plotdata <- data.frame(x=x, y=y)
    
    # If no information on the colors and line types of the boxplots is provided
    # (usually that returned by "plotVarDensity") boxplots with black lines are generated:
    if (is.null(plotres))
      p <- ggplot(plotdata, aes(x=.data$y, y=.data$x)) + theme_bw() + geom_boxplot() 
    else {
      # If information on the colors and line types of the boxplots is provided
      # boxplots with the specified colors and line types are generated:
      classes_dens <- plotres$allclasses
      colorsvec_dens <- plotres$colorsvec
      linetypesvec_dens <- plotres$linetypesvec
	  shapesvec <- plotres$shapesvec
      
      classtab <- table(y)
      
      classes_present <- names(classtab[classtab >= 1])
      
      colorsvec <- linetypesvec <- rep("", length(classes_present))
      colorsvec[classes_present %in% classes_dens] <- colorsvec_dens
      linetypesvec[classes_present %in% classes_dens] <- linetypesvec_dens
      
      # Classes for which no colors or line types are provided are depicted
      # in grey:
      colorsvec[colorsvec==""] <- "grey"
      linetypesvec[linetypesvec==""] <- "solid"
      
      ysub <- unique(plotdata$y)
      xsub <- sapply(ysub, function(x) median(plotdata$x[plotdata$y==x]))
      plotdatamedian <- data.frame(x=xsub, y=ysub)
      
      p <- ggplot(plotdata, aes(x=.data$y, y=.data$x, linetype=.data$y)) + theme_bw() + geom_boxplot(aes(color=.data$y)) +
	    geom_point(data=plotdatamedian, aes(shape = .data$y, color=.data$y), size=3) + 
        scale_color_manual(values=colorsvec) +
        scale_linetype_manual(values=linetypesvec) + scale_shape_manual(values = shapesvec) + theme(legend.position = "none")
    
	  shapes2 <-  c(1, 2, 5, 1, 2, 5)
	
	  for (i in 1:6) {
	    p <- p + geom_point(data=plotdatamedian[plotdatamedian$y==ysub[i],,drop=FALSE], aes(shape = .data$y), pch=shapes2[i], size=2.5)
	  }
	
	}
	
	
    
  }
  
  # Create a density plot for a factor covariate:
  
  if (inherits(x, "ordered") || inherits(x, "factor")) {
    
    if (inherits(x, "factor"))
      warning("The plot is likely not meaningful because the variable is an unordered factor.")
    
    x_levels <- levels(x)[levels(x) %in% unique(x)]
    
    # For plotting, the factor variable is transformed to a continuous variable:
    x <- as.numeric(x)
    
    plotdata <- data.frame(x=x, y=y)
    
    x_unique_sorted <- sort(unique(x))
    
    if (is.null(plotres))
      p <- ggplot(plotdata, aes(x=.data$y, y=.data$x)) + theme_bw() + geom_boxplot() +
      scale_y_continuous(breaks=x_unique_sorted, labels=x_levels)
    else {
      classes_dens <- plotres$allclasses
      colorsvec_dens <- plotres$colorsvec
      linetypesvec_dens <- plotres$linetypesvec
      
      classtab <- table(y)
      
      classes_present <- names(classtab[classtab >= 1])
      
      colorsvec <- linetypesvec <- rep("", length(classes_present))
      colorsvec[classes_present %in% classes_dens] <- colorsvec_dens
      linetypesvec[classes_present %in% classes_dens] <- linetypesvec_dens
      
      colorsvec[colorsvec==""] <- "grey"
      linetypesvec[linetypesvec==""] <- "solid"
      
      p <- ggplot(plotdata, aes(x=.data$y, y=.data$x, color=.data$y, linetype=.data$y)) + theme_bw() + geom_boxplot() +
        # The labels of the categories of the covariate are added to the x-axis:
        scale_y_continuous(breaks=x_unique_sorted, labels=x_levels) +
        scale_color_manual(values=colorsvec) +
        scale_linetype_manual(values=linetypesvec) + theme(legend.position = "none")
    }
    
  }
  
  # Add labels to the plot if provided:
  
  if (x_label=="")
    p <- p + theme(axis.title.x=element_blank())
  else
    p <- p + ylab(x_label)
  
  if (y_label=="")
    p <- p + xlab("class")
  else
    p <- p + xlab(y_label)
  
  if (plot_title!="")
    p <- p + ggtitle(plot_title)
  
  p
  
}
