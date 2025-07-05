#' Boxplot of OTU Relative Abundance by Group
#'
#' Generates a ggplot2-based boxplot showing the relative abundance of a specified OTU (Operational Taxonomic Unit)
#' across different sample groups.
#'
#' @param plot.otu Character vector of OTU names to plot. Must match column names of \code{count.data}.
#' @param count.data A numeric matrix or data frame of OTU counts (samples in rows, OTUs in columns).
#' @param groups A factor vector indicating group membership for each sample (must be same length as number of rows in \code{count.data}).
#'
#' @return A \code{ggplot} object displaying a boxplot of the relative abundance of the selected OTU across groups.
#'
#' @details
#' The function converts count data to relative abundance, selects the OTU specified in \code{plot.otu},
#' and creates a boxplot of relative abundance grouped by \code{groups}. If \code{groups} is not a factor,
#' the function will return an error. The mean of each group is marked with a red dot.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot stat_summary labs theme_minimal theme element_text element_rect
#' @export
#'
#' @examples
#' # Example 1: Basic usage
#' set.seed(123)
#' count.data <- matrix(rpois(300, lambda = 10), nrow = 30, ncol = 10)
#' colnames(count.data) <- paste0("OTU", 1:10)
#' groups <- factor(rep(c("Control", "Treatment"), each = 15))
#' otuboxplot("OTU3", count.data, groups)
#'
#' # Example 2: Create two binary indicators and recombine into a factor group
#' set.seed(456)
#' n <- 40
#' count.data <- matrix(rpois(n * 8, lambda = 12), nrow = n, ncol = 8)
#' colnames(count.data) <- paste0("OTU", 1:8)
#' group1 <- rbinom(n, 1, 0.5)
#' group2 <- rbinom(n, 1, 0.5)
#' group <- factor(paste0("Y", group1, "_C", group2))
#' otuboxplot("OTU5", count.data, group)
otuboxplot = function(plot.otu, count.data, groups){

  if (!inherits(groups, "factor")) {
    stop("\n The sample groups indicator needs to be factors!")
  }

  mylabels = levels(groups)
  nGroups = length(mylabels)

  count.data = as.matrix(count.data)
  otu.names = colnames(count.data)
  plot.otu.id = match(otu.names, plot.otu)
  if(sum(1*(!is.na(plot.otu.id))) == 0){
    stop("\n The requested OTU does not exit in OTU table OR the OTU names cannot be matched!")
  }

  rel_mat = count.data[, plot.otu, drop = FALSE]/rowSums(count.data)
  plotdata = data.frame(
    val   = as.vector(rel_mat),
    group = rep(groups, each = ncol(rel_mat))
  )

  pcommon = ggplot(data = plotdata, aes(x = group, y = val, fill = group)) +
    geom_boxplot(position = position_dodge(width = 1)) +
    stat_summary(aes(group = group), fun = mean, geom = "point", size = 3,
                 color = "red", fill = "red", position = position_dodge(width = 1)) +
    labs(
      title = plot.otu,
      x = NULL,
      y = "Relative Abundance",
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )

  return(pcommon)

}
