#-------------------------------------------------
#' Get legend
#'
#' This function takes omics data set and performs normalization and QC plotting to inspect data. Can be expanded with additional normalization methods.
#' 
#' @param plot1 a ggplot2 plot with a legend
#' 
#' @export
g_legend <- function(plot1){
  tmp <- ggplot_gtable(ggplot_build(plot1))
  leg <- which(sapply(tmp$grobs, function(x) x$name)=="guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}