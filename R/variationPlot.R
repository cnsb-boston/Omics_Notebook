#-------------------------------------------------
#' Variation Plot and Filter
#'
#' This function generates plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed, or name for the Omics set.
#' @param outputpath output file path for plots
#'
#' @return a vector of most variable features by MAD
#' 
#' @examples
#' 
#' @export
variationPlot <- function(eset, type, outputpath=output_plots_path, 
                                percent_choice=0.1) {
  emat<-exprs(eset)
  MEAN <- apply(emat, 1, mean)
  STDEV <- apply(emat, 1, sd)
  MAD <- apply(emat, 1, mad)
  MED <- apply(emat, 1, median)
  top_hits <- names(MAD)[order(MAD, decreasing=TRUE)[1:(percent_choice*length(MAD))]]
  emat_top <- emat[top_hits[top_hits!=""],]
  
  output_filename<-file.path(outputpath, paste("variation_",type,".pdf", sep=""))
  pdf(output_filename)
  plot(MEAN, STDEV, pch=".", cex=4, main=paste("mean vs. stdev: ", type, sep=""))
  plot(MED, MAD, pch=19, cex=0.5, log="", main=paste("median vs MAD: ",type, sep=""))
  points(MED[top_hits], MAD[top_hits], pch=19, cex=0.5, col="red")
  legend("topright", pch=20, col=c("black", "red"), legend=c("all points", paste("filtered top ", percent_choice*100, "%", sep="")) )
  dev.off()
  
  try({
  output_filename<-file.path(outputpath, paste("corrplot_",type,".pdf", sep=""))
  pdf(output_filename)
  corrplot::corrplot(cor(emat_top), type="upper", order="hclust", tl.col="black")
  dev.off()
  })
  try({
  output_filename<-file.path(outputpath, paste("pairs_",type,".jpeg", sep=""))
  jpeg(output_filename)
  pairs(emat, pch=".")
  dev.off()
  })
  return (top_hits[top_hits!=""])
}
