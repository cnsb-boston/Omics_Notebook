#-------------------------------------------------
#' Draw PCA
#'
#' This function generates PCA plots based on data
#' 
#' @param eset an ExpressionSet object with omcis data
#' @param type Type of data to be processed, or name for the Omics set.
#' @param outputpath output file path for plots
#'
#' @return a pca plot
#' 
#' @examples
#'     drawPCA(eset)
#' 
#' @export
drawPCA <- function(eset, x_axis="PC1", y_axis="PC2", type, outputpath=output_plots_path, show_sample_names=TRUE) {
  
  PC_data <- prcomp(t(exprs(eset)))
  percent_variance <- summary(PC_data)$importance["Proportion of Variance",]*100
  
  make_colors<-TRUE;
  if ("ColorsHex" %in% colnames(pData(eset)) ){ if(checkColor(pData(eset)[,"ColorsHex"])){
    colors_dots <- c(unique(as.character(pData(eset)[,"ColorsHex"])))
    names(colors_dots) <- unique(pData(eset)[,"Group"])
    make_colors<-FALSE;
  } }
  
  pca_graph <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis]))+ 
    { if( "Group2" %in% colnames(pData(eset)) ) {geom_point(aes(shape=pData(eset)$Group2,colour=pData(eset)$Group));} else {geom_point(aes(colour=pData(eset)$Group)); }} + 
    theme_bw() + theme(legend.title=element_blank()) + 
    { if(show_sample_names==TRUE) geom_text_repel(aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis], label=rownames(PC_data$x)), colour="black"); } +
    { if (!make_colors) scale_color_manual(values=colors_dots);  } +
    labs(title=paste("PCA - ",type,sep=""),
       x=paste(x_axis, sprintf(" (%2.0f%%)", percent_variance[x_axis]), sep=""),
       y=paste(y_axis, sprintf(" (%2.0f%%)", percent_variance[y_axis]), sep=""));
  try({
  pca_graph2 <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,"PC2"], y=PC_data$x[,"PC3"]))+ 
    { if( "Group2" %in% colnames(pData(eset)) ) {geom_point(aes(shape=pData(eset)$Group2,colour=pData(eset)$Group));} else {geom_point(aes(colour=pData(eset)$Group)); }} + 
    theme_bw() + theme(legend.title=element_blank()) + 
    { if(show_sample_names==TRUE) geom_text_repel(aes(x=PC_data$x[,"PC2"], y=PC_data$x[,"PC3"], label=rownames(PC_data$x)), colour="black"); } +
    { if (!make_colors) scale_color_manual(values=colors_dots);  } +
    labs(title=paste("PCA - ",type,sep=""),
       x=paste("PC2", sprintf(" (%2.0f%%)", percent_variance["PC2"]), sep=""),
       y=paste("PC3", sprintf(" (%2.0f%%)", percent_variance["PC3"]), sep=""));
  })
  pca_loading_graph <- ggplot(data=as.data.frame(PC_data$rotation), aes(x=(PC_data$rotation[,x_axis]), y=PC_data$rotation[,y_axis])) +
    geom_text(label=rownames(PC_data$rotation), colour="black", size=2) + 
    labs(x=paste(x_axis, sep=""), y=paste(y_axis, sep=""),
         title=paste("PC Factor Loadings\n ",type, sep=""))

  output_filename <- file.path(outputpath, paste("PCAplots_",type,".pdf",sep=""))
  pdf(output_filename, width=4, height=3)
  print(pca_graph)
  try({print(pca_graph2) })
  print(pca_loading_graph)
  plot(PC_data, type = 'l', main=paste("PCs vs. Variance",type,sep=""))
  dev.off()

  return(pca_graph)
}


