#-------------------------------------------------
#' Intensity Normalization
#'
#' This function takes omics data set and performs normalization and QC plotting to inspect data. Can be expanded with additional normalization methods.
#' 
#' @param eset an ExpressionSet object with omics data
#' @param annotate An annotation data frame
#' @param type Type of data to be processed, name for the Omics data set.
#' @param norm one of the supported normalization methods (Currently "none", "quantile", or "loess").
#' @param outputpath output file path for plots
#'
#' @return an expression set object, with normalized data
#' 
#' @examples
#'   eset <- intensityNorm(eset, norm=norm_method)}
#' 
#' @export
intensityNorm <- function(eset, norm, type, outputpath=output_plots_path,
                          annotate=annot, cutoff=0, data_format ){
  
  col_palette <- rainbow(length(levels(as.factor(pData(eset)$Group))))
  annotCol <- col_palette[as.factor(pData(eset)$Group)]
  try({if("Color_Groups" %in% colnames(pData(eset))) {
    if( checkColor(pData(eset)[,"Color_Groups"]) ){
      annotCol <- pData(eset)[,"Color_Groups"]
    }
  } })

  if( data_format!="Metabolites" ) { eset <- eset[ rowSums(exprs(eset)>0)>=cutoff*ncol(exprs(eset)), ]; }
  
  eset_matrix <- exprs(eset)
  df <- data.frame(reshape2::melt(eset_matrix, id.vars = NULL));
  colnames(df) <- c("Feature","Sample", "Intensity");
  plot1 <- ggplot(df, aes(x = Sample, y = Intensity)) + geom_boxplot(fill=annotCol) +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) +  
    labs(title=paste(type, " log2 Intensity: Before Normalization", sep=''));
  plot2 <- ggplot(df, aes(x=Intensity, colour=Sample)) + geom_density()+ theme_bw() + 
    labs(title=paste(type, " log2 Intensity: Before Normalization", sep=''));
  
  if(norm=='quantile'){
    eset_matrix_norm <- as.matrix(normalize.quantiles(eset_matrix))
    dimnames(eset_matrix_norm) <- dimnames(eset_matrix)
  } else if (norm=='loess'){
    eset_matrix_norm <- as.matrix(normalizeCyclicLoess(eset_matrix, method='pairs'))
  } else if (norm=='median'){
    eset_matrix_norm <- as.matrix(normalizeMedianValues(eset_matrix))
  } else if (norm=='z transform'){
    eset_matrix_norm <- as.matrix(scale(eset_matrix));
  } else {  eset_matrix_norm<-eset_matrix }
  
  if (norm != "none"){
    df <- data.frame(reshape2::melt(eset_matrix_norm, id.vars = NULL));
    colnames(df) <- c("Feature","Sample", "Intensity");
    plot3 <- ggplot(df, aes(x = Sample, y = Intensity)) + geom_boxplot(fill=annotCol) +
      theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
      labs(title=paste(type, " log2 Intensity: After ",
                       toupper(norm)," Normalization", sep=''));
    plot4 <- ggplot(df, aes(x=Intensity, colour=Sample)) + geom_density()+ theme_bw() + 
      labs(title=paste(type, " log2 Intensity: After ",
                       toupper(norm)," Normalization", sep=''));
  }
  
  output_filename<-file.path(outputpath, paste("boxplot_histogram_",type,".pdf", sep=''));
  pdf(output_filename);
  print(plot1);
  if(norm!="none") { print(plot3);}
  print(plot2);
  if(norm!="none") { print(plot4); }
  dev.off();
  
  if(norm!="none") { 
    grid.arrange(plot1+labs(title="Raw Data"), plot3+labs(title="Normalized"),
                 plot2+theme(legend.position="none")+labs(title="Raw Data"), 
                 plot4+theme(legend.position="none")+labs(title="Normalized"),
                 top=paste(type, ":", sep=""), ncol=4);
  } else {
    grid.arrange(plot1+labs(title="Intensity"),
                 plot2+theme(legend.position="none")+labs(title="Intensity"), 
                 top=paste(type, ":", sep=""), ncol=4);
    
  }
  
  output_filename<-file.path(outputpath, paste("MAplots_",type,".pdf", sep=''))
  pdf(output_filename)
  for (i in 1:ncol(eset)) {
    limma::plotMA(eset_matrix, array=i, main=paste(type, "MA Plot",i,sep=" "))
    if(norm!="none"){
      limma::plotMA(eset_matrix_norm, array=i, 
                  main=paste(type,toupper(norm),"Normalized MA Plot",i,sep=" "))
    }
  }
  dev.off();
  
  mn <- apply(eset_matrix_norm, 1, median);
  rle <- data.frame(sweep(eset_matrix_norm, MARGIN=1, STATS=mn, FUN='-'));
  output_filename <- file.path(outputpath, paste("RLE_",type,".pdf", sep=''));
  pdf(output_filename, width=11, height=8.5);
  par(mar=c(10,4,4,2)+0.1)
  boxplot(rle, main="RLE (Relative Log Expression)\nShould be centered on 0 (blue line)",
          names=colnames(eset_matrix_norm), las=2)
  lines(x=c(0,ncol(eset_matrix_norm)+1), y=rep(0,2), col="blue", lty=2)
  lines(x=c(0,ncol(eset_matrix_norm)+1), y=rep(0.1,2), col="red", lty=2)
  par(mar=c(5,4,4,2)+0.1)
  dev.off();
  
  eset_norm <- ExpressionSet(assayData=eset_matrix_norm);
  pData(eset_norm) <- pData(eset);
  rownames(pData(eset_norm)) <- colnames(eset_matrix_norm);
  fData(eset_norm) <- fData(eset);
  
  return (eset_norm);
}

