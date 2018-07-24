#-------------------------------------------------
#' QC and Normalization
#'
#' This function takes proteomics results and performs normalization and QC plotting to inspect data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param annotate An annotation data frame
#' @param type Type of data to be processed
#' @param norm one of the supported normalization methods
#' @param outputpath output file path for plots
#'
#' @return an expression set object, with normalized data
#' 
#' @examples
#'   if(isthereProteo) {eset_global <- QCandNorm(eset_global, norm=norm_method)}
#' 
#' @export
QCandNorm <- function(eset, norm, type, outputpath=output_plots_path,
                      annotate=annot ){
  
  col_palette <- rainbow(length(levels(pData(eset)$Group)))
  annotCol <- col_palette[pData(eset)$Group]
  
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
  }
  else if (norm=='loess'){
    eset_matrix_norm <- as.matrix(normalizeCyclicLoess(eset_matrix, method='pairs'))
  }
  else {  eset_matrix_norm<-eset_matrix }
  
  df <- data.frame(reshape2::melt(eset_matrix_norm, id.vars = NULL));
  colnames(df) <- c("Feature","Sample", "Intensity");
  plot3 <- ggplot(df, aes(x = Sample, y = Intensity)) + geom_boxplot(fill=annotCol) +
    theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    labs(title=paste(type, " log2 Intensity: After ",
                     toupper(norm)," Normalization", sep=''));
  plot4 <- ggplot(df, aes(x=Intensity, colour=Sample)) + geom_density()+ theme_bw() + 
    labs(title=paste(type, " log2 Intensity: After ",
                     toupper(norm)," Normalization", sep=''));
  
  output_filename<-file.path(outputpath, paste("boxplot_histogram_",type,".pdf", sep=''));
  pdf(output_filename);
  print(plot1);
  print(plot3);
  print(plot2);
  print(plot4);
  dev.off();
  
  grid.arrange(plot1+labs(title="Raw Data"), plot3+labs(title="Normalized"),
               plot2+theme(legend.position="none")+labs(title="Raw Data"), 
               plot4+theme(legend.position="none")+labs(title="Normalized"),
               top=paste(type, ":", sep=""), ncol=4);
  
  output_filename<-file.path(outputpath, paste("MAplots_",type,".pdf", sep=''))
  pdf(output_filename)
  for (i in 1:ncol(eset)) {
    limma::plotMA(eset_matrix, array=i, main=paste(type, "MA Plot",i,sep=" "))
    limma::plotMA(eset_matrix_norm, array=i, 
                  main=paste(type,toupper(norm),"Normalized MA Plot",i,sep=" "))
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

#-------------------------------------------------
#' Plot PCA
#'
#' This function generates PCA plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return a pca plot
#' 
#' @examples
#'     if(isthereProteo) { plotProteo <- plot.pca(eset_global) }
#' 
#' @export
plot.pca <- function(eset, x_axis="PC1", y_axis="PC2", type, outputpath=output_plots_path) {
  
  PC_data <- prcomp(t(exprs(eset)))
  percent_variance <- summary(PC_data)$importance["Proportion of Variance",]*100
  
  pca_graph <- ggplot(data=as.data.frame(PC_data$x), aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis], colour=pData(eset)$Group))+ 
    geom_point() + theme_bw() + theme(legend.title=element_blank()) + 
    labs(title=paste("Principal Component Analysis \nAcross All Features in All Samples - ",type,sep=""),
         x=paste(x_axis, sprintf(" (%2.0f%%)", percent_variance[x_axis]), sep=""),
         y=paste(y_axis, sprintf(" (%2.0f%%)", percent_variance[y_axis]), sep="")) +
    geom_text_repel(aes(x=PC_data$x[,x_axis], y=PC_data$x[,y_axis], label=rownames(PC_data$x)), colour="black")
  
  pca_loading_graph <- ggplot(data=as.data.frame(PC_data$rotation), aes(x=(PC_data$rotation[,x_axis]), y=PC_data$rotation[,y_axis])) +
    geom_text(label=rownames(PC_data$rotation), colour="black") + 
    labs(x=paste(x_axis, sep=""), y=paste(y_axis, sep=""),
         title=paste("Principal Component Factor Loadings\nAcross all Features and Samples - ",type, sep=""))
  
  output_filename <- file.path(outputpath, paste("PCAplots_",type,".pdf",sep=""))
  pdf(output_filename)
  print(pca_graph)
  print(pca_loading_graph)
  plot(PC_data, type = 'l', main=paste("Variances vs. Principal Components - ",type,sep=""))
  dev.off()
  print(pca_graph)
  
  return(pca_graph)
}

#-------------------------------------------------
#' Variation Plot and Filter
#'
#' This function generates plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return a matrix of most variable features by MAD
#' 
#' @examples
#' 
#' @export
variationPlotFilter <- function(eset, type, outputpath=output_plots_path, 
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
  
  output_filename<-file.path(outputpath, paste("variation_",type,".jpeg", sep=""))
  jpeg(output_filename)
  pairs(emat, pch=".")
  dev.off()
  return (emat_top)
}
#-------------------------------------------------
#' Draw the Heatmaps
#'
#' This function generates heatmap plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return a matrix of most variable features by MAD
#' 
#' @examples
#' 
#' @export
drawthemaps <- function(eset, emat_top=FALSE, type, outputpath=output_plots_path,
                        mapcolor=map_color){
  
  annotLab <- data.frame(Group = factor(pData(eset)$Group));
  annotCol <- list(Group = rainbow(length(levels(factor(pData(eset)$Group)))))
  names(annotCol$Group) <- levels(factor(annotLab$Group))
  ha_column <- HeatmapAnnotation(df=annotLab, col=annotCol)
  
  if(mapcolor=="viridis"){mapcolor <- viridis(11); maponeway <- viridis(11);
  } else {mapcolor <- (rev(brewer.pal(11, mapcolor)));
          maponeway <- rev(brewer.pal(9, "Blues")) }
  
  output_filename <- file.path(outputpath, paste("heatmaps_",type,".pdf", sep=''));
  pdf(output_filename);
  # all features, z-score
  emat_sel <- t(scale(t(exprs(eset)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  draw(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
               column_title=paste(type, ": All features, row z score", sep='') ))
  
  # all features, log2 intensity
  emat_sel <- exprs(eset)
  draw(Heatmap(matrix=emat_sel, col=maponeway, name="", top_annotation=ha_column,show_row_names=FALSE,
               column_title=paste(type, ": All proteins, log2 Intensity", sep='') ))
  
  # variation filter, z score
  if ( class(emat_top) =="matrix" ){
    emat_sel <- emat_top
    emat_sel <- t(scale(t(emat_sel))) # Z-score across rows
    emat_sel[emat_sel < -2] <- -2
    emat_sel[emat_sel > 2] <- 2
    draw(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
                 column_title=paste(type, ": Highest variation, row z score", sep='') ))
  }
  
  tmp<-dev.off();
}
#-------------------------------------------------
#' Make Interactive Heatmap
#'
#' This function generates an interactive heatmap
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#' @mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#'
#' @return
#' 
#' @examples
#' 
#' @export
makeInteractiveHM <- function(eset, type, outputpath=output_plots_path,mapcolor=map_color ){
  
  annotCol <- c("red", "green", "blue", "yellow", "green", "purple",
                "brown", "black", "grey", "orange", "white", "light green", 
                "pink", "light blue" );
  sampleCols <- annotCol[1:length(levels(pData(eset)$Group))][pData(eset)$Group];
  
  if(mapcolor=="viridis"){mapcolor <- viridis(11)
  } else {mapcolor <- colorRampPalette(rev(brewer.pal(11, mapcolor))) }
  
  #   emat.sel <- exprs(eset[rownames(eset) %in% rownames(limmaSig),])
  output_filename <- file.path(outputpath, paste("heatmap_all_",type,".html", sep=''));
  emat_sel <- exprs(eset)
  emat_sel <- t(scale(t(emat_sel))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  invisible(capture.output(tmp<-heatmaply(emat_sel, scale='none', margins=c(100,100,40,20),
                                          col_side_colors=pData(eset)$Group,
                                          colors=mapcolor,file=output_filename )));
  rm(tmp)
}
#-------------------------------------------------
#' Save Files
#'
#' This function saves txt files for subsequent use
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param outputpath output file path
#' @param workingdir specifies the working directory
#' @param saveTheEset specifies whether or not to save the eset as an RDS object
#'
#' @return
#' 
#' @examples
#' 
#' @export
saveTheFiles <- function(eset, type, outputpath=output_files_path, workingdir=working_dir,
                         saveTheEset=!newcontrastonly){
  if(saveTheEset){
    output_filename <- file.path(outputpath, paste("Data_eset_",type,".RDS", sep=""));             
    saveRDS(eset, file=output_filename); 
  }
  
  output_filename <- file.path(outputpath,paste("Expression_matrix_",type,".txt", sep=''));
  write.table(x=cbind(fData(eset)[,c("Gene", "Protein.names")], exprs(eset)), 
              file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
  
  output_filename <- file.path(outputpath,paste("Uniprot_IDs",type,".txt", sep=''));
  write.table(x=fData(eset)[,"Protein"], 
              file=output_filename, sep='\t',row.names=FALSE, col.names=FALSE, quote=FALSE);
}
#-------------------------------------------------
#' Write data to sheets
#'
#' This function outputs summary data to an excel sheet to share with collaborators
#' 
#' @param wb an openxlsx workbook object
#' @param eset an ExpressionSet object with proteomics data
#' @param type Type of data to be processed
#' @param data_format format of data to be processed
#' @param mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#'
#' @return wb with summary data formatted
#' 
#' @examples
#' 
#' @export
writeDataToSheets <- function(wb, eset, data_format, mapcolor=map_color, type){
  
  annotLab <- data.frame(Group = factor(pData(eset)$Group));
  annotCol <- list(Group = gsub('.{2}$','', rainbow(length(levels(pData(eset)$Group)))))
  sampleCols <- annotCol$Group[1:length(levels(pData(eset)$Group))][pData(eset)$Group];
  mapcolor <- rev(brewer.pal(7, "RdYlBu"))
  
  stName <- type
  addWorksheet(wb=wb, sheetName=stName)
  
  links <- fData(eset)$Link
  names(links) <- fData(eset)$Protein
  class(links) <- "hyperlink"
  
  emat_sel <- t(scale(t(exprs(eset)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  formatted_table <- data.frame(cbind(fData(eset)[,c("Gene", "Protein.names", "Protein")],emat_sel, 
                                      fData(eset)[,c("Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                                                     "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID")],
                                      exprs(eset) ))
  names <- make.unique(c("Gene", "Protein", "Uniprot", colnames(eset), "Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                         "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID", colnames(eset) ))
  
  # Add phospho site probabilities and data
  if(data_format=="Sites..MQ.") { 
    formatted_table <- data.frame(cbind(formatted_table, fData(eset)[,"Localization.prob"],
                                        fData(eset)[,grep("Probabilities", colnames(fData(eset)))],
                                        paste(fData(eset)[,"Amino.acid"],fData(eset)[,"Position"], sep=""),
                                        fData(eset)[,"Sequence.window"] ) )
    names <- c(names, "Localization Probability","Site Probabilities", "Amino Acid", "Peptide Sequence")
  }
  
  colnames(formatted_table) <- names
  
  writeDataTable(wb=wb, sheet=stName, x=formatted_table, xy=c("A",2), keepNA=FALSE, tableStyle="TableStyleLight1")
  writeData(wb=wb, sheet=stName, x=links, xy=c("C",3))
  
  conditionalFormatting(wb=wb, sheet=stName, type="colourScale", cols=4:(3+ncol(eset)), rows=3:(2+nrow(eset)), style=mapcolor[c(1,4,7)])
  
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"), rows=2, cols=c+3)
  }
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"), rows=2, cols=(c+3+9+ncol(eset)) )
  }
  
  mergeCells(wb=wb, sheet=stName, rows=1, cols=4:(3+ncol(eset)) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity, Row Z Score", xy=c(4,1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=4, stack=TRUE)
  
  mergeCells(wb=wb, sheet=stName, rows=1, cols=(13+ncol(eset)):(12+ncol(eset)+ncol(eset)) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity", xy=c((13+ncol(eset)),1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=(13+ncol(eset)), stack=TRUE)
  
  freezePane(wb=wb, sheet=stName, firstActiveRow=3, firstActiveCol=2)
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", fontColour='black'), rows=1:2, cols=1:(3+ncol(eset)+9+ncol(eset)),gridExpand=TRUE, stack=TRUE)   
  
  # Don't treat gene names as dates
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(numFmt="TEXT"), rows=3:(2+nrow(eset)), cols=1, gridExpand=TRUE )
  
  setRowHeights(wb=wb, sheet=stName, rows=2, heights=100)
  setColWidths(wb=wb, sheet=stName, cols=1:(3+ncol(eset)+9+ncol(eset)), widths=c(16,50,16,rep(4, ncol(eset)),50,50,50, 50,50,50, 8,8,8, rep(4, ncol(eset)) ))
}