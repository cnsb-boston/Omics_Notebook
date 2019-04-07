#-------------------------------------------------
#' Draw the Heatmaps
#'
#' This function generates heatmap plots based on data
#' 
#' @param eset an ExpressionSet object with omics data
#' @param type Type of data to be processed, or a name for the Omics data.
#' @param outputpath output file path for plots
#' @param limmaSig is list of significant genes
#'
#' @return a matrix of most variable features by MAD
#' 
#' @examples
#' 
#' @export
drawHeatmaps <- function(eset, emat_top=FALSE, type, outputpath=output_plots_path,
                        mapcolor=map_color, subset=FALSE, k_clust=0, 
                        outputcontrastpath=output_contrast_path, limmaSig=FALSE, 
                        show_row_names=FALSE, cluster_samples=TRUE){
  # Annotation Colors
  annotLab <- data.frame(Group = factor(pData(eset)$Group, levels=unique(pData(eset)$Group)));
  annotCol <- list(Group = rainbow(length(levels(factor(pData(eset)$Group)))))
  try({ if("ColorsHex" %in% colnames(pData(eset))) {
    if( checkColor(pData(eset)[,"ColorsHex"]) ){
      annotCol <- list( Group=unique(pData(eset)[,"ColorsHex"])[1:length(levels(factor(pData(eset)$Group)))] )
    }
  } })
  names(annotCol$Group) <- levels(factor(annotLab$Group))
  
  if("Group2" %in% colnames(pData(eset))){
    annotLab[,"Group2"] <- factor(pData(eset)$Group2, levels=unique(pData(eset)$Group2))
    annotCol <- append(annotCol, list(Group2=gray.colors(length(levels(factor(pData(eset)$Group2)))) ))
    names(annotCol$Group2) <- levels(factor(annotLab$Group2))
  }
  
  ha_column <- HeatmapAnnotation(df=annotLab, col=annotCol)
  
  if(mapcolor=="viridis"){mapcolor <- viridis(11); maponeway <- viridis(11);
  } else {mapcolor <- (rev(brewer.pal(11, mapcolor)));
          maponeway <- rev(brewer.pal(9, "Blues")) }
  
  # Standard Heatmaps
  if(class(limmaSig)=="logical"){
    output_filename <- file.path(outputpath, paste("heatmaps_",type,".pdf", sep=''));
    pdf(output_filename);
    
    # all features, z-score
    emat_sel <- na.omit(t(scale(t(exprs(eset))))) # Z-score across rows
    emat_sel[emat_sel < -2] <- -2
    emat_sel[emat_sel > 2] <- 2
    ht1 <- Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
                   cluster_columns=cluster_samples,
                   row_names_gp=gpar(fontsize=4),
                   column_title=paste(type, ": \nAll features, row z score", sep='') )
    print(ht1)
    #Optional k clustering
    if(k_clust !=0){
      kclus <- kmeans(emat_sel, 3);
      split <- paste0("Cluster ", kclus$cluster)
      ht1 <-Heatmap(matrix=(emat_sel), col=mapcolor, name="", top_annotation=ha_column,
                    cluster_columns=cluster_samples,
                    show_row_names=show_row_names,row_names_gp=gpar(fontsize=4), 
                    split=split,#km=k_clust,
                    column_title=paste(type, ": \nAll features, row z score", sep='') )
      print(ht1)
    }
    # Correlation
    #draw(Heatmap(matrix=cor(emat_sel), col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
    #             column_title=paste(type, ": Correlation, row z score", sep='') ))
  
    # Subset features only
    if( class(subset)!="logical" ){
      emat_sel <- exprs(eset)[subset,]
      emat_sel <- t(scale(t(emat_sel))) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      print(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,
                    cluster_columns=cluster_samples,
                    row_names_gp=gpar(fontsize=4),
                    column_title=paste(type, ": \nSelect genes, row z score", sep='') ) )
    }
  
    # variation filter, z score
    if ( class(emat_top) !="logical" ){
      emat_sel <- exprs(eset)[emat_top,]
      emat_sel <- na.omit(t(scale(t(emat_sel)))) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      print( Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
                     cluster_columns=cluster_samples,
                     column_title=paste(type, ": \nHighest variation, row z score", sep='') ) )
    }
  
    # all features, log2 intensity
    emat_sel <- na.omit(exprs(eset))
    print(Heatmap(matrix=emat_sel, col=maponeway, name="", top_annotation=ha_column,show_row_names=FALSE,
                  cluster_columns=cluster_samples,
                  column_title=paste(type, ": \nAll features, log2 Intensity", sep='') ))
  
    tmp<-dev.off();
    return(ht1);
    
  } else if(class(limmaSig)!="logical"){
  # Differential Heatmap
    output_filename <- file.path(outputcontrastpath, paste("heatmaps_",type,".pdf", sep=''));
    pdf(output_filename);

    # limma differential expression, z score
    emat_sel <- exprs(eset[rownames(eset) %in% limmaSig,])
    emat_sel <- na.omit(t(scale(t(emat_sel)))) # Z-score across rows
    emat_sel[emat_sel < -2] <- -2
    emat_sel[emat_sel > 2] <- 2 
    print(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,
                  cluster_columns=cluster_samples,
                  show_row_names=show_row_names,
                  row_names_gp=gpar(fontsize=4),
                  column_title=paste(type, ": \nDifferential Expression, row z score", sep='') ))
    tmp<-dev.off();
  }
}

