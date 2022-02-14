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
#' @return a heatmap object
#' 
#' @examples
#' 
#' @import Biobase
#' @export
drawHeatmaps <- function(eset, emat_top=FALSE, type, title_add="", 
                         outputpath=output_plots_path,
                         mapcolor=map_color, subset=FALSE, k_clust=0, 
                         outputcontrastpath=output_contrast_path, limmaSig=FALSE, 
                         show_row_names=NULL, cluster_samples=TRUE){
  
  eset <- eset[,order(pData(eset)$Group)]
  
  # Annotation Colors
  annotLab <- data.frame(Group = factor(pData(eset)$Group, levels=unique(pData(eset)$Group)));
  annotCol <- list(Group = grDevices::rainbow(length(levels(factor(pData(eset)$Group)))))
  try({ if("ColorsHex" %in% colnames(pData(eset))) {
    if( checkColor(pData(eset)[,"ColorsHex"]) ){
      annotCol <- list( Group=unique(pData(eset)[,"ColorsHex"])[1:length(levels(factor(pData(eset)$Group)))] )
    }
  } })
  names(annotCol$Group) <- levels(factor(annotLab$Group))
  
  if("Group2" %in% colnames(pData(eset))){
    annotLab[,"Group2"] <- factor(pData(eset)$Group2, levels=unique(pData(eset)$Group2))
    annotCol <- append(annotCol, list(Group2=grDevices::gray.colors(length(levels(factor(pData(eset)$Group2)))) ))
    names(annotCol$Group2) <- levels(factor(annotLab$Group2))
  }
  
  ComplexHeatmap::ht_opt("message" = FALSE)
  ha_column <- ComplexHeatmap::HeatmapAnnotation(df=annotLab, col=annotCol)
  
  if(mapcolor=="viridis"){mapcolor <- viridisLite::viridis(11); maponeway <- viridisLite::viridis(11);
  } else {mapcolor <- (rev(RColorBrewer::brewer.pal(11, mapcolor)));
          maponeway <- rev(RColorBrewer::brewer.pal(9, "Blues")) }

  isBigMap=function(emat){
    ncol(emat)>50 || nrow(emat)>100
  }

  getHeatmap=function(`matrix`,name="Z-score",show_row_names.=show_row_names,column_title,
                      cluster_columns=cluster_samples,split=NULL,title_names="",col=mapcolor) {
    ComplexHeatmap::Heatmap (`matrix`=`matrix`, col=col, name=name, top_annotation=ha_column,
            show_row_names=if(is.null(show_row_names.)) !isBigMap(`matrix`) else show_row_names.,
            row_labels=substring(rownames(`matrix`),0,50),
            cluster_columns=cluster_columns, use_raster=isBigMap(`matrix`), split=split,
            row_names_gp=grid::gpar(fontsize=4),
            column_title=paste(type, ": ", title_names, "\n", column_title, sep='') )
  }
  
  # Standard Heatmaps
  if(class(limmaSig)=="logical"){
    output_filename <- file.path(outputpath, paste("heatmaps_",type,".pdf", sep=''));
    pdf(output_filename);
    
    # all features, z-score
    emat_sel <- stats::na.omit(t(scale(t(exprs(eset))))) # Z-score across rows
    emat_sel[emat_sel < -2] <- -2
    emat_sel[emat_sel > 2] <- 2
    ht1 <- getHeatmap(matrix=emat_sel, show_row_names=FALSE, cluster_columns=cluster_samples,
                      column_title="All features, row z score")
    print(ht1)
    
    print(getHeatmap(matrix=emat_sel, show_row_names=FALSE, cluster_columns=FALSE,
                     column_title="All features, row z score"))
    
    #Optional k clustering
    if(k_clust !=0){
      kclus <- stats::kmeans(emat_sel, 3);
      split <- paste0("Cluster ", kclus$cluster)
      ht1 <-getHeatmap(matrix=emat_sel, cluster_columns=cluster_samples,
                       split=split,#km=k_clust,
                       column_title="All features, row z score")
      print(ht1)
    }
    # Correlation
    #draw(Heatmap(matrix=cor(emat_sel), col=mapcolor, name="Cor", top_annotation=ha_column,show_row_names=FALSE,
    #             column_title=paste(type, ": Correlation, row z score", sep='') ))
  
    # Subset features only
    if( class(subset)!="logical" ){ 
      for( k in 1:length(subset) ){ suppressWarnings({ try({
        emat_sel <- exprs(eset)[subset[[k]],]
        emat_sel <- na.omit(t(scale(t(emat_sel)))) # Z-score across rows
        emat_sel[emat_sel < -2] <- -2
        emat_sel[emat_sel > 2] <- 2
        print(getHeatmap(matrix=emat_sel, cluster_columns=TRUE,
                         column_title="Subset, row z score", title_names=names(subset)[k]))
        print(getHeatmap(matrix=emat_sel, cluster_columns=FALSE,
                         column_title="Subset, row z score", title_names=names(subset)[k]))
        # print(Heatmap(matrix=emat_sel, col=mapcolor, name="Z-score", top_annotation=ha_column,
        #               cluster_columns=TRUE,show_row_names=show_row_names,cluster_rows=F,
        #               row_names_gp=grid::gpar(fontsize=4),
        #               column_title=paste(type, ": ",names(subset)[k],"\n Subset, row z score", sep='') ) )
        # print(Heatmap(matrix=emat_sel, col=mapcolor, name="Z-score", top_annotation=ha_column,
        #               cluster_columns=FALSE,show_row_names=show_row_names,cluster_rows=F,
        #               row_names_gp=grid::gpar(fontsize=4),
        #               column_title=paste(type, ": ",names(subset)[k],"\n Subset, row z score", sep='') ) )
        
      }, silent=TRUE) }) }
    }
  
    # variation filter, z score
    if ( class(emat_top) !="logical" ){
      emat_sel <- exprs(eset)[emat_top,]
      emat_sel <- stats::na.omit(t(scale(t(emat_sel)))) # Z-score across rows
      emat_sel[emat_sel < -2] <- -2
      emat_sel[emat_sel > 2] <- 2
      print( getHeatmap(matrix=emat_sel, show_row_names=FALSE, cluster_columns=TRUE,
                        column_title="Highest variation, row z score"))
      print( getHeatmap(matrix=emat_sel, show_row_names=FALSE, cluster_columns=FALSE,
                        column_title="Highest variation, row z score"))
    }
  
    # all features, log2 intensity
    emat_sel <- stats::na.omit(exprs(eset))
    print(getHeatmap(matrix=emat_sel, col=maponeway, name="Value", show_row_names=FALSE,
                     cluster_columns=TRUE, column_title="All features, log2 Value"))
    print(getHeatmap(matrix=emat_sel, col=maponeway, name="Value", show_row_names=FALSE,
                     cluster_columns=FALSE, column_title="All features, log2 Value"))
  
    tmp<-dev.off();
    return(ht1);
    
  } else if(class(limmaSig)!="logical"){
  # Differential Heatmap
    output_filename <- file.path(outputcontrastpath, paste(type,"_heatmaps",".pdf", sep=''));
    pdf(output_filename);

    # limma differential expression, z score
    emat_sel <- exprs(eset[rownames(eset) %in% limmaSig,])
    emat_sel <- stats::na.omit(t(scale(t(emat_sel)))) # Z-score across rows
    emat_sel[emat_sel < -2] <- -2
    emat_sel[emat_sel > 2] <- 2 
    print(getHeatmap(matrix=emat_sel, cluster_columns=TRUE,
          title_names=title_add, column_title="Differential Features, z score"))
    print(getHeatmap(matrix=emat_sel, cluster_columns=FALSE,
          title_names=title_add, column_title="Differential Features, z score"))
    tmp<-dev.off();
  }
}

