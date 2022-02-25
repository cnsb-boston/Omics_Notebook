#-------------------------------------------------
#' Make Interactive Heatmap
#'
#' This function generates an interactive heatmap
#' 
#' @param eset an ExpressionSet object with omics data
#' @param type Type of data to be processed, or a name for the Omics set.
#' @param outputpath output file path for plots
#' @param mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#'
#' @return
#' 
#' @examples
#' 
#' @import Biobase
#' @export
interactiveHeatmap <- function(eset, type, outputpath=output_plots_path,mapcolor=map_color ){
  
  annotCol <- c("red", "green", "blue", "yellow", "green", "purple",
                "brown", "black", "grey", "orange", "white", "light green", 
                "pink", "light blue" );
  sampleCols <- annotCol[1:length(levels(pData(eset)$Group))][pData(eset)$Group];
  
  if(mapcolor=="viridis"){mapcolor <- viridisLite::viridis(11)
  } else {mapcolor <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, mapcolor))) }
  
  #   emat.sel <- exprs(eset[rownames(eset) %in% rownames(limmaSig),])
  output_filename <- file.path(outputpath, paste("heatmap_all_",type,".html", sep=''));
  emat_sel <- exprs(eset)
  emat_sel <- stats::na.omit(t(scale(t(emat_sel)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  suppressWarnings(invisible(utils::capture.output(tmp<-heatmaply::heatmaply(emat_sel, scale='none', margins=c(100,100,40,20),
                                          col_side_colors=pData(eset)$Group,
                                          colors=mapcolor,file=output_filename ))) );
  rm(tmp)
}
