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
  
  if(mapcolor=="viridis"){mapcolor <- viridisLite::viridis(3)
  } else {mapcolor <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(3, mapcolor))) }
  
  #   emat.sel <- exprs(eset[rownames(eset) %in% rownames(limmaSig),])
  emat_sel <- exprs(eset)
  emat_sel <- stats::na.omit(t(scale(t(emat_sel)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2

  jlist=list(rows=jsonlite::unbox(nrow(emat_sel)),
             columns=jsonlite::unbox(ncol(emat_sel)),
             seriesArrays=list(as.matrix(emat_sel)),
             seriesNames=type,
             seriesDataTypes="Float32",
             rowMetadataModel=list(vectors=list(list(name=jsonlite::unbox("id"),array=rownames(emat_sel)))),
             columnMetadataModel=list(vectors=list(list(name=jsonlite::unbox("id"),array=colnames(emat_sel)),
                                                   list(name=jsonlite::unbox("group"),array=pData(eset)$Group))))

  slist=function(...){
    pars=list(...)
    ret=list()
    for(i in 1:length(pars)){
      ret[[names(pars)[i]]]=jsonlite::unbox(pars[[i]])
    }
    ret
  }

  colanno=list(slist(field="id",display="text"),slist(field="group",display="color"))

  colorscheme=list(null=slist(min=0,max=1,missingcolor="#c0c0c0",scalingMode=0,stepped=F,transformValues=0))
  colorscheme$null$fractions=c(0,.5,1)
  # remove alpha channel, which breaks Morpheus rendering
  colorscheme$null$colors=sub("^(#[0-9a-fA-F]{6}).*","\\1",mapcolor)

  output_filename <- file.path(outputpath, paste("heatmap_all_",type,".json", sep=''));
  j = jsonlite::toJSON(list(dataset=jlist,
                            columns=colanno,
                            name=jsonlite::unbox(type),
                            colorScheme=list(valueToColorScheme=colorscheme)),matrix="rowmajor", pretty=2)
  writeLines(j,output_filename)

  output_filename
}
