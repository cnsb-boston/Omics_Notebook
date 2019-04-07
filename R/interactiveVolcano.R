#-------------------------------------------------
#' Make Interactive Volcano
#'
#' This function makes a volcano plot
#' 
#' @param eset is an expression set object
#' @param fit is the output of lmFit or eBayes
#' @param limmaSig is the output of topTable
#' @param dt is the result of decideTests()
#' @param type Type of data to be processed, or a name for the Omics analysis
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#' 
#' @export
interactiveVolcano <- function(eset, fit=NULL, limmaSig, dt=NULL, type, outputcontrastpath=output_contrast_path, col ){
  
  annotCol <- c("red", "green", "blue", "yellow", "green", "purple",
                "brown", "black", "grey", "orange", "white", "light green", 
                "pink", "light blue" );
  sampleCols <- annotCol[1:length(levels(pData(eset)$Group))][pData(eset)$Group];
  
  #   if(!newcontrast){
  #     glMDSPlot(exprs(eset), groups=pData(eset)$Group, html=paste("MDS-Plot_",type,sep=""),
  #            folder="InteractivePlots", path=outputpath, launch=FALSE  );
  #   }
  #   glMDPlot(fit, groups=pData(eset)$Group, status=dt, counts=exprs(eset),
  #            anno=fit$genes[,c("Gene", "Protein", "Protein.names","Uniprot")],
  #            xlab="Normalized log2 Intensity",
  #            display.columns=c("Gene", "Protein", "Protein.names", "Uniprot"),
  #            side.main="Gene", side.ylab="Normalized Intensity",
  #            main=colnames(dt), sample.cols=sampleCols,
  #            html=paste("MD-Plot_",type,sep=""),
  #            folder="InteractivePlots", path=outputcontrastpath, launch=FALSE );

  annot_columns <- "feature_identifier";
  if("Gene" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Gene"); }
  if("Protein" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Protein"); }
  if("Protein.names" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Protein.names"); }
  if("Uniprot" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "Uniprot"); }
  if("mz" %in% colnames(fData(eset))) { annot_columns <- c(annot_columns, "mz"); }
  if(length(annot_columns)==1) { annot_columns <- c(annot_columns, "feature_identifier"); }
  
  if(class(fit)!="NULL"){
    glXYPlot(x=fit$coef[,col], y=fit$lod[,col], xlab="logFC", status=dt[,col], counts=exprs(eset),
           groups=pData(eset)$Group, ylab="logOdds", #sample.cols=pData(eset)$Group,
           side.ylab="Normalized Intensity",
           anno=fData(eset)[,annot_columns], side.main='feature_identifier',
           html=paste("Volcano-Plot_",type, sep=""),
           folder="InteractivePlots", path=outputcontrastpath, launch=FALSE  );
  } else {
    glXYPlot(x=limmaSig$Mean, y=limmaSig$logFC, xlab="Mean",ylab="log2FC", #status=dt[,col], 
             counts=exprs(eset),
             groups=pData(eset)$Group, #sample.cols=pData(eset)$Group,
             side.ylab="Normalized Intensity",
             anno=fData(eset)[,annot_columns], side.main='feature_identifier',
             html=paste("Volcano-Plot_",type, sep=""),
             folder="InteractivePlots", path=outputcontrastpath, launch=FALSE  );
  }
  
}


