#-------------------------------------------------
#' Make Volcano
#'
#' This function makes a volcano plot
#' 
#' @param dat is the result of topTable()
#' @param dt is the result of decideTests()
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#'  
#' 
#' @export
makeVolcano <- function(dat, type, outputpath=output_contrast_path){ 
  threshold_fc_high <- sort(dat$logFC, decreasing=TRUE)[12]
  threshold_fc_low <- sort(dat$logFC, decreasing=FALSE)[12]
  threshold_p_low <- sort(dat$adj.P.Val, decreasing=FALSE)[20]
  
  plot1 <- ggplot(data=data.frame(dat)) + geom_point(aes(x=logFC, y=-log10(adj.P.Val),colour=(dat$adj.P.Val<0.05))) +
    scale_colour_manual(name="adj.P.Val<0.05", values=setNames(c("red", "grey"), c(T,F))) +
    labs(title=paste(type,"Volcano Plot", sep=" ")) + theme_bw() + 
    geom_text_repel(data=data.frame(subset(dat, logFC>threshold_fc_high)), aes(x=logFC, y=-log10(adj.P.Val), label=Gene)) + 
    geom_text_repel(data=data.frame(subset(dat, logFC<threshold_fc_low)), aes(x=logFC, y=-log10(adj.P.Val), label=Gene)) + 
    geom_text_repel(data=data.frame(subset(dat, adj.P.Val<threshold_p_low)), aes(x=logFC, y=-log10(adj.P.Val), label=Gene))
  
  output_filename <- file.path(outputpath, paste("volcano_",type,".pdf", sep=''));
  pdf(output_filename);
  plot(plot1)
  tmp<-dev.off();
}

#-------------------------------------------------
#' Make Interactive Volcano
#'
#' This function makes a volcano plot
#' 
#' @param eset is an expression set object
#' @param dt is the result of decideTests()
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#'  
#' 
#' @export
makeInteractiveV <- function(eset, fit, limmaSig, dt, type, outputcontrastpath=output_contrast_path, col ){
  
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
  
  glXYPlot(x=fit$coef[,col], y=fit$lod[,col], xlab="logFC", status=dt[,col], counts=exprs(eset),
           groups=pData(eset)$Group, ylab="logOdds", #sample.cols=pData(eset)$Group,
           side.ylab="Normalized Intensity",
           anno=fit$genes[,c("Gene", "Protein", "Protein.names", "Uniprot")], side.main='Gene',
           html=paste("Volcano-Plot_",type, sep=""),
           folder="InteractivePlots", path=outputcontrastpath, launch=FALSE  );
}

#-------------------------------------------------
#' Wrapper fnction for EnrichR
#'
#' This runs enrichment analysis with enrichr
#' 
#' @param genes, list of genes
#' @param type Type of data to be processed
#' @param search_data list of enrichr databases to search
#' @param outputpath output file path for plots
#'
#' @return 
#' 
#' @examples
#' 
#' @export
runEnrichR <- function (genes, type, search_dat=search_databases, 
                        outputpath=output_contrast_path) {
  
  genes_up <- genes[genes[,"logFC"]>0,"Gene"]
  genes_down <- genes[genes[,"logFC"]<0,"Gene"]
  
  enriched_up <- enrichr(unique(genes_up), databases = search_dat);
  enriched_up <- bind_rows(enriched_up, .id="databases")
  enriched_up <- enriched_up[order(enriched_up[,"Adjusted.P.value"]),]
  
  enriched_down <- enrichr(unique(genes_down), databases = search_dat);
  enriched_down <- bind_rows(enriched_down, .id="databases")
  enriched_down <- enriched_down[order(enriched_down[,"Adjusted.P.value"]),]
  
  output_filename <- file.path(outputpath,paste("enrichr_results_", type,"_up.txt", sep='')); 
  write.table(enriched_up, file=output_filename, sep="\t");
  output_filename <- file.path(outputpath,paste("enrichr_results_", type,"_down.txt", sep='')); 
  write.table(enriched_down, file=output_filename, sep="\t");
}

#-------------------------------------------------
#' Draw the Heatmaps - Differential
#'
#' This function generates heatmap plots based on data
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param limmaRes the results of topTable()
#' @param dt the result of decideTests
#' @param type Type of data to be processed
#' @param outputpath output file path for plots
#'
#' @return a matrix of most variable features by MAD
#' 
#' @examples
#' 
#' @export
drawthemapsDiff <- function(eset, limmaSig, type, outputcontrastpath=output_contrast_path,
                            mapcolor=map_color){
  
  annotLab <- data.frame(Group = factor(pData(eset)$Group));
  annotCol <- list(Group = rainbow(length(levels(factor(pData(eset)$Group)))));
  names(annotCol$Group) <- levels(annotLab$Group)
  ha_column <- HeatmapAnnotation(df=annotLab, col=annotCol)
  
  if(mapcolor=="viridis"){mapcolor <- viridis(11); maponeway <- viridis(11);
  } else {mapcolor <- (rev(brewer.pal(11, mapcolor))) }
  
  output_filename <- file.path(outputcontrastpath, paste("heatmaps_",type,".pdf", sep=''));
  pdf(output_filename);
  
  # limma differential expression, z score
  emat_sel <- exprs(eset[rownames(eset) %in% rownames(limmaSig),])
  emat_sel <- t(scale(t(emat_sel))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2 
  draw(Heatmap(matrix=emat_sel, col=mapcolor, name="", top_annotation=ha_column,show_row_names=FALSE,
               column_title=paste(type, ": Differential Expression, row z score", sep='') ))
  tmp<-dev.off();
}

#-------------------------------------------------
#' Save Files Differential
#'
#' This function saves txt files relevant to the differential analysis for subsequent use
#' 
#' @param eset an ExpressionSet object with proteomics data
#' @param limmaRes topTable output
#' @param type Type of data to be processed
#' @param contrast_name name of the contrast
#' @param outputcontrastpath output file path
#'
#' @return
#' 
#' @examples
#' 
#' @export
saveTheFilesDiff <- function(eset, limmaRes, type, contrast_name, outputcontrastpath=output_contrast_path){
  output_filename <- file.path(outputcontrastpath,paste("GSEA_",type, gsub("-","_",contrast_name),".rnk", sep=''));
  ranked <- cbind(limmaRes[,"Gene"],sign(limmaRes[,"logFC"]) * -log10(limmaRes[,"adj.P.Val"]))
  colnames(ranked)<-c("GeneName", "rank")
  colnames(ranked)<-c("GeneName", "rank")
  ranked <- ranked[order(as.numeric(ranked[,"rank"]), decreasing=TRUE),]
  ranked <- ranked[ranked[,"GeneName"]!="", ]
  ranked <- ranked[!duplicated(ranked[,"GeneName"]),]
  write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
}


#-------------------------------------------------
#' Write data to sheets Differential
#'
#' This function outputs summary data to an excel sheet to share with collaborators
#' 
#' @param wb an openxlsx workbook object
#' @param eset an ExpressionSet object with proteomics data
#' @param limmaFit a lmFit object
#' @param type Type of data to be processed
#' @param data_format format of data to be processed
#' @param mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#' @param contrastnames specifies names of the contrasts/coefficients in the lmFit object
#'
#' @return wb with summary data formatted
#' 
#' @examples
#' 
#' @export
#' 
writeDataToSheetsDiff <- function(wb, eset, limmaFit, data_format, mapcolor=map_color, type, contrastnames=contrastgroups){
  
  # Create table for each coefficient
  DiffList <- vector("list", (ncol(limmaFit$coefficients)) )
  for(i in 1:length(DiffList) ) {
    DiffList[[i]] = topTable(limmaFit, adjust.method="BH", n=Inf, sort.by='p', coef=i)
  }
  
  # Match the row order to the first contrast
  eset <- eset[match(rownames(DiffList[[1]]),rownames(eset)),]
  
  # Make colors for sample names
  annotLab <- data.frame(Group = factor(pData(eset)$Group));
  annotCol <- list(Group = gsub('.{2}$','', rainbow(length(levels(factor(pData(eset)$Group)))) ) )
  sampleCols <- annotCol$Group[1:length(levels(factor(pData(eset)$Group)))][factor(pData(eset)$Group)];
  mapcolor <- rev(brewer.pal(7, "RdYlBu"))
  
  # Create new sheet to add to the workbook
  stName <- as.character(type)
  addWorksheet(wb=wb, sheetName=stName)
  
  # Format Uniprot hyperlinks
  links <- fData(eset)$Link
  names(links) <- fData(eset)$Protein
  class(links) <- "hyperlink"
  
  # Make row z-score values for "heatmap"
  emat_sel <- t(scale(t(exprs(eset)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  
  # Create the data table
  formatted_table<-""
  formatted_table <- cbind(fData(eset)[,c("Gene", "Protein.names", "Protein")],emat_sel,stringsAsFactors=FALSE )
  for (i in 1:length(DiffList)) {
    DiffList[[i]] <- DiffList[[i]][match(rownames(DiffList[[1]]),rownames(DiffList[[i]])),]
    formatted_table <- cbind(formatted_table, DiffList[[i]][,c("P.Value","adj.P.Val","logFC")], stringsAsFactors=FALSE)
  }
  formatted_table<- data.frame(cbind(formatted_table, fData(eset)[,c("Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                                                     "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID")],
                                      exprs(eset), stringsAsFactors=FALSE ), stringsAsFactors=FALSE)
  # Create column names
  names <- make.unique(toupper(c("Gene", "Protein", "Uniprot", colnames(eset), rep(c("P.Value","adj.P.Val","logFC"),length(DiffList)),
                         "Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
                         "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID", colnames(eset) )))
  
  # Add phospho site probabilities and data
  if("Sites..MQ." %in% data_format) { 
    formatted_table <- data.frame(cbind(formatted_table, fData(eset)[,"Localization.prob"],
                                        fData(eset)[,grep("Probabilities", colnames(fData(eset)))],
                                        paste(fData(eset)[,"Amino.acid"],fData(eset)[,"Position"], sep=""),
                                        fData(eset)[,"Sequence.window"] ) )
    names <- c(names, "Localization Probability","Site Probabilities", "Amino Acid", "Peptide Sequence")
  }
  
  colnames(formatted_table) <-  names
  
  # Write data table to sheet
  writeDataTable(wb=wb, sheet=stName, x=formatted_table, xy=c("A",2), keepNA=FALSE, tableStyle="TableStyleLight1")
  writeData(wb=wb, sheet=stName, x=links, xy=c("C",3))

  # Add heatmap color
  conditionalFormatting(wb=wb, sheet=stName, type="colourScale", cols=4:(3+ncol(eset)), rows=3:(2+nrow(eset)), style=mapcolor[c(1,4,7)])

  # Add color bars for fold change
  for(i in 1:length(DiffList)){
    conditionalFormatting(wb=wb, sheet=stName, cols= ((3*i)+ncol(exprs(eset))+3), rows=3:(3+nrow(eset)), type="colourScale", style=mapcolor[c(1,4,7)])
  }
  
  
  # Rotate text for sample names
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"),
             rows=2, cols=c+3)
  }
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"),
             rows=2, cols=(c+3+(3*length(DiffList))+9+ncol(eset)) )
  }
  
  # Merge cells and add Intensity title
  mergeCells(wb=wb, sheet=stName, rows=1, cols=4:(3+ncol(eset)) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity, Row Z Score", xy=c(4,1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=4, stack=TRUE)

  # Add contrast titles
  for(i in 1:length(DiffList)){
    startCol<- 3+ ncol(eset) + 1 + (3*(i-1))
    mergeCells(wb=wb, sheet=stName, rows=1, cols=startCol:(startCol+2) )
    writeData(wb=wb, sheet=stName, x=colnames(limmaFit$coefficients)[i], xy=c(startCol,1))
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=startCol, stack=TRUE)
  }
  
  # Add final intensity title
  mergeCells(wb=wb, sheet=stName, rows=1, cols=(13+(4*length(DiffList))+ncol(eset)):(12+(4*length(DiffList))+(2*ncol(eset))) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity", xy=c((11+(4*length(DiffList))+ncol(eset)),1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=(13+(4*length(DiffList))+ncol(eset)), stack=TRUE)
  
  # Freeze columns/rows and bold column titles
  freezePane(wb=wb, sheet=stName, firstActiveRow=3, firstActiveCol=2)
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", fontColour='black'),
           rows=1:2, cols=1:length(names),gridExpand=TRUE, stack=TRUE)   
  
  # Don't treat gene names as dates
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(numFmt="TEXT"), rows=3:(2+nrow(eset)), cols=1, gridExpand=TRUE )
  
  # Set row heights and column widths
  setRowHeights(wb=wb, sheet=stName, rows=2, heights=100)
  setColWidths(wb=wb, sheet=stName, cols=1:(3+ncol(eset)+9+(3*length(DiffList))+ncol(eset)), widths=c(16,50,16,rep(4, ncol(eset)), rep(c(8,8,8),length(DiffList)), 50,50,50, 50,50,50, 8,8,8, rep(4, ncol(eset)) ))
}

#writeDataToSheetsDiff <- function(wb, eset, limmaRes, mapcolor=map_color, type, contrastnames=contrastgroups){
#  
#  eset <- eset[match(rownames(limmaRes),rownames(eset)),]
#  
#  annotLab <- data.frame(Group = factor(pData(eset)$Group));
#  annotCol <- list(Group = gsub('.{2}$','', rainbow(length(levels(pData(eset)$Group)))))
#  sampleCols <- annotCol$Group[1:length(levels(pData(eset)$Group))][pData(eset)$Group];
#  mapcolor <- rev(brewer.pal(7, "RdYlBu"))
#  
#  stName <- as.character(type)
#  addWorksheet(wb=wb, sheetName=stName)
#  
#  links <- fData(eset)$Link
#  names(links) <- fData(eset)$Protein
#  class(links) <- "hyperlink"
  
#  emat_sel <- t(scale(t(exprs(eset)))) # Z-score across rows
#  emat_sel[emat_sel < -2] <- -2
#  emat_sel[emat_sel > 2] <- 2
#  formatted_table <- data.frame(cbind(fData(eset)[,c("Gene", "Protein.names", "Protein")],emat_sel,
#                                      (limmaRes[,"adj.P.Val"]),(limmaRes[,"logFC"]),(limmaRes[,"AveExpr"]),(limmaRes[,"logFC"]),
#                                      fData(eset)[,c("Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
#                                                     "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID")],
#                                      exprs(eset), stringsAsFactors=FALSE ), stringsAsFactors=FALSE)
#  names <- make.unique(c("Gene", "Protein", "Uniprot", colnames(eset), "adj.P.Val", "logFC", "AveExpr", "Log Fold Change",
#                         "Uniprot_Function", "Uniprot_Cellular_Location", "Uniprot_Disease",
#                         "GO_biological_process", "GO_molecular_function", "GO_cellular_component", "GO_ID", "ReactomeID", "KEGG_ID", colnames(eset) ))
#  if(type=="phospho") { formatted_table <- data.frame(cbind(formatted_table, fData(eset)[,"Phospho..STY..Probabilities"]))
#                        names <- c(names, "Phospho..STY..Probabilities")}
#  colnames(formatted_table) <- names
#  
#  writeDataTable(wb=wb, sheet=stName, x=formatted_table, xy=c("A",2), keepNA=FALSE, tableStyle="TableStyleLight1")
#  writeData(wb=wb, sheet=stName, x=links, xy=c("C",3))
#  
#  conditionalFormatting(wb=wb, sheet=stName, type="colourScale", cols=4:(3+ncol(eset)), rows=3:(2+nrow(eset)), style=mapcolor[c(1,4,7)])
#  conditionalFormatting(wb=wb, sheet=stName, cols= (4+ncol(exprs(eset))+3), rows=3:(3+nrow(limmaRes)), type="databar", style=c("blue", "red"), showvalue=FALSE)
#  
#  for (c in 1:ncol(eset)){
#    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"), rows=2, cols=c+3)
#  }
#  for (c in 1:ncol(eset)){
#    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"), rows=2, cols=(c+3+4+9+ncol(eset)) )
#  }
#  
#  mergeCells(wb=wb, sheet=stName, rows=1, cols=4:(3+ncol(eset)) )
#  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity, Row Z Score", xy=c(4,1))
#  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=4, stack=TRUE)
#  
#  mergeCells(wb=wb, sheet=stName, rows=1, cols=(17+ncol(eset)):(16+ncol(eset)+ncol(eset)) )
#  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity", xy=c((17+ncol(eset)),1))
#  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=(17+ncol(eset)), stack=TRUE)
#  
#  freezePane(wb=wb, sheet=stName, firstActiveRow=3, firstActiveCol=2)
#  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", fontColour='black'),
#           rows=1:2, cols=1:(3+ncol(eset)+4+9+ncol(eset)),gridExpand=TRUE, stack=TRUE)   
#  
#  setRowHeights(wb=wb, sheet=stName, rows=2, heights=100)
#  setColWidths(wb=wb, sheet=stName, cols=1:(3+ncol(eset)+9+4+ncol(eset)), widths=c(16,50,16,rep(4, ncol(eset)), 8,8,8,16, 50,50,50, 50,50,50, 8,8,8, rep(4, ncol(eset)) ))
#}
