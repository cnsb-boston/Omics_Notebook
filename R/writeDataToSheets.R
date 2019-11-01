#-------------------------------------------------
#' Write data to sheets
#'
#' This function outputs summary data to an excel sheet to share with collaborators
#' 
#' @param wb an openxlsx workbook object
#' @param eset an ExpressionSet object with omics data
#' @param type Type of data to be processed, or a name for the Omics data
#' @param data_format format of data to be processed. See Omics_Notebook example for options.
#' @param mapcolor specifies color scale to use for heatmap: "viridis" "RdBu" "RdYlBu"
#' @param limmaFit a lmFit object or FALSE if no differential analysis
#'
#' @return wb with summary data formatted
#' 
#' @examples
#' 
#' @export
writeDataToSheets <- function(wb, eset, limmaFit=NULL, data_format, mapcolor=map_color, type,
                              coef_index=NULL, time_index=NULL, contrast_strings=NULL){
  
  eset <- eset[,order(pData(eset)$Group)]
  
  # Make colors for sample names
  annotLab <- data.frame(Group = factor(pData(eset)$Group));
  annotCol <- list(Group = gsub('.{2}$','', rainbow(length(levels(factor(pData(eset)$Group)))) ) )
  sampleCols <- annotCol$Group[1:length(levels(factor(pData(eset)$Group)))][factor(pData(eset)$Group)];
  mapcolor <- rev(brewer.pal(7, "RdYlBu"))
  
  if (class(limmaFit)!="NULL"){
    # Create table for each coefficient
    if(class(coef_index)!="NULL"){
      DiffList <- vector("list", length(coef_index) )
    } else {
      DiffList <- vector("list", (ncol(limmaFit$coefficients)) )
      coef_index <- 1:(ncol(limmaFit$coefficients))
    }
    if(class(coef_index)!="NULL" & class(time_index)!="NULL"){
      if(time_index>0){
        DiffList_T <- topTable(limmaFit, adjust.method="BH", n=Inf, coef=coef_index)
      } else {
        DiffList_T <- NULL
      }
    }
    if(length(DiffList)>1){
      DiffList_F <- topTable(limmaFit, adjust.method="BH", n=Inf)
    }
    for(i in 1:length(coef_index) ) {
      DiffList[[i]] = topTable(limmaFit, adjust.method="BH", n=Inf, sort.by='p', coef=coef_index[i])[,c("P.Value","adj.P.Val","logFC")]
    }
    # Match the row order to the first contrast
    eset <- eset[match(rownames(DiffList[[1]]),rownames(eset)),]
  } else { DiffList <- vector("list", 0 ) }
  
  # Create new sheet to add to the workbook
  stName <- as.character(type)
  addWorksheet(wb=wb, sheetName=stName)
  
  # Format Uniprot hyperlinks
  if ("Link" %in% colnames(fData(eset)) ){
    links <- fData(eset)$Link
    names(links) <- fData(eset)$Protein
    class(links) <- "hyperlink"
  }
  
  # Make row z-score values for "heatmap"
  emat_sel <- t(scale(t(exprs(eset)))) # Z-score across rows
  emat_sel[emat_sel < -2] <- -2
  emat_sel[emat_sel > 2] <- 2
  
  # Create the data table
  formatted_table<-""
  if ( "feature_identifier" %in% colnames(fData(eset)) ) {  
    formatted_table <- data.frame( fData(eset)[,c("feature_identifier")],emat_sel )
    names <- c("Feature", colnames(exprs(eset))  )
  } else if ( "Gene" %in% colnames(fData(eset)) ) {  
    formatted_table <- data.frame( fData(eset)[,c("Gene")],emat_sel )
    names <- c("Gene", colnames(exprs(eset))  )
  } else { 
    formatted_table <- data.frame(  rownames(exprs(eset)),emat_sel )
    names <- c("Feature", colnames(exprs(eset))  )
  }
  col_widths <- c(20, rep(4, ncol(eset)));
  
  if ( "Gene" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Gene")] );
    names <- c(names, "Gene");
    col_widths<-c(col_widths, 16);
  }
  if ( "Protein.names" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Protein.names")] );
    names <- c(names, "Protein");
    col_widths<-c(col_widths, 50);
  }
  if ( "Protein" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Protein")] );
    names <- c(names, "Uniprot"); 
    col_widths<-c(col_widths, 16);
  }
  
  if ( "mz" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("mz")] );
    names <- c(names, "MZ"); 
    col_widths<-c(col_widths, 16);
  }
  
  if ( "rt" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("rt")] );
    names <- c(names, "RT"); 
    col_widths<-c(col_widths, 16);
  }
  
  if ( "MaxFC" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("MaxFC")] );
    names <- c(names, "MaxFC"); 
    col_widths<-c(col_widths, 8);
  }
  
  col_index <- length(names)
  
  if(class(limmaFit)=="MArrayLM") { 
    for (i in 1:length(DiffList)) {
      DiffList[[i]] <- DiffList[[i]][match(rownames(DiffList[[1]]),rownames(DiffList[[i]])),];
      formatted_table <- data.frame(formatted_table, DiffList[[i]][,c("P.Value","adj.P.Val","logFC")] );
      names <- c(names, "P.Value","adj.P.Val","logFC");
      col_widths<-c(col_widths, 8,8,8);
    } 
    if(length(DiffList)>1){ try({
      DiffList_F <- DiffList_F[match(rownames(DiffList[[1]]),rownames(DiffList_F)),];
      formatted_table <- data.frame(formatted_table, DiffList_F[,"F"])
      names <- c(names, "F-Statistic");
      col_widths<-c(col_widths, 8);
    }) }
    try({ if(class(DiffList_T)!="NULL"){
      DiffList_T <- DiffList_T[match(rownames(DiffList[[1]]),rownames(DiffList_T)),];
      formatted_table <- data.frame(formatted_table, DiffList_T[,"F"])
      names <- c(names, "F-Statistic:Time Trajectory");
      col_widths<-c(col_widths, 8);
    } }, silent=T)
  }else if ( any( grepl("logfc_", colnames(fData(eset))) )){
    formatted_table <- data.frame(formatted_table, fData(eset)[,grepl("logfc_", colnames(fData(eset)) )]); 
    names <- c(names, colnames(fData(eset))[grepl("logfc_", colnames(fData(eset)) )])
    col_widths <- c(col_widths, rep(8, sum(grepl("logfc_", colnames(fData(eset)) ) )))
  }
  
  if ( any( grepl("mummichogID_", colnames(fData(eset))) )){
    formatted_table <- data.frame(formatted_table, fData(eset)[,grepl("mummichogID_", colnames(fData(eset)) )]); 
    names <- c(names, colnames(fData(eset))[grepl("mummichogID_", colnames(fData(eset)) )])
    col_widths <- c(col_widths, rep(8, sum(grepl("mummichogID_", colnames(fData(eset)) ) )))
  }
  
  if ( "Uniprot_Function" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Uniprot_Function")]); 
    names <- c(names,"Uniprot_Function"); col_widths<-c(col_widths, 50);
  }
  if ( "Uniprot_Cellular_Location" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Uniprot_Cellular_Location")] ); 
    names <- c(names,"Uniprot_Cellular_Location"); col_widths<-c(col_widths, 50);
  }
  if ( "Uniprot_Disease" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Uniprot_Disease")] );
    names <- c(names,"Uniprot_Disease"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_biological_process" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_biological_process")] ); 
    names <- c(names,"GO_biological_process"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_molecular_function" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_molecular_function")] ); 
    names <- c(names,"GO_molecular_function"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_cellular_component" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_cellular_component")] ); 
    names <- c(names,"GO_cellular_component"); col_widths<-c(col_widths, 50);
  }
  if ( "GO_ID" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("GO_ID")] ); 
    names <- c(names,"GO_ID"); col_widths<-c(col_widths, 8);
  }
  if ( "ReactomeID" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("ReactomeID")] ); 
    names <- c(names,"ReactomeID"); col_widths<-c(col_widths, 8);
  } 
  if ( "KEGG_ID" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("KEGG_ID")] ); 
    names <- c(names,"KEGG_ID"); col_widths<-c(col_widths, 8);
  }
  if ( "identifier" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("identifier")] ); 
    names <- c(names,"identifier"); col_widths<-c(col_widths, 16);
  }
  
  if ( "Sequence" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Sequence")] ); 
    names <- c(names,"Sequence"); col_widths<-c(col_widths, 16);
  }
  
  if ( "Proteins" %in% colnames(fData(eset)) ){ 
    formatted_table <- data.frame(formatted_table, fData(eset)[,c("Proteins")] ); 
    names <- c(names,"Proteins"); col_widths<-c(col_widths, 16);
  }
  
  # Add phospho site probabilities and data
  try({ 
    if(grepl("Sites", data_format)) { 
      formatted_table <- data.frame(formatted_table, fData(eset)[,"Localization.prob"],
                                    fData(eset)[,grep("Probabilities", colnames(fData(eset)))],
                                    paste(fData(eset)[,"Amino.acid"],fData(eset)[,"Position"], sep=""),
                                    fData(eset)[,"Sequence.window"] ) 
      names <- c(names, "Localization Probability","Site Probabilities", "Amino Acid", "Peptide Sequence")
      col_widths <- c(col_widths, 16, 16, 16, 16)
    }
  })

  formatted_table<- data.frame(formatted_table, exprs(eset) );  
  names <- c(names,colnames(eset) );
  col_widths <- c(col_widths, rep(4, ncol(eset)));
  
  
  formatted_table <- data.frame(formatted_table, stringsAsFactors=FALSE)
  colnames(formatted_table) <-  make.unique(toupper(names));
  
  # Write data table to sheet
  writeDataTable(wb=wb, sheet=stName, x=formatted_table, xy=c("A",2), keepNA=FALSE, tableStyle="TableStyleLight1")
  if ("Link" %in% colnames(fData(eset)) ){ writeData(wb=wb, sheet=stName, x=links, xy=c( 1+(ncol(eset)+3 ) ,3)) }
  
  # Add heatmap color
  conditionalFormatting(wb=wb, sheet=stName, type="colourScale", cols=2:(1+ncol(eset)), rows=3:(2+nrow(eset)), style=mapcolor[c(1,4,7)])

  # Add color bars for fold change
  if(length(DiffList)>0) { i<-1; #for(i in 1:length(DiffList)){
    conditionalFormatting(wb=wb, sheet=stName, cols= ( 3+col_index+(3*(i-1))) , rows=3:(3+nrow(limmaFit)), type="databar", style=c("blue", "red"), showvalue=FALSE)
  } #}
  
  # Rotate text for sample names
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"),
             rows=2, cols=c+1)
  }
  for (c in 1:ncol(eset)){
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(fgFill=sampleCols[c], textRotation=90, halign="center", valign="top"),
             rows=2, cols=(c + ( length(names) - ncol(eset) )  ) )
  }

  # Merge cells and add Intensity title
  mergeCells(wb=wb, sheet=stName, rows=1, cols=2:(1+ncol(eset)) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity, Row Z Score", xy=c(2,1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=2, stack=TRUE)
  
  # Add final intensity title
  mergeCells(wb=wb, sheet=stName, rows=1, cols=( length(names) : ( length(names) - ncol(eset) + 1 )) )
  writeData(wb=wb, sheet=stName, x="Normalized Log2 Intensity", xy=c(( length(names) - ncol(eset) + 1 ),1))
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=( length(names) - ncol(eset) + 1 ), stack=TRUE)
  
  # Add contrast titles
  if(class(contrast_strings)=="NULL"){
    contrast_strings <- colnames(limmaFit$coefficients)
  }
  if(length(DiffList)>0) { for(i in 1:length(DiffList)){
    startCol<- col_index + 1 + (3*(i-1))
    mergeCells(wb=wb, sheet=stName, rows=1, cols=startCol:(startCol+2) )
    writeData(wb=wb, sheet=stName, x=contrast_strings[i], xy=c(startCol,1))
    addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", halign="center"), rows=1, cols=startCol, stack=TRUE)
  } }
  
  # Freeze columns/rows and bold column titles
  freezePane(wb=wb, sheet=stName, firstActiveRow=3, firstActiveCol=2)
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(textDecoration="bold", fontColour='black'),
           rows=1:2, cols=1:length(names),gridExpand=TRUE, stack=TRUE)   
  
  # Don't treat gene names as dates
  addStyle(wb=wb, sheet=stName, style=openxlsx::createStyle(numFmt="TEXT"), rows=3:(2+nrow(eset)), cols=1, gridExpand=TRUE )
  
  # Set row heights and column widths
  setRowHeights(wb=wb, sheet=stName, rows=2, heights=100)
  setColWidths(wb=wb, sheet=stName, cols=1:ncol(formatted_table), widths= col_widths )
}

