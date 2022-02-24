#-------------------------------------------------
#' Save Files
#'
#' This function saves txt files for subsequent use
#' 
#' @param eset an ExpressionSet object with omics data
#' @param type Type of data to be processed, or a name for the Omics set
#' @param outputpath output file path
#' @param limmaRes topTable output
#' @param contrast_name name of the differential contrast
#' @param outputcontrastpath output file path for differntial
#'
#' 
#' @examples
#' 
#' @export
saveFiles <- function(data, type, outputpath=output_files_path, subset=FALSE, saveRDS=TRUE,
                      limmaRes=FALSE, contrast_name=FALSE, outputcontrastpath=output_contrast_path_files){
  # Save eSet RDS  
  if( saveRDS ){
    output_filename <- file.path(outputpath, paste("Data_",type,".RDS", sep=""));
    saveRDS(data, file=output_filename); 
  }
  
  # Save expression matrix
  select_cols<- c("Gene", "feature_identifier", "Protein", "mz","rt","Adduct", "identifier", "Metabolite.name", "KEGG")
  select_cols<- select_cols[which(select_cols %in% colnames(fData(data[["eSet"]])) )]
  output_table <- cbind(fData(data[["eSet"]])[,select_cols],exprs(data[["eSet"]]))
  output_filename <- file.path(outputpath,paste("Expression_matrix_",type,".txt", sep=''));
  write.table(x= output_table, file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
  
  if( class(subset)!="logical" ){ 
    for( k in 1:length(subset) ){ suppressWarnings({ try({
      output_table <-output_table[subset[[k]],]
      output_filename <- file.path(outputpath,paste("Expression_matrix_",type,names(subset)[k],".txt", sep=''));
      write.table(x= output_table, file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
    }) }) }  }
  
  
  # If differential analysis, save ranked list
  if ( class(limmaRes) != "logical") {
    # Save expression matrix
    output_filename <- file.path(outputcontrastpath,paste("DE_matrix_",type,"_", gsub("-","_",contrast_name),".txt", sep=''));
    select_cols<- c("Gene", "feature_identifier", "Protein", "mz","rt","Adduct", "identifier", "Metabolite.name", "KEGG",
                    "P.Value", "adj.P.Val", "logFC", "AveExpr","F","B")
    select_cols<- select_cols[which(select_cols %in% colnames(limmaRes))]
    write.table(x= limmaRes[,select_cols], file=output_filename, sep='\t',row.names=F, col.names=TRUE, quote=FALSE);
    
    if("Gene"%in% colnames(limmaRes)){
      # ranked list with direction
      output_filename <- file.path(outputcontrastpath,paste("GSEA_",type,"_", gsub("-","_",contrast_name),".rnk", sep=''));
      if( "P.Value" %in% colnames(limmaRes) ){
        ranked <- cbind(limmaRes[,"Gene"],sign(limmaRes[,"logFC"]) * -log10(limmaRes[,"P.Value"]))
      } else {
        ranked <- cbind(limmaRes[,"Gene"],limmaRes[,"logFC"] )
      }
      colnames(ranked)<-c("GeneName", "rank")
      ranked <- ranked[ranked[,"GeneName"]!="", ]
      ranked <- ranked[order(abs(as.numeric(ranked[,"rank"])), decreasing=TRUE),]
      ranked <- ranked[!duplicated(ranked[,"GeneName"]),]
      ranked <- ranked[order(as.numeric(ranked[,"rank"]), decreasing=TRUE),]
      write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
      
      # ranked list without direction
      output_filename <- file.path(outputcontrastpath,paste("GSEA_",type,"_", gsub("-","_",contrast_name),"_AbsVal.rnk", sep=''));
      ranked[,"rank"] <- abs(as.numeric(ranked[,"rank"]))
      write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
      
      
      # File for network analysis
      output_filename <- file.path(outputcontrastpath,paste("Network_",type, gsub("-","_",contrast_name),".txt", sep=''));
      if( "adj.P.Val" %in% colnames(limmaRes) ) {
        ranked <- cbind(limmaRes[,"Gene"],limmaRes["P.Value"],sign(limmaRes[,"logFC"]) * -log10(limmaRes[,"adj.P.Val"]))
        colnames(ranked)<-c("GeneName", "pValue", "rank")
        ranked <- ranked[order(as.numeric(ranked[,"pValue"]), decreasing=FALSE),]
      } else {
        ranked <- cbind(limmaRes[,"Gene"],limmaRes[,"logFC"] )
        colnames(ranked)<-c("GeneName", "logFC")
        ranked <- ranked[order(as.numeric(ranked[,"logFC"]), decreasing=FALSE),]
      }
      ranked <- ranked[ranked[,"GeneName"]!="", ]
      ranked <- ranked[!duplicated(ranked[,"GeneName"]),]
      ranked[,"DataSet"] <- type;
      write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
    }
    if( "mz"%in% colnames(limmaRes) & "rt"%in% colnames(limmaRes) & "P.Value" %in% colnames(limmaRes) & "t" %in% colnames(limmaRes)){
      output_filename <- file.path(outputcontrastpath,paste("Mummichog_",type, gsub("-","_",contrast_name),".txt", sep=''));
      ranked <- limmaRes[,c("mz", "rt", "P.Value", "t" ) ]
      colnames(ranked)<-c("mz", "rtime", "p-value", "t-score")
      write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);  
      
      # File for network analysis
      output_filename <- file.path(outputcontrastpath,paste("Network_",type, gsub("-","_",contrast_name),".txt", sep=''));
      ranked <- limmaRes[,c("feature_identifier","mz", "rt", "P.Value", "t" ) ]
      colnames(ranked)<-c("feature_identifier","mz", "rtime", "p-value", "t-score")
      ranked[,"DataSet"] <- type;
      ranked[,"Type"] <- "mz"
      ranked[,"Display_Name"] <- paste("mz", round(ranked[,"mz"], 4), sep="")
      write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);
      
      # MZ feature to HMDB interactions
      output_filename <- file.path(outputcontrastpath,paste("Network_",type, gsub("-","_",contrast_name),"_interactions.txt", sep=''));
      
      #ranked <- data.frame(Metabolite= gsub(" ","",gsub("\\)","",gsub("c\\(","",unlist(strsplit(limmaRes[1,"identifier"], ","))))) )
      #ranked[,"feature_identifier"] <- limmaRes[1,"feature_identifier"]
      #ranked[,"Edge_Type"] <- "mz_met"
      ranked=data.frame()
      for( i in 1:nrow(limmaRes) ){ try({
        if( nrow(data.frame(Metabolite= gsub(" ","",gsub("\\)","",gsub("c\\(","",unlist(strsplit(as.character(limmaRes[i,"identifier"]), ","))))) ))>0 ){
        toAdd <- data.frame(Metabolite= gsub(" ","",gsub("\\)","",gsub("c\\(","",unlist(strsplit(as.character(limmaRes[i,"identifier"]), ","))))) )
        toAdd[,"feature_identifier"] <- limmaRes[i,"feature_identifier"]
        toAdd[,"Edge_Type"] <- "mz_met"
        toAdd[,"Display_Name"] <- paste("mz", round(limmaRes[i,"mz"], 4), sep="")
        ranked <- rbind(ranked, toAdd)
      } }) }      
      write.table(x=ranked,file=output_filename, sep='\t',row.names=FALSE, col.names=TRUE, quote=FALSE);      
    }
  }
}

